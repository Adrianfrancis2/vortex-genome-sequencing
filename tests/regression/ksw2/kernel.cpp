#include <vx_spawn2.h>
#include "common.h"

// Both kernels implement the recurrence and boundaries of ref/gotoh.c (KSW2's
// ksw_extz): a gap of length l costs gapo + l*gape, out-of-band cells are -inf,
// and ties for the max go to the smallest anti-diagonal t+q, then the smallest t.

static inline int32_t max2(int32_t a, int32_t b) { return a > b ? a : b; }

struct scoring_t {
  int32_t match, mismatch, gapo, gape, qe, w;
};

// H on row 0 at column k (k query bases skipped), as in gotoh.c.
static inline int32_t top_h(const scoring_t& sc, int32_t k) {
  return k == 0 ? 0 : (k <= sc.w ? -(sc.gapo + sc.gape * k) : KSW2_NEG_INF);
}

// H on column 0 at row k (k target bases skipped), as in gotoh.c.
static inline int32_t left_h(const scoring_t& sc, int32_t k) {
  return k - 1 <= sc.w ? -(sc.gapo + sc.gape * k) : KSW2_NEG_INF;
}

static inline int32_t iabs(int32_t x) { return x < 0 ? -x : x; }

// One warp per pair. Anti-diagonal r holds the cells t + q = r, stored by t.
static void align_diag(const kernel_arg_t* arg, const scoring_t& sc_in) {
  auto pairs   = reinterpret_cast<const pair_t*>(arg->pairs_addr);
  auto seq     = reinterpret_cast<const uint8_t*>(arg->seq_addr);
  auto scratch = reinterpret_cast<int32_t*>(arg->scratch_addr);
  auto results = reinterpret_cast<result_t*>(arg->result_addr);

  uint32_t p = blockIdx.x;
  if (p >= arg->num_pairs)
    return;
  int32_t lane = threadIdx.x;
  int32_t nt   = blockDim.x;

  const pair_t pr = pairs[p];
  const uint8_t* qs = seq + pr.q_off;
  const uint8_t* ts = seq + pr.t_off;
  int32_t qlen = pr.qlen, tlen = pr.tlen;

  scoring_t sc = sc_in;
  if (sc.w < 0)
    sc.w = tlen > qlen ? tlen : qlen;

  int32_t* s  = scratch + pr.scratch_off;
  int32_t* H2 = s;            // anti-diagonal r-2
  int32_t* H1 = s + tlen;     // r-1
  int32_t* H0 = s + 2 * tlen; // r
  int32_t* E1 = s + 3 * tlen;
  int32_t* E0 = s + 4 * tlen;
  int32_t* F1 = s + 5 * tlen;
  int32_t* F0 = s + 6 * tlen;
  int32_t* red = s + 7 * tlen;

  int32_t best_h = 0, best_t = -1, best_q = -1;

  for (int32_t r = 0; r <= qlen + tlen - 2; ++r) {
    int32_t st = r - qlen + 1 > 0 ? r - qlen + 1 : 0;
    int32_t en = r < tlen - 1 ? r : tlen - 1;
    for (int32_t t = st + lane; t <= en; t += nt) {
      int32_t q = r - t;
      if (iabs(t - q) > sc.w) {
        H0[t] = E0[t] = F0[t] = KSW2_NEG_INF;
        if (t == tlen - 1 && q == qlen - 1)
          results[p].score = KSW2_NEG_INF;
        continue;
      }
      int32_t diag_h = t == 0 ? top_h(sc, q)
                     : q == 0 ? left_h(sc, t)
                     : H2[t - 1];
      int32_t up_h = t == 0 ? top_h(sc, q + 1) : H1[t - 1];
      int32_t up_e = t == 0 ? KSW2_NEG_INF : E1[t - 1];
      int32_t lf_h = q == 0 ? left_h(sc, t + 1) : H1[t];
      int32_t lf_f = q == 0 ? KSW2_NEG_INF : F1[t];

      int32_t e = max2(up_h - sc.qe, up_e - sc.gape);
      int32_t f = max2(lf_h - sc.qe, lf_f - sc.gape);
      int32_t sub = ts[t] == qs[q] ? sc.match : -sc.mismatch;
      int32_t h = max2(diag_h + sub, max2(e, f));
      H0[t] = h;
      E0[t] = e;
      F0[t] = f;

      // A lane visits r, then t, in increasing order, so keeping the first
      // of equal scores applies the tie rule.
      if (h > best_h) {
        best_h = h;
        best_t = t;
        best_q = q;
      }
      if (t == tlen - 1 && q == qlen - 1)
        results[p].score = h;
    }
    __syncthreads();
    int32_t* tmp = H2;
    H2 = H1; H1 = H0; H0 = tmp;
    tmp = E1; E1 = E0; E0 = tmp;
    tmp = F1; F1 = F0; F0 = tmp;
  }

  red[3 * lane + 0] = best_h;
  red[3 * lane + 1] = best_t;
  red[3 * lane + 2] = best_q;
  __syncthreads();

  if (lane == 0) {
    int32_t mh = 0, mt = -1, mq = -1;
    for (int32_t l = 0; l < nt; ++l) {
      int32_t h = red[3 * l], t = red[3 * l + 1], q = red[3 * l + 2];
      if (h <= 0)
        continue;
      if (h > mh || (h == mh && (t + q < mt + mq || (t + q == mt + mq && t < mt)))) {
        mh = h;
        mt = t;
        mq = q;
      }
    }
    results[p].max   = mh;
    results[p].max_t = mt;
    results[p].max_q = mq;
  }
}

// One thread per pair, rows of the target over columns of the query; this is
// gotoh.c with its row buffers in global memory.
static void align_thread(const kernel_arg_t* arg, const scoring_t& sc_in) {
  auto pairs   = reinterpret_cast<const pair_t*>(arg->pairs_addr);
  auto seq     = reinterpret_cast<const uint8_t*>(arg->seq_addr);
  auto scratch = reinterpret_cast<int32_t*>(arg->scratch_addr);
  auto results = reinterpret_cast<result_t*>(arg->result_addr);

  uint32_t p = blockIdx.x * blockDim.x + threadIdx.x;
  if (p >= arg->num_pairs)
    return;

  const pair_t pr = pairs[p];
  const uint8_t* qs = seq + pr.q_off;
  const uint8_t* ts = seq + pr.t_off;
  int32_t qlen = pr.qlen, tlen = pr.tlen;

  scoring_t sc = sc_in;
  if (sc.w < 0)
    sc.w = tlen > qlen ? tlen : qlen;

  int32_t* Hp = scratch + pr.scratch_off;
  int32_t* Ep = Hp + (qlen + 1);
  int32_t* Hc = Ep + (qlen + 1);
  int32_t* Ec = Hc + (qlen + 1);

  for (int32_t J = 0; J <= qlen; ++J) {
    Hp[J] = top_h(sc, J);
    Ep[J] = KSW2_NEG_INF;
  }

  int32_t mh = 0, mt = -1, mq = -1;
  for (int32_t I = 1; I <= tlen; ++I) {
    Hc[0] = left_h(sc, I);
    Ec[0] = KSW2_NEG_INF;
    int32_t f = KSW2_NEG_INF;
    uint8_t tb = ts[I - 1];
    for (int32_t J = 1; J <= qlen; ++J) {
      if (iabs(I - J) > sc.w) {
        Hc[J] = Ec[J] = f = KSW2_NEG_INF;
        continue;
      }
      int32_t sub = tb == qs[J - 1] ? sc.match : -sc.mismatch;
      int32_t e = max2(Hp[J] - sc.qe, Ep[J] - sc.gape);
      f = max2(Hc[J - 1] - sc.qe, f - sc.gape);
      int32_t h = max2(Hp[J - 1] + sub, max2(e, f));
      Hc[J] = h;
      Ec[J] = e;

      int32_t t = I - 1, q = J - 1;
      if (h > mh) {
        mh = h; mt = t; mq = q;
      } else if (h == mh && mh > 0
                 && (t + q < mt + mq || (t + q == mt + mq && t < mt))) {
        mt = t; mq = q;
      }
    }
    int32_t* tmp;
    tmp = Hp; Hp = Hc; Hc = tmp;
    tmp = Ep; Ep = Ec; Ec = tmp;
  }

  results[p].score = iabs(tlen - qlen) <= sc.w ? Hp[qlen] : KSW2_NEG_INF;
  results[p].max   = mh;
  results[p].max_t = mt;
  results[p].max_q = mq;
}

__kernel void kernel_main(kernel_arg_t* __UNIFORM__ arg) {
  scoring_t sc;
  sc.match    = arg->match;
  sc.mismatch = arg->mismatch;
  sc.gapo     = arg->gapo;
  sc.gape     = arg->gape;
  sc.qe       = arg->gapo + arg->gape;
  sc.w        = arg->w;
  if (arg->mode == KSW2_MODE_DIAG)
    align_diag(arg, sc);
  else
    align_thread(arg, sc);
}
