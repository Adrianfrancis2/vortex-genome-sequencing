#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <vector>
#include <vortex2.h>
#include "common.h"

#ifndef KSW2_DATA_DIR
#define KSW2_DATA_DIR "."
#endif

#define RT_CHECK(_expr)                                         \
   do {                                                         \
     int _ret = _expr;                                          \
     if (0 == _ret)                                             \
       break;                                                   \
     printf("Error: '%s' returned %d!\n", #_expr, (int)_ret);  \
     cleanup();                                                 \
     exit(-1);                                                  \
   } while (false)

// Scoring fixed for the project: minimap2's default single-affine values.
static const int32_t MATCH = 2, MISMATCH = 4, GAPO = 4, GAPE = 2;

struct host_pair_t {
  int id;
  std::string set;
  std::string q, t;
};

const char* kernel_file = "kernel.vxbin";
std::string pairs_path  = KSW2_DATA_DIR "/pairs.txt";
std::string golden_path = KSW2_DATA_DIR "/golden.txt";
std::string set_prefix;
uint32_t mode = KSW2_MODE_DIAG;
int32_t band = -1;
int max_pairs = -1;

vx_device_h device = nullptr;
vx_buffer_h pairs_buffer = nullptr;
vx_buffer_h seq_buffer = nullptr;
vx_buffer_h scratch_buffer = nullptr;
vx_buffer_h result_buffer = nullptr;
vx_queue_h  queue = nullptr;
vx_module_h module_ = nullptr;
vx_kernel_h kernel = nullptr;
kernel_arg_t kernel_arg = {};

static void show_usage() {
  std::cout << "KSW2 extension alignment on Vortex." << std::endl;
  std::cout << "Usage: [-m diag|thread] [-w band] [-s set_prefix] [-n max_pairs]"
               " [-p pairs.txt] [-g golden.txt] [-k kernel] [-h]" << std::endl;
  std::cout << "  -w applies a band; -g must then name a golden file made with the same band" << std::endl;
}

static void parse_args(int argc, char **argv) {
  int c;
  while ((c = getopt(argc, argv, "m:w:s:n:p:g:k:h")) != -1) {
    switch (c) {
    case 'm':
      if (std::string(optarg) == "diag") mode = KSW2_MODE_DIAG;
      else if (std::string(optarg) == "thread") mode = KSW2_MODE_THREAD;
      else { show_usage(); exit(-1); }
      break;
    case 'w': band = atoi(optarg); break;
    case 's': set_prefix = optarg; break;
    case 'n': max_pairs = atoi(optarg); break;
    case 'p': pairs_path = optarg; break;
    case 'g': golden_path = optarg; break;
    case 'k': kernel_file = optarg; break;
    case 'h': show_usage(); exit(0);
    default: show_usage(); exit(-1);
    }
  }
}

void cleanup() {
  if (device) {
    if (pairs_buffer)   vx_buffer_release(pairs_buffer);
    if (seq_buffer)     vx_buffer_release(seq_buffer);
    if (scratch_buffer) vx_buffer_release(scratch_buffer);
    if (result_buffer)  vx_buffer_release(result_buffer);
    if (kernel)  vx_kernel_release(kernel);
    if (module_) vx_module_release(module_);
    if (queue)   vx_queue_release(queue);
    vx_device_dump_perf(device, stdout);
    vx_device_release(device);
  }
}

static bool encode(const std::string& s, std::vector<uint8_t>& out) {
  for (char ch : s) {
    switch (ch) {
    case 'A': out.push_back(0); break;
    case 'C': out.push_back(1); break;
    case 'G': out.push_back(2); break;
    case 'T': out.push_back(3); break;
    default: return false;
    }
  }
  return true;
}

static std::vector<host_pair_t> load_pairs(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Error: cannot open " << path << std::endl;
    exit(-1);
  }
  std::vector<host_pair_t> pairs;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    host_pair_t p;
    size_t qlen, tlen;
    if (!(ls >> p.id >> p.set >> qlen >> tlen >> p.q >> p.t)
        || p.q.size() != qlen || p.t.size() != tlen) {
      std::cerr << "Error: bad line in " << path << ": " << line << std::endl;
      exit(-1);
    }
    if (!set_prefix.empty() && p.set.compare(0, set_prefix.size(), set_prefix) != 0)
      continue;
    pairs.push_back(p);
    if (max_pairs >= 0 && (int)pairs.size() >= max_pairs)
      break;
  }
  return pairs;
}

static std::unordered_map<int, result_t> load_golden(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Error: cannot open " << path << std::endl;
    exit(-1);
  }
  std::unordered_map<int, result_t> golden;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    int id;
    result_t r;
    if (ls >> id >> r.score >> r.max >> r.max_q >> r.max_t)
      golden[id] = r;
  }
  return golden;
}

// Cells inside the band, which is what the kernels compute.
static uint64_t band_cells(int32_t qlen, int32_t tlen, int32_t w) {
  if (w < 0)
    return (uint64_t)qlen * tlen;
  uint64_t n = 0;
  for (int32_t t = 0; t < tlen; ++t) {
    int32_t lo = t - w > 0 ? t - w : 0;
    int32_t hi = t + w < qlen - 1 ? t + w : qlen - 1;
    if (hi >= lo)
      n += hi - lo + 1;
  }
  return n;
}

int main(int argc, char *argv[]) {
  parse_args(argc, argv);

  auto hpairs = load_pairs(pairs_path);
  auto golden = load_golden(golden_path);
  uint32_t num_pairs = hpairs.size();
  if (num_pairs == 0) {
    std::cerr << "Error: no pairs selected" << std::endl;
    return -1;
  }

  std::cout << "open device connection" << std::endl;
  RT_CHECK(vx_device_open(0, &device));

  vx_queue_info_t qi = { sizeof(qi), nullptr, VX_QUEUE_PRIORITY_NORMAL, 0 };
  RT_CHECK(vx_queue_create(device, &qi, &queue));

  uint64_t num_threads;
  RT_CHECK(vx_device_query(device, VX_CAPS_NUM_THREADS, &num_threads));

  // pack sequences and per-pair descriptors
  std::vector<pair_t> h_pairs(num_pairs);
  std::vector<uint8_t> h_seq;
  uint64_t scratch_words = 0, cells = 0;
  for (uint32_t i = 0; i < num_pairs; ++i) {
    auto& hp = hpairs[i];
    auto& dp = h_pairs[i];
    dp.qlen = hp.q.size();
    dp.tlen = hp.t.size();
    dp.q_off = h_seq.size();
    dp.t_off = h_seq.size() + dp.qlen;
    if (!encode(hp.q, h_seq) || !encode(hp.t, h_seq)) {
      std::cerr << "Error: pair " << hp.id << " has a non-ACGT base" << std::endl;
      return -1;
    }
    dp.scratch_off = scratch_words;
    scratch_words += ksw2_scratch_words(dp.qlen, dp.tlen, num_threads);
    cells += band_cells(dp.qlen, dp.tlen, band);
  }
  while (h_seq.size() % 4)
    h_seq.push_back(0);

  std::cout << "ksw2: " << num_pairs << " pairs, mode="
            << (mode == KSW2_MODE_DIAG ? "diag" : "thread")
            << ", w=" << (band < 0 ? std::string("full") : std::to_string(band))
            << ", threads/warp=" << num_threads << std::endl;
  std::cout << "cells=" << cells << std::endl;

  uint64_t pairs_size   = num_pairs * sizeof(pair_t);
  uint64_t seq_size     = h_seq.size();
  uint64_t scratch_size = scratch_words * sizeof(int32_t);
  uint64_t result_size  = num_pairs * sizeof(result_t);

  std::cout << "allocate device memory" << std::endl;
  RT_CHECK(vx_buffer_create(device, pairs_size, VX_MEM_READ, &pairs_buffer));
  RT_CHECK(vx_buffer_address(pairs_buffer, &kernel_arg.pairs_addr));
  RT_CHECK(vx_buffer_create(device, seq_size, VX_MEM_READ, &seq_buffer));
  RT_CHECK(vx_buffer_address(seq_buffer, &kernel_arg.seq_addr));
  RT_CHECK(vx_buffer_create(device, scratch_size, VX_MEM_READ_WRITE, &scratch_buffer));
  RT_CHECK(vx_buffer_address(scratch_buffer, &kernel_arg.scratch_addr));
  RT_CHECK(vx_buffer_create(device, result_size, VX_MEM_WRITE, &result_buffer));
  RT_CHECK(vx_buffer_address(result_buffer, &kernel_arg.result_addr));

  std::cout << "upload source buffers" << std::endl;
  RT_CHECK(vx_enqueue_write(queue, pairs_buffer, 0, h_pairs.data(), pairs_size, 0, nullptr, nullptr));
  RT_CHECK(vx_enqueue_write(queue, seq_buffer, 0, h_seq.data(), seq_size, 0, nullptr, nullptr));

  std::cout << "load kernel module" << std::endl;
  RT_CHECK(vx_module_load_file(device, kernel_file, &module_));
  RT_CHECK(vx_module_get_kernel(module_, "main", &kernel));

  kernel_arg.num_pairs = num_pairs;
  kernel_arg.mode      = mode;
  kernel_arg.w         = band;
  kernel_arg.match     = MATCH;
  kernel_arg.mismatch  = MISMATCH;
  kernel_arg.gapo      = GAPO;
  kernel_arg.gape      = GAPE;

  // diag: a block is one warp and owns one pair; thread: one pair per thread
  uint32_t block_dim = num_threads;
  uint32_t grid_dim  = mode == KSW2_MODE_DIAG ? num_pairs
                                              : (num_pairs + block_dim - 1) / block_dim;

  std::cout << "launch kernel" << std::endl;
  vx_event_h launch_ev = nullptr, read_ev = nullptr;
  {
    vx_launch_info_t li = {};
    li.struct_size  = sizeof(li);
    li.kernel       = kernel;
    li.args_host    = &kernel_arg;
    li.args_size    = sizeof(kernel_arg);
    li.ndim         = 1;
    li.grid_dim[0]  = grid_dim;
    li.block_dim[0] = block_dim;
    RT_CHECK(vx_enqueue_launch(queue, &li, 0, nullptr, &launch_ev));
  }

  std::cout << "download results" << std::endl;
  std::vector<result_t> h_results(num_pairs);
  RT_CHECK(vx_enqueue_read(queue, h_results.data(), result_buffer, 0, result_size, 1, &launch_ev, &read_ev));

  std::cout << "wait for completion" << std::endl;
  RT_CHECK(vx_event_wait_value(read_ev, 1, VX_TIMEOUT_INFINITE));
  vx_event_release(read_ev);
  vx_event_release(launch_ev);

  std::cout << "verify result" << std::endl;
  int errors = 0, missing = 0;
  for (uint32_t i = 0; i < num_pairs; ++i) {
    auto& hp = hpairs[i];
    auto it = golden.find(hp.id);
    if (it == golden.end()) {
      ++missing;
      continue;
    }
    const result_t& g = it->second;
    const result_t& r = h_results[i];
    if (r.score != g.score || r.max != g.max || r.max_q != g.max_q || r.max_t != g.max_t) {
      if (errors < 20) {
        printf("*** pair %d (%s): got score=%d max=%d at q=%d t=%d, golden score=%d max=%d at q=%d t=%d\n",
               hp.id, hp.set.c_str(), r.score, r.max, r.max_q, r.max_t,
               g.score, g.max, g.max_q, g.max_t);
      }
      ++errors;
    }
  }

  std::cout << "cleanup" << std::endl;
  cleanup();

  if (missing != 0) {
    std::cout << missing << " pairs have no golden entry" << std::endl;
    std::cout << "FAILED!" << std::endl;
    return 1;
  }
  if (errors != 0) {
    std::cout << "Found " << std::dec << errors << " mismatching pairs!" << std::endl;
    std::cout << "FAILED!" << std::endl;
    return 1;
  }

  std::cout << "PASSED! " << num_pairs << "/" << num_pairs << " pairs match golden" << std::endl;
  return 0;
}
