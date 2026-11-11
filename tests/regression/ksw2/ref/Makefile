# Host-only KSW2 reference; not part of the Vortex regression build.
KSW2_DIR := ../../../../third_party/ksw2

CC      ?= gcc
CFLAGS  := -O2 -Wall -msse4.1 -I$(KSW2_DIR)

SRCS := ksw2_ref.c gotoh.c $(KSW2_DIR)/ksw2_extz.c $(KSW2_DIR)/ksw2_extz2_sse.c

all: ksw2_ref

ksw2_ref: $(SRCS) gotoh.h
	$(CC) $(CFLAGS) $(SRCS) -o $@

clean:
	rm -f ksw2_ref

.PHONY: all clean
