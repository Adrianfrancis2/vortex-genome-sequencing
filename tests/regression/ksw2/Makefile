ROOT_DIR := $(realpath ../../..)
include $(ROOT_DIR)/config.mk

PROJECT := ksw2

SRC_DIR := $(VORTEX_HOME)/tests/regression/$(PROJECT)

SRCS := $(SRC_DIR)/main.cpp

VX_SRCS := $(SRC_DIR)/kernel.cpp

CXXFLAGS += -DKSW2_DATA_DIR='"$(SRC_DIR)/ref/data"'

# The regression suite runs the edge cases only; the full set takes minutes.
OPTS ?= -s edge

KERNEL_LIB := vortex2

include ../common.mk
