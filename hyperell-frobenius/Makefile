# Build `euler`: Frobenius matrices for one curve at many primes, in parallel.
#
# Requirements: a C++17 compiler, NTL and GMP. If NTL is not under
# /usr/local, either copy config.mk.example to config.mk and edit it, or
# pass the settings on the command line:
#
#   make NTL_INC=$HOME/ntl/include NTL_LIB=$HOME/ntl/src NTL_LINK=$HOME/ntl/src/ntl.a
#
# Threads need NTL built with NTL_THREADS and NTL_THREAD_BOOST; without them
# use one thread (`euler 1`), which still batches.

BUILD    ?= build
CXX      ?= g++
CC       ?= cc
PYTHON   ?= python3
CXXFLAGS ?= -O2 -std=c++17
RFOREST_CFLAGS ?= -O3 -std=gnu11 -fomit-frame-pointer -funroll-loops
NTL_INC  ?= /usr/local/include
NTL_LIB  ?= /usr/local/lib
NTL_LINK ?= -lntl
NTL_RPATH ?= -Wl,-rpath,"$(NTL_LIB)"

GMP_INC  ?= $(NTL_INC)
GMP_LIB  ?= $(NTL_LIB)
CONFIG_MK ?= config.mk

-include $(CONFIG_MK)

OBJS := $(BUILD)/even_degree.o $(BUILD)/average.o
RFOREST_SRC := $(filter-out rforest/test_rforest.c,$(wildcard rforest/*.c rforest/fft62/*.c))
RFOREST_OBJS := $(patsubst rforest/%.c,$(BUILD)/rforest/%.o,$(RFOREST_SRC))
RFOREST_HEADERS := $(wildcard rforest/*.h rforest/fft62/*.h)

all: $(BUILD)/euler

$(BUILD):
	@mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.cpp | $(BUILD)
	$(CXX) -I"$(NTL_INC)" -I"$(GMP_INC)" -Isrc $(CPPFLAGS) $(CXXFLAGS) -pthread -c $< -o $@

$(BUILD)/average.o: src/average.h src/even_degree.h rforest/rforest.h
$(BUILD)/even_degree.o: src/even_degree.h

$(BUILD)/rforest/%.o: rforest/%.c $(RFOREST_HEADERS)
	@mkdir -p "$(@D)"
	$(CC) -I"$(GMP_INC)" -Irforest $(CPPFLAGS) $(RFOREST_CFLAGS) -c $< -o $@

$(BUILD)/euler: euler.cpp $(OBJS) $(RFOREST_OBJS) src/average.h | $(BUILD)
	$(CXX) -I"$(NTL_INC)" -I"$(GMP_INC)" -Isrc $(CPPFLAGS) $(CXXFLAGS) -pthread euler.cpp $(OBJS) $(RFOREST_OBJS) \
	    $(LDFLAGS) -L"$(NTL_LIB)" -L"$(GMP_LIB)" \
	    $(NTL_RPATH) -Wl,-rpath,"$(GMP_LIB)" $(NTL_LINK) -lgmpxx -lgmp -lm -lpthread -o $@

# Check full odd/even matrices and worker counts without importing Sage.
check: $(BUILD)/euler
	$(PYTHON) ../hyperell_regulator/frobenius.py --check --binary $(BUILD)/euler

clean:
	rm -rf $(BUILD)

.PHONY: all check clean
