# Compiler and flags
CXX = g++
MPICXX = mpic++
CXXFLAGS = -Wall -Wextra -O2
OPENMPFLAGS = -fopenmp

# Source directory
SRC_DIR = .
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
BINS := $(patsubst $(SRC_DIR)/%.cpp,%,$(SRCS))

.PHONY: all clean

all: $(BINS)

# Build rule
%: $(SRC_DIR)/%.cpp
	@if echo "$@" | grep -q "openmp"; then \
		$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -o $@ $<; \
	elif echo "$@" | grep -q "mpi"; then \
		$(MPICXX) $(CXXFLAGS) -o $@ $<; \
	else \
		$(CXX) $(CXXFLAGS) -o $@ $<; \
	fi

clean:
	rm -f $(BINS) *.o *~ core
