#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
MPI         = no
OPENMP      = yes
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
PAPI        = no
BENCHMARK   = no
PRECISION   = double

#===============================================================================
# Source Code List
#===============================================================================

case = models/c5g7/c5g7-cmfd.cpp

source = \
../src/Cell.cpp \
../src/Cmfd.cpp \
../src/CPUSolver.cpp \
../src/ExpEvaluator.cpp \
../src/Geometry.cpp \
../src/LocalCoords.cpp \
../src/log.cpp \
../src/Material.cpp \
../src/Point.cpp \
../src/PolarQuad.cpp \
../src/Solver.cpp \
../src/Surface.cpp \
../src/Timer.cpp \
../src/Track.cpp \
../src/TrackGenerator.cpp \
../src/Universe.cpp

cases = \
models/gradients/one-directional/one-directional-gradient.cpp \
models/gradients/two-directional/two-directional-gradient.cpp \
models/homogeneous/homogeneous-one-group.cpp \
models/c5g7/c5g7.cpp \
models/c5g7/c5g7-cmfd.cpp

#===============================================================================
# Sets Flags
#===============================================================================

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = g++
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icpc
  source += ../src/VectorizedSolver.cpp
endif

# IBM XL Compiler
ifeq ($(COMPILER),ibm)
  CC = mpic++
endif

# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  CC = bgxlc++
endif

# Standard Flags
CFLAGS := -std=c++11

# Linker Flags
LDFLAGS = -lm

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
endif

# Precision Flags
ifeq ($(PRECISION),single)
  CFLAGS += -DFP_PRECISION=float
  CFLAGS += -DSINGLE
endif
ifeq ($(PRECISION),double)
  CFLAGS += -DFP_PRECISION=double
  CFLAGS += -DDOUBLE
endif

# Python include flag
CFLAGS += -I/usr/include/python2.7

# Numpy include flag
CFLAGS += -I/usr/lib/python2.7/dist-packages/numpy/core/include 

# Vector Flags
CFLAGS += -DVEC_LENGTH=8
CFLAGS += -DVEC_ALIGNMENT=16

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
ifeq ($(COMPILER),gnu)
  CFLAGS += -Ofast -ffast-math -ftree-vectorize -msse2 -DGNU
endif
ifeq ($(COMPILER),intel)
  CFLAGS += -O3 -xhost -ansi-alias -no-prec-div -DINTEL
endif
ifeq ($(COMPILER),ibm)
  CFLAGS += -O5 -qhot -qsimd=auto -qalias=ansi:allptrs -qarch=qp -DIBM
endif
endif

# PAPI source (you may need to provide -I and -L pointing
# to PAPI depending on your installation
ifeq ($(PAPI),yes)
  CFLAGS += -DPAPI
  LDFLAGS += -lpapi
  OPENMP = yes
endif

# MPI
ifeq ($(MPI),yes)
  CC = mpic++
  CFLAGS += -DMPI
endif

# OpenMP
ifeq ($(OPENMP),yes)
ifeq ($(COMPILER), gnu)
  CFLAGS += -fopenmp
endif
ifeq ($(COMPILER), intel)
  CFLAGS += -openmp
endif
ifeq ($(COMPILER), ibm)
  CFLAGS += -qsmp
endif
  CFLAGS += -DOPENMP
endif

#===============================================================================
# Program Name and Object Files
#===============================================================================

programs = $(cases:.cpp=)
program = $(case:.cpp=)

obj = $(source:.cpp=.o)
headers = $(source:.cpp=.h)

#===============================================================================
# Targets to Build
#===============================================================================

$(program): $(obj) $(addsuffix .o, $(program)) $(headers)
	$(CC) $(CFLAGS) $(obj) $(addsuffix .o, $@) -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

all: $(programs)

clean:
	rm -rf $(programs) $(obj) $(addsuffix .o, $(programs))

edit:
	vim -p $(cases)

run:
	./$(program)
