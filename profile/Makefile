#===============================================================================
# User Options
#===============================================================================

COMPILER    = ibm
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

#case = homogeneous/homogeneous-one-group.cpp
case = c5g7/c5g7-cmfd.cpp

source = \
Cell.cpp \
Cmfd.cpp \
CPUSolver.cpp \
ExpEvaluator.cpp \
Geometry.cpp \
LocalCoords.cpp \
log.cpp \
Material.cpp \
Point.cpp \
PolarQuad.cpp \
Solver.cpp \
Surface.cpp \
Timer.cpp \
Track.cpp \
TrackGenerator.cpp \
Universe.cpp

cases = \
gradients/one-directional/one-directional-gradient.cpp \
gradients/two-directional/two-directional-gradient.cpp \
homogeneous/homogeneous-one-group.cpp \
c5g7/c5g7.cpp \
c5g7/c5g7-cmfd.cpp

#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
ifeq ($(COMPILER),ibm)
  CFLAGS += -std=gnu++0x
else
  CFLAGS += -std=c++11
endif

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = g++
  CFLAGS += -DGNU
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icpc
  source += VectorizedSolver.cpp
  CFLAGS += -DINTEL
  CFLAGS += -I/opt/intel/composer_xe_2013_sp1.3.174/mkl/include
  LDFLAGS += -mkl
endif

# IBM XL Compiler
ifeq ($(COMPILER),ibm)
  CC = mpic++
  CFLAGS += -DIBM
endif

# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  CC = bgxlc++
endif

# Clang Compiler
ifeq ($(COMPILER),clang)
  CC = clang
  CFLAGS += -DCLANG 
  CFLAGS += -Qunused-arguments 
  CFLAGS += -Wno-deprecated-register -Wno-parentheses-equality
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  ifeq ($(COMPILER),intel)
  CFLAGS += -vec-report
  endif
  ifeq ($(OPENMP),yes)
  CFLAGS += -openmp-report
  endif
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

# Vector Flags
CFLAGS += -DVEC_LENGTH=8
CFLAGS += -DVEC_ALIGNMENT=16

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
ifeq ($(COMPILER),gnu)
  CFLAGS += -Ofast -ffast-math -ftree-vectorize -msse2
endif
ifeq ($(COMPILER),intel)
  CFLAGS += -O3 -xhost -fast -ansi-alias -no-prec-div -ipo
endif
# XL compiler not compatible at the moment, use gcc flags
ifeq ($(COMPILER),ibm)
  CFLAGS += -O3 -ffast-math -ftree-vectorize
#  CFLAGS += -O5 -qhot -qsimd=auto -qalias=ansi:allptrs -qarch=qp
endif
ifeq ($(COMPILER),clang)
  CFLAGS += -O3 -ffast-math -fvectorize
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
  #CFLAGS += -qsmp
  CFLAGS += -fopenmp
endif
ifeq ($(COMPILER), clang)
  CFLAGS += -fopenmp
endif
  CFLAGS += -DOPENMP
endif

#===============================================================================
# Program Name and Object Files
#===============================================================================

obj := $(source:.cpp=.o)

source := $(addprefix ../src/, $(source))
headers = $(source:.cpp=.h)

obj := $(addprefix obj/, $(obj))
cases := $(addprefix models/, $(cases))
case := $(addprefix models/, $(case))

programs = $(cases:.cpp=)
program = $(case:.cpp=)

#===============================================================================
# Targets to Build
#===============================================================================

default: folder $(addsuffix .o, $(program)) $(program)

$(programs): $(obj)  $(headers)
	$(CC) $(CFLAGS) $(obj) $(addsuffix .o, $@) -o $@ $(LDFLAGS)

obj/%.o: ../src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

models/%.o: models/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

folder:
	mkdir -p obj

all: folder $(addsuffix .o, $(programs)) $(programs) 

clean:
	rm -rf $(programs) $(obj) $(addsuffix .o, $(programs))

edit:
	vim -p $(cases)

run:
	./$(program)
