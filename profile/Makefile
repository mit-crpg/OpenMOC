#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPENMP      = yes
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
PRECISION   = single

#===============================================================================
# Source Code List
#===============================================================================

case = c5g7/c5g7-cmfd.cpp

source = \
Cell.cpp \
Cmfd.cpp \
CPUSolver.cpp \
ExpEvaluator.cpp \
Geometry.cpp \
LocalCoords.cpp \
linalg.cpp \
log.cpp \
Material.cpp \
Matrix.cpp \
MOCKernel.cpp \
Point.cpp \
Quadrature.cpp \
Solver.cpp \
Surface.cpp \
Timer.cpp \
Track.cpp \
TrackGenerator.cpp \
TrackTraversingAlgorithms.cpp \
TraverseTracks.cpp \
Universe.cpp \
Vector.cpp

cases = \
gradients/one-directional/one-directional-gradient.cpp \
gradients/two-directional/two-directional-gradient.cpp \
homogeneous/homogeneous-one-group.cpp \
c5g7/c5g7.cpp \
c5g7/c5g7-cmfd.cpp

#===============================================================================
# Sets Flags
#===============================================================================

CFLAGS = -std=gnu++0x

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

# Clang Compiler
ifeq ($(COMPILER),clang)
  CC = clang++
  CFLAGS += -DCLANG
endif

# BG/Q gcc Cross-Compiler
# Note: make sure to set +mpiwrapper-gcc in your .soft file,
# as we only support the gcc cross-compiler wrapper
ifeq ($(COMPILER),bluegene)
  CC = mpic++
  CFLAGS += -DGNU
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
  CFLAGS += -O3 -ffast-math -fpic
endif
ifeq ($(COMPILER),intel)
  CFLAGS += -O3 -xhost -fast -ansi-alias -no-prec-div -ipo
endif
ifeq ($(COMPILER),bluegene)
  CFLAGS += -O3 -ffast-math -fpic
endif
ifeq ($(COMPILER),clang)
  CFLAGS += -O3 -ffast-math -fvectorize -fpic
endif
endif

# OpenMP
ifeq ($(OPENMP),yes)
ifeq ($(COMPILER), gnu)
  CFLAGS += -fopenmp
endif
ifeq ($(COMPILER), bluegene)
  CFLAGS += -fopenmp
endif
ifeq ($(COMPILER), intel)
  CFLAGS += -openmp
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

obj/%.o: ../src/%.cpp Makefile
	$(CC) $(CFLAGS) -c $< -o $@

models/%.o: models/%.cpp Makefile
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
