#===============================================================================
# User Options
#===============================================================================

COMPILER    = mpi
OPENMP      = yes
OPTIMIZE    = no
DEBUG       = yes
PROFILE     = no
PRECISION   = single
CMFD_PRECISION = single

#===============================================================================
# Source Code List
#===============================================================================

#case = load-geometry/load-full-core.cpp
#case = load-geometry/load-single-assembly.cpp
#case = load-geometry/load-tc-study.cpp
#case = load-geometry/load-few-assm.cpp
#case = load-geometry/load-2D-full-core.cpp
#case = load-geometry/trunc-core-2D.cpp
#case = load-geometry/load-corner.cpp
#case = load-geometry/load-assembly-lattice.cpp

#case = load-geometry/load-geometry.cpp
case = simple-lattice/simple-lattice-3d.cpp
#case = homogeneous/homogeneous.cpp

source = \
Cell.cpp \
Cmfd.cpp \
CPULSSolver.cpp \
CPUSolver.cpp \
ExpEvaluator.cpp \
Geometry.cpp \
LocalCoords.cpp \
linalg.cpp \
log.cpp \
Material.cpp \
Matrix.cpp \
Mesh.cpp \
MOCKernel.cpp \
Point.cpp \
Progress.cpp \
Quadrature.cpp \
Solver.cpp \
Surface.cpp \
Timer.cpp \
Track.cpp \
Track3D.cpp \
TrackGenerator.cpp \
TrackGenerator3D.cpp \
TrackTraversingAlgorithms.cpp \
TraverseSegments.cpp \
Universe.cpp \
Vector.cpp

cases = \
Takeda/Takeda-unrodded.cpp \
pin-cell/pin-cell-3d.cpp \
fixed-source/pin-cell-fixed-3d.cpp \
simple-lattice/simple-lattice-3d.cpp \
gradients/one-directional/one-directional-gradient.cpp \
gradients/two-directional/two-directional-gradient.cpp \
c5g7/c5g7-3d.cpp \
c5g7/c5g7-3d-cmfd.cpp \
c5g7/c5g7-rodded-B.cpp \
c5g7/c5g7-rodded-B-2x2.cpp \
c5g7/c5g7-ws.cpp \
homogeneous/homogeneous.cpp \
single-assembly/quarter-c5g7-assembly.cpp \
single-assembly/single-c5g7-assembly.cpp \
load-geometry/load-geometry.cpp \
load-geometry/load-single-assembly.cpp \
load-geometry/load-few-assm.cpp \
load-geometry/load-2D-full-core.cpp \
load-geometry/trunc-core-2D.cpp \
load-geometry/load-corner.cpp \
load-geometry/load-assembly-lattice.cpp \
load-geometry/load-tc-study.cpp \
load-geometry/load-full-core.cpp

#===============================================================================
# Sets Flags
#===============================================================================

CFLAGS = -std=gnu++0x

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  CFLAGS += -DGNU
  #CFLAGS += -march=native
  LDFLAGS += -lm
  LDFLAGS += -lstdc++
endif

# Regular gcc Compiler
ifeq ($(COMPILER),mpi)
  CC = mpicc
  CFLAGS += -DMPIx
  #CFLAGS += -march=native
  LDFLAGS += -lm
  LDFLAGS += -lstdc++
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
  CFLAGS += -O0
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
  CFLAGS += -O0
  CFLAGS += -fno-omit-frame-pointer
endif

# Precision Flags
ifeq ($(PRECISION),single)
  CFLAGS += -DFP_PRECISION=float
  CFLAGS += -DSINGLE
endif
ifeq ($(PRECISION),double)
  CFLAGS += -DFP_PRECISION=double
#  CFLAGS += -DDOUBLE
endif

ifeq ($(NEW_PRECISION),single)
  CFLAGS += -DNEW_PRECISION=float
endif
ifeq ($(NEW_PRECISION),double)
  CFLAGS += -DNEW_PRECISION=double
endif
ifeq ($(CMFD_PRECISION),single)
  CFLAGS += -DCMFD_PRECISION=float
endif
ifeq ($(CMFD_PRECISION),double)
  CFLAGS += -DCMFD_PRECISION=double
endif

# Vector Flags
CFLAGS += -DVEC_LENGTH=8
CFLAGS += -DVEC_ALIGNMENT=16

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
ifeq ($(COMPILER),gnu)
  CFLAGS += -O3 -ffast-math -fpic
endif
ifeq ($(COMPILER),mpi)
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
ifeq ($(COMPILER), mpi)
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
	rm -rf $(program) $(obj) $(addsuffix .o, $(programs))

edit:
	vim -p $(case) $(cases)

run:
	./$(program)
	#mpirun -np 2 ./$(program)
