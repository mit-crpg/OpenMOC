#===============================================================================
# User Options
#===============================================================================

COMPILER       ?= mpi
OPENMP         ?= yes
OPTIMIZE       ?= yes
DEBUG          ?= no
PROFILE        ?= no
INFO           ?= no
PRECISION      ?= single
CMFD_PRECISION ?= single

#===============================================================================
# Source Code List
#===============================================================================

#case = models/load-geometry/load-full-core.cpp
#case = models/load-geometry/load-single-assembly.cpp
#case = models/load-geometry/load-tc-study.cpp
#case = models/load-geometry/load-few-assm.cpp
#case = models/load-geometry/load-2D-full-core.cpp
#case = models/load-geometry/trunc-core-2D.cpp
#case = models/load-geometry/load-corner.cpp
#case = models/load-geometry/load-assembly-lattice.cpp

#case = models/load-geometry/load-geometry.cpp
#case = models/simple-lattice/simple-lattice-3d.cpp
#case = models/homogeneous/homogeneous.cpp
case ?= models/run_time_standard/run_time_standard.cpp
#case ?= models/assembly_1/load-assembly.cpp
#case ?= models/core/load-core.cpp

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
Region.cpp \
RunTime.cpp \
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
load-geometry/load-2D-full-core.cpp \
load-geometry/load-full-core.cpp

#===============================================================================
# Sets Flags
#===============================================================================

CFLAGS = -std=gnu++0x

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  CFLAGS += -DGNU
  LDFLAGS += -lm
  LDFLAGS += -lstdc++
endif

# MPI wrapped compiler
ifeq ($(COMPILER),mpi)
  CC = mpic++
  CFLAGS += -DMPIx
  LDFLAGS += -lm
  LDFLAGS += -lstdc++
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icpc
  #source += VectorizedSolver.cpp
  CFLAGS += -DINTEL
  #CFLAGS += -I/your/path/to/mkl/include
  #LDFLAGS += -mkl
endif

# MPI wrapped intel Compiler
ifeq ($(COMPILER),impi)
  CC = mpiicpc
  CFLAGS += -DINTEL
  CFLAGS += -DMPIx
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
  CFLAGS += -fno-omit-frame-pointer
else
  ifeq ($(COMPILER),gnu)
    CFLAGS += -Wno-unused-result
  endif
  ifeq ($(COMPILER),mpi)
    CFLAGS += -Wno-unused-result
  endif
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  # pprof flags
  CFLAGS += -pg
  CFLAGS += -O0
  CFLAGS += -fno-omit-frame-pointer
  # gprof and vtune flags
  #CFLAGS += -g
endif

# Precision Flags
ifeq ($(PRECISION),single)
  CFLAGS += -DFP_PRECISION=float
  CFLAGS += -DSINGLE
endif
ifeq ($(PRECISION),double)
  CFLAGS += -DFP_PRECISION=double
endif

ifeq ($(CMFD_PRECISION),single)
  CFLAGS += -DCMFD_PRECISION=float -DLINALG_TOL=1E-7
endif
ifeq ($(CMFD_PRECISION),double)
  CFLAGS += -DCMFD_PRECISION=double -DLINALG_TOL=1E-15
endif

# Vector Flags, take into account strip mining with VEC_LENGTH sized vectors
CFLAGS += -DVEC_ALIGNMENT=64
ifeq ($(PRECISION),single)
  CFLAGS += -DVEC_LENGTH=8
endif
ifeq ($(PRECISION),double)
  CFLAGS += -DVEC_LENGTH=4
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  ifeq ($(COMPILER),gnu)
    CFLAGS += -O3 -ffast-math -fpic -march=native #-flto
  endif
  ifeq ($(COMPILER),mpi)
    CFLAGS += -O3 -ffast-math -fpic -march=native #-flto
  endif
  ifeq ($(COMPILER),intel)
    CFLAGS += -O3 -xhost -ansi-alias -no-prec-div #-ipo
  endif
  ifeq ($(COMPILER),impi)
    CFLAGS += -O3 -xhost -ansi-alias -no-prec-div #-ipo
  endif
  ifeq ($(COMPILER),bluegene)
    CFLAGS += -O3 -ffast-math -fpic
  endif
  ifeq ($(COMPILER),clang)
    CFLAGS += -O3 -ffast-math -fvectorize -fpic
  endif
  # Set the number of groups at compile time
  #CFLAGS += -DNGROUPS=70
  # Set the source to linear source at compile time
  #CFLAGS += -DLINEARSOURCE
  # Set the dimension to 3D at compile time
  #CFLAGS += -DTHREED
  # Only vacuum boundary conditions, avoid double store of track fluxes
  #CFLAGS += -DONLYVACUUMBC
endif

# Optimization report flags
ifeq ($(INFO),yes)
  ifeq ($(COMPILER),intel)
    CFLAGS += -qopt-report n=5 -qopt-report-phase=all
  endif
  ifeq ($(COMPILER),impi)
    CFLAGS += -qopt-report n=5 -qopt-report-phase=all
  endif
  ifeq ($(COMPILER),gnu)
    CFLAGS += -fopt-info-all
  endif
  ifeq ($(COMPILER),mpi)
    CFLAGS += -fopt-info-all
  endif
  ifeq ($(OPENMP),yes)
    CFLAGS += -openmp-report
  endif
endif

# OpenMP flags
ifeq ($(OPENMP),yes)
  CFLAGS += -fopenmp
  CFLAGS += -DOPENMP
endif

CFLAGS += -I../src/
#===============================================================================
# Program Name and Object Files
#===============================================================================

obj := $(source:.cpp=.o)

source := $(addprefix ../src/, $(source))
headers = $(source:.cpp=.h)

obj := $(addprefix obj/, $(obj))
cases := $(addprefix models/, $(cases))
cases += $(case)

programs = $(cases:.cpp=)
program = $(case:.cpp=)

#===============================================================================
# Targets to Build
#===============================================================================

default: folder $(addsuffix .o, $(program)) $(program)

$(programs): $(obj)  $(headers) $(addsuffix .o, $(program))
	$(CC) $(CFLAGS) $(obj) $(addsuffix .o, $@) -o $@ $(LDFLAGS)

obj/%.o: ../src/%.cpp Makefile
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp Makefile
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
