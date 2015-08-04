%define DOCSTRING 
"A method of characteristics code for nuclear reactor physics calculations."
%enddef

%module(docstring=DOCSTRING) openmoc

/* Import rules for Python/C++ transferrable memory ownership */
%include thisown.i

%{
  #define SWIG_FILE_WITH_INIT
  #include <cstddef>
  #include "../src/constants.h"
  #include "../src/Cell.h"
  #include "../src/Geometry.h"
  #include "../src/LocalCoords.h"
  #include "../src/log.h"
  #include "../src/Material.h"
  #include "../src/Point.h"
  #include "../src/PolarQuad.h"
  #include "../src/Solver.h"
  #include "../src/CPUSolver.h"
  #include "../src/boundary_type.h"
  #include "../src/Surface.h"
  #include "../src/Timer.h"
  #include "../src/Track.h"
  #include "../src/TrackGenerator.h"
  #include "../src/Universe.h"
  #include "../src/Cmfd.h"
  #include "../src/Vector.h"
  #include "../src/Matrix.h"
  #include "../src/linalg.h"

  #ifdef ICPC
  #include "../src/VectorizedSolver.h"
  #endif

  #define printf PySys_WriteStdout

  /* Exception helpers */
  static int swig_c_error_num = 0;
  static char swig_c_err_msg[1024];

  const char* err_occurred(void) {
    if (swig_c_error_num) {
      swig_c_error_num = 0;
      return (const char*)swig_c_err_msg;
    }
    return NULL;
  }

  void set_err(const char *msg) {
    swig_c_error_num = 1;
    strncpy(swig_c_err_msg, msg, 1024);
  }

%}

%warnfilter(506) log_printf(logLevel level, const char *format, ...);
%warnfilter(511) swig::SwigPyIterator;
%warnfilter(511) Cell::setFill;

/* Methods for SWIG to ignore in generating Python API */
%ignore setFSRCentroid(int fsr, Point* centroid);
%ignore setFSRKeysMap(std::unordered_map<std::size_t, fsr_data>* FSR_keys_map);
%ignore setFSRsToKeys(std::vector<std::size_t>* FSRs_to_keys);
%ignore setFSRsToMaterialIDs(std::vector<int>* FSRs_to_material_IDs);
%ignore setFSRKeysMap(ParallelHashMap<std::size_t, fsr_data*>* FSR_keys_map);
%ignore initializeFSRVectors();

/* Instruct SWIG to ignore methods used in getting CSR Matrix format and Vector
 * attributes. These attributes should be used internally only by the Matrix and
 * Vector class methods and linear algrebra (linalg.h/linalg.cpp) methods. */
%ignore getArray();
%ignore getA();
%ignore getIA();
%ignore getJA();
%ignore getDiag();

%exception {
  try {
    $function
  } catch (const std::exception &e) {
      SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

/* Routines to allow parent classes to be cast to subclasses from Python */
%include casting.i

/* Routines which convert std::map return types to Python dictionaries */
%include map_to_dict.i

/* If the user uses the --no-numpy flag, then NumPy typemaps will not be used
 * and the NumPy C API will not be embedded in the source code. The NumPy
 * typemaps are used to allow users to pass NumPy arrays to/from the C++ source
 * code. This does not work well on BGQ, and some users may prefer not to embed
 * this into their code, so if --no-numpy is passed in we use SWIG typemaps to
 * allow users to pass in arrays of data as Python lists. */
#ifdef NO_NUMPY
%include swig_typemaps.i
#else
%include numpy_typemaps.i
#endif

%include <exception.i>
%include ../src/constants.h
%include ../src/Cell.h
%include ../src/Geometry.h
%include ../src/LocalCoords.h
%include ../src/log.h
%include ../src/Material.h
%include ../src/Point.h
%include ../src/PolarQuad.h
%include ../src/Solver.h
%include ../src/CPUSolver.h
%include ../src/boundary_type.h
%include ../src/Surface.h
%include ../src/Timer.h
%include ../src/Track.h
%include ../src/TrackGenerator.h
%include ../src/Universe.h
%include ../src/Cmfd.h
%include ../src/Vector.h
%include ../src/Matrix.h
%include ../src/linalg.h

#ifdef ICPC
%include ../src/VectorizedSolver.h
#endif

#define printf PySys_WriteStdout

#ifdef DOUBLE
typedef double FP_PRECISION;
#else
typedef float FP_PRECISION;
#endif
