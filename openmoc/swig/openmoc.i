%define DOCSTRING
"A method of characteristics code for nuclear reactor physics calculations."
%enddef

/* Include docstrings generated from Doxygen and doxy2swig.py */
%include docstring.i

%module(docstring=DOCSTRING, moduleimport="import _openmoc") openmoc

/* Import rules for Python/C++ transferrable memory ownership */
%include thisown.i

%{
  #define SWIG_FILE_WITH_INIT
  #include <cstddef>
  #include "../../src/boundary_type.h"
  #include "../../src/Cell.h"
  #include "../../src/Cmfd.h"
  #include "../../src/constants.h"
  #include "../../src/ExpEvaluator.h"
  #include "../../src/Geometry.h"
  #include "../../src/linalg.h"
  #include "../../src/log.h"
  #include "../../src/Material.h"
  #include "../../src/Matrix.h"
  #include "../../src/Mesh.h"
  #include "../../src/LocalCoords.h"
  #include "../../src/Point.h"
  #include "../../src/Progress.h"
  #include "../../src/Quadrature.h"
  #include "../../src/Region.h"
  #include "../../src/RunTime.h"
  #include "../../src/segmentation_type.h"
  #include "../../src/Solver.h"
  #include "../../src/CPUSolver.h"
  #include "../../src/CPULSSolver.h"
  #include "../../src/Surface.h"
  #include "../../src/Timer.h"
  #include "../../src/Track.h"
  #include "../../src/Track3D.h"
  #include "../../src/TrackGenerator.h"
  #include "../../src/TrackGenerator3D.h"
  #include "../../src/TraverseSegments.h"
  #include "../../src/TrackTraversingAlgorithms.h"
  #include "../../src/Universe.h"
  #include "../../src/Vector.h"

  #ifdef MPIx
  #include <mpi.h>
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
%warnfilter(511) std::vector;

/* Methods for SWIG to ignore in generating Python API */
%ignore setFSRCentroid(int fsr, Point* centroid);
%ignore printBGQMemory();
%ignore setFSRKeysMap(std::unordered_map<std::size_t, fsr_data>* FSR_keys_map);
%ignore setFSRsToKeys(std::vector<std::size_t>* FSRs_to_keys);
%ignore setFSRsToMaterialIDs(std::vector<int>* FSRs_to_material_IDs);
%ignore setFSRKeysMap(ParallelHashMap<std::size_t, fsr_data*>* FSR_keys_map);
%ignore initializeFSRVectors();
%ignore twiddleRead(int* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(bool* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(char* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(universeType* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(cellType* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(surfaceType* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(boundaryType* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(double* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore twiddleRead(long* ptr, size_t size, size_t nmemb, FILE* stream);
%ignore setRuntimeParameters(RuntimeParameters &RP, int argc, char *argv[]); 

/* Instruct SWIG to ignore methods used in getting CSR Matrix format and Vector
 * attributes. These attributes should be used internally only by the Matrix and
 * Vector class methods and linear algrebra (linalg.h/linalg.cpp) methods. */
%ignore Matrix::getArray();
%ignore Matrix::getA();
%ignore Matrix::getIA();
%ignore Matrix::getJA();
%ignore Matrix::getDiag();

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

/* Routines which pass / return NumPy arrays to / from C++ routine **/
%include numpy_typemaps.i

/* Routines which pass / return pointers to / from C++ routine **/
%include argument_typemaps.i

/* Include standard vector library for SWIG */
%include "std_vector.i"

/* Include standard string library for SWIG */
%include "std_string.i"

namespace std {
  %template(DoubleVector) vector<double>;
  %template(FloatVector) vector<float>;
  %template(IntVector) vector<int>;
  %template(LongVector) vector<long>;
  %template(Array) vector< vector<int> >;
  %template(DoubleArray) vector< vector<double> >;
}

/* Include the MPI library */
#ifdef MPIx
%include "mpi4py.i"
%mpi4py_typemap(Comm, MPI_Comm);
#endif

%include <exception.i>
%include ../../src/boundary_type.h
%include ../../src/Cell.h
%include ../../src/Cmfd.h
%include ../../src/constants.h
%include ../../src/ExpEvaluator.h
%include ../../src/exponentials.h
%include ../../src/Geometry.h
%include ../../src/linalg.h
%include ../../src/log.h
%include ../../src/Material.h
%include ../../src/Matrix.h
%include ../../src/Mesh.h
%include ../../src/LocalCoords.h
%include ../../src/Point.h
%include ../../src/Progress.h
%include ../../src/Quadrature.h
%include ../../src/Region.h
%include ../../src/segmentation_type.h
%include ../../src/RunTime.h
%include ../../src/Solver.h
%include ../../src/CPUSolver.h
%include ../../src/CPULSSolver.h
%include ../../src/Surface.h
%include ../../src/Timer.h
%include ../../src/Track.h
%include ../../src/Track3D.h
%include ../../src/TrackGenerator.h
%include ../../src/TrackGenerator3D.h
%include ../../src/TraverseSegments.h
%include ../../src/TrackTraversingAlgorithms.h
%include ../../src/Universe.h
%include ../../src/Vector.h

/* Use the toString routines for printing objects in Python */
%extend Cell{ std::string __repr__() { return self->toString(); } };
%extend Cmfd{ std::string __repr__() { return self->toString(); } };
%extend LocalCoords{ std::string __repr__() { return self->toString(); } };
%extend Material{ std::string __repr__() { return self->toString(); } };
%extend Point{ std::string __repr__() { return self->toString(); } };
%extend Quadrature{ std::string __repr__() { return self->toString(); } };
%extend Surface{ std::string __repr__() { return self->toString(); } };
%extend Universe{ std::string __repr__() { return self->toString(); } };

#define printf PySys_WriteStdout
