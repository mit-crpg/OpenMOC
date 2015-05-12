%module openmoc_cuda

%{
  #define SWIG_FILE_WITH_INIT

  #define PySys_WriteStdout printf

  #include <cstddef>
  #include "../../src/constants.h"
  #include "../../src/Solver.h"
  #include "../../src/accel/cuda/GPUSolver.h"
  #include "../../src/accel/cuda/GPUQuery.h"
  #include "../../src/accel/cuda/clone.h"

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

%exception {
  try {
    $function
  } catch (const std::exception &e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

#ifdef NO_NUMPY
#else
%include "../numpy.i"

%init %{
  import_array();
%}

/* The typemap used to match the method signature for the Solver's
 * computeFSRFissionRates method for the data processing routines in
 * openmoc.process */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* fission_rates, int num_FSRs)}


#endif

%include <exception.i>
%include <std_map.i>
%include ../../src/constants.h
%include ../../src/Solver.h
%include ../../src/accel/cuda/GPUSolver.h
%include ../../src/accel/cuda/GPUQuery.h
%include ../../src/accel/cuda/clone.h

#define PySys_WriteStdout printf

typedef float FP_PRECISION;
