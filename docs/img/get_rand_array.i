%module get_rand_array

%{
  #define SWIG_FILE_WITH_INIT
  #include "get_rand_array.h"
%}

/* Include the NumPy typemaps library */
%include "numpy.i"

%init %{
  import_array();
%}

/* Typemap for get_rand_array(double* output_array, int length) C/C++ routine */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* output_array, int length)};

%include "get_rand_array.h"
