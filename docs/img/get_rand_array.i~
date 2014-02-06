%module get_rand_array

%{
  #define SWIG_FILE_WITH_INIT
  #include "get_rand_array.h"
%}

%include "numpy.i"

%init %{
  import_array();
%}

%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* output_array, int length)};

%include "get_rand_array.h"
