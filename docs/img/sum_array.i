%module sum_array

%{
  #define SWIG_FILE_WITH_INIT
  #include "sum_array.h"
%}

/* Include the NumPy typemaps library */
%include "numpy.i"

%init %{
  import_array();
%}

/* Typemap for the sum_list(double* input_list, int length) C/C++ routine */
%apply (double* IN_ARRAY1, int DIM1) {(double* input_array, int length)};

%include "sum_array.h"
