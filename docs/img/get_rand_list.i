%module get_rand_list
%{
  #define SWIG_FILE_WITH_INIT
  #include "get_rand_list.h"
%}

%include "std_vector.i"

/* SWIG template for get_rand_list(int length) C++ routine */
namespace std {
   %template(DoubleVector) vector<double>;
}

%include "get_rand_list.h"


