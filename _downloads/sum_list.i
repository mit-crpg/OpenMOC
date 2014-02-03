%module sum_list
%{
  #define SWIG_FILE_WITH_INIT
  #include "sum_list.h"
%}

/* Include the SWIG typemaps library */
%include typemaps.i

/* Typemap for the sum_list(double* input_list, int length) C/C++ routine */
%typemap(in) (double* input_list, int length) {

  /* Check that the input is a Python list data structure */
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of values\n");
    return NULL;
  }

  /* Set the second parameter to the length of the Python list input */
  $2 = PySequence_Length($input);

  /* Allocate memory to convert the list into a C/C++ array */
  $1 = (double*) malloc($2 * sizeof(double));

  /* Loop over the values in the list */
  for (int i = 0; i < $2; i++) {
    
    /* Extract the value from the list at this location */
    PyObject *o = PySequence_GetItem($input,i);

    /* If the value is a number, cast it as an int and set the
     * input array value */
    if (PyNumber_Check(o)) {
      $1[i] = (double) PyFloat_AsDouble(o);
    }
    else {
      free($1);
      PyErr_SetString(PyExc_ValueError,"Expected a list of numbers\n");
      return NULL;
    }
  }
}

%include "sum_list.h"
