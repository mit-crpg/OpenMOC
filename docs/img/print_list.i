%module print_list
%{
  #define SWIG_FILE_WITH_INIT
  #include "print_list.h"
%}

/* Include the SWIG typemaps library */
%include typemaps.i

/* Typemap for the print(int length, double* list) C/C++ routine */
%typemap(in) (int length, double* list) {

  /* Check that the input is a Python list data structure */
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of values\n");
    return NULL;
  }

  /* Set the first parameter to the length of the Python list input */
  $1 = PySequence_Length($input);

  /* Allocate memory to convert the list into a C/C++ array */
  $2 = (double*) malloc($1 * sizeof(double));

  /* Loop over the values in the list */
  for (int i = 0; i < $1; i++) {
    
    /* Extract the value from the list at this location */
    PyObject *o = PySequence_GetItem($input,i);

    /* If the value is a number, cast it as an int and set the
     * input array value */
    if (PyNumber_Check(o)) {
      $2[i] = (double) PyFloat_AsDouble(o);
    }
    else {
      free($2);
      PyErr_SetString(PyExc_ValueError,"Expected a list of numbers\n");
      return NULL;
    }
  }
}

%include "print_list.h"
