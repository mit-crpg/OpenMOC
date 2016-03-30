/** Helper typemaps to enable array-based user input to be defined with Python
 *  lists when the --no-numpy build system flag is used at compile time. */

%module swig_typemaps

/* Typemap for the Material::set_____XS(double* xs, int num_groups)
 * method - allows users to pass in a Python list of cross-sections
 * for each energy group */
%typemap(in) (double* xs, int num_groups) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of values "
                    "for the cross-section array");
    return NULL;
  }

  $2 = PySequence_Length($input);  // num_groups
  $1 = (double*) malloc($2 * sizeof(double));  // cross-section array

  /* Loop over x */
  for (int i = 0; i < $2; i++) {

    /* Extract the value from the list at this location */
    PyObject *o = PySequence_GetItem($input,i);

    /* If value is a number, cast it as an int and set the input array value */
    if (PyNumber_Check(o)) {
      $1[i] = (double) PyFloat_AsDouble(o);
    }
    else {
      free($1);
      PyErr_SetString(PyExc_ValueError,"Expected a list of numbers "
                      "for cross-section values\n");
      return NULL;
    }
  }
}


/* Typemap for the Material::set_____XS(float* xs, int num_groups)
 * method - allows users to pass in a Python list of cross-sections
 * for each energy group */
%typemap(in) (float* xs, int num_groups) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of values "
                    "for the cross-section array");
    return NULL;
  }

  $2 = PySequence_Length($input);  // num_groups
  $1 = (float*) malloc($2 * sizeof(float));  // cross-section array

  /* Loop over x */
  for (int i = 0; i < $2; i++) {

    /* Extract the value from the list at this location */
    PyObject *o = PySequence_GetItem($input,i);

    /* If value is a number, cast it as an int and set input array value */
    if (PyNumber_Check(o)) {
      $1[i] = (float) PyFloat_AsFloat(o);
    }
    else {
      free($1);
      PyErr_SetString(PyExc_ValueError,"Expected a list of numbers "
                      "for cross-section values\n");
      return NULL;
    }
  }
}


/* Typemap for the Cmfd::setGroupStructure
 * (int* group_indices, int length_group_indices)
 * method - allows users to pass in a Python list of group indices
 * for each CMFD energy group */
%typemap(in) (int* group_indices, int length_group_indices) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of values "
                    "for the group indices array");
    return NULL;
  }

  $2 = PySequence_Length($input);  // length_group_indices
  $1 = (int*) malloc($2 * sizeof(int));  // group indices array

  /* Loop over x */
  for (int i = 0; i < $2; i++) {

    /* Extract the value from the list at this location */
    PyObject *o = PySequence_GetItem($input,i);

    /* If value is a number, cast it as an int and set the input array value */
    if (PyNumber_Check(o)) {
      $1[i] = (int) PyInt_AS_LONG(o);
    }
    else {
      free($1);
      PyErr_SetString(PyExc_ValueError,"Expected a list of numbers "
                      "for group indices values\n");
      return NULL;
    }
  }
}
