%module openmoc_bgq_double

%{
  #define SWIG_FILE_WITH_INIT
  #include <cstddef>
  #include "../../../src/Cell.h"
  #include "../../../src/Geometry.h"
  #include "../../../src/LocalCoords.h"
  #include "../../../src/log.h"
  #include "../../../src/Material.h"
  #include "../../../src/Point.h"
  #include "../../../src/PolarQuad.h"
  #include "../../../src/Solver.h"
  #include "../../../src/CPUSolver.h"
  #include "../../../src/Surface.h"
  #include "../../../src/Timer.h"
  #include "../../../src/Track.h"
  #include "../../../src/TrackGenerator.h"
  #include "../../../src/Universe.h"
  #include "../../../src/Cmfd.h"

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

%exception {
  try {
    $function
  } catch (const std::exception &e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

/* C++ casting helper method for openmoc.process computePinPowers routine */
%inline %{
  CellFill* castCellToCellFill(Cell* cell) {
    return dynamic_cast<CellFill*>(cell);
  }

  CellBasic* castCellToCellBasic(Cell* cell) {
    return dynamic_cast<CellBasic*>(cell);
  }

  Lattice* castUniverseToLattice(Universe* universe) {
    return dynamic_cast<Lattice*>(universe);
  }

  Universe* castLatticeToUniverse(Lattice* lattice) {
    return dynamic_cast<Universe*>(lattice);
  }

%}


/* If the user uses the --no-numpy flag, then NumPy typemaps will not be used
 * and the NumPy C API will not be embedded in the source code. The NumPy
 * typemaps are used to allow users to pass NumPy arrays to/from the C++ source
 * code. This does not work well on BGQ, and some users may prefer not to embed
 * this into their code, so if --no-numpy is passed in we use SWIG typemaps to
 * allow users to pass in arrays of data as Python lists. */
#ifdef NO_NUMPY

%include typemaps.i


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

    /* If the value is a number, cast it as an int and set the
     * input array value */
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

    /* If the value is a number, cast it as an int and set the
     * input array value */
    if (PyNumber_Check(o)) {
      $1[i] = (float) PyFloat_AsFloat(o);
    }
    else {
      free($1);
      PyErr_SetString(PyExc_ValueError,"Expected a list of numbers "
                      "as universe IDs when constructing Lattice cells\n");
      return NULL;
    }
  }
}


/* Typemap for Lattice::setLatticeCells(int num_x, int num_y, int* universes)
 * method - allows users to pass in Python lists of Universe IDs for each
 * lattice cell */
%typemap(in) (int num_x, int num_y, int* universes) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of integers "
                    "for the Lattice cells");
    return NULL;
  }

  $1 = PySequence_Length($input);  // num_x
  $2 = PySequence_Length(PyList_GetItem($input,0)); // num_y
  $3 = (int*) malloc(($1 * $2) * sizeof(int));  // universes

  /* Loop over x */
  for (int i = 0; i < $2; i++) {

    /* Get the inner list in the nested list for the lattice */
    PyObject* outer_list = PySequence_GetItem($input,i);

    /* Check that the length of this list is the same as the length
     * of the first list */
    if (PySequence_Length(outer_list) != $2) {
      PyErr_SetString(PyExc_ValueError, "Size mismatch. Expected $1 x $2 "
                      "elements for lattice\n");
      return NULL;
    }

    /* Loop over y */
    for (int j =0; j < $1; j++) {

      /* Extract the value from the list at this location */
      PyObject *o = PySequence_GetItem(outer_list,j);

      /* If the value is a number, cast it as an int and set the
       * input array value */
      if (PyNumber_Check(o)) {
        $3[i*$1 + j] = (int) PyInt_AS_LONG(o);
      }
      else {
        free($3);
        PyErr_SetString(PyExc_ValueError,"Expected a list of numbers as"
                        "Universe IDs when constructing Lattice cells\n");
        return NULL;
      }
    }
  }
}


/* If the user did not pass in the --no-numpy flag, then NumPy typemaps will be
 * used and the NumPy C API will be embedded in the source code. This will allow
 * users to pass arrays of data to/from the C++ source code (ie, setting group
 * cross-section values or extracting the scalar flux). */
#else

%include "../../numpy.i"


%init %{
  import_array();
%}

/* The typemape used to match the method signature for the
 * Lattice::setLatticeCells setter method. This allows users to set the lattice
 * cells (universe IDs) using a 2D NumPy array */
%apply (int DIM1, int DIM2, int* IN_ARRAY2) {(int num_x, int num_y, int* universes)}

/* The typemap used to match the method signature for the Material
 * cross-section setter methods. This allows users to set the cross-sections
 * using NumPy arrays */
%apply (double* IN_ARRAY1, int DIM1) {(double* xs, int num_groups)}

/* The typemape used to match the method signature for the TrackGenerator's
 * getter methods for track start and end coordinates for the plotting
 * routines in openmoc.plotter */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* coords, int num_tracks)}

/* The typemape used to match the method signature for the TrackGenerator's
 * getter methods for track segment start and end coordinates for the plotting
 * routines in openmoc.plotter */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* coords, int num_segments)}

/* The typemap used to match the method signature for the Solver's
 * computeFSRFissionRates method for the data processing routines in
 * openmoc.process */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* fission_rates, int num_FSRs)}

/* The typemap used to match the method signature for the Universe's
 * getCellIds method for the data processing routines in openmoc.process */
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* cell_ids, int num_cells)}


#endif


%include <exception.i>
%include ../../../src/Cell.h
%include ../../../src/Geometry.h
%include ../../../src/LocalCoords.h
%include ../../../src/log.h
%include ../../../src/Material.h
%include ../../../src/Point.h
%include ../../../src/PolarQuad.h
%include ../../../src/Solver.h
%include ../../../src/CPUSolver.h
%include ../../../src/Surface.h
%include ../../../src/Timer.h
%include ../../../src/Track.h
%include ../../../src/TrackGenerator.h
%include ../../../src/Universe.h
%include ../../../src/Cmfd.h

typedef double FP_PRECISION;
