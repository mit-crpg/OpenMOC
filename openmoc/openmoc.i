%define DOCSTRING 
"A method of characteristics code for nuclear reactor physics calculations."
%enddef

%module(docstring=DOCSTRING) openmoc

%{
  #define SWIG_FILE_WITH_INIT
  #include <cstddef>
  #include "../src/Cell.h"
  #include "../src/Geometry.h"
  #include "../src/boundary_type.h"
  #include "../src/LocalCoords.h"
  #include "../src/log.h"
  #include "../src/Material.h"
  #include "../src/Point.h"
  #include "../src/PolarQuad.h"
  #include "../src/Solver.h"
  #include "../src/CPUSolver.h"
  #include "../src/Surface.h"
  #include "../src/Timer.h"
  #include "../src/Track.h"
  #include "../src/TrackGenerator.h"
  #include "../src/Universe.h"
  #include "../src/Cmfd.h"

  #ifdef ICPC
  #include "../src/VectorizedSolver.h"
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

%exception {
  try {
    $function
  } catch (const std::exception &e) {
      SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

/* C++ casting helper method for openmoc.process computePinPowers
 * routine and the OpenCG compatibility module */
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

  Plane* castSurfaceToPlane(Surface* plane) {
    return dynamic_cast<Plane*>(plane);
  }

  XPlane* castSurfaceToXPlane(Surface* xplane) {
    return dynamic_cast<XPlane*>(xplane);
  }

  YPlane* castSurfaceToYPlane(Surface* yplane) {
    return dynamic_cast<YPlane*>(yplane);
  }

  ZPlane* castSurfaceToZPlane(Surface* zplane) {
    return dynamic_cast<ZPlane*>(zplane);
  }

  Circle* castSurfaceToCircle(Surface* circle) {
    return dynamic_cast<Circle*>(circle);
  }

%}


/* If the user uses the --no-numpy flag, then NumPy typemaps will not be used
 * and the NumPy C API will not be embedded in the source code. The NumPy
 * typemaps are used to allow users to pass NumPy arrays to/from the C++ source
 * code. This does not work well on BGQ, and some users may prefer not to embed
 * this into their code, so if --no-numpy is passed in we use SWIG typemaps to
 * allow users to pass in arrays of data as Python lists. */
#ifdef NO_NUMPY


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



/* If the user did not pass in the --no-numpy flag, then NumPy typemaps will be
 * used and the NumPy C API will be embedded in the source code. This will allow
 * users to pass arrays of data to/from the C++ source code (ie, setting group
 * cross-section values or extracting the scalar flux). */
#else

%include "numpy.i"

%init %{
  import_array();
%}

/* The typemap used to match the method signature for the
 * Lattice::setLatticeCells setter method. This allows users to set the lattice
 * cells (universe IDs) using a 2D NumPy array */
%apply (int DIM1, int DIM2, int* IN_ARRAY2) {(int num_x, int num_y, int* universes)}

/* The typemap used to match the method signature for the
 * Cmfd::setGroupStructure method. This allows users to set the CMFD group 
 * structure using a NumPy array */
%apply (int* IN_ARRAY1, int DIM1) {(int* group_indices, int length_group_indices)}

/* The typemap used to match the method signature for the Material
 * cross-section setter methods. This allows users to set the cross-sections
 * using NumPy arrays */
%apply (double* IN_ARRAY1, int DIM1) {(double* xs, int num_groups)}

/* The typemap used to match the method signature for the TrackGenerator's
 * getter methods for track start and end coordinates for the plotting
 * routines in openmoc.plotter */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* coords, int num_tracks)}

/* The typemap used to match the method signature for the TrackGenerator's
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

/* The typemap used to match the method signature for the 
 * PolarQuad::setSinThetas method. This allows users to set the polar angle 
 * quadrature sine thetas using a NumPy array */
%apply (double* IN_ARRAY1, int DIM1) {(double* sin_thetas, int num_polar)}

/* The typemap used to match the method signature for the 
 * PolarQuad::setWeights method. This allows users to set the polar angle 
 * quadrature weights using a NumPy array */
%apply (double* IN_ARRAY1, int DIM1) {(double* weights, int num_polar)}

#endif


/* Typemap for all methods which return a std::map<int, Cell*>. This includes
 * the Geometry::getAllCells(), Universe::getAllCells(), etc. These methods
 * are particularly useful for OpenCG compatibility. */
%include <std_map.i>
%clear std::map<int, Cell*>;
%typemap(out) std::map<int, Cell*> {

  $result = PyDict_New();
  int size = $1.size();

  std::map<int, Cell*>::iterator iter;
  Cell* cell;
  int cell_id;

  for (iter = $1.begin(); iter != $1.end(); ++iter) {
    cell_id = iter->first;
    cell = iter->second;
    PyObject* value =
         SWIG_NewPointerObj(SWIG_as_voidptr(cell), $descriptor(Cell*), 0);
    PyDict_SetItem($result, PyInt_FromLong(cell_id), value);
  }
}


/* Typemap for all methods which return a std::map<int, surface_halfspace>.
 * This includes the Cell::getSurfaces() method, which is useful for OpenCG
 * compatibility. */
%include <std_map.i>
%clear std::map<int, surface_halfspace>;
%typemap(out) std::map<int, surface_halfspace> {

  $result = PyDict_New();
  int size = $1.size();

  std::map<int, surface_halfspace>::iterator iter;
  surface_halfspace* surf;
  int surf_id;

  for (iter = $1.begin(); iter != $1.end(); ++iter) {
    surf_id = iter->first;
    surf = &iter->second;
    PyObject* value =
         SWIG_NewPointerObj(SWIG_as_voidptr(surf),
                            $descriptor(surface_halfspace*), 0);
    PyDict_SetItem($result, PyInt_FromLong(surf_id), value);
  }
}


/* Typemap for all methods which return a std::map<int, Material*>.
 * This includes the Geometry::getAllMaterials() method, which is useful 
 * for OpenCG compatibility. */
%include <std_map.i>
%clear std::map<int, Material*>;
%typemap(out) std::map<int, Material*> {

  $result = PyDict_New();
  int size = $1.size();

  std::map<int, Material*>::iterator iter;
  Material* mat;
  int mat_id;

  for (iter = $1.begin(); iter != $1.end(); ++iter) {
    mat_id = iter->first;
    mat = iter->second;
    PyObject* value =
         SWIG_NewPointerObj(SWIG_as_voidptr(mat), $descriptor(Material*), 0);
    PyDict_SetItem($result, PyInt_FromLong(mat_id), value);
  }
}


/* Typemap for all methods which return a std::map<int, Universe*>.
 * This includes the Lattice::getUniqueUniverses() method which is ueseful for
 * OpenCG compatibility. */
%include <std_map.i>
%clear std::map<int, Universe*>;
%typemap(out) std::map<int, Universe*> {

  $result = PyDict_New();
  int size = $1.size();

  std::map<int, Universe*>::iterator iter;
  Universe* univ;
  int univ_id;

  for (iter = $1.begin(); iter != $1.end(); ++iter) {
    univ_id = iter->first;
    univ = iter->second;
    PyObject* value =
         SWIG_NewPointerObj(SWIG_as_voidptr(univ), $descriptor(Universe*), 0);
    PyDict_SetItem($result, PyInt_FromLong(univ_id), value);
  }
}


/* Typemap for Lattice::setUniverses(int num_x, int num_y, Universe** universes)
 * method - allows users to pass in Python lists of Universes for each
 * lattice cell */
%typemap(in) (int num_x, int num_y, Universe** universes) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of integers "
                    "for the Lattice cells");
    return NULL;
  }

  $1 = PySequence_Length($input);  // num_x
  $2 = PySequence_Length(PyList_GetItem($input,0)); // num_y
  $3 = (Universe**) malloc(($1 * $2) * sizeof(Universe*)); // universes

  /* Loop over x */
  for (int i = 0; i < $1; i++) {

    /* Get the inner list in the nested list for the lattice */
    PyObject* outer_list = PyList_GetItem($input,i);

    /* Check that the length of this list is the same as the length
     * of the first list */
    if (PySequence_Length(outer_list) != $2) {
      PyErr_SetString(PyExc_ValueError, "Size mismatch in Universes "
                      "list for Lattice which must be a 2D list of lists");
      return NULL;
    }

    /* Loop over y */
    for (int j =0; j < $2; j++) {
      /* Extract the value from the list at this location and convert
       * SWIG wrapper to pointer to underlying C++ class instance */
      PyObject* o = PyList_GetItem(outer_list, j);
      void *p1 = 0;
      SWIG_ConvertPtr(o, &p1, SWIGTYPE_p_Universe, 0 | 0);
      $3[i*$2+j] = (Universe*) p1;
    }
  }
}


%include <exception.i>
%include <std_map.i>
%include ../src/Cell.h
%include ../src/Geometry.h
%include ../src/boundary_type.h
%include ../src/LocalCoords.h
%include ../src/log.h
%include ../src/Material.h
%include ../src/Point.h
%include ../src/PolarQuad.h
%include ../src/Solver.h
%include ../src/CPUSolver.h
%include ../src/Surface.h
%include ../src/Timer.h
%include ../src/Track.h
%include ../src/TrackGenerator.h
%include ../src/Universe.h
%include ../src/Cmfd.h

#ifdef ICPC
%include "../src/VectorizedSolver.h"
#endif

#define printf PySys_WriteStdout

#ifdef DOUBLE
typedef double FP_PRECISION;
#else
typedef float FP_PRECISION;
#endif
