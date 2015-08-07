/** Routines and typemaps to convert C++ std::map return types to Python dict */

%module map_to_dict

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
    for (int j = 0; j < $2; j++) {
      /* Extract the value from the list at this location and convert
       * SWIG wrapper to pointer to underlying C++ class instance */
      PyObject* o = PyList_GetItem(outer_list, j);
      void *p1 = 0;
      SWIG_ConvertPtr(o, &p1, SWIGTYPE_p_Universe, 0 | 0);
      $3[i*$2+j] = (Universe*) p1;
    }
  }
}
