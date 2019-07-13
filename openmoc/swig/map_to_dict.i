/** Routines and typemaps to convert C++ std::map return types to Python dict */

%module map_to_dict

/* Typemap for all methods which return a std::map<int, Cell*>. This includes
 * the Geometry::getAllCells(), Universe::getAllCells(), etc. These methods
 * are particularly useful for OpenMC compatibility. */
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


/* Typemap for all methods which return a std::map<int, Surface*>. This
 * includes the Geometry::getAllSurfaces() method, which is useful for 
 * OpenMC compatibility. */
%clear std::map<int, Surface*>;
%typemap(out) std::map<int, Surface*> {

  $result = PyDict_New();
  int size = $1.size();

  std::map<int, Surface*>::iterator iter;
  Surface* surf;
  int surf_id;

  for (iter = $1.begin(); iter != $1.end(); ++iter) {
    surf_id = iter->first;
    surf = iter->second;
    PyObject* value =
         SWIG_NewPointerObj(SWIG_as_voidptr(surf),
                            $descriptor(Surface*), 0);
    PyDict_SetItem($result, PyInt_FromLong(surf_id), value);
  }
}


/* Typemap for all methods which return a std::map<int, Halfspace*>.
 * This includes the Cell::getSurfaces() method, which is useful for OpenMC
 * compatibility. */
%clear std::map<int, Halfspace*>;
%typemap(out) std::map<int, Halfspace*> {
   $result = PyDict_New();
  int size = $1.size();
   std::map<int, Halfspace*>::iterator iter;
  Halfspace* surf;
  int surf_id;
   for (iter = $1.begin(); iter != $1.end(); ++iter) {
    surf_id = iter->first;
    surf = iter->second;
    PyObject* value =
         SWIG_NewPointerObj(SWIG_as_voidptr(surf),
                            $descriptor(Halfspace*), 0);
    PyDict_SetItem($result, PyInt_FromLong(surf_id), value);
  }
}



/* Typemap for all methods which return a std::map<int, Material*>.
 * This includes the Geometry::getAllMaterials() method, which is useful
 * for OpenMC compatibility. */
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
 * OpenMC compatibility. */
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


/* Typemap for Lattice::setUniverses(int num_z, int num_y, int num_x,
 *                                     Universe** universes)
 * method - allows users to pass in Python lists of Universes for each
 * lattice cell */
%typemap(in) (int num_z, int num_y, int num_x, Universe** universes) {

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a Python list of integers "
                    "for the Lattice cells");
    return NULL;
  }

  $1 = PySequence_Length($input);  // num_z
  $2 = PySequence_Length(PyList_GetItem($input,0)); // num_y
  $3 = PySequence_Length(PyList_GetItem(PyList_GetItem($input,0), 0)); // num_x
  $4 = (Universe**) malloc(($1 * $2 * $3) * sizeof(Universe*)); // universes

  /* Loop over the xy-planes */
  for (int k = 0; k < $1; k++) {

    /* Get the 2D list of universes in the k-th xy-plane */
    PyObject* outer_outer_list = PyList_GetItem($input,k);

    /* Loop over y */
    for (int j = 0; j < $2; j++) {

      /* Get the list of universes in the j-th row of the k-th xy-plane */
      PyObject* outer_list = PyList_GetItem(outer_outer_list, j);

      /* Check that the number of universes in the j-th row of the k-th xy-plane
       * is the same as the number of universes in the 1st row of the 1st
       * xy-plane */
      if (PySequence_Length(outer_list) != $3) {
        PyErr_SetString(PyExc_ValueError, "Size mismatch in dimensions of 3D "
                        "list of Universes in input to Lattice:setUniverses"
                        " method");
        return NULL;
      }

      /* Loop over universes in j-th row of the k-th xy-plane */
      for (int i =0; i < $3; i++) {
        /* Extract the value from the list at this location and convert
         * SWIG wrapper to pointer to underlying C++ class instance */
        PyObject* o = PyList_GetItem(outer_list, i);
        void *p1 = 0;
        SWIG_ConvertPtr(o, &p1, SWIGTYPE_p_Universe, 0 | 0);
        $4[k*($2*$3) + j*$3 + i] = (Universe*) p1;
      }
    }
  }
}
