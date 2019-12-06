/** Helper typemaps to enable array-based user input to be defined with
 *  NumPy arrays for the default build configuration of OpenMOC */

%module numpy_typemaps

%include "numpy.i"

%init %{
  import_array();
%}

/* The typemap used to match the method signature for the
 * Lattice::setLatticeCells setter method. This allows users to set the lattice
 * cells (universe IDs) using a 2D NumPy array */
%apply (int DIM1, int DIM2, int* IN_ARRAY2) {(int num_x, int num_y, int* universes)}

/* The typemap used to match the method signature for the Material
 * cross-section setter methods. This allows users to set the cross-sections
 * using NumPy arrays */
%apply (double* IN_ARRAY1, int DIM1) {(double* xs, int num_groups)}

/* The typemap used to match the method signature for the Cell rotation
 * angle setter method. This allows users to set the rotation angles
 * using NumPy arrays */
%apply (double* IN_ARRAY1, int DIM1) {(double* rotation, int num_axes)}

/* The typemap used to match the method signature for the Cell translation
 * setter method. This allows users to set translations using NumPy arrays */
%apply (double* IN_ARRAY1, int DIM1) {(double* translation, int num_axes)}

/* The typemap used to match the method signature for the Cell's
 * getter method for rotations used by the OpenMC compatibility module. */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* rotations, int num_axes)}

/* The typemap used to match the method signature for the Cell's
 * getter method for translations used by the OpenMC compatibility module. */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* translations, int num_axes)}

/* The typemap used to match the method signature for the TrackGenerator's
 * getter methods for track start and end coordinates for the plotting
 * routines in openmoc.plotter */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* coords, long num_tracks)}

/* The typemap used to match the method signature for the TrackGenerator's
 * getter methods for track segment start and end coordinates for the plotting
 * routines in openmoc.plotter */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* coords, long num_segments)}

/* The typemap used to match the method signature for the Solver's
 * computeFSRFissionRates method for the data processing routines in
 * openmoc.process */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* fission_rates, long num_FSRs)}

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

/* The typemap used to match the method signature for Geometry::loadSPHFactors */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* sph_factors, int num_domains_groups),
                                           (double* sph_to_domain_ids, int num_sph_domains)}

/* The typemap used to match the method signature for Solver::getFluxes */
%apply (FP_PRECISION* ARGOUT_ARRAY1, int DIM1) {(FP_PRECISION* out_fluxes, int num_fluxes)}

/* The typemap used to match the method signature for Solver::setFluxes */
%apply (FP_PRECISION* INPLACE_ARRAY1, int DIM1) {(FP_PRECISION* in_fluxes, int num_fluxes)}

/* The typemap used to match the method signature for Mesh::getFormattedReactionRates */
%typemap(out) std::vector<std::vector<std::vector<FP_PRECISION> > >& 
{
  for(int i = 0; i < $1->size(); ++i) {
    for(int j = 0; j < $1->data()[i].size(); ++j) {
      int subLength = $1->data()[i]->data()[j].size();
      npy_doublep dims[] = { subLength };
      PyObject* temp = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, $1->data()[i]->data()[j].data());
      $result = SWIG_Python_AppendOutput($result, temp);
    }
  }
}

/* The typemap to match the method signature for RuntimeParameters::setRuntimeParameters */
//NOTE Make sure to use UTF8 encoding with Python strings
%typemap(in) (int argc, char *argv[]) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $1 = PyList_Size($input);
  $2 = (char **) malloc(($1+1)*sizeof(char *));
  for (i = 0; i < $1; i++) {
    PyObject *s = PyList_GetItem($input,i);
    //To know type of data: std::cout << pytype_string(s) << std::endl;
    if (!PyString_Check(s)) {
        free($2);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[i] = PyString_AsString(s);
  }
  $2[i] = 0;
}
