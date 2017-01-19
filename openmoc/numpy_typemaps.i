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
 * getter method for rotations used by the OpenCG compatibility module. */
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double* rotations, int num_axes)}

/* The typemap used to match the method signature for the Cell's
 * getter method for translations used by the OpenCG compatibility module. */
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

/* The typemap used to match the method signature for Solver::getFluxes */
%apply (FP_PRECISION* ARGOUT_ARRAY1, int DIM1) {(FP_PRECISION* out_fluxes, int num_fluxes)}

/* The typemap used to match the method signature for Solver::setFluxes */
%apply (FP_PRECISION* INPLACE_ARRAY1, int DIM1) {(FP_PRECISION* in_fluxes, int num_fluxes)}
