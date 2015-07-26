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
