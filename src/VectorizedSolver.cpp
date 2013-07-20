#include "VectorizedSolver.h"


/**
 * @brief Constructor initializes empty arrays for source, flux, etc.
 * @details The construcor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator, and uses this to initialie empty arrays for the 
 *          flat source regions, boundary angular fluxes, scalar flatsourcergion
 *          fluxes, flatsourceregion sources and flatsourceregion powers. The 
 *          constructor initalizes the number of threads to a default of 1.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the trackgenerator
 */
VectorizedSolver::VectorizedSolver(Geometry* geometry, 
				   TrackGenerator* track_generator) :
    CPUSolver(geometry, track_generator) { 
  
    _thread_taus = NULL;
    _thread_exponentials = NULL;

    if (geometry != NULL)
        setGeometry(geometry);

    if (track_generator != NULL)
        setTrackGenerator(track_generator);

    vmlSetMode(VML_EP);
}


/**
 * @brief Destructor deletes arrays of boundary angular flux for all tracks,
 *        scalar flux and source for each flat source region.
 */
VectorizedSolver::~VectorizedSolver() {

    if (_boundary_flux != NULL) {
        _mm_free(_boundary_flux);
        _boundary_flux = NULL;
    }

    if (_boundary_leakage != NULL) {
        _mm_free(_boundary_leakage);
        _boundary_leakage = NULL;
    }

    if (_scalar_flux != NULL) {
        _mm_free(_scalar_flux);
	_scalar_flux = NULL;
    }

    if (_thread_fsr_flux != NULL) {
        _mm_free(_thread_fsr_flux);
        _thread_fsr_flux = NULL;
    }

    if (_fission_sources != NULL) {
        _mm_free(_fission_sources);
	_fission_sources = NULL;
    }

    if (_scatter_sources != NULL) {
        _mm_free(_scatter_sources);
	_scatter_sources = NULL;
    }

    if (_source != NULL) {
        _mm_free(_source);
	_source = NULL;
    }

    if (_old_source != NULL) {
        _mm_free(_old_source);
	_old_source = NULL;
    }

    if (_reduced_source != NULL) {
        _mm_free(_reduced_source);
	_reduced_source = NULL;
    }

    if (_thread_taus != NULL) {
        _mm_free(_thread_taus);
        _thread_taus = NULL;
    }

    if (_thread_exponentials != NULL) {
        _mm_free(_thread_exponentials);
        _thread_exponentials = NULL;
    }

}


/**
 * @brief Returns the number of vector lengths required to fit the number
 *        of energy groups.
 * @return The number of vector widths
 */
int VectorizedSolver::getNumVectorWidths() {
    return _num_vector_lengths;
}


/**
 * @brief Sets the geometry for the solver and aligns all material cross-section
 *        data for SIMD vector instructions.
 * @param geometry a pointer to the geometry
 */
void VectorizedSolver::setGeometry(Geometry* geometry) {

    CPUSolver::setGeometry(geometry);

    /* Compute the number of SIMD vector widths needed to fit energy groups */
    _num_vector_lengths = (_num_groups / VEC_LENGTH) + 1;

    /* Reset the number of energy groups by rounding up for the number
     * of vector widths needed to accomodate the energy groups */
    _num_groups = _num_vector_lengths * VEC_LENGTH;

    _polar_times_groups = _num_groups * _num_polar;

    std::map<short int, Material*> materials = geometry->getMaterials();
    std::map<short int, Material*>::iterator iter;

    /* Iterate over each material and replace it's xs with a new one 
     * array that is a multiple of VEC_LENGTH long */
    for (iter=materials.begin(); iter != materials.end(); ++iter)
        (*iter).second->alignData();
}


/**
 * @brief Allocates memory for the exponential prefactor table
 */
void VectorizedSolver::precomputePrefactors() {

    CPUSolver::precomputePrefactors();

    /* Deallocates memory for the exponentials if it was allocated for a
     * previous simulation */
    if (_thread_exponentials != NULL)
        _mm_free(_thread_exponentials);

    /* Allocates memory for an array of exponential values for each thread
     * - this is not used by default, but can be to allow for vectorized
     * evaluation of the exponentials. Unfortunately this does not appear
     * to give any performance boost. */
    int size = _num_threads * _polar_times_groups * sizeof(FP_PRECISION);
    _thread_exponentials = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
}



/**
 * @brief Allocates memory for track boundary angular fluxes and 
 *        flat source region scalar fluxes and leakages.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void VectorizedSolver::initializeFluxArrays() {
   
    /* Delete old flux arrays if they exist */
    if (_boundary_flux != NULL)
        _mm_free(_boundary_flux);

    if (_boundary_leakage != NULL)
        _mm_free(_boundary_leakage);

    if (_scalar_flux != NULL)
        _mm_free(_scalar_flux);

    if (_thread_fsr_flux != NULL)
        _mm_free(_thread_fsr_flux);

    if (_thread_taus != NULL)
        _mm_free(_thread_taus);

    int size;

    /* Allocate aligned memory for all flux arrays */
    try{

        size = 2 * _tot_num_tracks * _num_groups * _num_polar;
	size *= sizeof(FP_PRECISION);
	_boundary_flux = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_boundary_leakage = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);

	size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
	_scalar_flux = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);

	size = _num_threads * _num_groups * sizeof(FP_PRECISION);
	_thread_fsr_flux = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);

	size = _num_threads * _polar_times_groups * sizeof(FP_PRECISION);
	_thread_taus = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flat source region source arrays.
 * @details Deletes memory for old source arrays if they were allocated from
 *          previous simulation.
 */
void VectorizedSolver::initializeSourceArrays() {

    /* Delete old sources arrays if they exist */
    if (_fission_sources != NULL)
        _mm_free(_fission_sources);

    if (_scatter_sources != NULL)
        _mm_free(_scatter_sources);

    if (_source != NULL)
        _mm_free(_source);

    if (_old_source != NULL)
        _mm_free(_old_source);

    if (_reduced_source != NULL)
        _mm_free(_reduced_source);

    if (_source_residuals != NULL)
        _mm_free(_source_residuals);

    int size;

    /* Allocate aligned memory for all source arrays */
    try{
        size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
	_fission_sources = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_scatter_sources = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_source = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_old_source = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_reduced_source = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_source_residuals = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source region sources array. Backtrace:%s", e.what());
    }
}


/**
 * @brief Normalizes all flat source region scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times \f$ \nu \f$).
 */
void VectorizedSolver::normalizeFluxes() {

    double* nu_sigma_f;
    FP_PRECISION volume;
    FP_PRECISION tot_fission_source;
    FP_PRECISION norm_factor;

    /* Compute total fission source for each region, energy group */
    #pragma omp parallel for private(volume, nu_sigma_f) \
      reduction(+:tot_fission_source) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get pointers to important data structures */
	nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
	volume = _FSR_volumes[r];

	/* Loop over energy group vector lengths */
	for (int v=0; v < _num_vector_lengths; v++) {

            /* Loop over each energy group within this vector */
            #pragma simd vectorlength(VEC_LENGTH)
            for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) {
	        _fission_sources(r,e) = nu_sigma_f[e] * _scalar_flux(r,e);
	        _fission_sources(r,e) *= volume;
            }
        }
    }

    /* Compute the total fission source */
    int size = _num_FSRs * _num_groups;
    #ifdef SINGLE
    tot_fission_source = cblas_sasum(size, _fission_sources, 1);
    #else
    tot_fission_source = cblas_dasum(size, _fission_sources, 1);
    #endif

    /* Compute the normalization factor */
    norm_factor = 1.0 / tot_fission_source;

    log_printf(DEBUG, "tot fiss src = %f, Normalization factor = %f", 
               tot_fission_source, norm_factor);

    /* Normalize the flat source region scalar fluxes */
    #ifdef SINGLE
    cblas_sscal(size, norm_factor, _scalar_flux, 1);
    #else
    cblas_dscal(size, norm_factor, _scalar_flux, 1);
    #endif

    /* Normalize the boundary flux */
    size = 2 * _tot_num_tracks * _num_polar * _num_groups;

    #ifdef SINGLE
    cblas_sscal(size, norm_factor, _boundary_flux, 1);
    #else
    cblas_dscal(size, norm_factor, _boundary_flux, 1);
    #endif

    return;
}


/**
 * @brief Computes the total source (fission and scattering) in each flat 
 *        source region.
 * @details This method computes the total source in each region based on
 *          this iteration's current approximation to the scalar flux. A
 *          residual for the source with respect to the source compute on
 *          the previous iteration is computed and returned. The residual
 *          is determined as follows:
 *          /f$ res = \sqrt{\frac{\displaystyle\sum \displaystyle\sum 
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{# FSRs}}} \f$
 *
 * @return the residual between this source and the previous source
 */
FP_PRECISION VectorizedSolver::computeFSRSources() {

    FP_PRECISION scatter_source;
    FP_PRECISION fission_source;
    double* nu_sigma_f;
    double* sigma_s;
    double* sigma_t;
    double* chi;
    Material* material;

    FP_PRECISION source_residual = 0.0;

    /* For all regions, find the source */
    #pragma omp parallel for private(material, nu_sigma_f, chi,	\
      sigma_s, sigma_t, fission_source, scatter_source) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        material = _FSR_materials[r];
	nu_sigma_f = material->getNuSigmaF();
	chi = material->getChi();
	sigma_s = material->getSigmaS();
        sigma_t = material->getSigmaT();

	for (int v=0; v < _num_vector_lengths; v++) {

	    /* Compute fission source for each group */
            #pragma simd vectorlength(VEC_LENGTH)
            for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
	        _fission_sources(r,e) = _scalar_flux(r,e) * nu_sigma_f[e];
        }

        #ifdef SINGLE
        fission_source = cblas_sasum(_num_groups, &_fission_sources(r,0), 1);
        #else
        fission_source = cblas_dasum(_num_groups, &_fission_sources(r,0), 1);
        #endif

	/* Compute total scattering source for group G */
        for (int G=0; G < _num_groups; G++) {
            scatter_source = 0;

	    for (int v=0; v < _num_vector_lengths; v++) {

                #pragma simd vectorlength(VEC_LENGTH)
                for (int g=v*VEC_LENGTH; g < (v+1)*VEC_LENGTH; g++)
		    _scatter_sources(r,g) = sigma_s[G*_num_groups+g]*_scalar_flux(r,g);
            }

            #ifdef SINGLE
	    scatter_source = cblas_sasum(_num_groups, &_scatter_sources(r,0), 1);
            #else
	    scatter_source = cblas_dasum(_num_groups, &_scatter_sources(r,0), 1);
            #endif

	    /* Set the total source for region r in group G */
	    _source(r,G) = ((1.0 / _k_eff) * fission_source *
                           chi[G] + scatter_source) * ONE_OVER_FOUR_PI;

	    _reduced_source(r,G) = _source(r,G) / sigma_t[G];

	    /* Compute the norm of residual of the source in the region, group */
	    if (fabs(_source(r,G)) > 1E-10)
	        _source_residuals(r,G) = pow((_source(r,G) - _old_source(r,G)) 
					    / _source(r,G), 2);
	    
	    /* Update the old source */
	    _old_source(r,G) = _source(r,G);
        }
    }

    /* Sum up the residuals from each group and in each region */
    #ifdef SINGLE
    source_residual = cblas_sasum(_num_FSRs * _num_groups, _source_residuals, 1);
    #else
    source_residual = cblas_dasum(_num_FSRs * _num_groups, _source_residuals, 1);
    #endif

    source_residual = sqrt(source_residual / _num_FSRs);

    return source_residual;
}



/**
 * @brief Add the source term contribution in the transport equation to 
 *        the flat source region scalar flux
 */
void VectorizedSolver::addSourceToScalarFlux() {

    FP_PRECISION volume;
    double* sigma_t;

    /* Add in source term and normalize flux to volume for each region */
    /* Loop over flat source regions, energy groups */
    #pragma omp parallel for private(volume, sigma_t) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        volume = _FSR_volumes[r];
	sigma_t = _FSR_materials[r]->getSigmaT();

	/* Loop over each energy group vector length */
	for (int v=0; v < _num_vector_lengths; v++) {

	    /* Loop over energy groups within this vector */
            #pragma simd vectorlength(VEC_LENGTH)
            for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) {
                _scalar_flux(r,e) *= 0.5;
		_scalar_flux(r,e) = FOUR_PI * _reduced_source(r,e) + 
		  (_scalar_flux(r,e) / (sigma_t[e] * volume));
	    }
        }
    }
    
    return;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from the total fission and absorption rates.
 * @details This method computes the current approximation to the 
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum \displaystyle\sum \nu
 *                        \Sigma_f \Phi V}{\displaystyle\sum 
 *                        \displaystyle\sum \Sigma_a \Phi V} \f$
 */
void VectorizedSolver::computeKeff() {

    Material* material;
    double* sigma_a;
    double* nu_sigma_f;
    FP_PRECISION volume;

    double tot_abs = 0.0;
    double tot_fission = 0.0;

    FP_PRECISION* absorption_rates = new FP_PRECISION[_num_FSRs*_num_groups];
    FP_PRECISION* fission_rates = new FP_PRECISION[_num_FSRs*_num_groups];

    /* Loop over all flat source regions and compute the volume-weighted
     * fission and absorption rates */
    #pragma omp parallel for private(volume, material, \
      sigma_a, nu_sigma_f) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        volume = _FSR_volumes[r];
	material = _FSR_materials[r];
	sigma_a = material->getSigmaA();
	nu_sigma_f = material->getNuSigmaF();

	/* Loop over each energy group vector length */
	for (int v=0; v < _num_vector_lengths; v++) {

	    /* Loop over energy groups within this vector */
            #pragma simd vectorlength(VEC_LENGTH)
	    for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) {
	        absorption_rates[r*_num_groups+e] = sigma_a[e] * _scalar_flux(r,e);
		fission_rates[r*_num_groups+e] = nu_sigma_f[e] * _scalar_flux(r,e);
		absorption_rates[r*_num_groups+e] *= volume;
		fission_rates[r*_num_groups+e] *= volume;
	    }
        }
    }

    /* Reduce absorption and fission rates across FSRs, energy groups */
    int size = _num_FSRs * _num_groups;

    #ifdef SINGLE
    tot_abs = cblas_sasum(size, absorption_rates, 1);
    tot_fission = cblas_sasum(size, fission_rates, 1);
    #else
    tot_abs = cblas_dasum(size, absorption_rates ,1);
    tot_fission = cblas_dasum(size, fission_rates, 1);
    #endif

    /** Reduce leakage array across tracks, energy groups, polar angles */
    size = 2 * _tot_num_tracks * _polar_times_groups;

    #ifdef SINGLE
    _leakage = cblas_sasum(size, _boundary_leakage, 1) * 0.5;
    #else
    _leakage = cblas_dasum(size, _boundary_leakage, 1) * 0.5;
    #endif

    _k_eff = tot_fission / (tot_abs + _leakage);

    log_printf(DEBUG, "abs = %f, fission = %f, leakage = %f, "
	       "k_eff = %f", tot_abs, tot_fission, _leakage, _k_eff);

    delete [] absorption_rates;
    delete [] fission_rates;

    return;
}



/**
 * @brief Computes the contribution to the flat source region scalar flux
 *        from a single track segment.
 * @details This method integrates the angular flux for a track segment across
 *        energy groups and polar angles, and tallies it into the flat
 *        source region scalar flux, and updates the track's angular flux.
 * @param curr_segment a pointer to the segment of interest
 * @param track_flux a pointer to the track's angular flux
 * @param fsr_flux a pointer to the temporary flat source region flux buffer
 */
void VectorizedSolver::scalarFluxTally(segment* curr_segment,
   	                               FP_PRECISION* track_flux,
	                               FP_PRECISION* fsr_flux){

    int tid = omp_get_thread_num();
    int fsr_id = curr_segment->_region_id;
    FP_PRECISION length = curr_segment->_length;
    double* sigma_t = curr_segment->_material->getSigmaT();

    /* The average flux along this segment in the flat source region */
    FP_PRECISION psibar;
    FP_PRECISION* exponentials = &_thread_exponentials[tid*_polar_times_groups];

    computeExponentials(curr_segment, exponentials);

    /* Set the flat source region flux buffer to zero */
    memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

    /* Tally the flux contribution from segment to FSR's scalar flux */
    /* Loop over polar angles */
    for (int p=0; p < _num_polar; p++){

        /* Loop over each energy group vector length */
        for (int v=0; v < _num_vector_lengths; v++) {

	    /* Loop over energy groups within this vector */
            #pragma simd vectorlength(VEC_LENGTH) private(psibar)
            for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) {
	        psibar = (track_flux(p,e) - _reduced_source(fsr_id,e)) * 
		          exponentials(p,e);
	        fsr_flux[e] += psibar * _polar_weights[p];
		track_flux(p,e) -= psibar;
	    }
	}
    }

    /* Atomically increment the FSR scalar flux from the temporary array */
    omp_set_lock(&_FSR_locks[fsr_id]);
    {
        #ifdef SINGLE
        vsAdd(_num_groups, &_scalar_flux(fsr_id,0), fsr_flux, 
	      &_scalar_flux(fsr_id,0));
        #else
        vdAdd(_num_groups, &_scalar_flux(fsr_id,0), fsr_flux, 
	      &_scalar_flux(fsr_id,0));
        #endif
    }
    omp_unset_lock(&_FSR_locks[fsr_id]);

    return;
}


/**
 * @brief Computes an array of the exponential in the transport equation,
 *        \f$ exp(-\frac{\Sigma_t * l}{sin(\theta)}) \f$, for each energy group
 *        and polar angle for a given segment.
 * @param curr_segment pointer to the segment of interest
 * @param exponentials the array to store the exponential values
 */
void VectorizedSolver::computeExponentials(segment* curr_segment, 
	                                   FP_PRECISION* exponentials) {

    FP_PRECISION length = curr_segment->_length;
    double* sigma_t = curr_segment->_material->getSigmaT();

    /* Evaluate the exponentials using the lookup table - linear interpolation */
    if (_interpolate_exponential) {
        FP_PRECISION tau;
        int index;

	for (int e=0; e < _num_groups; e++) {
	    for (int p=0; p < _num_polar; p++) {
	        tau = sigma_t[e] * length;
		index = int(tau * _inverse_prefactor_spacing) * 
		        _two_times_num_polar;
		exponentials(p,e) = (1. - 
				     (_prefactor_array[index+2 * p] * tau + 
				      _prefactor_array[index + 2 * p +1]));
	    }
        }
    }

    /* Evalute the exponentials using the intrinsic exp function */
    else {

	int tid = omp_get_thread_num();

        FP_PRECISION* sinthetas = _quad->getSinThetas();
	FP_PRECISION* taus = &_thread_taus[tid*_polar_times_groups];

	/* Initialize the tau argument for the exponentials */
	for (int p=0; p < _num_polar; p++) {

	    for (int v=0; v < _num_vector_lengths; v++) {
	        
                #pragma simd vectorlength(VEC_LENGTH)
	        for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) 
		    taus(p,e) = -sigma_t[e] * length / sinthetas[p];
	    }
	}

        /* Evaluate the negative of the exponentials using Intel's MKL */
        #ifdef SINGLE
	vsExp(_polar_times_groups, taus, exponentials);
        #else
	vdExp(_polar_times_groups, taus, exponentials);
        #endif

	/* Compute one minus the exponentials */
	for (int p=0; p < _num_polar; p++) {

	    for (int v=0; v < _num_vector_lengths; v++) {
	        
                #pragma simd vectorlength(VEC_LENGTH)
	        for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) 
		    exponentials(p,e) = 1.0 - exponentials(p,e);
	    }
	}
    }
}


/**
 * @brief Updates the boundary flux for a track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the track is given to the reflecting track. For vacuum
 *          boundary conditions, the outgoing flux tallied as leakage.
 * @param track_id the ID number for the track of interest
 * @param direction the track direction (forward - true, reverse - false)
 * @param track_flux a pointer to the track's outgoing angular flux
 */
void VectorizedSolver::transferBoundaryFlux(int track_id, 
					    bool direction,
					    FP_PRECISION* track_flux) {
    int start;
    bool bc;
    FP_PRECISION* track_leakage;
    int track_out_id;

    /* Extract boundary conditions for this track and the pointer to the 
     * outgoing reflective track, and index into the leakage array */

    /* For the "forward" direction */
    if (direction) {
        start = _tracks[track_id]->isReflOut() * _polar_times_groups;
        track_leakage = &_boundary_leakage(track_id,0);
        track_out_id = _tracks[track_id]->getTrackOut()->getUid();
        bc = _tracks[track_id]->getBCOut();
    }

    /* For the "reverse" direction */
    else {
        start = _tracks[track_id]->isReflIn() * _polar_times_groups;
        track_leakage = &_boundary_leakage(track_id,_polar_times_groups);
        track_out_id = _tracks[track_id]->getTrackIn()->getUid();
        bc = _tracks[track_id]->getBCIn();
    }

    FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,0,0,start);

    /* Loop over polar angles and energy groups */
    for (int p=0; p < _num_polar; p++) {

        /* Loop over each energy group vector length */
        for (int v=0; v < _num_vector_lengths; v++) {

	    /* Loop over energy groups within this vector */
            #pragma simd vectorlength(VEC_LENGTH)
	    for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) {
	        track_out_flux(p,e) = track_flux(p,e) * bc;
		track_leakage(p,e) = track_flux(p,e) * 
		                     _polar_weights[p] * (!bc);
	    }
	}
    }
}
