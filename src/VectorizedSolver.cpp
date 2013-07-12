#include "CPUSolver.h"


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
CPUSolver::CPUSolver(Geometry* geom, TrackGenerator* track_generator) :
  Solver(geom, track_generator) {

    setNumThreads(1);
}


/**
 * @brief Destructor deletes arrays of boundary angular flux for all tracks,
 *        scalar flux and source for each flatsourceregion.
 */
CPUSolver::~CPUSolver() { 

  if (_FSR_locks != NULL)
      delete [] _FSR_locks;

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

  if (_old_scalar_flux != NULL) {
      _mm_free(_old_scalar_flux);
      _old_scalar_flux = NULL;
  }

  if (_old_scalar_flux != NULL) {
      _mm_free(_old_scalar_flux);
      _old_scalar_flux = NULL;
  }

  if (_fission_source != NULL) {
      _mm_free(_fission_source);
      _fission_source = NULL;
  }

  if (_source != NULL) {
      _mm_free(_source);
      _source = NULL;
  }

  if (_old_source != NULL) {
      _mm_free(_old_source);
      _old_source = NULL;
  }

  if (_ratios != NULL) {
      _mm_free(_ratios);
      _ratios = NULL;
  }
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int CPUSolver::getNumThreads() {
    return _num_threads;
}


/**
 * @brief Returns the number of energy groups divided by the vector width
 * @details The vector width is defined by VEC_LENGTH and is used to 
 *          for alignment of data structures for SIMD vector instructions.
 * @return the number energy groups divided by the ector width
 */
int CPUSolver::getNumGroupVectorWidths() {
    return _num_groups_vec;
}


FP_PRECISION CPUSolver::getFSRScalarFlux(int fsr_id, int energy_group) {

    /* Error checking */
    if (fsr_id >= _num_FSRs)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		 "in enery group %d since the solver only contains FSR with "
		   "IDs greater than or equal to %d", 
		   fsr_id, energy_group, _num_FSRs-1);
    if (fsr_id < 0)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		  "in energy group %d since FSRs do not have negative IDs", 
		  fsr_id, energy_group);
    if (energy_group-1 >= _num_groups)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		   "in energy group %d since the solver only has %d energy "
		   "groups", fsr_id, energy_group, _num_groups);
    if (energy_group <= 0)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		 "in energy group %d since energy groups are greater than 1",
		 fsr_id, energy_group);

    return _scalar_flux(fsr_id,energy_group);
}

/**
 * @brief Return a 2D array indexed by flatsourceregion IDs and energy groups 
 *        which contains the corresponding fluxes for each flatsourceregion.
 * @return a 2D array of flatsourceregion scalar fluxes
 */
FP_PRECISION* CPUSolver::getFSRScalarFluxes() {
    if (_scalar_flux == NULL)
        log_printf(ERROR, "Unable to returns the Solver's scalar flux array "
		 "since it has not yet been allocated in memory");

    return _scalar_flux;
}


/**
 * @brief Return an array indexed by flatsourceregion IDs with the
 *        corresponding flatsourceregion power.
 * @return an array of flatsourceregion powers
 */
FP_PRECISION* CPUSolver::getFSRPowers() {
    if (_FSRs_to_powers == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR power array "
		 "since it has not yet been allocated in memory");

    return _FSRs_to_powers;
}


/**
 * @brief Return an array indexed by flatsourceregion IDs with the
 *        corresponding pin cell power.
 * @return an array of flatsourceregion pin powers
 */
FP_PRECISION* CPUSolver::getFSRPinPowers() {
    if (_FSRs_to_pin_powers == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR pin power array "
		 "since it has not yet been allocated in memory");

    return _FSRs_to_pin_powers;
}



/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
 * @details This method sets the number of threads to be no larger than the
 *          number of complementary angle pairs - reflecting angles which form
 *          cyclical tracks - since that is how OpenMOC parallizes loops.
 * @param num_threads the number of threads
 */
void CPUSolver::setNumThreads(int num_threads) {
    if (num_threads <= 0)
        log_printf(ERROR, "Unable to set the number of threads for the Solver "
		   "to %d since it is less than or equal to 0", num_threads);

    _num_threads = num_threads;

    /* Set the number of threads for OpenMP */
    omp_set_num_threads(_num_threads);
}


void CPUSolver::setGeometry(Geometry* geometry) {
    Solver::setGeometry(geometry);
    _num_groups_vec = ceil(_num_groups / VEC_LENGTH);
}


/**
 * @brief Allocates memory for track boundary angular fluxes and 
 *        flatsourceregion scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void CPUSolver::initializeFluxArrays() {
   
    /* Delete old flux arrays if they exist */
    if (_boundary_flux != NULL)
        _mm_free(_boundary_flux);

    if (_boundary_leakage != NULL)
        _mm_free(_boundary_leakage);

    if (_scalar_flux != NULL)
        _mm_free(_scalar_flux);

    if (_old_scalar_flux != NULL)
        _mm_free(_old_scalar_flux);

    int size;

    /* Allocate aligned memory for all flux arrays */
    try{

        size = 2 * _tot_num_tracks * _polar_times_groups 
	        * sizeof(FP_PRECISION);
	_boundary_flux = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_boundary_leakage = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);

	size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
	_scalar_flux = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_old_scalar_flux = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flatsourceregion source arrays.
 * @details Deletes memory for old source arrays if they were allocated from
 *          previous simulation.
 */
void CPUSolver::initializeSourceArrays() {

    /* Delete old sources arrays if they exist */
    if (_fission_source != NULL)
        _mm_free(_fission_source);

    if (_source != NULL)
        _mm_free(_source);

    if (_old_source != NULL)
        _mm_free(_old_source);

    if (_ratios != NULL)
        _mm_free(_ratios);

    int size;

    /* Allocate aligned memory for all source arrays */
    try{
        size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
	_fission_source = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_source = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_old_source = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
	_ratios = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source region sources array. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flatsourceregion power arrays.
 * @details Deletes memory for power arrays if they were allocated from
 *          previous simulation.
 */
void CPUSolver::initializePowerArrays() {

    /* Delete old power arrays if they exist */
    if (_FSRs_to_powers != NULL)
        delete [] _FSRs_to_powers;

    if (_FSRs_to_pin_powers != NULL)
        delete [] _FSRs_to_pin_powers;

    /* Allocate memory for FSR power and pin power arrays */
    try{
	_FSRs_to_powers = new FP_PRECISION[_num_FSRs];
	_FSRs_to_pin_powers = new FP_PRECISION[_num_FSRs];
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's FSR "
		   "power arrays. Backtrace:%s", e.what());
    }
}


/**
 * @brief Creates a polar quadrature object for the solver.
 */
void CPUSolver::initializePolarQuadrature() {
    /* Deletes the old quadrature if one existed */
    if (_quad != NULL)
        delete _quad;

    _quad = new Quadrature(_quadrature_type, _num_polar);
    _polar_times_groups = _num_groups * _num_polar;
}

 
/**
 * @brief Pre-computes exponential pre-factors for each segment of each track 
 *        for each polar angle. 
 * @details This method will generate a hashmap which contains values of the 
 *          pre-factor for specific segment lengths (the keys into the hashmap).
 */
void CPUSolver::precomputePrefactors() {

    /* Build exponential prefactor array based on table look up with linear 
     * interpolation */
    log_printf(INFO, "Building exponential prefactor hashtable...");

    FP_PRECISION azim_weight;

    _polar_weights = new FP_PRECISION[_num_azim*_num_polar];

    /* Precompute the total azimuthal weight for tracks at each polar angle */
    #pragma omp parallel for private (azim_weight)
    for (int i=0; i < _num_azim; i++) {
        azim_weight = _azim_weights[i];

        for (int p=0; p < _num_polar; p++)
	    _polar_weights(i,p) = azim_weight*_quad->getMultiple(p)*FOUR_PI;
    }

    /* Set size of prefactor array */
    int num_array_values = 10 * sqrt(1. / (8. * _source_convergence_thresh));
    _prefactor_spacing = 10. / num_array_values;
    _prefactor_array_size = _two_times_num_polar * num_array_values;
    _prefactor_max_index = _prefactor_array_size - _two_times_num_polar - 1.;
    
    log_printf(DEBUG, "Prefactor array size: %i, max index: %i",
	       _prefactor_array_size, _prefactor_max_index);

    /* allocate arrays */
    _prefactor_array = new FP_PRECISION[_prefactor_array_size];

    FP_PRECISION expon;
    FP_PRECISION intercept;
    FP_PRECISION slope;

    /* Create prefactor array */
    for (int i=0; i < num_array_values; i ++){
        for (int p=0; p < _num_polar; p++){
	    expon = exp(- (i * _prefactor_spacing) / _quad->getSinTheta(p));
	    slope = - expon / _quad->getSinTheta(p);
	    intercept = expon * (1 + (i * _prefactor_spacing) /
				 _quad->getSinTheta(p));
	    _prefactor_array[_two_times_num_polar * i + 2 * p] = slope;
	    _prefactor_array[_two_times_num_polar * i + 2 * p + 1] = intercept;
	}
    }

    /* Compute the reciprocal of the prefactor spacing */
    _inverse_prefactor_spacing = 1.0 / _prefactor_spacing;

    return;
}


/**
 * @brief Initializes each of the flatsourceregion objects inside the solver's
 *        array of flatsourceregions. 
 * @details This method assigns each flatsourceregion a unique, monotonically
 *          increasing ID, sets the material for each flatsourceregion, and 
 *          assigns a volume based on the cumulative length of all of the 
 *          segments inside the flatsourceregion.
 */
void CPUSolver::initializeFSRs() {

    log_printf(INFO, "Initializing flat source regions...");

    /* Delete old FSRs array if it exists */
    if (_FSR_volumes != NULL)
        delete [] _FSR_volumes;

    if (_FSR_materials != NULL)
        delete [] _FSR_materials;

    _FSR_volumes = new FP_PRECISION[_num_FSRs];
    _FSR_materials = new Material*[_num_FSRs];
    _FSR_locks = new omp_lock_t[_num_FSRs];

    std::vector<segment*> segments;
    std::vector<segment*>::iterator iter;
    FP_PRECISION volume;
    CellBasic* cell;
    Material* material;
    Universe* univ_zero = _geometry->getUniverse(0);

    /* Initialize the FSR volumes to zero */
    memset(_FSR_volumes, FP_PRECISION(0.), _num_FSRs*sizeof(FP_PRECISION));

    /* Set each FSR's "volume" by accumulating the total length of all tracks
     * inside the FSR. Loop over azimuthal angle, track and segment. 
     * Note: this code region cannot be parallelized without a mutex lock
     * on FSR volume due to race conditions. */
    for (int i=0; i < _tot_num_tracks; i++) {
        
        int azim_index = _tracks[i]->getAzimAngleIndex();
        segments = _tracks[i]->getSegments();
	
	for (iter=segments.begin(); iter != segments.end(); ++iter) {
	    volume = (*iter)->_length * _azim_weights[azim_index];
	    _FSR_volumes[(*iter)->_region_id] += volume;
	}
    }

    /* Loop over all FSRs to extract FSR material pointers */
    #pragma omp parallel for private(cell, material)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get the cell corresponding to this FSR from the geometry */
        cell = static_cast<CellBasic*>(_geometry->findCell(univ_zero, r));

	/* Get the cell's material and assign it to the FSR */
	material = _geometry->getMaterial(cell->getMaterial());
	_FSR_materials[r] = material;

	log_printf(DEBUG, "FSR id = %d has cell id = %d and material id = %d "
                  "and volume = %f", r, cell->getId(), 
                   _FSR_materials[r]->getUid(), _FSR_volumes[r]);
    }

    /* Loop over all FSRs to initialize OpenMP locks */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++)
        omp_init_lock(&_FSR_locks[r]);

    return;
}


/**
 * @brief Zero each track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
void CPUSolver::zeroTrackFluxes() {

    /* Loop over azimuthal angle, track, polar angle, energy group
     * and set each track's incoming and outgoing flux to zero */
    #pragma omp parallel for
    for (int i=0; i < _tot_num_tracks; i++) {
        for (int pe2=0; pe2 < 2*_polar_times_groups; pe2++)
    	    _boundary_flux(i,pe2) = 0.0;
    }

    return;
}


/**
 * @brief Set the scalar flux for each energy group inside each 
 *        flatsourceregion to a constant value.
 * @param value the value to assign to each flat source region flux
 */
void CPUSolver::flattenFSRFluxes(FP_PRECISION value) {

    /* Loop over all FSRs and energy groups */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++) {
	    _scalar_flux(r,e) = value;
	    _old_scalar_flux(r,e) = value;
         }
     }

    return;
}


/**
 * @brief Set the source for each energy group inside each flatsourceregion
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
void CPUSolver::flattenFSRSources(FP_PRECISION value) {

    /* Loop over all FSRs and energy groups */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++) {
            _source(r,e) = value;
    	    _old_source(r,e) = value;
        }
    }

    return;
}


/**
 * @brief Normalizes all flatsourceregion scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times nu).
 */
void CPUSolver::normalizeFluxes() {

    double* nu_sigma_f;
    FP_PRECISION volume;
    FP_PRECISION tot_fission_source;
    FP_PRECISION norm_factor;

    memset(_fission_source, 0, _num_FSRs * _num_groups);

    /* Compute total fission source for each region, energy group */
    #pragma omp parallel for private(volume, nu_sigma_f) \
      reduction(+:tot_fission_source)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get pointers to important data structures */
	nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
	volume = _FSR_volumes[r];

	for (int e=0; e < _num_groups; e++)
	    _fission_source(r,e) = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
    }

    /* Compute the total fission source */
    tot_fission_source = pairwise_sum<FP_PRECISION>(_fission_source, 
						    _num_FSRs*_num_groups);
    /* Normalize scalar fluxes in each region */
    norm_factor = 1.0 / tot_fission_source;

    log_printf(DEBUG, "Normalization factor = %f", norm_factor);

    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++)
	    _scalar_flux(r,e) *= norm_factor;
    }

    /* Normalize angular boundary fluxes for each track */
    #pragma omp parallel for
    for (int i=0; i < _track_generator->getNumTracks(); i++) {
        for (int pe2=0; pe2 < 2*_polar_times_groups; pe2++)
	    _boundary_flux(i,pe2) *= norm_factor;
    }

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
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{# FSRs}} \f$
 *
 * @return the residual between this source and the previous source
 */
FP_PRECISION CPUSolver::computeFSRSources() {

    FP_PRECISION scatter_source;
    FP_PRECISION fission_source;
    double* nu_sigma_f;
    double* sigma_s;
    double* sigma_t;
    double* chi;
    Material* material;

    FP_PRECISION* source_residuals = new FP_PRECISION[_num_groups*_num_FSRs];
    FP_PRECISION source_residual = 0.0;

    /* For all regions, find the source */
    //TODO: This can be parallelized! Need to privatize some variable and
    //      reduce the residual at the end
    #pragma omp parallel for private(material, nu_sigma_f, chi, \
      sigma_s, sigma_t, fission_source, scatter_source)
    for (int r=0; r < _num_FSRs; r++) {

        FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];
        FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];

        material = _FSR_materials[r];
	nu_sigma_f = material->getNuSigmaF();
	chi = material->getChi();
	sigma_s = material->getSigmaS();
        sigma_t = material->getSigmaT();

	/* Compute fission source for each group */
	for (int e=0; e < _num_groups; e++)
	    fission_sources[e] = _scalar_flux(r,e) * nu_sigma_f[e];

	fission_source = pairwise_sum<FP_PRECISION>(fission_sources, 
						    _num_groups);
	
	/* Compute total scattering source for group G */
        for (int G=0; G < _num_groups; G++) {
            scatter_source = 0;

	    for (int g=0; g < _num_groups; g++)
                scatter_sources[g] = sigma_s[G*_num_groups+g] * _scalar_flux(r,g);

	    scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources, 
                                                        _num_groups);

	    /* Set the total source for region r in group G */
	    _source(r,G) = ((1.0 / _k_eff) * fission_source *
                           chi[G] + scatter_source) * ONE_OVER_FOUR_PI;
	

	    _ratios(r,G) = _source(r,G) / sigma_t[G];

	    /* Compute the norm of residual of the source in this region, group */
	    if (fabs(_source(r,G)) > 1E-10)
	      source_residuals(r,G) = pow((_source(r,G) - _old_source(r,G)) 
                                      / _source(r,G), 2);
	    
	    /* Update the old source */
	    _old_source(r,G) = _source(r,G);
        }

	delete [] scatter_sources;
	delete [] fission_sources;
    }

    /* Sum up the residuals from each group and in each region */
    source_residual = pairwise_sum<FP_PRECISION>(source_residuals, 
                                                 _num_FSRs*_num_groups);
    source_residual = sqrt(source_residual / _num_FSRs);

    delete [] source_residuals;

    return source_residual;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from the total fission and absorption rates.
 * @details This method computes the current approximation to the 
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum \displaystyle\sum \nu
 *                        \Sigma_f \Phi V}{\displaystyle\sum 
 *                        \displaystyle\sum \Sigma_a \Phi V} \f$
 */
void CPUSolver::computeKeff() {

    Material* material;
    double* sigma_a;
    double* nu_sigma_f;
    FP_PRECISION volume;

    double tot_abs = 0.0;
    double tot_fission = 0.0;

    FP_PRECISION* absorption_rates = new FP_PRECISION[_num_FSRs*_num_groups];
    FP_PRECISION* fission_rates = new FP_PRECISION[_num_FSRs*_num_groups];

    #pragma omp parallel for private(volume, material, sigma_a, nu_sigma_f)
    for (int r=0; r < _num_FSRs; r++) {

        volume = _FSR_volumes[r];
	material = _FSR_materials[r];
	sigma_a = material->getSigmaA();
	nu_sigma_f = material->getNuSigmaF();

	for (int e=0; e < _num_groups; e++) {
            absorption_rates[r*_num_groups+e] = sigma_a[e] * _scalar_flux(r,e) * volume;
	    fission_rates[r*_num_groups+e] = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
        }
    }

    /* Reduce absorptoin and fission rates across FSRs, energy groups */
    tot_abs = pairwise_sum<FP_PRECISION>(absorption_rates, _num_FSRs*_num_groups);
    tot_fission = pairwise_sum<FP_PRECISION>(fission_rates, _num_FSRs*_num_groups);

    /** Reduce leakage array across tracks, energy groups, polar angles */
    _leakage = pairwise_sum<FP_PRECISION>(_boundary_leakage, 
					  2*_tot_num_tracks*_polar_times_groups);
    _leakage *= 0.5;
    

    _k_eff = tot_fission / (tot_abs + _leakage);

    //    printf("abs = %1.15f, fiss = %1.15f, leak = %1.15f, keff = %1.15f\n", 
    //	   tot_abs, tot_fission, _leakage, _k_eff);

    log_printf(DEBUG, "tot_abs = %f, tot_fission = %f, leakage = %f, "
	       "k_eff = %f", tot_abs, tot_fission, _leakage, _k_eff);

    delete [] absorption_rates;
    delete [] fission_rates;

    return;
}


/**
 * This method performs on or more fixed source iterations by integrating
 * the flux along each track and updating the boundary fluxes for the
 * corresponding output track, while updating the scalar flux in each
 * flat source region
 * @param max_iterations the maximum number of iterations allowed
 */
void CPUSolver::transportSweep() {

    Track* curr_track;
    int num_segments;
    segment* curr_segment;    
    FP_PRECISION* track_flux;

    log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

    /* Initialize flux in each region to zero */
    flattenFSRFluxes(0.0);

    /* Loop over azimuthal angle halfspaces */
    for (int i=0; i < 2; i++) {

        int min = i * (_tot_num_tracks / 2);
	int max = (i + 1) * (_tot_num_tracks / 2);
	
	/* Loop over each thread within this azimuthal angle halfspace */
        #pragma omp parallel for private(curr_track, num_segments, \
	  curr_segment, track_flux)
	for (int track_id=min; track_id < max; track_id++) {

	    /* TODO: Allocate this up front */
	    int size = _num_FSRs * sizeof(FP_PRECISION);
	    FP_PRECISION* fsr_flux = (FP_PRECISION*)_mm_malloc(size,
							       VEC_ALIGNMENT);
	    /* Initialize local pointers to important data structures */	
	    curr_track = _tracks[track_id];
	    num_segments = curr_track->getNumSegments();
	    track_flux = &_boundary_flux(track_id,0);

	    /* Loop over each segment in forward direction */
	    for (int s=0; s < num_segments; s++) {
	        curr_segment = curr_track->getSegment(s);
		scalarFluxTally(curr_segment, track_flux, fsr_flux);
	    }

	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, true, track_flux);
	    
	    /* Loop over each segment in reverse direction */
	    track_flux += _polar_times_groups;
	    
	    for (int s=num_segments-1; s > -1; s--) {
	        curr_segment = curr_track->getSegment(s);
		scalarFluxTally(curr_segment, track_flux, fsr_flux);
	    }
	    
	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, false, track_flux);
	}
    }

    return;
}


void CPUSolver::scalarFluxTally(segment* curr_segment,
                                FP_PRECISION* track_flux,
                                FP_PRECISION* fsr_flux){

    FP_PRECISION delta;
    FP_PRECISION sigma_t_l;
    int index;

    int fsr_id = curr_segment->_region_id;
    FP_PRECISION length = curr_segment->_length;
    double* sigma_t = curr_segment->_material->getSigmaT();

    /* Loop over energy groups */
    for (int e=0; e < _num_groups; e++) {
        fsr_flux[e] = 0.;
	sigma_t_l = sigma_t[e] * length;
	index = prefactorindex(sigma_t_l);

	/* Loop over polar angles */
	#pragma novector
	for (int p=0; p < _num_polar; p++){
	    delta = (track_flux(p,e) - 
	    _ratios(fsr_id,e)) * 
	      prefactor(index,p,sigma_t_l);
	    fsr_flux[e] += delta * _polar_weights[p];
	    track_flux(p,e) -= delta;
	}
    }

    /* Atomically increment the FSR scalar flux from the temporary array */
    omp_set_lock(&_FSR_locks[fsr_id]);
    {
        for (int e=0; e < _num_groups; e++) {
	    _scalar_flux(fsr_id,e) += fsr_flux[e];
	}
    }
    omp_unset_lock(&_FSR_locks[fsr_id]);

    return;
}


void CPUSolver::transferBoundaryFlux(int track_id, bool direction,
				     FP_PRECISION* track_flux) {
    int start;
    bool bc;
    FP_PRECISION* track_leakage;
    int track_out_id;

    if (direction) {
        start = _tracks[track_id]->isReflOut() * _polar_times_groups;
        bc = _tracks[track_id]->getBCOut();
        track_leakage = &_boundary_leakage(track_id,0);
        track_out_id = _tracks[track_id]->getTrackOut()->getUid();
    }
    else {
        start = _tracks[track_id]->isReflIn() * _polar_times_groups;
        bc = _tracks[track_id]->getBCIn();
        track_leakage = &_boundary_leakage(track_id,_polar_times_groups);
        track_out_id = _tracks[track_id]->getTrackIn()->getUid();
    }

    FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,start);

    #pragma novector
    for (int p=0; p < _num_polar; p++) {
        #pragma novector
        for (int e=0; e < _num_groups; e++) {
	    track_out_flux(p,e) = track_flux(p,e) * bc;
	    track_leakage(p,e) = track_flux(p,e) * 
	      _polar_weights[p] * (!bc);
	}
    }
}


void CPUSolver::addSourceToScalarFlux() {

    FP_PRECISION volume;
    double* sigma_t;

    /* Add in source term and normalize flux to volume for each region */
    /* Loop over flat source regions, energy groups */
    #pragma omp parallel for private(volume, sigma_t)
    for (int r=0; r < _num_FSRs; r++) {
        volume = _FSR_volumes[r];
	sigma_t = _FSR_materials[r]->getSigmaT();

	for (int e=0; e < _num_groups; e++) {
            _scalar_flux(r,e) *= 0.5;
	    _scalar_flux(r,e) = FOUR_PI * _ratios(r,e) + 
	      (_scalar_flux(r,e) / (sigma_t[e] * volume));
        }
    }
    
    return;
}

 


/**
 * @brief Compute the fission rates in each flatsourceregion and stores them 
 *        in an array indexed by flatsourceregion ID.
 */
void CPUSolver::computePinPowers() {

    log_printf(INFO, "Computing FSR pin powers...");

    double* sigma_f;
    FP_PRECISION tot_pin_power = 0.;
    FP_PRECISION avg_pin_power = 0.;
    FP_PRECISION num_nonzero_pins = 0.;
    FP_PRECISION curr_pin_power = 0.;
    FP_PRECISION prev_pin_power = 0.;

    /* Loop over all FSRs and compute the fission rate*/
    #pragma omp parallel for private (sigma_f)
    for (int r=0; r < _num_FSRs; r++) {
        sigma_f = _FSR_materials[r]->getSigmaF();

        for (int e=0; e < _num_groups; e++)
	    _FSRs_to_powers[r] += sigma_f[e] * _scalar_flux(r,e);
    }

    /* Compute the pin powers by adding up the powers of FSRs in each
     * lattice cell, saving lattice cell powers to files, and saving the
     * pin power corresponding to each FSR id in FSR_to_pin_powers */
    _geometry->computePinPowers(_FSRs_to_powers, _FSRs_to_pin_powers);


    /* Compute the total power based by accumulating the power of each unique
     * pin with a nonzero power */
    for (int r=0; r < _num_FSRs; r++) {
        curr_pin_power = _FSRs_to_pin_powers[r];

	/* If this pin power is unique and nozero (doesn't match the previous
	 * pin's power), then tally it */
	if (curr_pin_power > 0. && curr_pin_power != prev_pin_power) {
	    tot_pin_power += curr_pin_power;
	    num_nonzero_pins++;
	    prev_pin_power = curr_pin_power;
	}
    }

    /* Compute the average pin power */
    avg_pin_power = tot_pin_power / num_nonzero_pins;

    /* Normalize each pin power to the average non-zero pin power */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++)
        _FSRs_to_pin_powers[r] /= avg_pin_power;

    return;
}
