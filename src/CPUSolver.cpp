#include "CPUSolver.h"


/**
 * @brief Constructor initializes array pointers for tracks and materials.
 * @details The constructor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator. The constructor initalizes the number of threads to a 
 *          default of 1.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the trackgenerator
 */
CPUSolver::CPUSolver(Geometry* geom, TrackGenerator* track_generator) :
    Solver(geom, track_generator) {

    setNumThreads(1);

    _FSR_locks = NULL;
    _thread_fsr_flux = NULL;

    _interpolate_exponent = true;
}


/**
 * @brief Destructor deletes array for OpenMP atomic locks for scalar flux
 *        updates, and calls Solver subclass destructor to deletes arrays
 *        for fluxes and sources.
 */
CPUSolver::~CPUSolver() { 

    if (_FSR_locks != NULL)
        delete [] _FSR_locks;

    if (_thread_fsr_flux != NULL)
        delete [] _thread_fsr_flux;

    if (_exponentials != NULL)
        delete [] _exponentials;
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int CPUSolver::getNumThreads() {
    return _num_threads;
}


/**
 * @brief Returns the scalar flux for some energy group for a flat source region
 * @param fsr_id the ID for the FSR of interest
 * @param energy_group the energy group of interest
 */
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
 * @brief Return an array indexed by flat source region IDs and energy groups 
 *        which contains the corresponding fluxes for each flat source region.
 * @return an array of flat source region scalar fluxes
 */
FP_PRECISION* CPUSolver::getFSRScalarFluxes() {

    if (_scalar_flux == NULL)
        log_printf(ERROR, "Unable to returns the Solver's scalar flux array "
		 "since it has not yet been allocated in memory");

    return _scalar_flux;
}


/**
 * @brief Return an array indexed by flat source region IDs with the
 *        corresponding flat source region power.
 * @return an array of flat source region powers
 */
FP_PRECISION* CPUSolver::getFSRPowers() {

    if (_FSRs_to_powers == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR power array "
		 "since it has not yet been allocated in memory");

    return _FSRs_to_powers;
}


/**
 * @brief Return an array indexed by flat source region IDs with the
 *        corresponding pin cell power.
 * @return an array of flat source region pin powers
 */
FP_PRECISION* CPUSolver::getFSRPinPowers() {

    if (_FSRs_to_pin_powers == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR pin power array "
		 "since it has not yet been allocated in memory");

    return _FSRs_to_pin_powers;
}



/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
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


/**
 * @brief Sets the solver to use linear interpolation to compute the exponential
 *        in the transport equation
 */
void CPUSolver::useExponentialInterpolation() {
    _interpolate_exponent = true;
}


/**
 * @brief Sets the solver to use the exponential intrinsic function to 
 *        compute the exponential in the transport equation
 */
void CPUSolver::useExponentialIntrinsic() {
    _interpolate_exponent = false;
}


/**
 * @brief Allocates memory for track boundary angular fluxes and leakages
 *        flat source region scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void CPUSolver::initializeFluxArrays() {
   
    /* Delete old flux arrays if they exist */
    if (_boundary_flux != NULL)
        delete [] _boundary_flux;

    if (_boundary_leakage != NULL)
        delete [] _boundary_leakage;

    if (_scalar_flux != NULL)
        delete [] _scalar_flux;

    if (_thread_fsr_flux != NULL)
        delete [] _thread_fsr_flux;

    int size;

    /* Allocate memory for the flux and leakage arrays */
    try{

        size = 2 * _tot_num_tracks * _polar_times_groups;
	_boundary_flux = new FP_PRECISION[size];
	_boundary_leakage = new FP_PRECISION[size];

	/* Allocate an array for the scalar flux */
	size = _num_FSRs * _num_groups;
	_scalar_flux = new FP_PRECISION[size];

	/* Allocate a thread local local memory buffer for FSR scalar flux */
	size = _num_groups * _num_threads * sizeof(FP_PRECISION);
	_thread_fsr_flux = new FP_PRECISION[size];
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
void CPUSolver::initializeSourceArrays() {

    /* Delete old sources arrays if they exist */
    if (_fission_source != NULL)
        delete [] _fission_source;

    if (_source != NULL)
        delete [] _source;

    if (_old_source != NULL)
        delete [] _old_source;

    if (_ratios != NULL)
        delete [] _ratios;

    int size;

    /* Allocate memory for all source arrays */
    try{
        size = _num_FSRs * _num_groups;
	_fission_source = new FP_PRECISION[size];
	_source = new FP_PRECISION[size];
	_old_source = new FP_PRECISION[size];
	_ratios = new FP_PRECISION[size];
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source region sources array. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flat source region power arrays.
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
 *          prefactor for specific segment lengths (the keys into the hashmap).
 */
void CPUSolver::precomputePrefactors() {

    /* Build exponential prefactor array based on table look up with linear 
     * interpolation */
    log_printf(INFO, "Building exponential prefactor hashtable...");

    FP_PRECISION azim_weight;

    _polar_weights = new FP_PRECISION[_num_azim*_num_polar];
    _exponentials = new FP_PRECISION[_num_threads * _polar_times_groups];

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
 * @brief Initializes each of the flat source region objects inside the solver's
 *        array of flatsourceregions. 
 * @details This method assigns each flat source region a unique, monotonically
 *          increasing ID, sets the material for each flat source region, and 
 *          assigns a volume based on the cumulative length of all of the 
 *          segments inside the flat source region.
 */
void CPUSolver::initializeFSRs() {

    log_printf(INFO, "Initializing flat source regions...");

    /* Delete old FSR arrayS if they exist */
    if (_FSR_volumes != NULL)
        delete [] _FSR_volumes;

    if (_FSR_materials != NULL)
        delete [] _FSR_materials;

    _FSR_volumes = new FP_PRECISION[_num_FSRs];
    _FSR_materials = new Material*[_num_FSRs];
    _FSR_locks = new omp_lock_t[_num_FSRs];

    int num_segments;
    segment* curr_segment;
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
	num_segments = _tracks[i]->getNumSegments();

	for (int s=0; s < num_segments; s++) {
            curr_segment = _tracks[i]->getSegment(s);
	    volume = curr_segment->_length * _azim_weights[azim_index];
	    _FSR_volumes[curr_segment->_region_id] += volume;
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
    int size = 2 * _tot_num_tracks * _polar_times_groups * sizeof(FP_PRECISION);
    memset(_boundary_flux, FP_PRECISION(0.), size);
    return;
}


/**
 * @brief Set the scalar flux for each energy group inside each flat source 
 *        region to a constant value.
 * @param value the value to assign to each flat source region flux
 */
void CPUSolver::flattenFSRFluxes(FP_PRECISION value) {
    int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
    memset(_scalar_flux, FP_PRECISION(value), size);
    return;
}


/**
 * @brief Set the source for each energy group inside each flat source region
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
void CPUSolver::flattenFSRSources(FP_PRECISION value) {
    int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
    memset(_source, value, size);
    memset(_old_source, value, size);
    return;
}


/**
 * @brief Normalizes all flat source region scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times $\nu$).
 */
void CPUSolver::normalizeFluxes() {

    double* nu_sigma_f;
    FP_PRECISION volume;
    FP_PRECISION tot_fission_source;
    FP_PRECISION norm_factor;

    memset(_fission_source, 0, _num_FSRs * _num_groups);

    /* Compute total fission source for each region, energy group */
    #pragma omp parallel for private(volume, nu_sigma_f)	\
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

    log_printf(DEBUG, "tot fiss src = %f, Normalization factor = %f", 
               tot_fission_source, norm_factor);

    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++)
	    _scalar_flux(r,e) *= norm_factor;
    }

    /* Normalize angular boundary fluxes for each track */
    #pragma omp parallel for
    for (int i=0; i < _tot_num_tracks; i++) {
        for (int j=0; j < 2; j++) {
	    for (int p=0; p < _num_polar; p++) {
	        for (int e=0; e < _num_groups; e++) {
		    _boundary_flux(i,j,p,e) *= norm_factor;
		}
	    }
	}
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
    #pragma omp parallel for private(material, nu_sigma_f, chi,	\
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
                scatter_sources[g] = sigma_s[G*_num_groups+g]*_scalar_flux(r,g);

	    scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources, 
                                                        _num_groups);

	    /* Set the total source for region r in group G */
	    _source(r,G) = ((1.0 / _k_eff) * fission_source *
                           chi[G] + scatter_source) * ONE_OVER_FOUR_PI;

	    _ratios(r,G) = _source(r,G) / sigma_t[G];

	    /* Compute the norm of residual of the source in the region, group */
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

    /* Loop over all flat source regions and compute the volume-weighted
     * fission and absorption rates */
    #pragma omp parallel for private(volume, material, sigma_a, nu_sigma_f)
    for (int r=0; r < _num_FSRs; r++) {

        volume = _FSR_volumes[r];
	material = _FSR_materials[r];
	sigma_a = material->getSigmaA();
	nu_sigma_f = material->getNuSigmaF();

	for (int e=0; e < _num_groups; e++) {
            absorption_rates[r*_num_groups+e] = sigma_a[e] * _scalar_flux(r,e);
	    fission_rates[r*_num_groups+e] = nu_sigma_f[e] * _scalar_flux(r,e);
	    absorption_rates[r*_num_groups+e] *= volume;
	    fission_rates[r*_num_groups+e] *= volume;
        }
    }

    /* Reduce absorptoin and fission rates across FSRs, energy groups */
    int size = _num_FSRs * _num_groups;
    tot_abs = pairwise_sum<FP_PRECISION>(absorption_rates, size);
    tot_fission = pairwise_sum<FP_PRECISION>(fission_rates, size);

    /** Reduce leakage array across tracks, energy groups, polar angles */
    size = 2 * _tot_num_tracks * _polar_times_groups;
    _leakage = pairwise_sum<FP_PRECISION>(_boundary_leakage, size) * 0.5;

    _k_eff = tot_fission / (tot_abs + _leakage);

    log_printf(DEBUG, "tot_abs = %f, tot_fission = %f, leakage = %f, "
	       "k_eff = %f", tot_abs, tot_fission, _leakage, _k_eff);

    delete [] absorption_rates;
    delete [] fission_rates;

    return;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles, 
 *        tracks, segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the 
 *          boundary fluxes for the corresponding output track, while updating 
 *          the scalar flux in each flat source region
 */
void CPUSolver::transportSweep() {

    int tid;
    Track* curr_track;
    int num_segments;
    segment* curr_segment;    
    FP_PRECISION* track_flux;

    /* Normalize angular boundary fluxes for each track */
    #pragma omp parallel for
    for (int i=0; i < _tot_num_tracks; i++) {
        for (int j=0; j < 2; j++) {
	    for (int p=0; p < _num_polar; p++) {
	        for (int e=0; e < _num_groups; e++) {
		    if (isnan(_boundary_flux(i,j,p,e)))
		        printf("_boundary_flux = %f\n", _boundary_flux(i,j,p,e));
		}
	    }
	}
    }



    log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

    /* Initialize flux in each region to zero */
    flattenFSRFluxes(0.0);

    /* Loop over azimuthal angle halfspaces */
    for (int i=0; i < 2; i++) {

        /* Compute the minimum and maximum track IDs corresponding to 
         * this azimuthal angular halfspace */
        int min = i * (_tot_num_tracks / 2);
	int max = (i + 1) * (_tot_num_tracks / 2);
	
	/* Loop over each thread within this azimuthal angle halfspace */
	#pragma omp parallel for private(curr_track, num_segments, \
	  curr_segment, track_flux, tid )
	for (int track_id=min; track_id < max; track_id++) {

	    tid = omp_get_thread_num();

	    /* Initialize local pointers to important data structures */	
	    curr_track = _tracks[track_id];
	    num_segments = curr_track->getNumSegments();
	    track_flux = &_boundary_flux(track_id,0,0,0);

	    /* Loop over each segment in forward direction */
	    for (int s=0; s < num_segments; s++) {
	        curr_segment = curr_track->getSegment(s);
		scalarFluxTally(curr_segment, track_flux, 
	                        &_thread_fsr_flux(tid));
	    }

	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, true, track_flux);
	    
	    /* Loop over each segment in reverse direction */
	    track_flux += _polar_times_groups;
	    
	    for (int s=num_segments-1; s > -1; s--) {
	        curr_segment = curr_track->getSegment(s);
		scalarFluxTally(curr_segment, track_flux, 
				&_thread_fsr_flux(tid));
	    }
	    
	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, false, track_flux);
	}
    }

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
void CPUSolver::scalarFluxTally(segment* curr_segment,
   	                        FP_PRECISION* track_flux,
	                        FP_PRECISION* fsr_flux){

    int tid = omp_get_thread_num();
    int fsr_id = curr_segment->_region_id;

    /* The average flux along this segment in the flat source region */
    FP_PRECISION psibar;

    FP_PRECISION* exponentials = &_exponentials[tid * _polar_times_groups];
    computeExponentials(curr_segment, exponentials);

    /* Set the flat source region flux buffer to zero */
    memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

    /* Loop over energy groups */
    for (int e=0; e < _num_groups; e++) {

	/* Loop over polar angles */
	for (int p=0; p < _num_polar; p++){
	    psibar = (track_flux(p,e) - _ratios(fsr_id,e)) * exponentials(p,e);
	    fsr_flux[e] += psibar * _polar_weights[p];
	    track_flux(p,e) -= psibar;
	}
    }

    /* Atomically increment the FSR scalar flux from the temporary array */
    omp_set_lock(&_FSR_locks[fsr_id]);
    {
        for (int e=0; e < _num_groups; e++)
	    _scalar_flux(fsr_id,e) += fsr_flux[e];
    }
    omp_unset_lock(&_FSR_locks[fsr_id]);

    return;
}


void CPUSolver::computeExponentials(segment* curr_segment, 
				    FP_PRECISION* exponentials) {

    FP_PRECISION length = curr_segment->_length;
    double* sigma_t = curr_segment->_material->getSigmaT();

    if (_interpolate_exponent) {
        FP_PRECISION tau;
        int index;

        for (int e=0; e < _num_groups; e++) {

            tau = sigma_t[e] * length;
	    index = prefactorindex(tau);

	    for (int p=0; p < _num_polar; p++)
	        exponentials(p,e) = prefactor(index,p,tau);
        }
    }
    else {

        FP_PRECISION* sinthetas = _quad->getSinThetas();

	for (int e=0; e < _num_groups; e++) {

            for (int p=0; p < _num_polar; p++)
	        exponentials(p,e) = 1.0 - exp(-sigma_t[e] * length / sinthetas[p]);
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
void CPUSolver::transferBoundaryFlux(int track_id, bool direction,
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
        bc = _tracks[track_id]->getBCOut();
        track_leakage = &_boundary_leakage(track_id,0);
        track_out_id = _tracks[track_id]->getTrackOut()->getUid();
    }

    /* For the "reverse" direction */
    else {
        start = _tracks[track_id]->isReflIn() * _polar_times_groups;
        bc = _tracks[track_id]->getBCIn();
        track_leakage = &_boundary_leakage(track_id,_polar_times_groups);
        track_out_id = _tracks[track_id]->getTrackIn()->getUid();
    }

    FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,0,0,start);

    /* Loop over polar angles and energy groups */
    for (int p=0; p < _num_polar; p++) {

        for (int e=0; e < _num_groups; e++) {
	    track_out_flux(p,e) = track_flux(p,e) * bc;
	    track_leakage(p,e) = track_flux(p,e) * 
	      _polar_weights[p] * (!bc);
	}
    }
}


/**
 * @brief Add the source term contribution in the transport equation to 
 *        the flat source region scalar flux
 */
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
 * @brief Compute the fission rates in each flat source region and stores them 
 *        in an array indexed by flat source region ID.
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
