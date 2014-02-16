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
CPUSolver::CPUSolver(Geometry* geometry, TrackGenerator* track_generator, 
  Cmfd* cmfd) : Solver(geometry, track_generator, cmfd) {

    setNumThreads(1);

    _FSR_locks = NULL;
    _mesh_surface_locks = NULL;
    _thread_fsr_flux = NULL;
}


/**
 * @brief Destructor deletes array for OpenMP atomic locks for scalar flux
 *        updates, and calls Solver parent class destructor to deletes arrays
 *        for fluxes and sources.
 */
CPUSolver::~CPUSolver() {

    if (_FSR_locks != NULL)
        delete [] _FSR_locks;

    if (_mesh_surface_locks != NULL)
        delete [] _mesh_surface_locks;

    if (_thread_fsr_flux != NULL)
        delete [] _thread_fsr_flux;

    if (_surface_currents != NULL)
        delete [] _surface_currents;
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
 * @return the flat source region scalar flux
 */
FP_PRECISION CPUSolver::getFSRScalarFlux(int fsr_id, int energy_group) {

    /* Error checking */
    if (fsr_id >= _num_FSRs)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
                "in enery group %d since the solver only contains FSR with "
                "IDs less than or equal to %d", 
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

    return _scalar_flux(fsr_id,energy_group-1);
}


/**
 * @brief Returns the source for some energy group for a flat source region
 * @param fsr_id the ID for the FSR of interest
 * @param energy_group the energy group of interest
 * @return the flat source region source
 */
FP_PRECISION CPUSolver::getFSRSource(int fsr_id, int energy_group) {

    /* Error checking */
    if (fsr_id >= _num_FSRs)
        log_printf(ERROR, "Unable to return a source for FSR id = %d "
                "in enery group %d since the solver only contains FSR with "
                "IDs less than or equal to %d", 
                fsr_id, energy_group, _num_FSRs-1);

    if (fsr_id < 0)
        log_printf(ERROR, "Unable to return a source for FSR id = %d "
                "in energy group %d since FSRs do not have negative IDs", 
                fsr_id, energy_group);

    if (energy_group-1 >= _num_groups)
        log_printf(ERROR, "Unable to return a source for FSR id = %d "
                "in energy group %d since the solver only has %d energy "
                "groups", fsr_id, energy_group, _num_groups);

    if (energy_group <= 0)
        log_printf(ERROR, "Unable to return a source for FSR id = %d "
                "in energy group %d since energy groups are greater than 1",
                fsr_id, energy_group);

    return _source(fsr_id,energy_group-1);
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
 * @brief Return an array indexed by surface IDs and energy groups
 *        which contains the corresponding surface currents
 *        for each mesh cell surface.
 * @return an array of mesh cell surface currents
 */
double* CPUSolver::getSurfaceCurrents() {

    if (_surface_currents == NULL)
        log_printf(ERROR, "Unable to returns the Solver's surface currents "
                "array since it has not yet been allocated in memory");

    return _surface_currents;
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
        size = _num_groups * _num_threads;
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
    if (_fission_sources != NULL)
        delete [] _fission_sources;

    if (_scatter_sources != NULL)
        delete [] _scatter_sources;

    if (_source != NULL)
        delete [] _source;

    if (_old_source != NULL)
        delete [] _old_source;

    if (_reduced_source != NULL)
        delete [] _reduced_source;

    if (_source_residuals != NULL)
        delete [] _source_residuals;

    int size;

    /* Allocate memory for all source arrays */
    try{
        size = _num_FSRs * _num_groups;
        _fission_sources = new FP_PRECISION[size];
        _source = new FP_PRECISION[size];
        _old_source = new FP_PRECISION[size];
        _reduced_source = new FP_PRECISION[size];

        size = _num_threads * _num_groups;
        _scatter_sources = new FP_PRECISION[size];

        size = _num_FSRs;
        _source_residuals = new FP_PRECISION[size];

    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
                "source region sources array. Backtrace:%s", e.what());
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

    /* Precompute the total azimuthal weight for tracks at each polar angle */
    #pragma omp parallel for private(azim_weight) schedule(guided)
    for (int i=0; i < _num_azim; i++) {

        azim_weight = _azim_weights[i];

        for (int p=0; p < _num_polar; p++)
            _polar_weights(i,p) = azim_weight*_quad->getMultiple(p)*FOUR_PI;
    }

    /* Set size of prefactor array */
    int num_array_values = 10 * sqrt(1. / (8. * _source_convergence_thresh * 1e-2));
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
 * @brief Initializes the volumes and material arrays for each flat source 
 *        region. 
 * @details This method assigns each flat source region a unique, monotonically
 *          increasing ID, sets the material for each flat source region, and 
 *          assigns a volume based on the cumulative length of all of the 
 *          segments inside the flat source region.
 */
void CPUSolver::initializeFSRs() {

    log_printf(INFO, "Initializing flat source regions...");

    /* Delete old FSR arrays if they exist */
    if (_FSR_volumes != NULL)
        delete [] _FSR_volumes;

    if (_FSR_materials != NULL)
        delete [] _FSR_materials;

    _FSR_volumes = (FP_PRECISION*)calloc(_num_FSRs, sizeof(FP_PRECISION));
    _FSR_materials = new Material*[_num_FSRs];
    _FSR_locks = new omp_lock_t[_num_FSRs];

    int num_segments;
    segment* curr_segment;
    segment* segments;
    FP_PRECISION volume;
    CellBasic* cell;
    Material* material;
    Universe* univ_zero = _geometry->getUniverse(0);

    /* Set each FSR's "volume" by accumulating the total length of all tracks
     * inside the FSR. Loop over azimuthal angle, track and segment. 
     * Note: this code region cannot be parallelized without a mutex lock
     * on FSR volume due to race conditions. */
    for (int i=0; i < _tot_num_tracks; i++) {
        
        int azim_index = _tracks[i]->getAzimAngleIndex();
        num_segments = _tracks[i]->getNumSegments();
        segments = _tracks[i]->getSegments();

        for (int s=0; s < num_segments; s++) {
            curr_segment = &segments[s];
            volume = curr_segment->_length * _azim_weights[azim_index];
            _FSR_volumes[curr_segment->_region_id] += volume;
        }
    }

    /* Loop over all FSRs to extract FSR material pointers */
    #pragma omp parallel for private(cell, material) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get the cell corresponding to this FSR from the geometry */
        cell = _geometry->findCellContainingFSR(r);

        /* Get the cell's material and assign it to the FSR */
        material = _geometry->getMaterial(cell->getMaterial());
        _FSR_materials[r] = material;

        log_printf(DEBUG, "FSR id = %d has cell id = %d and material id = %d "
                "and volume = %f", r, cell->getId(), 
                _FSR_materials[r]->getUid(), _FSR_volumes[r]);
    }

    /* Loop over all FSRs to initialize OpenMP locks */
    #pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_FSRs; r++)
        omp_init_lock(&_FSR_locks[r]);
    
    return;
}


/** 
 * @brief Initializes Cmfd prior to source iteration.
 * @details Instantiates a dummy Cmfd object if one was not assigned to
 *          the solver by the user and initializes FSRs, materials, fluxes
 *          and the mesh. This method intializes a global array for the
 *          surface currents.
 */
void CPUSolver::initializeCmfd() {

  /* Call parent class method */
  Solver::initializeCmfd();
  
  /* Delete old surface currents array it it exists */
  if (_surface_currents != NULL)
    delete [] _surface_currents;
  
  int size;

  /* Allocate memory for the surface currents array */
  try{

    /* Allocate an array for the surface currents */
    if (_cmfd->getMesh()->getCmfdOn()){ 
      size = _num_mesh_cells * _cmfd->getNumCmfdGroups() * 8;
      _surface_currents = new double[size];
    }

  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the solver's fluxes. "
	       "Backtrace:%s", e.what());
  }

  if (_cmfd->getMesh()->getCmfdOn()){ 

    /* Initialize an array of OpenMP locks for each Cmfd mesh surface */
    _mesh_surface_locks = new omp_lock_t[_cmfd->getMesh()->getNumCells() * 8];

    /* Loop over all mesh cells to initialize OpenMP locks */
    #pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_mesh_cells*8; r++)
        omp_init_lock(&_mesh_surface_locks[r]);
    }

  return;
}


/**
 * @brief Zero each track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
void CPUSolver::zeroTrackFluxes() {

    #pragma omp parallel for schedule(guided)
    for (int t=0; t < _tot_num_tracks; t++) {
        for (int d=0; d < 2; d++) {
            for (int p=0; p < _num_polar; p++) {
                for (int e=0; e < _num_groups; e++) {
                    _boundary_flux(t,d,p,e) = 0.0;
                }  
            }
        }
    }

  return;
}


/**
 * @brief Set the scalar flux for each energy group inside each flat source 
 *        region to a constant value.
 * @param value the value to assign to each flat source region flux
 */
void CPUSolver::flattenFSRFluxes(FP_PRECISION value) {

    #pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++)
            _scalar_flux(r,e) = value;
    }

    return;
}


 /**
  * @brief Set the surface currents for each energy group inside each
  *                         mesh cell to zero
  */
void CPUSolver::zeroSurfaceCurrents() {

    #pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_mesh_cells; r++) {
        for (int s=0; s < 8; s++) {
            for (int e=0; e < _num_groups; e++)
                _surface_currents(r*8+s,e) = 0.0;
        }
    }

    return;
}


/**
 * @brief Set the source for each energy group inside each flat source region
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
void CPUSolver::flattenFSRSources(FP_PRECISION value) {

    #pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++) {
            _source(r,e) = value;
            _old_source(r,e) = value;
        }
    }

    return;
}


/**
 * @brief Normalizes all flat source region scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times \f$ \nu \f$).
 */
void CPUSolver::normalizeFluxes() {

    FP_PRECISION* nu_sigma_f;
    FP_PRECISION volume;
    FP_PRECISION tot_fission_source;
    FP_PRECISION norm_factor;

    /* Compute total fission source for each region, energy group */
    #pragma omp parallel for private(volume, nu_sigma_f)        \
        reduction(+:tot_fission_source) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get pointers to important data structures */
        nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
        volume = _FSR_volumes[r];

        for (int e=0; e < _num_groups; e++)
            _fission_sources(r,e) = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
    }

    /* Compute the total fission source */
    tot_fission_source = pairwise_sum<FP_PRECISION>(_fission_sources, 
                                                    _num_FSRs*_num_groups);

    /* Normalize scalar fluxes in each region */
    norm_factor = 1.0 / tot_fission_source;

    log_printf(DEBUG, "tot fiss src = %f, Normalization factor = %f", 
            tot_fission_source, norm_factor);

    #pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++)
            _scalar_flux(r,e) *= norm_factor;
    }

    /* Normalize angular boundary fluxes for each track */
    #pragma omp parallel for schedule(guided)
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
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{\# FSRs}}} \f$
 *
 * @return the residual between this source and the previous source
 */
FP_PRECISION CPUSolver::computeFSRSources() {

    int tid;
    Material* material;
    FP_PRECISION scatter_source;
    FP_PRECISION fission_source;
    FP_PRECISION* nu_sigma_f;
    FP_PRECISION* sigma_s;
    FP_PRECISION* sigma_t;
    FP_PRECISION* chi;

    FP_PRECISION source_residual = 0.0;
    
    FP_PRECISION inverse_k_eff = 1.0 / _k_eff;

    /* For all regions, find the source */
    #pragma omp parallel for private(tid, material, nu_sigma_f, chi,        \
     sigma_s, sigma_t, fission_source, scatter_source) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        tid = omp_get_thread_num();
        material = _FSR_materials[r];
        nu_sigma_f = material->getNuSigmaF();
        chi = material->getChi();
        sigma_s = material->getSigmaS();
        sigma_t = material->getSigmaT();

        /* Initialize the source residual to zero */
        _source_residuals[r] = 0.;

        /* Compute fission source for each group */
        if (material->isFissionable()) {
            for (int e=0; e < _num_groups; e++)
                _fission_sources(r,e) = _scalar_flux(r,e) * nu_sigma_f[e];

            fission_source = pairwise_sum<FP_PRECISION>(&_fission_sources(r,0),
                                                        _num_groups);
            fission_source *= inverse_k_eff;
        }

        else
            fission_source = 0.0;         

        /* Compute total scattering source for group G */
        for (int G=0; G < _num_groups; G++) {
            scatter_source = 0;

            for (int g=0; g < _num_groups; g++)
                _scatter_sources(tid,g) = sigma_s[G*_num_groups+g] * 
                                           _scalar_flux(r,g);

            scatter_source=pairwise_sum<FP_PRECISION>(&_scatter_sources(tid,0),
                                                       _num_groups);

            /* Set the total source for region r in group G */
            _source(r,G) = (fission_source * chi[G] + scatter_source) * 
                           ONE_OVER_FOUR_PI;

            _reduced_source(r,G) = _source(r,G) / sigma_t[G];

            /* Compute the norm of residual of the source in the region */
            if (fabs(_source(r,G)) > 1E-10)
               _source_residuals[r] += pow((_source(r,G) - _old_source(r,G)) 
                                           / _source(r,G), 2);

            /* Update the old source */
            _old_source(r,G) = _source(r,G);
        }
    }

    /* Sum up the residuals from each region */
    source_residual = pairwise_sum<FP_PRECISION>(_source_residuals, _num_FSRs);
    source_residual = sqrt(source_residual / (_num_FSRs * _num_groups));

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

    int tid;
    Material* material;
    FP_PRECISION* sigma_a;
    FP_PRECISION* nu_sigma_f;
    FP_PRECISION volume;

    FP_PRECISION tot_abs = 0.0;
    FP_PRECISION tot_fission = 0.0;

    FP_PRECISION* FSR_rates = new FP_PRECISION[_num_FSRs];
    FP_PRECISION* group_rates = new FP_PRECISION[_num_threads * _num_groups];

    /* Loop over all flat source regions and compute the volume-weighted
     * absorption rates */
    #pragma omp parallel for private(tid, volume, \
        material, sigma_a) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        tid = omp_get_thread_num() * _num_groups;
        volume = _FSR_volumes[r];
        material = _FSR_materials[r];
        sigma_a = material->getSigmaA();
        
        for (int e=0; e < _num_groups; e++)
            group_rates[tid+e] = sigma_a[e] * _scalar_flux(r,e);

        FSR_rates[r]=pairwise_sum<FP_PRECISION>(&group_rates[tid], 
                                                _num_groups);
        FSR_rates[r] *= volume;
    }

    /* Reduce absorption rates across FSRs */
    tot_abs = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);


    /* Loop over all flat source regions and compute the volume-weighted
     * fission rates */
    #pragma omp parallel for private(tid, volume, \
        material, nu_sigma_f) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        tid = omp_get_thread_num() * _num_groups;
        volume = _FSR_volumes[r];
        material = _FSR_materials[r];
        nu_sigma_f = material->getNuSigmaF();
        
        for (int e=0; e < _num_groups; e++)
            group_rates[tid+e] = nu_sigma_f[e] * _scalar_flux(r,e);

        FSR_rates[r]=pairwise_sum<FP_PRECISION>(&group_rates[tid], 
                                                _num_groups);
        FSR_rates[r] *= volume;
    }

    /* Reduce fission rates across FSRs */
    tot_fission = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);
    
    /** Reduce leakage array across tracks, energy groups, polar angles */
    int size = 2 * _tot_num_tracks * _polar_times_groups;
    _leakage = pairwise_sum<FP_PRECISION>(_boundary_leakage, size) * 0.5;

    _k_eff = tot_fission / (tot_abs + _leakage);
    
    log_printf(DEBUG, "abs = %f, fission = %f, leakage = %f, "
                "k_eff = %f", tot_abs, tot_fission, _leakage, _k_eff);

    delete [] FSR_rates;
    delete [] group_rates;

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
    int min_track, max_track;
    Track* curr_track;
    int azim_index;
    int num_segments;
    segment* curr_segment;
    segment* segments;
    FP_PRECISION* track_flux;

    log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

    /* Initialize flux in each region to zero */
    flattenFSRFluxes(0.0);

    if (_cmfd->getMesh()->getCmfdOn())
      zeroSurfaceCurrents();

    /* Loop over azimuthal angle halfspaces */
    for (int i=0; i < 2; i++) {

        /* Compute the minimum and maximum track IDs corresponding to 
         * this azimuthal angular halfspace */
        min_track = i * (_tot_num_tracks / 2);
        max_track = (i + 1) * (_tot_num_tracks / 2);

        /* Loop over each thread within this azimuthal angle halfspace */
        #pragma omp parallel for private(curr_track, azim_index, \
            num_segments, curr_segment, segments, track_flux, \
            tid) schedule(guided)
        for (int track_id=min_track; track_id < max_track; track_id++) {

            tid = omp_get_thread_num();

            /* Initialize local pointers to important data structures */        
            curr_track = _tracks[track_id];
            azim_index = curr_track->getAzimAngleIndex();
            num_segments = curr_track->getNumSegments();
            segments = curr_track->getSegments();
            track_flux = &_boundary_flux(track_id,0,0,0);

            /* Loop over each segment in forward direction */
            for (int s=0; s < num_segments; s++) {
                curr_segment = &segments[s];
                scalarFluxTally(curr_segment, azim_index, track_flux, 
                                &_thread_fsr_flux(tid),true);
            }

            /* Transfer flux to outgoing track */
            transferBoundaryFlux(track_id, azim_index, true, track_flux);

            /* Loop over each segment in reverse direction */
            track_flux += _polar_times_groups;

            for (int s=num_segments-1; s > -1; s--) {
                curr_segment = &segments[s];
                scalarFluxTally(curr_segment, azim_index, track_flux, 
                                &_thread_fsr_flux(tid),false);
            }

            /* Transfer flux to outgoing track */
            transferBoundaryFlux(track_id, azim_index, false, track_flux);
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
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux a pointer to the track's angular flux
 * @param fsr_flux a pointer to the temporary flat source region flux buffer
 */
void CPUSolver::scalarFluxTally(segment* curr_segment,
                                int azim_index,
                                FP_PRECISION* track_flux,
                                FP_PRECISION* fsr_flux,
                                bool fwd){

    int tid = omp_get_thread_num();
    int fsr_id = curr_segment->_region_id;
    FP_PRECISION length = curr_segment->_length;
    FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();

    /* The change in angular flux along this segment in the FSR */
    FP_PRECISION deltapsi;
    FP_PRECISION exponential;

    /* Set the flat source region flux buffer to zero */
    memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

    /* Loop over energy groups */
    for (int e=0; e < _num_groups; e++) {

        /* Loop over polar angles */
        for (int p=0; p < _num_polar; p++){
            exponential = computeExponential(sigma_t[e], length, p);
            deltapsi = (track_flux(p,e)-_reduced_source(fsr_id,e))*exponential;
            fsr_flux[e] += deltapsi * _polar_weights[p];
            track_flux(p,e) -= deltapsi;
        }
    }

    if (_cmfd->getMesh()->getCmfdOn()){
        if (curr_segment->_mesh_surface_fwd != -1 && fwd){

                /* set polar angle * energy group to 0 */
                int pe = 0;

                /* Atomically increment the meshSurface current from the temporary array */
                omp_set_lock(&_mesh_surface_locks[curr_segment->_mesh_surface_fwd]);

                /* loop over energy groups */
                for (int e = 0; e < _num_groups; e++) {

                        /* loop over polar angles */
                        for (int p = 0; p < _num_polar; p++){

                                /* increment current (polar and azimuthal weighted flux, group)*/
                                _surface_currents(curr_segment->_mesh_surface_fwd,e) += track_flux(p,e)*_polar_weights[p]/2.0;

                                pe++;
                        }
                }

                /* Release mesh surface lock */
            omp_unset_lock(&_mesh_surface_locks[curr_segment->_mesh_surface_fwd]);

        }
        else if (curr_segment->_mesh_surface_bwd != -1 && !fwd){

                /* set polar angle * energy group to 0 */
                int pe = 0;

                /* Atomically increment the meshSurface current from the temporary array */
                omp_set_lock(&_mesh_surface_locks[curr_segment->_mesh_surface_bwd]);

                /* loop over energy groups */
                for (int e = 0; e < _num_groups; e++) {

                        /* loop over polar angles */
                        for (int p = 0; p < _num_polar; p++){

                                /* increment current (polar and azimuthal weighted flux, group)*/
                                _surface_currents(curr_segment->_mesh_surface_bwd,e) += track_flux(p,e)*_polar_weights[p]/2.0;

                                pe++;
                        }
                }

                /* Release mesh surface lock */
                omp_unset_lock(&_mesh_surface_locks[curr_segment->_mesh_surface_bwd]);
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


/**
 * @brief Computes the exponential term in the transport equation for a
 *        track segment.
 * @details This method computes \f$ 1 - exp(-l\Sigma^T_g/sin(\theta_p)) \f$ 
 *          for a segment with total group cross-section and for
 *          some polar angle.
 * @param sigma_t the total group cross-section at this energy
 * @param length the length of the line segment projected in the xy-plane
 * @param p the polar angle index
 * @return the evaluated exponential
 */
FP_PRECISION CPUSolver::computeExponential(FP_PRECISION sigma_t,
                                           FP_PRECISION length, int p) {

    FP_PRECISION exponential;
    FP_PRECISION tau = sigma_t * length;

    /* Evaluate the exponential using the lookup table - linear interpolation */
    if (_interpolate_exponential) {
        int index;

        index = round_to_int(tau * _inverse_prefactor_spacing);
        index *= _two_times_num_polar;
        exponential = (1. - (_prefactor_array[index+2 * p] * tau + 
                          _prefactor_array[index + 2 * p +1]));
    }

    /* Evalute the exponential using the intrinsic exp function */
    else {
        FP_PRECISION sintheta = _quad->getSinTheta(p);
        exponential = 1.0 - exp(- tau / sintheta);
    }

    return exponential;
}


/**
 * @brief Updates the boundary flux for a track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the track is given to the reflecting track. For vacuum
 *          boundary conditions, the outgoing flux tallied as leakage.
 * @param track_id the ID number for the track of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param direction the track direction (forward - true, reverse - false)
 * @param track_flux a pointer to the track's outgoing angular flux
 */
void CPUSolver::transferBoundaryFlux(int track_id, int azim_index,
                                     bool direction, FP_PRECISION* track_flux) {
    int start;
    int bc;
    FP_PRECISION* track_leakage;
    int track_out_id;

    /* Extract boundary conditions for this track and the pointer to the 
     * outgoing reflective track, and index into the leakage array */

    /* For the "forward" direction */
    if (direction) {
        start = _tracks[track_id]->isReflOut() * _polar_times_groups;
        bc = (int)_tracks[track_id]->getBCOut();
        track_leakage = &_boundary_leakage(track_id,0);
        track_out_id = _tracks[track_id]->getTrackOut()->getUid();
    }

    /* For the "reverse" direction */
    else {
        start = _tracks[track_id]->isReflIn() * _polar_times_groups;
        bc = (int)_tracks[track_id]->getBCIn();
        track_leakage = &_boundary_leakage(track_id,_polar_times_groups);
        track_out_id = _tracks[track_id]->getTrackIn()->getUid();
    }

    FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,0,0,start);

    /* Loop over polar angles and energy groups */
    for (int e=0; e < _num_groups; e++) {
        for (int p=0; p < _num_polar; p++) {
            track_out_flux(p,e) = track_flux(p,e) * bc;
            track_leakage(p,e) = track_flux(p,e) * 
              _polar_weights(azim_index,p) * (!bc);
        }
    }
}


/**
 * @brief Add the source term contribution in the transport equation to 
 *        the flat source region scalar flux
 */
void CPUSolver::addSourceToScalarFlux() {

    FP_PRECISION volume;
    FP_PRECISION* sigma_t;

    /* Add in source term and normalize flux to volume for each region */
    /* Loop over flat source regions, energy groups */
    #pragma omp parallel for private(volume, sigma_t) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

        volume = _FSR_volumes[r];
        sigma_t = _FSR_materials[r]->getSigmaT();

        for (int e=0; e < _num_groups; e++) {
            _scalar_flux(r,e) *= 0.5;
            _scalar_flux(r,e) = FOUR_PI * _reduced_source(r,e) + 
              (_scalar_flux(r,e) / (sigma_t[e] * volume));
        }
    }
    
    return;
}


/**
 * @brief Computes the volume-weighted, energy integrated fission rate in 
 *        each flat source region and stores them in an array indexed by 
 *        flat source region ID.
 * @param fission_rates an array to store the fission rates, passed in as a
 *        numpy array from Python
 * @param num_FSRs the number of FSRs passed in from Python
 */
void CPUSolver::computeFSRFissionRates(double* fission_rates, int num_FSRs) {

    log_printf(INFO, "Computing FSR fission rates...");

    FP_PRECISION* sigma_f;

    FP_PRECISION* scalar_flux = getFSRScalarFluxes();

    /* Loop over all FSRs and compute the volume-weighted fission rate */
    #pragma omp parallel for private (sigma_f) schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {
        sigma_f = _FSR_materials[r]->getSigmaF();

        for (int e=0; e < _num_groups; e++)
            fission_rates[r] += sigma_f[e] * _scalar_flux(r,e);
    }

    return;
}

