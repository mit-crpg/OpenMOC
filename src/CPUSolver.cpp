#include "CPUSolver.h"
#include <unordered_map>
#include <numeric>

/**
 * @brief Constructor initializes array pointers for Tracks and Materials.
 * @details The constructor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of OpenMP threads to a default of 1.
 * @param track_generator an optional pointer to the TrackGenerator
 */
CPUSolver::CPUSolver(TrackGenerator* track_generator)
  : Solver(track_generator) {

  setNumThreads(1);
  _FSR_locks = NULL;
  _source_type = "Flat";
#ifdef MPIx
  _track_message_size = 0;
  _MPI_requests = NULL;
  _MPI_sends = NULL;
  _MPI_receives = NULL;
  _neighbor_connections.clear();
#endif
}


/**
 * @brief Destructor deletes array for OpenMP mutual exclusion locks for
 *        FSR scalar flux updates, and calls Solver parent class destructor
 *        to deletes arrays for fluxes and sources.
 */
CPUSolver::~CPUSolver() {
#ifdef MPIx
  deleteMPIBuffers();
#endif
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int CPUSolver::getNumThreads() {
  return _num_threads;
}


/**
 * @brief Fills an array with the scalar fluxes.
 * @details This class method is a helper routine called by the OpenMOC
 *          Python "openmoc.krylov" module for Krylov subspace methods.
 *          Although this method appears to require two arguments, in
 *          reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_fluxes = num_groups * num_FSRs
 *          fluxes = solver.getFluxes(num_fluxes)
 * @endcode
 *
 * @param out_fluxes an array of FSR scalar fluxes in each energy group
 * @param num_fluxes the total number of FSR flux values
 */
void CPUSolver::getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) {

  if (num_fluxes != _NUM_GROUPS * _geometry->getNumTotalFSRs())
    log_printf(ERROR, "Unable to get FSR scalar fluxes since there are "
               "%d groups and %d FSRs which does not match the requested "
               "%d flux values", _NUM_GROUPS, _geometry->getNumTotalFSRs(),
               num_fluxes);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to get FSR scalar fluxes since they "
               "have not yet been allocated");

  /* Copy the fluxes into the input array */
  else {
#pragma omp parallel for schedule(static)
    for (long r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _NUM_GROUPS; e++)
        out_fluxes[r*_NUM_GROUPS+e] = _scalar_flux(r,e);
    }
  }
  /* Reduce domain data for domain decomposition */
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {

    /* Allocate buffer for communication */
    long num_total_FSRs = _geometry->getNumTotalFSRs();
    FP_PRECISION* temp_fluxes = new FP_PRECISION[num_total_FSRs*_NUM_GROUPS];

    int rank = 0;
    MPI_Comm comm = _geometry->getMPICart();
    MPI_Comm_rank(comm, &rank);
    for (long r=0; r < num_total_FSRs; r++) {

      /* Determine the domain and local FSR ID */
      long fsr_id = r;
      int domain = 0;
      _geometry->getLocalFSRId(r, fsr_id, domain);

      /* Set data if in the correct domain */
      if (domain == rank)
        for (int e=0; e < _NUM_GROUPS; e++)
          temp_fluxes[r*_NUM_GROUPS+e] = out_fluxes[fsr_id*_NUM_GROUPS+e];
      else
        for (int e=0; e < _NUM_GROUPS; e++)
          temp_fluxes[r*_NUM_GROUPS+e] = 0.0;
    }

    /* Determine the type of FP_PRECISION and communicate fluxes */
    MPI_Datatype flux_type;
    if (sizeof(FP_PRECISION) == 4)
      flux_type = MPI_FLOAT;
    else
      flux_type = MPI_DOUBLE;

    MPI_Allreduce(temp_fluxes, out_fluxes, num_total_FSRs*_NUM_GROUPS,
                  flux_type, MPI_SUM, comm);
    delete [] temp_fluxes;
  }
#endif
}


/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
 * @param num_threads the number of threads
 */
void CPUSolver::setNumThreads(int num_threads) {
  if (num_threads <= 0)
    log_printf(ERROR, "Unable to set the number of threads to %d "
               "since it is less than or equal to 0", num_threads);

#ifdef MPIx
  /* Check the MPI library has enough thread support */
  int provided;
  MPI_Query_thread(&provided);
  if (num_threads > 1 && provided < MPI_THREAD_SERIALIZED)
    log_printf(WARNING, "Not enough thread support level in the MPI library, "
               "re-compile with another library. Thread support level should"
               "be at least MPI_THREAD_SERIALIZED.");
#endif

  if (_track_generator != NULL && num_threads > 1)
    if ((_track_generator->getSegmentFormation() == OTF_STACKS ||
       _track_generator->getSegmentFormation() == OTF_TRACKS) &&
       _track_generator->getNumThreads() != num_threads)
      log_printf(WARNING, "The number of threads used in track generation (%d) "
                 "should match the number of threads used in the solver (%d) "
                 "for OTF ray-tracing methods, as threaded buffers are shared.",
                 _track_generator->getNumThreads(), num_threads);

  /* Set the number of threads for OpenMP */
  _num_threads = num_threads;
  omp_set_num_threads(_num_threads);
}


/**
 * @brief Set the flux array for use in transport sweep source calculations.
 * @details This is a helper method for the checkpoint restart capabilities,
 *          as well as the IRAMSolver in the openmoc.krylov submodule. This
 *          routine may be used as follows from within Python:
 *
 * @code
 *          fluxes = numpy.random.rand(num_FSRs * num_groups, dtype=np.float)
 *          solver.setFluxes(fluxes)
 * @endcode
 *
 *          NOTE: This routine stores a pointer to the fluxes for the Solver
 *          to use during transport sweeps and other calculations. Hence, the
 *          flux array pointer is shared between NumPy and the Solver.
 *
 * @param in_fluxes an array with the fluxes to use
 * @param num_fluxes the number of flux values (# groups x # FSRs)
 */
void CPUSolver::setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) {
  if (num_fluxes != _NUM_GROUPS * _num_FSRs)
    log_printf(ERROR, "Unable to set an array with %d flux values for %d "
               " groups and %d FSRs", num_fluxes, _NUM_GROUPS, _num_FSRs);

  /* Allocate array if flux arrays have not yet been initialized */
  if (_scalar_flux == NULL)
    initializeFluxArrays();

  /* Set the scalar flux array pointer to the array passed in from NumPy */
  _scalar_flux = in_fluxes;
  _user_fluxes = true;
}


/**
 * @brief Assign a fixed source for a flat source region and energy group.
 * @details Fixed sources should be scaled to reflect the fact that OpenMOC
 *          normalizes the scalar flux such that the total energy- and
 *          volume-integrated production rate sums to 1.0.
 * @param fsr_id the flat source region ID
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void CPUSolver::setFixedSourceByFSR(long fsr_id, int group,
                                    FP_PRECISION source) {

  Solver::setFixedSourceByFSR(fsr_id, group, source);

  /* Allocate the fixed sources array if not yet allocated */
  if (_fixed_sources == NULL) {
    long size = _num_FSRs * _NUM_GROUPS;
    _fixed_sources = new FP_PRECISION[size]();
  }

  /* Warn the user if a fixed source has already been assigned to this FSR */
  if (fabs(_fixed_sources(fsr_id,group-1)) > FLT_EPSILON)
    log_printf(WARNING, "Overriding fixed source %f in FSR ID=%d group %g, with"
               " %f", _fixed_sources(fsr_id,group-1), fsr_id, group, source);

  /* Store the fixed source for this FSR and energy group */
  _fixed_sources(fsr_id,group-1) = source;
}


/**
 * @brief Reset all fixed sources and fixed sources moments to 0.
 */
void CPUSolver::resetFixedSources() {

  /* Reset fixed source by FSR map */
  std::map< std::pair<int, int>, FP_PRECISION >::iterator fsr_iter;
  for (fsr_iter = _fix_src_FSR_map.begin();
       fsr_iter != _fix_src_FSR_map.end(); ++fsr_iter)
    fsr_iter->second = 0;

  /* Reset fixed source by cell map */
  std::map< std::pair<Cell*, int>, FP_PRECISION >::iterator cell_iter;
  for (cell_iter = _fix_src_cell_map.begin();
       cell_iter != _fix_src_cell_map.end(); ++cell_iter)
    cell_iter->second = 0;

  /* Reset fixed source by material map */
  std::map< std::pair<Material*, int>, FP_PRECISION >::iterator mat_iter;
  for (mat_iter = _fix_src_material_map.begin();
       mat_iter != _fix_src_material_map.end(); ++mat_iter)
    mat_iter->second = 0;

  /* Reset array of fixed sources */
  if (_fixed_sources != NULL)
    memset(_fixed_sources, 0, _num_FSRs * _NUM_GROUPS * sizeof(FP_PRECISION));
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 * @details This method allocates and initializes an array of OpenMP
 *          mutual exclusion locks for each FSR for use in the
 *          transport sweep algorithm.
 */
void CPUSolver::initializeFSRs() {

#ifdef LINEARSOURCE
  if (strcmp(_source_type.c_str(), "Flat") == 0)
    log_printf(ERROR, "OpenMOC was compiled for linear sources only. Remove "
               "-DLINEARSOURCE optimization flag to use a flat source solver");
#endif

  Solver::initializeFSRs();

  /* Get FSR locks from TrackGenerator */
  _FSR_locks = _track_generator->getFSRLocks();
}


/**
 * @brief Allocates memory for Track boundary angular flux and leakage
 *        and FSR scalar flux arrays.
 * @details Deletes memory for old flux arrays if they were allocated
 *          for a previous simulation.
 */
void CPUSolver::initializeFluxArrays() {

  /* Delete old flux arrays if they exist */
  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_start_flux != NULL)
    delete [] _start_flux;

  if (_boundary_leakage != NULL)
    delete [] _boundary_leakage;

  if (_scalar_flux != NULL)
    delete [] _scalar_flux;

  if (_old_scalar_flux != NULL)
    delete [] _old_scalar_flux;

  if (_stabilizing_flux != NULL)
    delete [] _stabilizing_flux;

#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    deleteMPIBuffers();
#endif

  long size;

  /* Allocate memory for the Track boundary fluxes and leakage arrays */
  try {
    size = 2 * _tot_num_tracks * _fluxes_per_track;
    long max_size = size;
#ifdef MPIX
    if (_geometry->isDomainDecomposed())
      MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX,
                    _geometry->getMPICart());
#endif
    double max_size_mb = (double) (2 * max_size * sizeof(float))
        / (double) (1e6);
#ifdef ONLYVACUUMBC
    max_size_mb /= 2;
#endif
    log_printf(NORMAL, "Max boundary angular flux storage per domain = %6.2f "
               "MB", max_size_mb);

    _boundary_flux = new float[size]();
#ifndef ONLYVACUUMBC
    _start_flux = new float[size]();
#endif

    /* Allocate memory for boundary leakage if necessary */
    if (!_keff_from_fission_rates) {
      _boundary_leakage = new float[_tot_num_tracks]();
    }

    /* Determine the size of arrays for the FSR scalar fluxes */
    size = _num_FSRs * _NUM_GROUPS;
    max_size = size;
#ifdef MPIX
    if (_geometry->isDomainDecomposed())
      MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX,
                    _geometry->getMPICart());
#endif

    /* Determine the amount of memory allocated */
    int num_flux_arrays = 2;
    if (_stabilize_transport)
      num_flux_arrays++;

    max_size_mb = (double) (num_flux_arrays * max_size * sizeof(FP_PRECISION))
        / (double) (1e6);
    log_printf(NORMAL, "Max scalar flux storage per domain = %6.2f MB",
               max_size_mb);

    /* Allocate scalar fluxes */
    _scalar_flux = new FP_PRECISION[size]();
    _old_scalar_flux = new FP_PRECISION[size]();

#ifdef ONLYVACUUMBC
    _track_flux_sent.resize(2);
    _track_flux_sent.at(0).resize(_tot_num_tracks, false);
    _track_flux_sent.at(1).resize(_tot_num_tracks, false);
#endif

    /* Allocate stabilizing flux vector if necessary */
    if (_stabilize_transport) {
      _stabilizing_flux = new FP_PRECISION[size]();
    }

#ifdef MPIx
    /* Allocate memory for angular flux exchanging buffers */
    if (_geometry->isDomainDecomposed())
      setupMPIBuffers();
#endif
  }

  catch (std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the fluxes");
  }
}


/**
 * @brief Allocates memory for FSR source arrays.
 * @details Deletes memory for old source arrays if they were allocated for a
 *          previous simulation.
 */
void CPUSolver::initializeSourceArrays() {

  /* Delete old sources arrays if they exist */
  if (_reduced_sources != NULL)
    delete [] _reduced_sources;
  if (_fixed_sources != NULL && !_fixed_sources_initialized)
    delete [] _fixed_sources;

  long size = _num_FSRs * _NUM_GROUPS;

  /* Allocate memory for all source arrays */
  _reduced_sources = new FP_PRECISION[size]();
  if (_fixed_sources_on && !_fixed_sources_initialized)
    _fixed_sources = new FP_PRECISION[size]();

  long max_size = size;
#ifdef MPIX
  if (_geometry->isDomainDecomposed())
    MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX,
                  _geometry->getMPICart());
#endif
  double max_size_mb = (double) (max_size * sizeof(FP_PRECISION))
        / (double) (1e6);
  if (_fixed_sources_on)
    max_size_mb *= 2;
  log_printf(NORMAL, "Max source storage per domain = %6.2f MB",
             max_size_mb);

  /* Populate fixed source array with any user-defined sources */
  if (_fixed_sources_on && !_fixed_sources_initialized)
    initializeFixedSources();
}


/**
 * @brief Populates array of fixed sources assigned by FSR.
 */
void CPUSolver::initializeFixedSources() {

  Solver::initializeFixedSources();

  long fsr_id;
  int group;
  std::pair<int, int> fsr_group_key;
  std::map< std::pair<int, int>, FP_PRECISION >::iterator fsr_iter;

  /* Populate fixed source array with any user-defined sources */
  for (fsr_iter = _fix_src_FSR_map.begin();
       fsr_iter != _fix_src_FSR_map.end(); ++fsr_iter) {

    /* Get the FSR with an assigned fixed source */
    fsr_group_key = fsr_iter->first;
    fsr_id = fsr_group_key.first;
    group = fsr_group_key.second;

    if (group <= 0 || group > _NUM_GROUPS)
      log_printf(ERROR,"Unable to use fixed source for group %d in "
                 "a %d energy group problem", group, _NUM_GROUPS);

    if (fsr_id < 0 || fsr_id >= _num_FSRs)
      log_printf(ERROR,"Unable to use fixed source for FSR %d with only "
                 "%d FSRs in the geometry", fsr_id, _num_FSRs);

    _fixed_sources(fsr_id, group-1) = _fix_src_FSR_map[fsr_group_key];
  }

  /* Remember initialization to avoid re-initializing unless it's necessary */
  _fixed_sources_initialized = true;
}


/**
 * @brief Zero each Track's boundary fluxes for each energy group
 *        (and polar angle in 2D) in the "forward" and "reverse" directions.
 */
void CPUSolver::zeroTrackFluxes() {

#pragma omp parallel for schedule(static)
  for (long t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int pe=0; pe < _fluxes_per_track; pe++) {
        _boundary_flux(t, d, pe) = 0.0;
#ifndef ONLYVACUUMBC
        _start_flux(t, d, pe) = 0.0;
#endif
      }
    }
  }
}


/**
 * @brief Copies values from the start flux into the boundary flux array
 *        for both the "forward" and "reverse" directions.
 */
void CPUSolver::copyBoundaryFluxes() {

#pragma omp parallel for schedule(static)
  for (long t=0; t < _tot_num_tracks; t++)
    for (int d=0; d < 2; d++)
      for (int pe=0; pe < _fluxes_per_track; pe++)
        _boundary_flux(t,d,pe) = _start_flux(t, d, pe);
}


/**
 * @brief Computes the total current impingent on boundary CMFD cells from
 *        starting angular fluxes
 */
//FIXME: make suitable for 2D
void CPUSolver::tallyStartingCurrents() {

#pragma omp parallel for schedule(static)
  for (long t=0; t < _tot_num_tracks; t++) {

    /* Get 3D Track data */
    TrackGenerator3D* track_generator_3D =
        dynamic_cast<TrackGenerator3D*>(_track_generator);
    if (track_generator_3D != NULL) {
      TrackStackIndexes tsi;
      Track3D track;
      track_generator_3D->getTSIByIndex(t, &tsi);
      track_generator_3D->getTrackOTF(&track, &tsi);

      /* Determine the first and last CMFD cells of each track */
      double azim = track.getPhi();
      double polar = track.getTheta();
      double delta_x = cos(azim) * sin(polar) * TINY_MOVE;
      double delta_y = sin(azim) * sin(polar) * TINY_MOVE;
      double delta_z = cos(polar) * TINY_MOVE;
      Point* start = track.getStart();
      Point* end = track.getEnd();

      /* Get the track weight */
      int azim_index = track.getAzimIndex();
      int polar_index = track.getPolarIndex();
      double weight = _quad->getWeightInline(azim_index, polar_index);

      /* Tally currents */
      _cmfd->tallyStartingCurrent(start, delta_x, delta_y, delta_z,
                                  &_start_flux(t, 0, 0), weight);
      _cmfd->tallyStartingCurrent(end, -delta_x, -delta_y, -delta_z,
                                  &_start_flux(t, 1, 0), weight);
    }
    else {
      log_printf(ERROR, "Starting currents not implemented yet for 2D MOC");
    }
  }
}


#ifdef MPIx
/**
 * @brief Buffers used to transfer angular flux information are initialized
 * @details Track connection book-keeping information is also saved for
 *          efficiency during angular flux packing.
 */
void CPUSolver::setupMPIBuffers() {

  /* Determine the size of the buffers */
  _track_message_size = _fluxes_per_track + 3;
  int message_length = TRACKS_PER_BUFFER * _track_message_size;

  /* Initialize MPI requests and status */
  if (_geometry->isDomainDecomposed()) {

    if (_send_buffers.size() > 0)
      deleteMPIBuffers();

    log_printf(NORMAL, "Setting up MPI Buffers for angular flux exchange...");

    /* Fill the hash map of send buffers */
    int idx = 0;
    for (int dx=-1; dx <= 1; dx++) {
      for (int dy=-1; dy <= 1; dy++) {
        for (int dz=-1; dz <= 1; dz++) {
          if (abs(dx) + abs(dy) == 1 ||
              (dx == 0 && dy == 0 && dz != 0)) {
            int domain = _geometry->getNeighborDomain(dx, dy, dz);
            if (domain != -1) {
              _neighbor_connections.insert({domain, idx});
              _neighbor_domains.push_back(domain);
              idx++;

              /* Inititalize vector that shows how filled send_buffers are */
              _send_buffers_index.push_back(0);
            }
          }
        }
      }
    }

    /* Estimate and print size of flux transfer buffers */
    int num_domains = _neighbor_domains.size();
    int size = 2 * message_length * num_domains * sizeof(float);
    int max_size;
    MPI_Allreduce(&size, &max_size, 1, MPI_INT, MPI_MAX,
                  _geometry->getMPICart());
    log_printf(INFO_ONCE, "Max track fluxes transfer buffer storage = %.2f MB",
               max_size / 1e6);

    /* Allocate track fluxes transfer buffers */
    _send_buffers.resize(num_domains);
    _receive_buffers.resize(num_domains);
    for (int i=0; i < num_domains; i++) {
#ifdef ONLYVACUUMBC
      /* Increase capacity because buffers will overflow and need a resize */
      _send_buffers.at(i).reserve(3*message_length);
      _receive_buffers.at(i).reserve(3*message_length);
#endif
      _send_buffers.at(i).resize(message_length);
      _receive_buffers.at(i).resize(message_length);
    }

    /* Setup Track communication information for all neighbor domains */
    _boundary_tracks.resize(num_domains);
    for (int i=0; i < num_domains; i++) {

      /* Initialize Track ID's to -1 */
      int start_idx = _fluxes_per_track + 1;
      for (int idx = start_idx; idx < message_length;
           idx += _track_message_size) {
        long* track_info_location =
             reinterpret_cast<long*>(&_send_buffers.at(i)[idx]);
        track_info_location[0] = -1;
        track_info_location =
             reinterpret_cast<long*>(&_receive_buffers.at(i)[idx]);
        track_info_location[0] = -1;
      }
    }

    /* Allocate vector of send/receive buffer sizes */
    _send_size.resize(num_domains, 0);
#ifdef ONLYVACUUMBC
    _receive_size.resize(num_domains, TRACKS_PER_BUFFER);
#endif

    /* Build array of Track connections */
    _track_connections.resize(2);
    _track_connections.at(0).resize(_tot_num_tracks);
    _track_connections.at(1).resize(_tot_num_tracks);

#ifdef ONLYVACUUMBC
    _domain_connections.resize(2);
    _domain_connections.at(0).resize(_tot_num_tracks);
    _domain_connections.at(1).resize(_tot_num_tracks);
#endif

    /* Determine how many Tracks communicate with each neighbor domain */
    log_printf(NORMAL, "Initializing Track connections accross domains...");
    std::vector<long> num_tracks;
    num_tracks. resize(num_domains, 0);

#pragma omp parallel for
    for (long t=0; t<_tot_num_tracks; t++) {

      Track* track;
      /* Get 3D Track data */
      if (_SOLVE_3D) {
        TrackStackIndexes tsi;
        track = new Track3D();
        TrackGenerator3D* track_generator_3D =
          dynamic_cast<TrackGenerator3D*>(_track_generator);
        track_generator_3D->getTSIByIndex(t, &tsi);
        track_generator_3D->getTrackOTF(dynamic_cast<Track3D*>(track), &tsi);
      }
      /* Get 2D Track data */
      else {
        Track** tracks = _track_generator->get2DTracksArray();
        track = tracks[t];
      }

      /* Save the index of the forward and backward connecting Tracks */
      _track_connections.at(0).at(t) = track->getTrackNextFwd();
      _track_connections.at(1).at(t) = track->getTrackNextBwd();

      /* Determine the indexes of connecting domains */
      int domains[2];
      domains[0] = track->getDomainFwd();
      domains[1] = track->getDomainBwd();
      bool interface[2];
      interface[0] = track->getBCFwd() == INTERFACE;
      interface[1] = track->getBCBwd() == INTERFACE;
      for (int d=0; d < 2; d++) {
        if (domains[d] != -1 && interface[d]) {
          int neighbor = _neighbor_connections.at(domains[d]);
#pragma omp atomic update
          num_tracks[neighbor]++;
        }
      }
      if (_SOLVE_3D)
        delete track;

    }

    /* Resize the buffers for the counted number of Tracks */
    for (int i=0; i < num_domains; i++) {
      _boundary_tracks.at(i).resize(num_tracks[i]);
      num_tracks[i] = 0;
    }

    /* Determine which Tracks communicate with each neighbor domain */
#ifndef ONLYVACUUMBC
#pragma omp parallel for
#endif
    for (long t=0; t<_tot_num_tracks; t++) {

      Track* track;
      /* Get 3D Track data */
      if (_SOLVE_3D) {
        TrackStackIndexes tsi;
        track = new Track3D();
        TrackGenerator3D* track_generator_3D =
          dynamic_cast<TrackGenerator3D*>(_track_generator);
        track_generator_3D->getTSIByIndex(t, &tsi);
        track_generator_3D->getTrackOTF(dynamic_cast<Track3D*>(track), &tsi);
      }
      /* Get 2D Track data */
      else {
        Track** tracks = _track_generator->get2DTracksArray();
        track = tracks[t];
      }

      /* Determine the indexes of connecting domains */
      int domains[2];
      domains[0] = track->getDomainFwd();
      domains[1] = track->getDomainBwd();
      bool interface[2];
      interface[0] = track->getBCFwd() == INTERFACE;
      interface[1] = track->getBCBwd() == INTERFACE;
      for (int d=0; d < 2; d++) {
        if (domains[d] != -1 && interface[d]) {
          int neighbor = _neighbor_connections.at(domains[d]);

          long slot;
#pragma omp atomic capture
          {
            slot = num_tracks[neighbor];
            num_tracks[neighbor]++;
          }
          _boundary_tracks.at(neighbor).at(slot) = 2*t + d;
#ifdef ONLYVACUUMBC
          //NOTE _boundary_tracks needs to be ordered if ONLYVACUUMBC is used
          _domain_connections.at(d).at(t) = domains[d];
#endif
        }
#ifdef ONLYVACUUMBC
        /* Keep a list of the tracks that start at vacuum BC */
        else {
#pragma omp critical
          {
            _tracks_from_vacuum.push_back(2*t + 1 - d);
          }
          //NOTE domains[d] can be set for a track ending on an interface
          // with vacuum and another track
          _domain_connections.at(d).at(t) = -1;
        }
#endif
      }
      if (_SOLVE_3D)
        delete track;
    }

    printLoadBalancingReport();
    log_printf(NORMAL, "Finished setting up MPI buffers...");

    /* Setup MPI communication bookkeeping */
    _MPI_requests = new MPI_Request[2*num_domains];
    _MPI_sends = new bool[num_domains];
    _MPI_receives = new bool[num_domains];
    for (int i=0; i < num_domains; i++) {
      _MPI_sends[i] = false;
      _MPI_receives[i] = false;
    }
  }
}


/**
 * @brief The arrays used to store angular flux information are deleted along
 *        with book-keeping information for track connections.
 */
void CPUSolver::deleteMPIBuffers() {
  for (int i=0; i < _send_buffers.size(); i++) {
    _send_buffers.at(i).clear();
  }
  _send_buffers.clear();

  for (int i=0; i < _receive_buffers.size(); i++) {
    _receive_buffers.at(i).clear();
  }
  _receive_buffers.clear();
  _neighbor_domains.clear();

  for (int i=0; i < _boundary_tracks.size(); i++)
    _boundary_tracks.at(i).clear();
  _boundary_tracks.clear();

  delete [] _MPI_requests;
  delete [] _MPI_sends;
  delete [] _MPI_receives;
}
#endif


#ifdef ONLYVACUUMBC
/**
 * @brief Resets the track angular fluxes at the vacuum boundary condition.
 */
void CPUSolver::resetBoundaryFluxes() {

  if (_geometry->isDomainDecomposed()) {
    log_printf(INFO, "Setting %ld incoming vaccum boundary fluxes to 0.",
               _tracks_from_vacuum.size());

#pragma omp parallel for
    for (long i=0; i < _tracks_from_vacuum.size(); i++) {

      long t_id = _tracks_from_vacuum.at(i) / 2;
      int dir = _tracks_from_vacuum.at(i) - 2 * t_id;
      for (int pe=0; pe < _fluxes_per_track; pe++)
        _boundary_flux(t_id, dir, pe) = 0.;
    }
  }
  else {
#pragma omp parallel for
    for (long t=0; t < _tot_num_tracks; t++)
      for (int d=0; d < 2; d++)
        for (int pe=0; pe < _fluxes_per_track; pe++)
          _boundary_flux(t,d,pe) = 0.;
  }
}
#endif


#ifdef MPIx
/**
 * @brief Prints out tracking information for cycles, traversing domain
 *        interfaces.
 * @details This function prints Track starting and ending points for a cycle
 *          that traverses the entire Geometry.
 * @param track_start The starting Track ID from which the cycle is followed
 * @param domain_start The domain for the starting Track
 * @param length The number of Tracks to follow across the cycle
 */
void CPUSolver::printCycle(long track_start, int domain_start, int length) {

  /* Initialize buffer for MPI communication */
  int message_size = sizeof(sendInfo);

  /* Initialize MPI requests and status */
  MPI_Comm MPI_cart = _geometry->getMPICart();
  int num_ranks;
  MPI_Comm_size(MPI_cart, &num_ranks);
  MPI_Status stat;
  MPI_Request request[num_ranks];

  int rank;
  MPI_Comm_rank(MPI_cart, &rank);

  /* Loop over all tracks and exchange fluxes */
  long curr_track = track_start;
  int curr_rank = domain_start;
  bool fwd = true;
  for (int t=0; t < length; t++) {

    /* Check if this rank is sending the Track */
    if (rank == curr_rank) {

      /* Get 3D Track data */
      TrackStackIndexes tsi;
      Track3D track;
      TrackGenerator3D* track_generator_3D =
        dynamic_cast<TrackGenerator3D*>(_track_generator);
      track_generator_3D->getTSIByIndex(curr_track, &tsi);
      track_generator_3D->getTrackOTF(&track, &tsi);

      /* Get connecting tracks */
      long connect;
      bool connect_fwd;
      Point* start;
      Point* end;
      int next_domain;
      if (fwd) {
        connect = track.getTrackPrdcFwd();
        connect_fwd = track.getNextFwdFwd();
        start = track.getStart();
        end = track.getEnd();
        next_domain = track.getDomainFwd();
      }
      else {
        connect = track.getTrackPrdcBwd();
        connect_fwd = track.getNextBwdFwd();
        start = track.getEnd();
        end = track.getStart();
        next_domain = track.getDomainBwd();
      }

      /* Write information */
      log_printf(NODAL, "Rank %d: Track (%f, %f, %f) -> (%f, %f, %f)", rank,
                 start->getX(), start->getY(), start->getZ(), end->getX(),
                 end->getY(), end->getZ());

      /* Check domain for reflected boundaries */
      if (next_domain == -1) {
        next_domain = curr_rank;
        if (fwd)
          connect = track.getTrackNextFwd();
        else
          connect = track.getTrackNextBwd();
      }

      /* Pack the information */
      sendInfo si;
      si.track_id = connect;
      si.domain = next_domain;
      si.fwd = connect_fwd;

      /* Send the information */
      for (int i=0; i < num_ranks; i++)
        if (i != rank)
          MPI_Isend(&si, message_size, MPI_BYTE, i, 0, MPI_cart, &request[i]);

      /* Copy information */
      curr_rank = next_domain;
      fwd = connect_fwd;
      curr_track = connect;

      /* Wait for sends to complete */
      bool complete = false;
      while (!complete) {
        complete = true;
        for (int i=0; i < num_ranks; i++) {
          if (i != rank) {
            int flag;
            MPI_Test(&request[i], &flag, &stat);
            if (flag == 0)
              complete = false;
          }
        }
      }
    }

    /* Receiving info */
    else {

      /* Create object to receive sent information */
      sendInfo si;

      /* Issue the receive from the current node */
      MPI_Irecv(&si, message_size, MPI_BYTE, curr_rank, 0, MPI_cart,
                &request[0]);

      /* Wait for receive to complete */
      bool complete = false;
      while (!complete) {
        complete = true;
        int flag;
        MPI_Test(&request[0], &flag, &stat);
        if (flag == 0)
          complete = false;
      }

      /* Copy information */
      curr_rank = si.domain;
      fwd = si.fwd;
      curr_track = si.track_id;
    }

    MPI_Barrier(MPI_cart);
  }

  /* Join MPI at the end of communication */
  MPI_Barrier(MPI_cart);
}


/**
 * @brief Angular flux transfer information is packed into buffers.
 * @details On each domain, angular flux and track connection information
 *          is packed into buffers. Each buffer pertains to a neighboring
 *          domain. This function proceeds packing buffers until for each
 *          neighboring domain either all the tracks have been packed or the
 *          associated buffer is full. This provided integer array contains
 *          the index of the last track handled for each neighboring domain.
 *          These numbers are updated at the end with the last track handled.
 * @arg packing_indexes index of last track sent for each neighbor domain
 */
void CPUSolver::packBuffers(std::vector<long> &packing_indexes) {

  /* Fill send buffers for every domain */
  int num_domains = packing_indexes.size();
#pragma omp parallel for num_threads(num_domains)
  for (int i=0; i < num_domains; i++) {

    /* Reset send buffers : start at beginning if the buffer has not been
       prefilled, else start after what has been prefilled */
    int start_idx = _send_buffers_index.at(i) * _track_message_size +
                    _fluxes_per_track + 1;
    int max_idx = _track_message_size * TRACKS_PER_BUFFER;
    for (int idx = start_idx; idx < max_idx; idx += _track_message_size) {
      long* track_info_location =
        reinterpret_cast<long*>(&_send_buffers.at(i)[idx]);
      track_info_location[0] = -1;
    }

    /* Fill send buffers with Track information :
       -start from last track packed (packing_indexes)
       -take into account the pre-filling
       (pre-filling only if ONLYVACUUMBC flag is used)
       -only fill to max_buffer_fill */
    int max_buffer_idx = _boundary_tracks.at(i).size() -
          packing_indexes.at(i);

    if (_send_buffers_index.at(i) + max_buffer_idx > TRACKS_PER_BUFFER)
      max_buffer_idx = TRACKS_PER_BUFFER - _send_buffers_index.at(i);

    /* Keep track of buffer size to avoid sending more fluxes than needed */
    _send_size.at(i) = std::max(_send_buffers_index.at(i) + max_buffer_idx,
         _send_buffers_index.at(i));

    for (int b=0; b < max_buffer_idx; b++) {

      long boundary_track_idx = packing_indexes.at(i) + b;
      long buffer_index = (_send_buffers_index.at(i)+b) * _track_message_size;

#ifdef ONLYVACUUMBC
      /* Exit loop if all fluxes have been sent */
      if (boundary_track_idx >= _boundary_tracks.at(i).size())
        break;
#endif

      /* Get 3D Track data */
      long boundary_track = _boundary_tracks.at(i).at(boundary_track_idx);
      long t = boundary_track / 2;
      int d = boundary_track - 2*t;
      long connect_track = _track_connections.at(d).at(t);

      /* Fill buffer with angular fluxes */
      for (int pe=0; pe < _fluxes_per_track; pe++)
        _send_buffers.at(i)[buffer_index + pe] = _boundary_flux(t,d,pe);

      /* Assign the connecting Track information */
      long idx = buffer_index + _fluxes_per_track;
      _send_buffers.at(i)[idx] = d;
      long* track_info_location =
        reinterpret_cast<long*>(&_send_buffers.at(i)[idx+1]);
      track_info_location[0] = connect_track;

 #ifdef ONLYVACUUMBC
      /* Invalidate track transfer if it has already been sent by prefilling */
      if (_track_flux_sent.at(d).at(t)) {
        track_info_location[0] = long(-2);

        /* Use freed-up spot in send_buffer and keep track of it */
        b--;
        packing_indexes.at(i)++;
      }
#endif
    }

    /* Record the next Track ID, reset index of fluxes in send_buffers */
    packing_indexes.at(i) += std::max(0, max_buffer_idx);
    _send_buffers_index.at(i) = 0;
  }
}


/**
 * @brief Transfers all angular fluxes at interfaces to their appropriate
 *        domain neighbors
 * @details The angular fluxes stored in the _boundary_flux array that
 *          intersect INTERFACE boundaries are transfered to their appropriate
 *          neighbor's _start_flux array at the periodic indexes.
 */
void CPUSolver::transferAllInterfaceFluxes() {

  /* Initialize MPI requests and status */
  MPI_Comm MPI_cart = _geometry->getMPICart();
  MPI_Status stat;

  /* Wait for all MPI Ranks to be done with sweeping */
  _timer->startTimer();
  MPI_Barrier(MPI_cart);
  _timer->stopTimer();
  _timer->recordSplit("Idle time");

  /* Initialize timer for total transfer cost */
  _timer->startTimer();

  /* Get rank of each process */
  int rank;
  MPI_Comm_rank(MPI_cart, &rank);

  /* Create bookkeeping vectors */
  std::vector<long> packing_indexes;

  /* Resize vectors to the number of domains */
  int num_domains = _neighbor_domains.size();
  packing_indexes.resize(num_domains, 0);

  /* Start communication rounds */
  int round_counter = -1;
  while (true) {

    round_counter++;

    /* Pack buffers with angular flux data */
    _timer->startTimer();
    packBuffers(packing_indexes);
    _timer->stopTimer();
    _timer->recordSplit("Packing time");

#ifndef ONLYVACUUMBC
    _timer->startTimer();

    /* Send and receive from all neighboring domains */
    bool communication_complete = true;

    for (int i=0; i < num_domains; i++) {

      /* Get the communicating neighbor domain */
      int domain = _neighbor_domains.at(i);

      /* Check if a send/receive needs to be created */
      long* first_track_idx =
        reinterpret_cast<long*>(&_send_buffers.at(i)[_fluxes_per_track+1]);
      long first_track = first_track_idx[0];

      if (first_track != -1) {

        /* Send outgoing flux */
        if (!_MPI_sends[i]) {
          MPI_Isend(&_send_buffers.at(i)[0], _track_message_size *
                    _send_size.at(i), MPI_FLOAT, domain, 1, MPI_cart,
                    &_MPI_requests[i*2]);
          _MPI_sends[i] = true;
        }
        else
          if (!_MPI_sends[i])
            _MPI_requests[i*2] = MPI_REQUEST_NULL;

        /* Receive incoming flux */
        if (!_MPI_receives[i]) {
          MPI_Irecv(&_receive_buffers.at(i)[0], _track_message_size *
                    TRACKS_PER_BUFFER, MPI_FLOAT, domain, 1, MPI_cart,
                    &_MPI_requests[i*2+1]);
          _MPI_receives[i] = true;
        }
        else
          if (!_MPI_receives[i])
            _MPI_requests[i*2+1] = MPI_REQUEST_NULL;

        /* Mark communication as ongoing */
        communication_complete = false;
      }
      else {
        if (!_MPI_sends[i])
          _MPI_requests[i*2] = MPI_REQUEST_NULL;
        if (!_MPI_receives[i])
          _MPI_requests[i*2+1] = MPI_REQUEST_NULL;
      }
    }

    /* Check if communication is done */
    if (communication_complete) {
      _timer->stopTimer();
      _timer->recordSplit("Communication time");
      break;
    }

    /* Block for communication round to complete */
    //FIXME Not necessary, buffers could be unpacked while waiting
    MPI_Waitall(2 * num_domains, _MPI_requests, MPI_STATUSES_IGNORE);
    _timer->stopTimer();
    _timer->recordSplit("Communication time");

    /* Reset status for next communication round and copy fluxes */
    _timer->startTimer();
    for (int i=0; i < num_domains; i++) {

      /* Reset send */
      _MPI_sends[i] = false;

      /* Copy angular fluxes if necessary */
      if (_MPI_receives[i]) {

        /* Get the buffer for the connecting domain */
        for (int t=0; t < TRACKS_PER_BUFFER; t++) {

          /* Get the Track ID */
          float* curr_track_buffer = &_receive_buffers.at(i)[
                                     t*_track_message_size];
          long* track_idx =
            reinterpret_cast<long*>(&curr_track_buffer[_fluxes_per_track+1]);
          long track_id = track_idx[0];

          /* Break out of loop once buffer is finished */
          if (track_id == -1)
            break;

          /* Check if the angular fluxes are active */
          if (track_id > -1) {
            int dir = curr_track_buffer[_fluxes_per_track];

            for (int pe=0; pe < _fluxes_per_track; pe++)
              _start_flux(track_id, dir, pe) = curr_track_buffer[pe];
          }
        }
      }

      /* Reset receive flag */
      _MPI_receives[i] = false;
    }

    _timer->stopTimer();
    _timer->recordSplit("Unpacking time");
  }

  /* Join MPI at the end of communication */
  MPI_Barrier(MPI_cart);
  _timer->stopTimer();
  _timer->recordSplit("Total transfer time");
}
#else
    /* In while(true) loop, timer started */
    /* Number of communication rounds is bounded */
    long max_boundary_tracks = 0;
    for (int i=0; i < num_domains; i++)
      max_boundary_tracks = std::max(max_boundary_tracks,
                                     long(_boundary_tracks.at(i).size()));
    bool active_communication =
         max_boundary_tracks > (round_counter * TRACKS_PER_BUFFER);

    MPI_Request _MPI_req[2*num_domains];

    /* Set size of received messages, adjust buffer if needed */
    _timer->startTimer();
    for (int i=0; i < num_domains; i++) {

      /* Size of received message, in number of tracks */
      _receive_size.at(i) = -1;
      int domain = _neighbor_domains.at(i);
      if (active_communication) {
        /* Communicate _send_buffers' sizes to adapt _receive_buffers' sizes */
        MPI_Isend(&_send_size.at(i), 1, MPI_INT, domain, 0, MPI_cart,
                  &_MPI_req[i*2]);
        MPI_Irecv(&_receive_size.at(i), 1, MPI_INT, domain, 0, MPI_cart,
                  &_MPI_req[i*2 + 1]);
      }
      else {
        _MPI_req[i*2] = MPI_REQUEST_NULL;
        _MPI_req[i*2 + 1] = MPI_REQUEST_NULL;
      }
    }

    /* Send and receive from all neighboring domains */
    bool communication_complete = true;

    /* Start all sends and receives when the buffers' sizes are known to
       reduce synchronization */
    bool all_transfers_started = false;
    while (!all_transfers_started) {
      all_transfers_started = true;

      for (int i=0; i < num_domains; i++) {

        /* Get the communicating neighbor domain */
        int domain = _neighbor_domains.at(i);

        /* Send/receive fluxes if there are fluxes to be sent, if the size of
           the message is known and if they haven't been sent already */
        if (active_communication) {

          /* Send outgoing flux */
          if (_send_size.at(i) > 0 && !_MPI_sends[i]) {
            MPI_Isend(&_send_buffers.at(i)[0], _track_message_size *
                      _send_size.at(i), MPI_FLOAT, domain, 1, MPI_cart,
                      &_MPI_requests[i*2]);
            _MPI_sends[i] = true;
          }
          else
            if (!_MPI_sends[i])
              _MPI_requests[i*2] = MPI_REQUEST_NULL;

          /* Receive incoming flux */
          if (_receive_size.at(i) > 0 && !_MPI_receives[i]) {

            /* Adjust receiving buffer if incoming message is too large */
            if (_num_iterations == 0)
              if (_receive_size.at(i) > _receive_buffers.at(i).size() /
                                        _track_message_size)
                _receive_buffers.at(i).resize(_receive_size.at(i) *
                                              _track_message_size);

            MPI_Irecv(&_receive_buffers.at(i)[0], _track_message_size *
                      _receive_size.at(i), MPI_FLOAT, domain, 1, MPI_cart,
                      &_MPI_requests[i*2+1]);
            _MPI_receives[i] = true;
          }
          else
            if (!_MPI_receives[i])
              _MPI_requests[i*2+1] = MPI_REQUEST_NULL;

          /* Mark communication as ongoing */
          communication_complete = false;
        }
        else {
          if (!_MPI_sends[i])
            _MPI_requests[i*2] = MPI_REQUEST_NULL;
          if (!_MPI_receives[i])
            _MPI_requests[i*2+1] = MPI_REQUEST_NULL;
        }
        /* Check that all MPI receive calls have been made */
        if (active_communication && !_MPI_receives[i] &&
            _receive_size.at(i) != 0) {
          all_transfers_started = false;
          int flag;
          if (_receive_size.at(i) == -1)
            MPI_Test(&_MPI_req[i*2 + 1], &flag, MPI_STATUSES_IGNORE);
        }
      }
    }

    /* Check if communication is done */
    if (communication_complete) {
      _timer->stopTimer();
      _timer->recordSplit("Communication time");
      break;
    }

    /* Block for communication round to complete */
    //FIXME Not necessary, buffers could be unpacked while waiting
    MPI_Waitall(2 * num_domains, _MPI_requests, MPI_STATUSES_IGNORE);
    _timer->stopTimer();
    _timer->recordSplit("Communication time");

    /* Reset status for next communication round and copy fluxes */
    _timer->startTimer();
    for (int i=0; i < num_domains; i++) {

      /* Reset send */
      _MPI_sends[i] = false;

      /* Copy angular fluxes if necessary */
      if (_MPI_receives[i]) {

        /* Get the buffer for the connecting domain */
        for (int t=0; t < _receive_size.at(i); t++) {

          /* Get the Track ID */
          float* curr_track_buffer = &_receive_buffers.at(i)[
                                     t*_track_message_size];
          long* track_idx =
            reinterpret_cast<long*>(&curr_track_buffer[_fluxes_per_track+1]);
          long track_id = track_idx[0];

          /* Break out of loop once buffer is finished */
          if (track_id == -1)
            break;

          /* Check if the angular fluxes are active */
          /* -2 : already transfered through pre-filling
           * -1 : padding of buffer */
          if (track_id > -1) {
            int dir = curr_track_buffer[_fluxes_per_track];

            /* Before copying an incoming flux over an unsent flux, save
             * the unsent flux in the send buffer (pre-filling) */

            int send_domain = _domain_connections.at(dir).at(track_id);
            long boundary_track = 2 * track_id + dir;

            /* Only save destination flux if it doesn't go to vacuum */
            if (send_domain >= 0) {
              int i_next = _neighbor_connections.at(send_domain);

              /* Check that all tracks to i_next haven't been sent */
              if (packing_indexes.at(i_next) <
                  _boundary_tracks.at(i_next).size()) {

                /* Check that the current track hasn't been sent already */
                if (boundary_track >= _boundary_tracks.at(i_next).at(
                    packing_indexes.at(i_next))) {

                  /* Keep track of the send buffer prefilling */
                  int buffer_index = _send_buffers_index.at(i_next);
                  buffer_index *= _track_message_size;
                  _send_buffers_index.at(i_next)++;

                  if (buffer_index >= _send_buffers.at(i_next).size()) {
                    log_printf(WARNING, "MPI angular flux communication buffer"
                               " from rank %d to %d overflowed. Buffer memory "
                               "increased dynamically.", rank, send_domain);
                    _send_buffers.at(i_next).resize(_send_buffers.at(
                          i_next).size() + TRACKS_PER_BUFFER *
                          _track_message_size);
                  }

                  /* Copy flux, direction and next track in send_buffer */
                  for (int pe=0; pe < _fluxes_per_track; pe++)
                    _send_buffers.at(i_next).at(buffer_index + pe) =
                         _boundary_flux(track_id, dir, pe);
                  _send_buffers.at(i_next).at(buffer_index + _fluxes_per_track) =
                       dir;
                  long* track_info_location =
                       reinterpret_cast<long*>(&_send_buffers.at(i_next).at(
                       buffer_index + _fluxes_per_track + 1));
                  track_info_location[0] = _track_connections.at(dir).at(
                       track_id);

                  /* Remember that track flux has been placed in send buffer
                   * to avoid sending a wrong track flux when packing buffer */
                  //NOTE Track fluxes are always communicated in the same order
                  if (_num_iterations == 0)
                    _track_flux_sent.at(dir).at(track_id) = true;
                }
              }
            }

            for (int pe=0; pe < _fluxes_per_track; pe++)
              _boundary_flux(track_id, dir, pe) = curr_track_buffer[pe];
          }
        }
      }

      /* Reset receive */
      _MPI_receives[i] = false;
    }

    _timer->stopTimer();
    _timer->recordSplit("Unpacking time");
  }

  /* Join MPI at the end of communication */
  MPI_Barrier(MPI_cart);
  _timer->stopTimer();
  _timer->recordSplit("Total transfer time");
}
#endif


/**
 * @brief A debugging tool used to check track links across domains
 * @details Domains are traversed in rank order. For each domain, all tracks
 *          are traversed and for each track, information is requested about
 *          the connecting track - specifically its angles, and its location.
 *          When the information is returned, it is checked by the requesting
 *          domain for consistency.
 *          //NOTE : This routine may only be called after the track fluxes
 *          have been exchanged.
 */
void CPUSolver::boundaryFluxChecker() {

  /* Get MPI information */
  MPI_Comm MPI_cart = _geometry->getMPICart();
  MPI_Request req;
  MPI_Status stat;
  int my_rank;
  MPI_Comm_rank(MPI_cart, &my_rank);
  int num_ranks;
  MPI_Comm_size(MPI_cart, &num_ranks);

  /* Loop over all domains for requesting information */
  int tester = 0;
  int new_tester = 0;
  while (tester < num_ranks) {

    if (tester == my_rank) {

      log_printf(NODAL, "Checking boundary fluxes for process %d", tester);

      /* Loop over all tracks */
      for (long t=0; t < _tot_num_tracks; t++) {

        /* Get the Track */
        TrackStackIndexes tsi;
        Track3D track;
        TrackGenerator3D* track_generator_3D =
          dynamic_cast<TrackGenerator3D*>(_track_generator);
        track_generator_3D->getTSIByIndex(t, &tsi);
        track_generator_3D->getTrackOTF(&track, &tsi);

        /* Check connection information on both directions */
        for (int dir=0; dir < 2; dir++) {

          boundaryType bc;
          if (dir == 0)
            bc = track.getBCFwd();
          else
            bc = track.getBCBwd();

          /* Check domain boundaries */
          if (bc == INTERFACE) {

            /* Get connection information */
            int dest;
            long connection[2];
            if (dir == 0) {
              dest = track.getDomainFwd();
              connection[0] = track.getTrackNextFwd();
            }
            else {
              dest = track.getDomainBwd();
              connection[0] = track.getTrackNextBwd();
            }
            connection[1] = dir;

            /* Check for a valid destination */
            if (dest == -1)
              log_printf(ERROR, "Track %d on domain %d has been found to have "
                         "an INTERFACE boundary but no connecting domain", t,
                         my_rank);

            /* Send a request for info */
            MPI_Send(connection, 2, MPI_LONG, dest, 0, MPI_cart);

            /* Receive infomation */
            int receive_size = _fluxes_per_track + 2 * 5;
            float buffer[receive_size];
            MPI_Recv(buffer, receive_size, MPI_FLOAT, dest, 0, MPI_cart, &stat);

            /* Unpack received information */
            float angular_fluxes[_fluxes_per_track];
            for (int i=0; i < _fluxes_per_track; i++)
              angular_fluxes[i] = buffer[i];

            double track_info[5];
            for (int i=0; i < 5; i++) {
              int idx = _fluxes_per_track + 2 * i;
              double* track_info_location =
                reinterpret_cast<double*>(&buffer[idx]);
              track_info[i] = track_info_location[0];
            }
            double x = track_info[0];
            double y = track_info[1];
            double z = track_info[2];
            double theta = track_info[3];
            double phi = track_info[4];

            /* Check received information */

            /* Get the connecting point */
            Point* point;
            if (dir == 0)
              point = track.getEnd();
            else
              point = track.getStart();

            /* Check position */
            if (fabs(point->getX() - x) > 1e-5 ||
                fabs(point->getY() - y) > 1e-5 ||
                fabs(point->getZ() - z) > 1e-5)
              log_printf(ERROR, "Track linking error: Track %d in domain %d "
                         "with connecting point (%f, %f, %f) does not connect "
                         "with \n Track %d in domain %d at point (%f, %f, %f)",
                         t, my_rank, point->getX(), point->getY(),
                         point->getZ(), connection[0], dest, x, y, z);

            /* Check double reflection */
            bool x_min = fabs(point->getX() - _geometry->getMinX()) < 1e-5;
            bool x_max = fabs(point->getX() - _geometry->getMaxX()) < 1e-5;
            bool x_bound = x_min || x_max;
            bool z_min = fabs(point->getZ() - _geometry->getMinZ()) < 1e-5;
            bool z_max = fabs(point->getZ() - _geometry->getMaxZ()) < 1e-5;
            bool z_bound = z_min || z_max;

            /* Forgive angle differences on double reflections */
            if (x_bound && z_bound) {
              phi = track.getPhi();
              theta = track.getTheta();
            }

            if (fabs(track.getPhi() - phi) > 1e-5 ||
                fabs(track.getTheta() - theta) > 1e-5)
              log_printf(ERROR, "Track linking error: Track %d in domain %d "
                         "with direction (%f, %f) does not match Track %d in "
                         " domain %d with direction (%f, %f)",
                         t, my_rank, track.getTheta(), track.getPhi(),
                         connection[0], dest, theta, phi);

            for (int pe=0; pe < _fluxes_per_track; pe++) {
              if (fabs(angular_fluxes[pe] - _boundary_flux(t, dir, pe))
                  > 1e-7) {
                std::string dir_string;
                if (dir == 0)
                  dir_string = "FWD";
                else
                  dir_string = "BWD";
                log_printf(ERROR, "Angular flux mismatch found on Track %d "
                           "in domain %d in %s direction at index %d. Boundary"
                           " angular flux at this location is %f but the "
                           "starting flux at connecting Track %d in domain %d "
                           "in the -- direction is %f", t, my_rank,
                           dir_string.c_str(), pe, _boundary_flux(t, dir, pe),
                           connection[0], dest, angular_fluxes[pe]);
              }
            }
          }

          /* Check on-node boundaries */
          else {

            /* Get the connecting Track */
            long connecting_idx;
            if (dir == 0)
              connecting_idx = track.getTrackNextFwd();
            else
              connecting_idx = track.getTrackNextBwd();

            TrackStackIndexes connecting_tsi;
            Track3D connecting_track;
            track_generator_3D->getTSIByIndex(connecting_idx, &connecting_tsi);
            track_generator_3D->getTrackOTF(&connecting_track,
                                            &connecting_tsi);

            /* Extract Track information */
            double x, y, z;
            bool connect_fwd;
            Point* point;
            if (dir == 0) {
              connect_fwd = track.getNextFwdFwd();
              point = track.getEnd();
            }
            else {
              connect_fwd = track.getNextBwdFwd();
              point = track.getStart();
            }
            if (connect_fwd) {
              x = connecting_track.getStart()->getX();
              y = connecting_track.getStart()->getY();
              z = connecting_track.getStart()->getZ();
            }
            else {
              x = connecting_track.getEnd()->getX();
              y = connecting_track.getEnd()->getY();
              z = connecting_track.getEnd()->getZ();
            }
            double phi = connecting_track.getPhi();
            double theta = connecting_track.getTheta();

            /* Check for a vacuum boundary condition */
            if (bc == VACUUM)
              memset(&_boundary_flux(t, dir, 0), 0, sizeof(float) *
                                                    _fluxes_per_track);

            /* Check angular fluxes */
            for (int pe=0; pe < _fluxes_per_track; pe++) {
              if (fabs(_start_flux(connecting_idx, !connect_fwd, pe)
                  - _boundary_flux(t, dir, pe)) > 1e-7) {
                std::string dir_string;
                std::string dir_conn_string;
                if (dir == 0)
                  dir_string = "FWD";
                else
                  dir_string = "BWD";
                if (connect_fwd)
                  dir_conn_string = "FWD";
                else
                  dir_conn_string = "BWD";
                log_printf(ERROR, "Angular flux mismatch found on Track %d "
                           "in domain %d in %s direction at index %d. Boundary"
                           " angular flux at this location is %f but the "
                           "starting flux at connecting Track %d in domain %d "
                           "in the %s direction is %f", t, my_rank, dir_string.c_str(),
                           pe, _boundary_flux(t, dir, pe), connecting_idx,
                           my_rank,  dir_conn_string.c_str(),
                           _start_flux(connecting_idx, !connect_fwd, pe));
              }
            }

            /* Test reflective boundaries */
            if (bc == REFLECTIVE) {

              /* Check that the reflecting Track has a different direction */
              if (fabs(phi - track.getPhi()) < 1e-5 &&
                  fabs(theta - track.getTheta()) < 1e-5)
                log_printf(ERROR, "Reflective boundary found on Track %d "
                           "with azimuthal angle %f and polar angle %f but "
                           "the reflective Track at index %d has the same "
                           "angles.", t, phi, theta, connecting_idx);

              /* Check that the reflecting Track shares the connecting point */
              if (fabs(point->getX() - x) > 1e-5 ||
                  fabs(point->getY() - y) > 1e-5 ||
                  fabs(point->getZ() - z) > 1e-5) {
                log_printf(ERROR, "Track linking error: Reflective Track %d "
                           "with connecting point (%f, %f, %f) does not "
                           "connect with Track %d at point (%f, %f, %f)",
                           t, point->getX(), point->getY(), point->getZ(),
                           connecting_idx, x, y, z);
              }
            }

            /* Test periodic boundaries */
            if (bc == PERIODIC) {

              /* Check that the periodic Track has the same direction */
              if (fabs(phi - track.getPhi()) < 1e-5 ||
                  fabs(theta - track.getTheta()) < 1e-5)
                log_printf(ERROR, "Periodic boundary found on Track %d "
                           "with azimuthal angle %f and polar angle %f but "
                           "the periodic Track at index %d has azimuthal "
                           " angle %f and polar angle %f", t, track.getPhi(),
                           track.getTheta(), connecting_idx, phi, theta);

              /* Check that the periodic Track does not share the same
                 connecting point */
              if (fabs(point->getX() - x) < 1e-5 &&
                  fabs(point->getY() - y) < 1e-5 &&
                  fabs(point->getZ() - z) < 1e-5)
                log_printf(ERROR, "Periodic boundary found on Track %d "
                           "at connecting point (%f, %f, %f) but the "
                           "connecting periodic Track at index %d has the "
                           "same connecting point", t, x, y, z,
                           connecting_idx);
            }
          }
        }
      }

      /* Broadcast new tester */
      tester++;
      long broadcast[2];
      broadcast[0] = -1;
      broadcast[1] = tester;
      for (int i=0; i < my_rank; i++) {
        MPI_Send(broadcast, 2, MPI_LONG, i, 0, MPI_cart);
      }
      for (int i = my_rank+1; i < num_ranks; i++) {
        MPI_Send(broadcast, 2, MPI_LONG, i, 0, MPI_cart);
      }
    }
    /* Responder */
    else {

      /* Look for messages */
      int message;
      MPI_Iprobe(tester, MPI_ANY_TAG, MPI_cart, &message, &stat);

      /* Check for an information request */
      if (message) {

        /* Receive the request for information */
        long connection[2];
        MPI_Recv(connection, 2, MPI_LONG, tester, 0, MPI_cart, &stat);

        /* Check for a broadcast of the new tester rank */
        if (connection[0] == -1) {
          tester = connection[1];
        }
        else {
          /* Handle an information request */

          /* Fill the requested information */
          long t = connection[0];
          int dir = connection[1];

          /* Get the Track */
          TrackStackIndexes tsi;
          Track3D track;
          TrackGenerator3D* track_generator_3D =
            dynamic_cast<TrackGenerator3D*>(_track_generator);
          track_generator_3D->getTSIByIndex(t, &tsi);
          track_generator_3D->getTrackOTF(&track, &tsi);

          /* Check that the track should be transfered */
          boundaryType bc_transfer;
          int dest_transfer;
          if (dir) {
            bc_transfer = track.getBCFwd();
            dest_transfer = track.getDomainFwd();
          }
          else {
            bc_transfer = track.getBCBwd();
            dest_transfer = track.getDomainBwd();
          }
          if (bc_transfer != INTERFACE)
            log_printf(NODAL, "Node %d requested track %ld that doesn't end at"
                              "an interface boundary condition.", tester, t);
          if (dest_transfer != tester)
            log_printf(NODAL, "Node %d requested track %ld that connects with "
                              "domain %d.", tester, t, dest_transfer);

          /* Fill the information */
          int send_size = _fluxes_per_track + 2 * 5;
          float buffer[send_size];
          for (int pe=0; pe < _fluxes_per_track; pe++)
            buffer[pe] = _start_flux(t, dir, pe);

          /* Get the connecting point */
          Point* point;
          if (dir == 0)
            point = track.getStart();
          else
            point = track.getEnd();

          /* Fill tracking data */
          double track_data[5];
          track_data[0] = point->getX();
          track_data[1] = point->getY();
          track_data[2] = point->getZ();
          track_data[3] = track.getTheta();
          track_data[4] = track.getPhi();

          for (int i=0; i < 5; i++) {
            int idx = _fluxes_per_track + 2 * i;
            double* track_info_location =
              reinterpret_cast<double*>(&buffer[idx]);
            track_info_location[0] = track_data[i];
          }

          /* Send the information */
          MPI_Send(buffer, send_size, MPI_FLOAT, tester, 0, MPI_cart);
        }
      }
    }
  }
  MPI_Barrier(MPI_cart);
  log_printf(NORMAL, "Passed boundary flux check");
}
#endif


/**
 * @brief Set the scalar flux for each FSR and energy group to some value.
 * @param value the value to assign to each FSR scalar flux
 */
void CPUSolver::flattenFSRFluxes(FP_PRECISION value) {

#pragma omp parallel for schedule(static)
  for (long r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _NUM_GROUPS; e++)
      _scalar_flux(r,e) = value;
  }
}


/**
 * @brief Set the scalar flux for each FSR to a chi spectrum
 */
void CPUSolver::flattenFSRFluxesChiSpectrum() {
  if (_chi_spectrum_material == NULL)
    log_printf(ERROR, "A flattening of the FSR fluxes for a chi spectrum was "
               "requested but no chi spectrum material was set.");

  FP_PRECISION* chi = _chi_spectrum_material->getChi();
#pragma omp parallel for schedule(static)
  for (long r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _NUM_GROUPS; e++)
      _scalar_flux(r,e) = chi[e];
  }
}


/**
 * @brief Stores the FSR scalar fluxes in the old scalar flux array.
 */
void CPUSolver::storeFSRFluxes() {

#pragma omp parallel for schedule(static)
  for (long r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _NUM_GROUPS; e++)
      _old_scalar_flux(r,e) = _scalar_flux(r,e);
  }
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
double CPUSolver::normalizeFluxes() {

  double* int_fission_sources = _regionwise_scratch;

  /* Compute total fission source for each FSR, energy group */
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    FP_PRECISION* group_fission_sources = _groupwise_scratch.at(tid);
#pragma omp for schedule(static)
    for (long r=0; r < _num_FSRs; r++) {

      /* Get pointers to important data structures */
      FP_PRECISION* nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
      FP_PRECISION volume = _FSR_volumes[r];

      for (int e=0; e < _NUM_GROUPS; e++)
        group_fission_sources[e] = nu_sigma_f[e] * _scalar_flux(r,e) * volume;

      int_fission_sources[r] = pairwise_sum<FP_PRECISION>(group_fission_sources,
                                                        _NUM_GROUPS);
    }
  }

  /* Compute the total fission source */
  double tot_fission_source = pairwise_sum<double>(int_fission_sources,
                                                         _num_FSRs);

  /* Get the total number of source regions */
  long total_num_FSRs = _num_FSRs;

#ifdef MPIx
  /* Reduce total fission rates across domains */
  if (_geometry->isDomainDecomposed()) {

    /* Get the communicator */
    MPI_Comm comm = _geometry->getMPICart();

    /* Reduce fission rates */
    double reduced_fission;
    MPI_Allreduce(&tot_fission_source, &reduced_fission, 1, MPI_DOUBLE,
                  MPI_SUM, comm);
    tot_fission_source = reduced_fission;

    /* Get total number of FSRs across all domains */
    MPI_Allreduce(&_num_FSRs, &total_num_FSRs, 1, MPI_LONG, MPI_SUM, comm);
  }
#endif

  /* Normalize scalar fluxes in each FSR */
  double norm_factor = total_num_FSRs / tot_fission_source;

  log_printf(DEBUG, "Tot. Fiss. Src. = %f, Norm. factor = %f",
             tot_fission_source, norm_factor);

#pragma omp parallel for schedule(static)
  for (long r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _NUM_GROUPS; e++)
      _scalar_flux(r, e) *= norm_factor;
  }

  /* Normalize angular boundary fluxes for each Track */
#pragma omp parallel for schedule(static)
  for (long idx=0; idx < 2 * _tot_num_tracks * _fluxes_per_track; idx++) {
#ifndef ONLYVACUUMBC
    _start_flux[idx] *= norm_factor;
#endif
    _boundary_flux[idx] *= norm_factor;
  }

  return norm_factor;
}


/**
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void CPUSolver::computeFSRSources(int iteration) {

  long num_negative_fsrs = 0;
  long num_negative_sources = 0;

  /* For all FSRs, find the source */
#pragma omp parallel for schedule(static)
  for (long r=0; r < _num_FSRs; r++) {

    Material* material = _FSR_materials[r];
    FP_PRECISION* sigma_s = material->getSigmaS();
    FP_PRECISION fiss_mat;
    FP_PRECISION fission_sources[_NUM_GROUPS];
    FP_PRECISION scatter_sources[_NUM_GROUPS];
    bool negative_source_in_fsr = false;

    /* Compute total (fission+scatter+fixed) source for group G */
    for (int G=0; G < _NUM_GROUPS; G++) {
      int first_idx = G * _NUM_GROUPS;
      fiss_mat = 0;
      for (int g=0; g < _NUM_GROUPS; g++) {
        if (material->isFissionable())
          fiss_mat = material->getFissionMatrixByGroup(g+1,G+1);
        scatter_sources[g] = sigma_s[first_idx+g] * _scalar_flux(r,g);
        fission_sources[g] = _scalar_flux(r,g) * fiss_mat;
      }
      double scatter_source =
          pairwise_sum<FP_PRECISION>(scatter_sources, _NUM_GROUPS);
      double fission_source = pairwise_sum<FP_PRECISION>(fission_sources,
                                                  _NUM_GROUPS);
      fission_source /= _k_eff;
      _reduced_sources(r,G) = fission_source;
      _reduced_sources(r,G) += scatter_source;
      if (_fixed_sources_on)
        _reduced_sources(r,G) += _fixed_sources(r,G);
      _reduced_sources(r,G) *= ONE_OVER_FOUR_PI;

      /* Correct negative sources to (near) zero */
      if (_reduced_sources(r,G) < 0.0) {
#pragma omp atomic update
        num_negative_sources++;
        negative_source_in_fsr = true;
        if (iteration < 30 && !_negative_fluxes_allowed)
          _reduced_sources(r,G) = FLUX_EPSILON;
      }
    }

    if (negative_source_in_fsr)
#pragma omp atomic update
      num_negative_fsrs++;
  }

  /* Tally the total number of negative source across the entire problem */
  long total_num_negative_sources = num_negative_sources;
  long total_num_negative_fsrs = num_negative_fsrs;
  int num_negative_source_domains = (num_negative_sources > 0);
  int total_num_negative_source_domains = num_negative_source_domains;
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    MPI_Allreduce(&num_negative_sources, &total_num_negative_sources, 1,
                  MPI_LONG, MPI_SUM, _geometry->getMPICart());
    MPI_Allreduce(&num_negative_fsrs, &total_num_negative_fsrs, 1,
                  MPI_LONG, MPI_SUM, _geometry->getMPICart());
    MPI_Allreduce(&num_negative_source_domains,
                  &total_num_negative_source_domains, 1,
                  MPI_INT, MPI_SUM, _geometry->getMPICart());
  }
#endif

  /* Report negative sources */
  if (total_num_negative_sources > 0 && !_negative_fluxes_allowed) {
    if (_geometry->isRootDomain()) {
      log_printf(WARNING, "Computed %ld negative sources in %ld fsrs on %d "
                 "domains", total_num_negative_sources,
                 total_num_negative_fsrs, total_num_negative_source_domains);
      if (iteration < 30)
        log_printf(WARNING, "Negative sources corrected to zero");
    }

    /* Output negative sources for debugging */
    if ((_print_negative_sources || get_log_level() == DEBUG) && _cmfd != NULL)
      printNegativeSources(_num_iterations, _cmfd->getNumX(), _cmfd->getNumY(),
                           _cmfd->getNumZ());
  }
}


/**
 * @brief Computes the total fission source in each FSR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeFSRFissionSources() {

#pragma omp parallel default(none)
  {
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION fiss_mat;
    FP_PRECISION fission_source;
    FP_PRECISION fission_sources[_NUM_GROUPS];

    /* Compute the total source for each FSR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      material = _FSR_materials[r];

      /* Compute fission source for group g */
      //NOTE use full fission matrix instead of chi because of transpose
      for (int g=0; g < _NUM_GROUPS; g++) {
        fiss_mat = 0;
        for (int g_prime=0; g_prime < _NUM_GROUPS; g_prime++) {
          if (material->isFissionable())
            fiss_mat = material->getFissionMatrixByGroup(g_prime+1,g+1);
          fission_sources[g_prime] = fiss_mat * _scalar_flux(r,g_prime);
        }

        fission_source = pairwise_sum<FP_PRECISION>(fission_sources,
                                                    _NUM_GROUPS);

        /* Compute total (fission) reduced source */
        _reduced_sources(r,g) = fission_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI;
      }
    }
  }
}


/**
 * @brief Computes the total scattering source in each FSR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeFSRScatterSources() {

#pragma omp parallel default(none)
  {
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s;
    FP_PRECISION scatter_source;
    FP_PRECISION scatter_sources[_NUM_GROUPS];

    /* Compute the total source for each FSR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      material = _FSR_materials[r];

      /* Compute scatter source for group g */
      for (int g=0; g < _NUM_GROUPS; g++) {
        for (int g_prime=0; g_prime < _NUM_GROUPS; g_prime++) {
          sigma_s = material->getSigmaSByGroup(g_prime+1,g+1);
          scatter_sources[g_prime] = sigma_s * _scalar_flux(r,g_prime);
        }

        scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources,
                                                    _NUM_GROUPS);

        /* Compute total (scatter) reduced source */
        _reduced_sources(r,g) = scatter_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI;
      }
    }
  }
}


/**
 * @brief Computes the residual between source/flux iterations.
 * @param res_type the type of residuals to compute
 *        (SCALAR_FLUX, FISSION_SOURCE, TOTAL_SOURCE)
 * @return the average residual in each FSR
 */
double CPUSolver::computeResidual(residualType res_type) {

  long norm;
  double residual;
  double* residuals = _regionwise_scratch;
  memset(residuals, 0, _num_FSRs * sizeof(double));

  FP_PRECISION* reference_flux = _old_scalar_flux;
  if (_calculate_residuals_by_reference)
    reference_flux = _reference_flux;

  if (res_type == SCALAR_FLUX) {

    norm = _num_FSRs;

#pragma omp parallel for schedule(static)
    for (long r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _NUM_GROUPS; e++)
        if (reference_flux(r,e) > 0.) {
          residuals[r] += pow((_scalar_flux(r,e) - reference_flux(r,e)) /
                              reference_flux(r,e), 2);
      }
    }
  }

  else if (res_type == FISSION_SOURCE) {

    norm = _num_fissionable_FSRs;

    double new_fission_source, old_fission_source;
    FP_PRECISION* nu_sigma_f;
    Material* material;

#pragma omp parallel for private(new_fission_source, old_fission_source,\
     material, nu_sigma_f) schedule(static)
    for (long r=0; r < _num_FSRs; r++) {
      new_fission_source = 0.;
      old_fission_source = 0.;
      material = _FSR_materials[r];

      if (material->isFissionable()) {
        nu_sigma_f = material->getNuSigmaF();

        for (int e=0; e < _NUM_GROUPS; e++) {
          new_fission_source += _scalar_flux(r,e) * nu_sigma_f[e];
          old_fission_source += reference_flux(r,e) * nu_sigma_f[e];
        }

        if (old_fission_source > 0.)
          residuals[r] = pow((new_fission_source -  old_fission_source) /
                              old_fission_source, 2);
      }
    }
  }

  else if (res_type == TOTAL_SOURCE) {

    norm = _num_FSRs;

    double new_total_source, old_total_source;
    double inverse_k_eff = 1.0 / _k_eff;
    FP_PRECISION* nu_sigma_f;
    Material* material;

#pragma omp parallel for private(new_total_source, old_total_source,\
     material, nu_sigma_f) schedule(static)
    for (long r=0; r < _num_FSRs; r++) {
      new_total_source = 0.;
      old_total_source = 0.;
      material = _FSR_materials[r];

      if (material->isFissionable()) {
        nu_sigma_f = material->getNuSigmaF();

        for (int e=0; e < _NUM_GROUPS; e++) {
          new_total_source += _scalar_flux(r,e) * nu_sigma_f[e];
          old_total_source += reference_flux(r,e) * nu_sigma_f[e];
        }

        new_total_source *= inverse_k_eff;
        old_total_source *= inverse_k_eff;
      }

      /* Compute total scattering source for group G */
      FP_PRECISION* sigma_s = material->getSigmaS();
      for (int G=0; G < _NUM_GROUPS; G++) {
        int first_idx = G * _NUM_GROUPS;
        for (int g=0; g < _NUM_GROUPS; g++) {
          new_total_source += sigma_s[first_idx+g] * _scalar_flux(r,g);
          old_total_source += sigma_s[first_idx+g] * reference_flux(r,g);
        }
      }

      if (old_total_source > 0.)
        residuals[r] = pow((new_total_source -  old_total_source) /
                            old_total_source, 2);
    }
  }

  /* Sum up the residuals from each FSR and normalize */
  residual = pairwise_sum<double>(residuals, _num_FSRs);

#ifdef MPIx
  /* Reduce residuals across domains */
  if (_geometry->isDomainDecomposed()) {

    /* Get the communicator */
    MPI_Comm comm = _geometry->getMPICart();

    /* Reduce residuals */
    double reduced_res;
    MPI_Allreduce(&residual, &reduced_res, 1, MPI_DOUBLE, MPI_SUM, comm);
    residual = reduced_res;

    /* Reduce normalization factors */
    long reduced_norm;
    MPI_Allreduce(&norm, &reduced_norm, 1, MPI_LONG, MPI_SUM, comm);
    norm = reduced_norm;
  }
#endif

  if (res_type == FISSION_SOURCE && norm == 0)
      log_printf(ERROR, "The Solver is unable to compute a "
                 "FISSION_SOURCE residual without fissionable FSRs");

  /* Error check residual componenets */
  if (residual < 0.0) {
    log_printf(WARNING, "MOC residual mean square error %6.4f less than zero",
               residual);
    residual = 0.0;
  }
  if (norm <= 0) {
    log_printf(WARNING, "MOC residual norm %d less than one", norm);
    norm = 1;
  }

  /* Compute RMS residual */
  residual = sqrt(residual / norm);
  return residual;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from successive fission sources.
 */
void CPUSolver::computeKeff() {

  double* FSR_rates = _regionwise_scratch;
  double rates[3];

  int num_rates = 1;
  if (!_keff_from_fission_rates)
    num_rates = 2;

  /* Loop over all FSRs and compute the volume-integrated total rates */
  for (int type=0; type < num_rates; type++) {
#pragma omp parallel for schedule(static)
    for (long r=0; r < _num_FSRs; r++) {

      int tid = omp_get_thread_num();
      FP_PRECISION* group_rates = _groupwise_scratch.at(tid);
      FP_PRECISION volume = _FSR_volumes[r];
      Material* material = _FSR_materials[r];

      /* Get cross section for desired rate */
      FP_PRECISION* sigma;
      if (type == 0)
        sigma = material->getNuSigmaF();
      else
        sigma = material->getSigmaA();

      for (int e=0; e < _NUM_GROUPS; e++)
        group_rates[e] = sigma[e] * _scalar_flux(r,e);

      FSR_rates[r] = pairwise_sum<FP_PRECISION>(group_rates, _NUM_GROUPS);
      FSR_rates[r] *= volume;
    }

    /* Reduce new fission rates across FSRs */
    rates[type] = pairwise_sum<double>(FSR_rates, _num_FSRs);
  }

  /* Compute total leakage rate */
  if (!_keff_from_fission_rates) {
    rates[2] = pairwise_sum<float>(_boundary_leakage, _tot_num_tracks);
    num_rates=3;
  }

  /* Get the total number of source regions */
  long total_num_FSRs = _num_FSRs;

#ifdef MPIx
  /* Reduce rates across domians */
  if (_geometry->isDomainDecomposed()) {

    /* Get the communicator */
    MPI_Comm comm = _geometry->getMPICart();

    /* Copy local rates */
    double local_rates[num_rates];
    for (int i=0; i < num_rates; i++)
      local_rates[i] = rates[i];

     /* Reduce computed rates */
    MPI_Allreduce(local_rates, rates, num_rates, MPI_DOUBLE, MPI_SUM, comm);

    /* Get total number of FSRs across all domains */
    MPI_Allreduce(&_num_FSRs, &total_num_FSRs, 1, MPI_LONG, MPI_SUM, comm);
  }
#endif
  if (!_keff_from_fission_rates)
    /* Compute k-eff from fission, absorption, and leakage rates */
    _k_eff = rates[0] / (rates[1] + rates[2]);
  else
    _k_eff *= rates[0] / total_num_FSRs;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void CPUSolver::transportSweep() {

  log_printf(DEBUG, "Transport sweep with %d OpenMP threads",
      _num_threads);

  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->zeroCurrents();

  /* Initialize flux in each FSR to zero */
  flattenFSRFluxes(0.0);

#ifndef ONLYVACUUMBC
  /* Copy starting flux to current flux */
  copyBoundaryFluxes();
#endif

  /* Tally the starting fluxes to boundaries */
  if (_cmfd != NULL)
    if (_cmfd->isSigmaTRebalanceOn())
      tallyStartingCurrents();

  /* Zero boundary leakage tally */
  if (_boundary_leakage != NULL)
    memset(_boundary_leakage, 0, _tot_num_tracks * sizeof(float));

  /* Tracks are traversed and the MOC equations from this CPUSolver are applied
     to all Tracks and corresponding segments */
  _timer->startTimer();
  if (_OTF_transport) {
    TransportSweepOTF sweep_tracks(_track_generator);
    sweep_tracks.setCPUSolver(this);
    sweep_tracks.execute();
  }
  else {
    TransportSweep sweep_tracks(this);
    sweep_tracks.execute();
  }

  /* Record sweep time */
  _timer->stopTimer();
  _timer->recordSplit("Transport Sweep");

#ifdef MPIx
  /* Transfer all interface fluxes after the transport sweep */
  if (_track_generator->getGeometry()->isDomainDecomposed())
    transferAllInterfaceFluxes();
#endif

#ifdef ONLYVACUUMBC
  resetBoundaryFluxes();
#endif
}


/**
 * @brief Computes the contribution to the FSR scalar flux from a Track segment.
 * @details This method integrates the angular flux for a Track segment across
 *          energy groups and polar angles, and tallies it into the FSR
 *          scalar flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index azimuthal angle index for this segment
 * @param fsr_flux buffer to store the contribution to the region's scalar flux
 * @param track_flux a pointer to the Track's angular flux
 */
void CPUSolver::tallyScalarFlux(segment* curr_segment,
                                int azim_index,
                                FP_PRECISION* __restrict__ fsr_flux,
                                float* track_flux) {

  long fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();

  if (_SOLVE_3D) {

    // The for loop is cut in chunks of size VEC_LENGTH (strip-mining) to ease
    // vectorization of the loop by the compiler
    // Determine number of SIMD vector groups
    const int num_vector_groups = _NUM_GROUPS / VEC_LENGTH;

    for (int v=0; v < num_vector_groups; v++) {
      int start_vector = v * VEC_LENGTH;

#pragma omp simd aligned(sigma_t, fsr_flux)
      for (int e=start_vector; e < start_vector + VEC_LENGTH; e++) {

        FP_PRECISION tau = sigma_t[e] * length;

        /* Compute the exponential */
        FP_PRECISION exponential;
        expF1_fractional(tau, &exponential);

        /* Compute attenuation and tally the contribution to the scalar flux */
        FP_PRECISION delta_psi = (tau * track_flux[e] - length *
                _reduced_sources(fsr_id, e)) * exponential;
        track_flux[e] -= delta_psi;
        fsr_flux[e] += delta_psi;
      }
    }

    // The rest of the loop is treated separately
#pragma omp simd aligned(sigma_t, fsr_flux)
    for (int e=num_vector_groups * VEC_LENGTH; e < _NUM_GROUPS; e++) {
      FP_PRECISION tau = sigma_t[e] * length;

      /* Compute the exponential */
      FP_PRECISION exponential;
      expF1_fractional(tau, &exponential);

      /* Compute attenuation and tally the contribution to the scalar flux */
      FP_PRECISION delta_psi = (tau * track_flux[e] - length *
              _reduced_sources(fsr_id, e)) * exponential;
      track_flux[e] -= delta_psi;
      fsr_flux[e] += delta_psi;
    }
  }
  else {
//FIXME: Implement strip mining for the 2D flat source solver
    ExpEvaluator* exp_evaluator = _exp_evaluators[azim_index][0];
    const int num_polar_2 = _num_polar / 2;

    /* Compute tau in advance to simplify attenuation loop */
    FP_PRECISION tau[_NUM_GROUPS * num_polar_2]
                 __attribute__ ((aligned(VEC_ALIGNMENT)));

#pragma omp simd aligned(tau)
    for (int pe=0; pe < num_polar_2 * _NUM_GROUPS; pe++)
      tau[pe] = sigma_t[pe % _NUM_GROUPS] * length;

    FP_PRECISION delta_psi[_NUM_GROUPS * num_polar_2]
                 __attribute__ ((aligned(VEC_ALIGNMENT)));

    /* Loop over polar angles and energy groups */
#pragma omp simd aligned(tau, delta_psi)
    for (int pe=0; pe < num_polar_2 * _NUM_GROUPS; pe++) {

      FP_PRECISION wgt = _quad->getWeightInline(azim_index,
                                                int(pe/_NUM_GROUPS));

      /* Compute the exponential */
      FP_PRECISION exponential = exp_evaluator->computeExponential(tau[pe],
                                                int(pe/_NUM_GROUPS));

      /* Compute attenuation of the track angular flux */
      delta_psi[pe] = (tau[pe] * track_flux[pe] - length *
                      _reduced_sources(fsr_id, pe%_NUM_GROUPS)) * exponential;

      track_flux[pe] -= delta_psi[pe];
      delta_psi[pe] *= wgt;
    }

    /* Tally to scalar flux buffer */
    //TODO Change loop to accept 'pe' indexing, and keep vectorized
    for (int p=0; p < num_polar_2; p++) {
#pragma omp simd aligned(fsr_flux)
      for (int e=0; e < _NUM_GROUPS; e++)
        fsr_flux[e] += delta_psi[p*_NUM_GROUPS + e];
    }
  }
}


/**
 * @brief Move the segment(s)' contributions to the scalar flux from the buffer
 * to the global scalar flux array.
 * @param fsr_id the id of the fsr
 * @param weight the quadrature weight (only for 3D ray tracing)
 * @param fsr_flux the buffer containing the segment(s)' contribution
 */
void CPUSolver::accumulateScalarFluxContribution(long fsr_id,
                                                 FP_PRECISION weight,
                                                 FP_PRECISION* __restrict__
                                                 fsr_flux) {

  // Atomically increment the FSR scalar flux from the temporary array
  omp_set_lock(&_FSR_locks[fsr_id]);

  // Add to global scalar flux vector
#pragma omp simd aligned(fsr_flux)
  for (int e=0; e < _NUM_GROUPS; e++)
    _scalar_flux(fsr_id,e) += weight * fsr_flux[e];

  omp_unset_lock(&_FSR_locks[fsr_id]);
#ifdef INTEL
#pragma omp flush
#endif

  /* Reset buffers */
  memset(fsr_flux, 0, _NUM_GROUPS * sizeof(FP_PRECISION));
}


/**
 * @brief Tallies the current contribution from this segment across the
 *        the appropriate CMFD mesh cell surface.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index the azimuthal index for this segment
 * @param polar_index the polar index for this segmenbt
 * @param track_flux a pointer to the Track's angular flux
 * @param fwd boolean indicating direction of integration along segment
 */
void CPUSolver::tallyCurrent(segment* curr_segment, int azim_index,
                             int polar_index, float* track_flux,
                             bool fwd) {

  /* Tally surface currents if CMFD is in use */
  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->tallyCurrent(curr_segment, track_flux, azim_index, polar_index, fwd);
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the Track is given to the reflecting Track. For vacuum
 *          boundary conditions, the outgoing flux is tallied as leakage.
 * @param track a pointer to the Track of interest
 * @param azim_index azimuthal angle index for this segment
 * @param polar_index polar angle index for this segment
 * @param direction the Track direction (forward - true, reverse - false)
 * @param track_flux a pointer to the Track's outgoing angular flux
 */
void CPUSolver::transferBoundaryFlux(Track* track,
                                     int azim_index, int polar_index,
                                     bool direction,
                                     float* track_flux) {

  /* Extract boundary conditions for this Track and the pointer to the
   * outgoing reflective Track, and index into the leakage array */
  boundaryType bc_out;
  long track_out_id;
  int start_out;

  /* For the "forward" direction */
  if (direction) {
    bc_out = track->getBCFwd();
    track_out_id = track->getTrackNextFwd();
    start_out = _fluxes_per_track * (!track->getNextFwdFwd());
  }

  /* For the "reverse" direction */
  else {
    bc_out = track->getBCBwd();
    track_out_id = track->getTrackNextBwd();
    start_out = _fluxes_per_track * (!track->getNextBwdFwd());
  }

  /* Determine if flux should be transferred */
  if (bc_out == REFLECTIVE || bc_out == PERIODIC) {
    float* track_out_flux = &_start_flux(track_out_id, 0, start_out);
    memcpy(track_out_flux, track_flux, _fluxes_per_track * sizeof(float));
  }
  /* For vacuum boundary conditions, losing the flux is enough */

  /* Tally leakage if applicable */
  if (!_keff_from_fission_rates) {
    if (bc_out == VACUUM) {
      long track_id = track->getUid();
      FP_PRECISION weight = _quad->getWeightInline(azim_index, polar_index);
      for (int pe=0; pe < _fluxes_per_track; pe++)
        _boundary_leakage[track_id] += weight * track_flux[pe];
    }
  }
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void CPUSolver::addSourceToScalarFlux() {

  FP_PRECISION volume;
  FP_PRECISION* sigma_t;
  long num_negative_fluxes = 0;

  /* Add in source term and normalize flux to volume for each FSR */
  /* Loop over FSRs, energy groups */
#pragma omp parallel for private(volume, sigma_t) schedule(static)
  for (long r=0; r < _num_FSRs; r++) {
    volume = _FSR_volumes[r];
    sigma_t = _FSR_materials[r]->getSigmaT();

    /* Handle zero volume source region case */
    if (volume < FLT_EPSILON)
      volume = 1e30;

    for (int e=0; e < _NUM_GROUPS; e++) {

      _scalar_flux(r, e) /= (sigma_t[e] * volume);
      _scalar_flux(r, e) += FOUR_PI * _reduced_sources(r, e) / sigma_t[e];

      if (_scalar_flux(r, e) < 0.0 && !_negative_fluxes_allowed) {
        _scalar_flux(r, e) = FLUX_EPSILON;
#pragma omp atomic update
        num_negative_fluxes++;
      }
    }
  }

  /* Tally the total number of negative fluxes across the entire problem */
  long total_num_negative_fluxes = num_negative_fluxes;
  int num_negative_flux_domains = (num_negative_fluxes > 0);
  int total_num_negative_flux_domains = num_negative_flux_domains;
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    MPI_Allreduce(&num_negative_fluxes, &total_num_negative_fluxes, 1,
                  MPI_LONG, MPI_SUM, _geometry->getMPICart());
    MPI_Allreduce(&num_negative_flux_domains,
                  &total_num_negative_flux_domains, 1,
                  MPI_INT, MPI_SUM, _geometry->getMPICart());
  }
#endif

  /* Report negative fluxes */
  if (total_num_negative_fluxes > 0  && !_negative_fluxes_allowed) {
    if (_geometry->isRootDomain()) {
      log_printf(WARNING, "Computed %ld negative fluxes on %d domains",
                 total_num_negative_fluxes, total_num_negative_flux_domains);
    }
  }
}


/**
 * @brief Computes the stabilizing flux for transport stabilization
 */
void CPUSolver::computeStabilizingFlux() {

  if (_stabilization_type == DIAGONAL) {

    /* Loop over all flat source regions */
#pragma omp parallel for schedule(static)
    for (long r=0; r < _num_FSRs; r++) {

      /* Extract the scattering matrix */
      FP_PRECISION* scattering_matrix = _FSR_materials[r]->getSigmaS();

      /* Extract total cross-sections */
      FP_PRECISION* sigma_t = _FSR_materials[r]->getSigmaT();

      for (int e=0; e < _NUM_GROUPS; e++) {

        /* Extract the in-scattering (diagonal) element */
        FP_PRECISION sigma_s = scattering_matrix[e*_NUM_GROUPS+e];

        /* For negative cross-sections, add the absolute value of the
           in-scattering rate to the stabilizing flux */
        if (sigma_s < 0.0)
          _stabilizing_flux(r, e) = -_scalar_flux(r,e) * _stabilization_factor
              * sigma_s / sigma_t[e];
      }
    }
  }
  else if (_stabilization_type == YAMAMOTO) {

    /* Treat each group */
#pragma omp parallel for schedule(static)
    for (int e=0; e < _NUM_GROUPS; e++) {

      /* Look for largest absolute scattering ratio */
      FP_PRECISION max_ratio = 0.0;
      for (long r=0; r < _num_FSRs; r++) {

        /* Extract the scattering value */
        FP_PRECISION scat = _FSR_materials[r]->getSigmaSByGroup(e+1, e+1);

        /* Extract total cross-sections */
        FP_PRECISION total = _FSR_materials[r]->getSigmaTByGroup(e+1);

        /* Determine scattering ratio */
        FP_PRECISION ratio = std::abs(scat / total);
        if (ratio > max_ratio)
          max_ratio = ratio;
      }
      max_ratio *= _stabilization_factor;
      for (long r=0; r < _num_FSRs; r++) {
        _stabilizing_flux(r, e) = _scalar_flux(r,e) * max_ratio;
      }
    }
  }
  else if (_stabilization_type == GLOBAL) {

    /* Get the multiplicative factor */
    FP_PRECISION mult_factor = 1.0 / _stabilization_factor - 1.0;

    /* Apply the global muliplicative factor */
#pragma omp parallel for schedule(static)
    for (long r=0; r < _num_FSRs; r++)
      for (int e=0; e < _NUM_GROUPS; e++)
        _stabilizing_flux(r, e) = mult_factor * _scalar_flux(r,e);
  }
}


/**
 * @brief Adjusts the scalar flux for transport stabilization
 */
void CPUSolver::stabilizeFlux() {

  if (_stabilization_type == DIAGONAL) {

    /* Loop over all flat source regions */
#pragma omp parallel for schedule(static)
    for (long r=0; r < _num_FSRs; r++) {

      /* Extract the scattering matrix */
      FP_PRECISION* scattering_matrix = _FSR_materials[r]->getSigmaS();

      /* Extract total cross-sections */
      FP_PRECISION* sigma_t = _FSR_materials[r]->getSigmaT();

      for (int e=0; e < _NUM_GROUPS; e++) {

        /* Extract the in-scattering (diagonal) element */
        FP_PRECISION sigma_s = scattering_matrix[e*_NUM_GROUPS+e];

        /* For negative cross-sections, add the stabilizing flux
           and divide by the diagonal matrix element used to form it so that
           no bias is introduced but the source iteration is stabilized */
        if (sigma_s < 0.0) {
          _scalar_flux(r, e) += _stabilizing_flux(r,e);
          _scalar_flux(r, e) /= (1.0 - _stabilization_factor * sigma_s /
                                 sigma_t[e]);
        }
      }
    }
  }
  else if (_stabilization_type == YAMAMOTO) {

    /* Treat each group */
#pragma omp parallel for schedule(static)
    for (int e=0; e < _NUM_GROUPS; e++) {

      /* Look for largest absolute scattering ratio */
      FP_PRECISION max_ratio = 0.0;
      for (long r=0; r < _num_FSRs; r++) {

        /* Extract the scattering value */
        FP_PRECISION scat = _FSR_materials[r]->getSigmaSByGroup(e+1, e+1);

        /* Extract total cross-sections */
        FP_PRECISION total = _FSR_materials[r]->getSigmaTByGroup(e+1);

        /* Determine scattering ratio */
        FP_PRECISION ratio = std::abs(scat / total);
        if (ratio > max_ratio)
          max_ratio = ratio;
      }
      max_ratio *= _stabilization_factor;
      for (long r=0; r < _num_FSRs; r++) {
        _scalar_flux(r, e) += _stabilizing_flux(r, e);
        _scalar_flux(r, e) /= (1 + max_ratio);
      }
    }
  }
  else if (_stabilization_type == GLOBAL) {

    /* Apply the damping factor */
#pragma omp parallel for schedule(static)
    for (long r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _NUM_GROUPS; e++) {
        _scalar_flux(r, e) += _stabilizing_flux(r, e);
        _scalar_flux(r, e) *= _stabilization_factor;
      }
    }
  }
}


/**
 * @brief Computes the volume-averaged, energy-integrated fission or nu-fission
 *        rate in each FSR and stores them in an array indexed by FSR ID.
 * @details This is a helper method for SWIG to allow users to retrieve
 *          FSR fission or nu-fission rates as a NumPy array. An example of
 *          how this method can be called from Python is as follows:
 *
 * @code
 *          num_FSRs = geometry.getNumFSRs()
 *          fission_rates = solver.computeFSRFissionRates(num_FSRs)
 * @endcode
 *
 * @param fission_rates an array to store the fission rates (implicitly
 *                      passed in as a NumPy array from Python)
 * @param num_FSRs the number of FSRs passed in from Python
 * @param nu whether return nu-fission (true) or fission (default, false) rates
 */
void CPUSolver::computeFSRFissionRates(double* fission_rates, long num_FSRs,
                                       bool nu) {

  if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to compute FSR fission rates since the "
               "source distribution has not been calculated");

  log_printf(INFO, "Computing FSR fission rates...");

  FP_PRECISION* sigma_f;
  FP_PRECISION vol;

  /* Initialize fission rates to zero */
  for (long r=0; r < _num_FSRs; r++)
    fission_rates[r] = 0.0;

  /* Loop over all FSRs and compute the volume-weighted fission rate */
#pragma omp parallel for private (sigma_f, vol) schedule(static)
  for (long r=0; r < _num_FSRs; r++) {
    if (nu) {
      sigma_f = _FSR_materials[r]->getNuSigmaF();
    }
    else {
      sigma_f = _FSR_materials[r]->getSigmaF();
    }
    vol = _FSR_volumes[r];

    for (int e=0; e < _NUM_GROUPS; e++)
      fission_rates[r] += sigma_f[e] * _scalar_flux(r,e) * vol;
  }

  /* Reduce domain data for domain decomposition */
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {

    /* Allocate buffer for communication */
    long num_total_FSRs = _geometry->getNumTotalFSRs();
    double* temp_fission_rates = new double[num_total_FSRs];
    for (int i=0; i < num_total_FSRs; i++)
      temp_fission_rates[i] = 0;

    int rank = 0;
    MPI_Comm comm = _geometry->getMPICart();
    MPI_Comm_rank(comm, &rank);
    for (long r=0; r < num_total_FSRs; r++) {

      /* Determine the domain and local FSR ID */
      long fsr_id = r;
      int domain = 0;
      _geometry->getLocalFSRId(r, fsr_id, domain);

      /* Set data if in the correct domain */
      if (domain == rank)
        temp_fission_rates[r] = fission_rates[fsr_id];
    }

    MPI_Allreduce(temp_fission_rates, fission_rates, num_total_FSRs,
                  MPI_DOUBLE, MPI_SUM, comm);
    delete [] temp_fission_rates;
  }
#endif
}


/**
 * @brief A function that prints a summary of the input parameters
 */
void CPUSolver::printInputParamsSummary() {

  /* Print general solver input parameters summary */
  Solver::printInputParamsSummary();

  /* Print threads used */
  log_printf(NORMAL, "Using %d threads", _num_threads);
}


#ifdef MPIx
/**
 * @brief A function that prints the repartition of integrations and tracks
 *        among domains and interfaces.
 */
void CPUSolver::printLoadBalancingReport() {

  /* Give a measure of the load imbalance for the sweep step (segments) */
  int num_ranks = 1;
  long num_segments = _track_generator->getNumSegments();
  long min_segments = num_segments, max_segments = num_segments,
       total_segments = num_segments;
  if (_geometry->isDomainDecomposed()) {
    MPI_Comm_size(_geometry->getMPICart(), &num_ranks);
    MPI_Reduce(&num_segments, &min_segments, 1, MPI_LONG, MPI_MIN, 0,
               _geometry->getMPICart());
    MPI_Reduce(&num_segments, &max_segments, 1, MPI_LONG, MPI_MAX, 0,
               _geometry->getMPICart());
    MPI_Reduce(&num_segments, &total_segments, 1, MPI_LONG, MPI_SUM, 0,
               _geometry->getMPICart());
  }
  FP_PRECISION mean_segments = float(total_segments) / num_ranks;
  log_printf(INFO_ONCE, "Min / max / mean number of segments in domains: "
             "%.1e / %.1e / %.1e", float(min_segments), float(max_segments),
             mean_segments);

  /* Give a measure of load imbalance for the communication phase */
  FP_PRECISION tracks_x = 0, tracks_y = 0, tracks_z = 0;
  int domain = _geometry->getNeighborDomain(0, 0, 1);
  if (domain != -1)
    tracks_z = _boundary_tracks.at(_neighbor_connections.at(domain)).size();

  domain = _geometry->getNeighborDomain(0, 1, 0);
  if (domain != -1)
    tracks_y = _boundary_tracks.at(_neighbor_connections.at(domain)).size();

  domain = _geometry->getNeighborDomain(1, 0, 0);
  if (domain != -1)
    tracks_x = _boundary_tracks.at(_neighbor_connections.at(domain)).size();

  long sum_border_tracks_200 = std::max(FP_PRECISION(1),
                                        tracks_x + tracks_y + tracks_z) / 100.;
  log_printf(INFO_ONCE, "Percentage of tracks exchanged in X/Y/Z direction: "
             "%.2f / %.2f / %.2f %", tracks_x / sum_border_tracks_200, tracks_y
             / sum_border_tracks_200, tracks_z / sum_border_tracks_200);
}
#endif


/**
 * @brief A function that prints the source region fluxes on a 2D mesh grid
 * @param dim1 coordinates of the mesh grid in the first direction
 * @param dim2 coordinates of the mesh grid in the second direction
 * @param offset The location of the mesh grid center on the perpendicular axis
 * @param plane 'xy', 'xz' or 'yz' the plane in which the mesh grid lies
 */
void CPUSolver::printFSRFluxes(std::vector<double> dim1,
                               std::vector<double> dim2,
                               double offset, const char* plane) {

  int rank = 0;
#ifdef MPIx
  MPI_Comm comm;
  if (_geometry->isDomainDecomposed()) {
    comm = _geometry->getMPICart();
    MPI_Comm_rank(comm, &rank);
  }

  MPI_Datatype precision;
  if (sizeof(FP_PRECISION) == 4)
    precision = MPI_FLOAT;
  else
    precision = MPI_DOUBLE;
#endif
  std::vector<long> fsr_ids = _geometry->getSpatialDataOnGrid(dim1, dim2,
                                                              offset, plane,
                                                              "fsr");
  std::vector<int> domain_contains_coords(fsr_ids.size());
  std::vector<int> num_contains_coords(fsr_ids.size());
#pragma omp parallel for
  for (long r=0; r < fsr_ids.size(); r++) {
    if (fsr_ids.at(r) != -1)
      domain_contains_coords.at(r) = 1;
    else
      domain_contains_coords.at(r) = 0;
  }

#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Allreduce(&domain_contains_coords[0], &num_contains_coords[0],
                  fsr_ids.size(), MPI_INT, MPI_SUM, comm);
#endif
  if (!_geometry->isDomainDecomposed())
    for (int i=0; i < fsr_ids.size(); i++)
      num_contains_coords[i] = domain_contains_coords[i];

  for (int e=0; e < _NUM_GROUPS; e++) {

    std::vector<FP_PRECISION> domain_fluxes(fsr_ids.size(), 0);
    std::vector<FP_PRECISION> total_fluxes(fsr_ids.size());

#pragma omp parallel for
    for (long r=0; r < fsr_ids.size(); r++) {
      if (domain_contains_coords.at(r) != 0)
        domain_fluxes.at(r) = getFlux(fsr_ids.at(r), e+1);
    }

#ifdef MPIx
    if (_geometry->isDomainDecomposed())
      MPI_Allreduce(&domain_fluxes[0], &total_fluxes[0],
                    fsr_ids.size(), precision, MPI_SUM, comm);
#endif
    if (!_geometry->isDomainDecomposed())
      for (int i=0; i < fsr_ids.size(); i++)
        total_fluxes[i] = domain_fluxes[i];

    if (rank == 0) {
      for (int i=0; i<dim1.size(); i++) {
        for (int j=0; j<dim2.size(); j++) {
          int r = i + j*dim1.size();
          double flux = total_fluxes.at(r) / num_contains_coords.at(r);
          log_printf(NORMAL, "(%d: %f, %d: %f) -> %f", i, dim1.at(i), j,
                     dim2.at(j), flux);
        }
      }
    }
  }
}


/**
 * @brief A function that prints fsr fluxes in xy plane at z=middle
 * @details Recommend deletion, since redundant printFluxes
 */
void CPUSolver::printFluxesTemp() {

  Universe* root = _geometry->getRootUniverse();

  int nx = 100;
  int ny = 100;
  int nz = 1;

  double x_min = root->getMinX() + 2*TINY_MOVE;
  double x_max = root->getMaxX() - 2*TINY_MOVE;
  double y_min = root->getMinY() + 2*TINY_MOVE;
  double y_max = root->getMaxY() - 2*TINY_MOVE;
  double z_min = root->getMinZ() + 2*TINY_MOVE;
  double z_max = root->getMaxZ() - 2*TINY_MOVE;

  std::vector<double> x(nx);
  std::vector<double> y(ny);
  for (int i=0; i < nx; i++)
    x.at(i) = x_min + i * (x_max - x_min) / nx;
  for (int j=0; j < ny; j++)
    y.at(j) = y_min + j * (y_max - y_min) / ny;

  double z_mid = (z_min + z_max) / 2 + TINY_MOVE;

  printFSRFluxes(x, y, z_mid, "xy");
}


/**
 * @brief A function that prints the number of FSRs with negative sources in
 *        the whole geometry subdivided by a 3D lattice. The number of negative
 *        sources per energy group is also printed out.
 * @param iteration the current iteration
 * @param num_x number of divisions in X direction
 * @param num_y number of divisions in Y direction
 * @param num_z number of divisions in Z direction
 */
void CPUSolver::printNegativeSources(int iteration, int num_x, int num_y,
                                     int num_z) {

  long long iter = iteration;
  std::string fname = "k_negative_sources_iter_";
  std::string iter_num = std::to_string(iter);
  fname += iter_num;
  std::ofstream out(fname);

  /* Create a lattice */
  Lattice lattice;
  lattice.setNumX(num_x);
  lattice.setNumY(num_y);
  lattice.setNumZ(num_z);

  /* Get the root universe */
  Universe* root_universe = _geometry->getRootUniverse();

  /* Determine the geometry widths in each direction */
  double width_x = (root_universe->getMaxX() - root_universe->getMinX())/num_x;
  double width_y = (root_universe->getMaxY() - root_universe->getMinY())/num_y;
  double width_z = (root_universe->getMaxZ() - root_universe->getMinZ())/num_z;

  /* Determine the center-point of the geometry */
  double offset_x = (root_universe->getMinX() + root_universe->getMaxX()) / 2;
  double offset_y = (root_universe->getMinY() + root_universe->getMaxY()) / 2;
  double offset_z = (root_universe->getMinZ() + root_universe->getMaxZ()) / 2;

  /* Create the Mesh lattice */
  lattice.setWidth(width_x, width_y, width_z);
  lattice.setOffset(offset_x, offset_y, offset_z);
  lattice.computeSizes();

  /* Create a group-wise negative source mapping */
  int by_group[_NUM_GROUPS];
  for (int e=0; e < _NUM_GROUPS; e++)
    by_group[e] = 0;

  int mapping[num_x*num_y*num_z];
  for (int i=0; i < num_x*num_y*num_z; i++)
    mapping[i] = 0;

  /* Loop over all flat source regions */
  for (long r=0; r < _num_FSRs; r++) {

    /* Determine the Mesh cell containing the FSR */
    Point* pt = _geometry->getFSRPoint(r);
    int lat_cell = lattice.getLatticeCell(pt);

    /* Determine the number of negative sources */
    for (int e=0; e < _NUM_GROUPS; e++) {
      if (_reduced_sources(r,e) < 10 * FLUX_EPSILON) {
        by_group[e]++;
        mapping[lat_cell]++;
      }
    }
  }

  /* If domain decomposed, do a reduction */
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    int size = num_x * num_y * num_z;
    int neg_src_send[size];
    for (int i=0; i < size; i++)
      neg_src_send[i] = mapping[i];
    MPI_Allreduce(neg_src_send, mapping, size, MPI_INT, MPI_SUM,
                  _geometry->getMPICart());

    int neg_src_grp_send[size];
    for (int e=0; e < _NUM_GROUPS; e++)
        neg_src_grp_send[e] = by_group[e];
    MPI_Allreduce(neg_src_grp_send, by_group, _NUM_GROUPS, MPI_INT, MPI_SUM,
                  _geometry->getMPICart());
  }
#endif


  /* Print negative source distribution to file */
  if (_geometry->isRootDomain()) {
    out << "[NORMAL]  Group-wise distribution of negative sources:"
        << std::endl;
    for (int e=0; e < _NUM_GROUPS; e++)
      out << "[NORMAL]  Group "  << e << ": " << by_group[e] << std::endl;
    out << "[NORMAL]  Spatial distribution of negative sources:" << std::endl;
    for (int z=0; z < num_z; z++) {
      out << " -------- z = " << z << " ----------" << std::endl;
      for (int y=0; y < num_y; y++) {
        for (int x=0; x < num_x; x++) {
          int ind = (z * num_y + y) * num_x + x;
          out << mapping[ind] << " ";
        }
        out << std::endl;
      }
    }
  }
  out.close();
}
