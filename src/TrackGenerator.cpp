#include "TrackGenerator.h"
#include "TrackTraversingAlgorithms.h"
#include <iomanip>

/**
 * @brief Constructor for the TrackGenerator assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param azim_spacing azimuthal track spacing (cm)
 */
TrackGenerator::TrackGenerator(Geometry* geometry, int num_azim,
                               double azim_spacing) {

  setGeometry(geometry);
  setNumThreads(1);
  setNumAzim(num_azim);
  setDesiredAzimSpacing(azim_spacing);
  _contains_2D_tracks = false;
  _contains_2D_segments = false;
  _quadrature = NULL;
  _z_coord = 0;
  _segment_formation = EXPLICIT_2D;
  _max_optical_length = std::numeric_limits<FP_PRECISION>::max();
  _max_num_segments = 0;
  _FSR_volumes = NULL;
  _dump_segments = true;
  _segments_centered = false;
  _FSR_locks = NULL;
  _tracks_2D_array = NULL;
  _tracks_per_azim = NULL;
  _timer = new Timer();
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator::~TrackGenerator() {

  if (_contains_2D_tracks) {

    /* Delete 2D tracks */
    for (int a=0; a < _num_azim/2; a++)
      delete [] _tracks_2D[a];
    delete [] _tracks_2D;

    delete[] _tracks_2D_array;

    /* Delete track laydown information */
    delete [] _num_x;
    delete [] _num_y;
  }

  if (_tracks_per_azim != NULL)
    delete [] _tracks_per_azim;

  if (_FSR_locks != NULL)
    delete [] _FSR_locks;

  if (_FSR_volumes != NULL)
    delete [] _FSR_volumes;

  delete _quadrature;
  delete _timer;
}


/**
 * @brief Return the number of tracks for each azimuthal angle.
 * @return the number of tracks for each azimuthal angle
 */
long* TrackGenerator::getTracksPerAzim() {
  return _tracks_per_azim;
}


/**
 * @brief Return the number of azimuthal angles in \f$ [0, 2\pi] \f$.
 * @return the number of azimuthal angles in \f$ 2\pi \f$
 */
int TrackGenerator::getNumAzim() {
  return _num_azim;
}


/**
 * @brief Return the track azimuthal spacing (cm).
 * @details This will return the user-specified track spacing and NOT the
 *          effective track spacing which is computed and used to generate
 *          cyclic tracks.
 * @return the track azimuthal spacing (cm)
 */
double TrackGenerator::getDesiredAzimSpacing() {
  return _azim_spacing;
}


/**
 * @brief Return the Geometry for this TrackGenerator if one has been set.
 * @return a pointer to the Geometry
 */
Geometry* TrackGenerator::getGeometry() {
  if (_geometry == NULL)
    log_printf(ERROR, "Unable to return the TrackGenerator's Geometry "
               "since it has not yet been set");

  return _geometry;
}


/**
 * @brief Return the array of FSR locks for atomic FSR operations.
 * @return an array of FSR locks
 */
omp_lock_t* TrackGenerator::getFSRLocks() {
  if (_FSR_locks == NULL)
    log_printf(ERROR, "Unable to return the TrackGenerator's FSR locks "
               "since they have not yet been created");

  return _FSR_locks;
}


/**
 * @brief Computes the volumes/areas of each Cell and Material in the Geometry.
 * @details This computes the volumes/areas of each Cell and Material in the
 *          Geometry from the length of track segments crossing each Cell. This
 *          routine assigns the volume and number of instances to each Cell and
 *          Material. This is a helper routine that is called after track
 *          segmentation is complete in TrackGenerator::generateTracks().
 */
void TrackGenerator::initializeVolumes() {

  if (!containsTracks())
    log_printf(ERROR, "Unable to initialize volumes since tracks "
               "have not yet been generated");

  Cell* cell;
  Material* material;
  int num_FSRs = _geometry->getNumFSRs();

  /* Initialize a volume buffer and compute volumes */
  initializeFSRVolumesBuffer();
  FP_PRECISION* fsr_volumes = getFSRVolumes();

  /* Reset cell and material volumes */
  for (int i=0; i < num_FSRs; i++) {
    cell = _geometry->findCellContainingFSR(i);
    cell->setVolume(0);

    material = cell->getFillMaterial();
    material->setVolume(0);
  }

  /* Compute volume and number of instances for each Cell and Material */
  for (int i=0; i < num_FSRs; i++) {
    cell = _geometry->findCellContainingFSR(i);
    cell->incrementVolume(fsr_volumes[i]);
    cell->incrementNumInstances();

    material = cell->getFillMaterial();
    material->incrementVolume(fsr_volumes[i]);
    material->incrementNumInstances();
  }
}


/**
 * @brief Initialize an array to contain the FSR volumes.
 */
void TrackGenerator::initializeFSRVolumesBuffer() {

  if (_FSR_volumes != NULL)
    delete[] _FSR_volumes;

#pragma omp critical
  {
    long num_FSRs = _geometry->getNumFSRs();
    _FSR_volumes = new FP_PRECISION[num_FSRs]();
  }
}


/**
 * @brief Return the array used to store the FSR volumes.
 * @return _FSR_volumes the FSR volumes array indexed by FSR ID
 */
FP_PRECISION* TrackGenerator::getFSRVolumesBuffer() {

  return _FSR_volumes;
}


/**
 * @brief Return the total number of Tracks across the Geometry.
 * @return the total number of Tracks
 */
long TrackGenerator::getNumTracks() {
  return getNum2DTracks();
}


/**
 * @brief Return the total number of 2D Tracks across the Geometry.
 * @return the total number of 2D Tracks
 */
long TrackGenerator::getNum2DTracks() {

  long num_2D_tracks = 0;

  for (int a=0; a < _num_azim/2; a++)
    num_2D_tracks += _num_x[a] + _num_y[a];

  return num_2D_tracks;
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
long TrackGenerator::getNumSegments() {
  return getNum2DSegments();
}


/**
 * @brief Return the total number of 2D Track segments across the Geometry.
 * @return the total number of 2D Track segments
 */
long TrackGenerator::getNum2DSegments() {

  if (!containsSegments())
    log_printf(ERROR, "Cannot get the number of 2D segments since they "
               "have not been generated.");

  long num_2D_segments = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      num_2D_segments += _tracks_2D[a][i].getNumSegments();
    }
  }

  return num_2D_segments;
}


/**
 * @brief Returns an array of the Track pointers by increasing UID
 * @details An array of pointers to all 2D Track objects in the Geometry is
 *          returned, arranged by increasing unique identifier (UID).
 * @return the array of Track pointers
 */
Track** TrackGenerator::get2DTracksArray() {

  if (!TrackGenerator::containsTracks())
    log_printf(ERROR, "Unable to return the 1D array of Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_2D_array;
}


/**
 * @brief Returns an array of the Track pointers by increasing UID.
 * @details Calls TrackGenerator::get2DTracksArray to return all 2D Tracks
 * @return the array of Track pointers
 */
Track** TrackGenerator::getTracksArray() {
  return get2DTracksArray();
}


/**
 * @brief Returns a 2D jagged array of the 2D Tracks.
 * @details The first index into the array is the azimuthal angle and the
 *          second index is the Track number.
 * @return the 2D jagged array of 2D Tracks
 */
Track** TrackGenerator::get2DTracks() {

  if (!TrackGenerator::containsTracks())
    log_printf(ERROR, "Unable to return the 3D ragged array of the 2D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_2D;
}


/**
 * @brief Calculates and returns the maximum optical length for any segment
 *        in the Geometry.
 * @details The _max_optical_length value is recomputed, updated, and returned.
 *          This value determines when segments must be split during ray
 *          tracing.
 * @return the maximum optical length of any segment in the Geometry
 */
FP_PRECISION TrackGenerator::getMaxOpticalLength() {
  MaxOpticalLength update_max_optical_length(this);
  update_max_optical_length.execute();
  return _max_optical_length;
}


/**
 * @brief Returns the maximum number of segments along a single track.
 * @details The TrackGenerator::countSegments routine must be called before
 *          this function will return a correct value
 * @return the maximum number of segments on a single track
 */
int TrackGenerator::getMaxNumSegments() {
  return _max_num_segments;
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int TrackGenerator::getNumThreads() {
  return _num_threads;
}


/**
 * @brief Returns the number of 2D Tracks in the x-direction for a given
 *        azimuthal angle index.
 * @param azim the azimuthal angle index
 * @return the number of 2D Tracks in the x-direction of the Geometry
 */
int TrackGenerator::getNumX(int azim) {
  return _num_x[azim];
}


/**
 * @brief Returns the number of 2D Tracks in the y-direction for a given
 *        azimuthal angle index.
 * @param azim the azimuthal angle index
 * @return the number of 2D Tracks in the y-direction of the Geometry
 */
int TrackGenerator::getNumY(int azim) {
  return _num_y[azim];
}


/**
 * @brief FSR volumes are copied to an array input by the user.
 * @param out_volumes The array to which FSR volumes are copied
 * @param num_fsrs The number of FSR volumes to copy. The first num_fsrs
 *        volumes stored in the FSR volumes array are copied.
 */
void TrackGenerator::exportFSRVolumes(double* out_volumes, int num_fsrs) {

  for (int i=0; i < num_fsrs; i++)
    out_volumes[i] = _FSR_volumes[i];

}


/**
 * @brief Computes and returns an array of volumes indexed by FSR.
 * @details Note: The memory is stored in the FSR volumes buffer of the
 *          TrackGenerator and is freed during deconstruction.
 * @return a pointer to the array of FSR volumes
 */
FP_PRECISION* TrackGenerator::getFSRVolumes() {

  /* Reset FSR volumes to zero */
  long num_FSRs = _geometry->getNumFSRs();
  if (_FSR_volumes != NULL)
    memset(_FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));

  /* Create volume calculator and calculate new FSR volumes */
  VolumeCalculator volume_calculator(this);
  volume_calculator.execute();

  /* Check to ensure all FSRs are crossed by at least one track */
  long num_zero_volume_fsrs = 0;
  for (long i=0; i < num_FSRs; i++) {
    if (fabs(_FSR_volumes[i]) < FLT_EPSILON) {
      num_zero_volume_fsrs++;
      log_printf(WARNING, "Zero volume calculated for FSR %d, point (%f, %f, %f)",
                 i, _geometry->getFSRPoint(i)->getX(),
                 _geometry->getFSRPoint(i)->getY(),
                 _geometry->getFSRPoint(i)->getZ());
    }
  }
  if (num_zero_volume_fsrs > 0)
    log_printf(NODAL, "Zero volume calculated in %ld FSR regions since "
               "no track traversed the FSRs. Use a finer track laydown to "
               "ensure every FSR is traversed.", num_zero_volume_fsrs);

  /* Print out total domains' volumes for debugging purposes */
  double total_volume = 0;
  for (long r=0; r < num_FSRs; r++)
    total_volume += _FSR_volumes[r];
  log_printf(DEBUG, "Total volume %f cm3", total_volume);

  return _FSR_volumes;
}


/**
 * @brief Returns the volume of an FSR.
 * @param fsr_id the ID for the FSR of interest
 * @return the FSR volume
 */
FP_PRECISION TrackGenerator::getFSRVolume(long fsr_id) {

  if (_FSR_volumes == NULL)
    log_printf(ERROR, "Unable to get the FSR volume since FSR volumes "
               "have not yet been generated");

  else if (fsr_id < 0 || fsr_id > _geometry->getNumFSRs())
    log_printf(ERROR, "Unable to get the volume for FSR %d since the FSR IDs "
               "lie in the range (0, %d)", fsr_id, _geometry->getNumFSRs());
  return _FSR_volumes[fsr_id];
}


/**
 * @brief Returns the z-coord of the radial plane used in 2D calculations.
 * @return the z-coord of the 2D calculation
 */
double TrackGenerator::getZCoord() {
  return _z_coord;
}


/**
 * @brief Returns the Quadrature object.
 * @return the Quadrature object
 */
Quadrature* TrackGenerator::getQuadrature() {
  return _quadrature;
}


/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
 * @param num_threads the number of threads
 */
void TrackGenerator::setNumThreads(int num_threads) {

  if (num_threads <= 0)
    log_printf(ERROR, "Unable to set the number of threads for the "
               "TrackGenerator to %d since it is less than or equal to 0"
               , num_threads);

#ifdef MPIx
  /* Check that the MPI library has enough thread support */
  int provided;
  MPI_Query_thread(&provided);
  if (num_threads > 1 && provided < MPI_THREAD_SERIALIZED)
    log_printf(WARNING, "Not enough thread support level in the MPI library, "
               "re-compile with another library. Thread support level should "
               "be at least MPI_THREAD_SERIALIZED.");
#endif

  _num_threads = num_threads;

  /* Set the number of threads for OpenMP */
  omp_set_num_threads(_num_threads);
  if (_geometry != NULL)
    _geometry->setNumThreads(num_threads);

  /* Print CPU assignments, useful for NUMA where by-socket is the preferred
   * CPU grouping */
  std::vector<int> cpus;
  cpus.reserve(num_threads);

#pragma omp parallel for schedule(static)
  for (int i=0; i<num_threads; i++) {
#pragma omp critical
    {
      cpus.push_back(sched_getcpu());
    }
  }
  sort(cpus.begin(), cpus.end());

  std::stringstream str_cpus;
  for (int i=0; i<cpus.size(); i++)
      str_cpus << cpus.at(i) << " ";

  /* Get rank of process */
  int rank = 0;
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Comm_rank(_geometry->getMPICart(), &rank);
#endif

  if (num_threads > 1)
    log_printf(NODAL, "CPUs on rank %d process: %s", rank,
               str_cpus.str().c_str());
}


/**
 * @brief Set the number of azimuthal angles in \f$ [0, 2\pi] \f$.
 * @param num_azim the number of azimuthal angles in \f$ 2\pi \f$
 */
void TrackGenerator::setNumAzim(int num_azim) {

  if (num_azim < 0)
    log_printf(ERROR, "Unable to set a negative number of azimuthal angles "
               "%d for the TrackGenerator.", num_azim);

  if (num_azim % 4 != 0)
    log_printf(ERROR, "Unable to set the number of azimuthal angles to %d for "
               "the TrackGenerator since it is not a multiple of 4", num_azim);

  _num_azim = num_azim;
  resetStatus();
}


/**
 * @brief Set the suggested azimuthal track spacing (cm).
 * @param spacing the suggested track azimuthal spacing
 */
void TrackGenerator::setDesiredAzimSpacing(double spacing) {
  if (spacing < 0)
    log_printf(ERROR, "Unable to set a negative track azimuthal spacing "
               "%f for the TrackGenerator.", spacing);

  _azim_spacing = spacing;
  resetStatus();
}


/**
 * @brief Set a pointer to the Geometry to use for track generation.
 * @param geometry a pointer to the Geometry
 */
void TrackGenerator::setGeometry(Geometry* geometry) {
  _geometry = geometry;
  _x_min = geometry->getMinX();
  _y_min = geometry->getMinY();
  _z_min = geometry->getMinZ();
  _x_max = geometry->getMaxX();
  _y_max = geometry->getMaxY();
  _z_max = geometry->getMaxZ();
  resetStatus();
}


/**
 * @brief Sets the z-coord of the radial plane used in 2D calculations.
 * @param z_coord the z-coord of the radial plane
 */
void TrackGenerator::setZCoord(double z_coord) {
  _z_coord = z_coord;

  /* Move the CMFD lattice near the plane of interest */
  if (_geometry->getCmfd() != NULL)
    _geometry->getCmfd()->getLattice()->getOffset()->setZ(z_coord - 0.5);
}


/**
 * @brief Sets the Quadrature used for integrating the MOC equations.
 * @param quadrature a pointer to the Quadrature object used in calculation
 */
void TrackGenerator::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
}


/**
 * @brief Returns whether or not the TrackGenerator contains Tracks
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains Tracks; false otherwise
 */
bool TrackGenerator::containsTracks() {
  return _contains_2D_tracks;
}


/**
 * @brief Returns whether or not the TrackGenerator contains segments
 *        for its current number of azimuthal angles, track spacing and
 *        geometry for its current segmentation type.
 * @return true if the TrackGenerator contains segments; false otherwise
 */
bool TrackGenerator::containsSegments() {
  return _contains_2D_segments;
}


/**
 * @brief Fills an array with the x,y,z coordinates for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords(num_tracks*4)
 * @endcode
 *
 * @param coords an array of coords of length 6 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator::retrieveTrackCoords(double* coords, long num_tracks) {
  retrieve2DTrackCoords(coords, num_tracks);
}


/**
 * @brief Fills an array with the x,y,z coordinates for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNum2DTracks()
 *          coords = track_generator.retrieve2DTrackCoords(num_tracks*4)
 * @endcode
 *
 * @param coords an array of coords of length 6 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator::retrieve2DTrackCoords(double* coords, long num_tracks) {

  if (num_tracks != NUM_VALUES_PER_RETRIEVED_TRACK * getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the Track coordinates since the "
               "TrackGenerator contains %d Tracks with %d coordinates but an "
               "array of length %d was input",
               getNum2DTracks(), NUM_VALUES_PER_RETRIEVED_TRACK *
               getNum2DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      coords[counter]   = _tracks_2D[a][i].getStart()->getX();
      coords[counter+1] = _tracks_2D[a][i].getStart()->getY();
      coords[counter+2] = _tracks_2D[a][i].getStart()->getZ();
      coords[counter+3] = _tracks_2D[a][i].getEnd()->getX();
      coords[counter+4] = _tracks_2D[a][i].getEnd()->getY();
      coords[counter+5] = _tracks_2D[a][i].getEnd()->getZ();
      //NOTE The Z-coordinate is constant

      counter += NUM_VALUES_PER_RETRIEVED_TRACK;
    }
  }
}


/**
 * @brief Fills an array with the x,y,z coordinates for each Track segment.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          segments. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_segments = track_generator.getNumSegments()
 *          coords = track_generator.retrieveSegmentCoords(num_segments*5)
 * @endcode
 *
 * @param coords an array of coords of length 7 times the number of segments
 * @param num_segments the total number of Track segments
 */
void TrackGenerator::retrieveSegmentCoords(double* coords, long num_segments) {
  retrieve2DSegmentCoords(coords, num_segments);
}


/**
 * @brief Fills an array with the x,y,z coordinates for each Track segment.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          segments. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_segments = track_generator.getNum2DSegments()
 *          coords = track_generator.retrieve2DSegmentCoords(num_segments*5)
 * @endcode
 *
 * @param coords an array of coords of length 7 times the number of segments
 * @param num_segments the total number of Track segments
 */
void TrackGenerator::retrieve2DSegmentCoords(double* coords, long num_segments) {

  if (num_segments != NUM_VALUES_PER_RETRIEVED_SEGMENT * getNum2DSegments())
    log_printf(ERROR, "Unable to retrieve the Track segment coordinates since "
               "the TrackGenerator contains %d segments with %d coordinates "
               "but an array of length %d was input",
               getNum2DSegments(), NUM_VALUES_PER_RETRIEVED_SEGMENT *
               getNum2DSegments(), num_segments);

  segment* curr_segment = NULL;
  double x0, x1, y0, y1, z0, z1;
  double phi;
  segment* segments;

  int counter = 0;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      x0    = _tracks_2D[a][i].getStart()->getX();
      y0    = _tracks_2D[a][i].getStart()->getY();
      z0    = _tracks_2D[a][i].getStart()->getZ();
      phi   = _tracks_2D[a][i].getPhi();

      segments = _tracks_2D[a][i].getSegments();

      for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++) {
        curr_segment = &segments[s];

        coords[counter] = curr_segment->_region_id;

        coords[counter+1] = x0;
        coords[counter+2] = y0;
        coords[counter+3] = z0;

        x1 = x0 + cos(phi) * curr_segment->_length;
        y1 = y0 + sin(phi) * curr_segment->_length;
        z1 = z0;

        coords[counter+4] = x1;
        coords[counter+5] = y1;
        coords[counter+6] = z1;

        x0 = x1;
        y0 = y1;
        z0 = z1;

        counter += NUM_VALUES_PER_RETRIEVED_SEGMENT;
      }
    }
  }
}


/**
 * @brief Checks the boundary conditions for all 2D surfaces for inconsistent
 *        periodic boundary conditions.
 */
void TrackGenerator::checkBoundaryConditions() {

  /* Extract the X and Y boundaries for whole Geometry */
  Universe* root_universe = _geometry->getRootUniverse();
  boundaryType min_x_bound = root_universe->getMinXBoundaryType();
  boundaryType max_x_bound = root_universe->getMaxXBoundaryType();
  boundaryType min_y_bound = root_universe->getMinYBoundaryType();
  boundaryType max_y_bound = root_universe->getMaxYBoundaryType();

  /* Check X and Y boundaries for consistency */
  if ((min_x_bound == PERIODIC && max_x_bound != PERIODIC) ||
      (min_x_bound != PERIODIC && max_x_bound == PERIODIC))
    log_printf(ERROR, "Cannot create tracks with only one x boundary"
               " set to PERIODIC");

  else if ((min_y_bound == PERIODIC && max_y_bound != PERIODIC) ||
           (min_y_bound != PERIODIC && max_y_bound == PERIODIC))
    log_printf(ERROR, "Cannot create tracks with only one y boundary"
               " set to PERIODIC");

  /* Check that there are no periodic boundaries if domain decomposed */
  if (_geometry->isDomainDecomposed())
    if (min_x_bound == PERIODIC || min_y_bound == PERIODIC)
      log_printf(ERROR, "Periodic boundaries are not supported for domain "
                 "decomposition");

  /* Check for correct track method if a PERIODIC bc is present */
  if (_geometry->isDomainDecomposed() || min_x_bound == PERIODIC ||
      min_y_bound == PERIODIC)
    _periodic = true;
  else
    _periodic = false;
}


/**
 * @brief Generates tracks for some number of azimuthal angles and track spacing
 * @details Computes the effective angles and track spacing. Computes the
 *          number of Tracks for each azimuthal angle, allocates memory for
 *          all Tracks at each angle and sets each Track's starting and ending
 *          Points, azimuthal angle, and azimuthal angle quadrature weight.
 */
void TrackGenerator::generateTracks() {

  /* Start recording track generation time */
  _timer->startTimer();

  /* Check for valid quadrature */
  if (_quadrature != NULL) {
    if (_quadrature->getNumAzimAngles() != _num_azim) {
      delete _quadrature;
      _quadrature = NULL;
    }
  }

  /* Generate Tracks, perform ray tracing across the geometry, and store
   * the data to a Track file */
  try {

    /* Create default quadrature set if user one has not been set */
    if (_quadrature == NULL)
      initializeDefaultQuadrature();

    /* Initialize the quadrature set */
    _quadrature->initialize();

    /* Check periodic BCs for symmetry */
    checkBoundaryConditions();

    /* Lay down Tracks accross the Geometry */
    if (_geometry == NULL)
      log_printf(ERROR, "Unable to lay down Tracks since no Geometry "
                 "has been set for the TrackGenerator");

    /* Initialize the Tracks */
    initializeTracks();

    /* Initialize the track file directory and read in tracks if they exist */
    //NOTE Useful for 2D simulations, currently broken
    //initializeTrackFileDirectory();

    /* If track file not present, generate segments */
    if (_use_input_file == false) {

      /* Segmentize the tracks */
      segmentize();
      //FIXME HERE dumpSegmentsToFile();
    }

    /* Allocate array of mutex locks for each FSR */
    long num_FSRs = _geometry->getNumFSRs();
    _FSR_locks = new omp_lock_t[num_FSRs];

    /* Print number of FSRs as soon as it's available */
    long total_num_FSRs = num_FSRs;
#ifdef MPIx
    if (_geometry->isDomainDecomposed()) {
      MPI_Comm MPI_cart = _geometry->getMPICart();
      MPI_Allreduce(&num_FSRs, &total_num_FSRs, 1, MPI_LONG,
                    MPI_SUM, MPI_cart);
    }
#endif
    log_printf(NORMAL, "Total number of FSRs %ld", total_num_FSRs);
    log_printf(DEBUG, "Number of FSRs in domain %ld", num_FSRs);

    /* Loop over all FSRs to initialize OpenMP locks */
#pragma omp parallel for schedule(guided)
    for (long r=0; r < num_FSRs; r++)
      omp_init_lock(&_FSR_locks[r]);

    /* Precompute the quadrature weights */
    _quadrature->precomputeWeights(_segment_formation != EXPLICIT_2D);

  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to allocate memory needed to generate "
               "Tracks. Backtrace:\n%s", e.what());
  }

  /* Stop recording track generation time and print */
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Barrier(_geometry->getMPICart());
#endif
  _timer->stopTimer();
  _timer->recordSplit("Track Generation Time");
}


/**
 * @brief Allocates a new Quadrature with the default Quadrature.
 * @details The default quadrature for 2D calculations is the TY quadrature
 */
void TrackGenerator::initializeDefaultQuadrature() {

  if (_quadrature != NULL)
    delete _quadrature;

  log_printf(NORMAL, "Initializing a default angular quadrature...");
  _quadrature = new TYPolarQuad();
  _quadrature->setNumAzimAngles(_num_azim);
  _quadrature->setNumPolarAngles(6);
}


/**
 * @brief Calculates the least common multiple of two numbers a and b
 * @param a first number
 * @param b second number (order does not matter)
 * @return the least common multiple of a and b
 */
double TrackGenerator::leastCommonMultiple(double a, double b) {

  bool _found = false;
  int lcm_a = 1;
  int lcm_b;
  double residual;

  /* For efficiency, make a the longer length */
  if (a < b) {
    double a_temp = a;
    a = b;
    b = a_temp;
  }

  while (!_found) {

    lcm_b = (int) round((lcm_a * a) / b);
    residual = fabs(lcm_a * a - lcm_b * b);

    if (residual < LCM_TOLERANCE)
      _found = true;
    else
      lcm_a++;
  }

  return lcm_a * a;
}


/**
 * @brief Returns the type of ray tracing used for segment formation.
 * @return the segmentation type
 */
segmentationType TrackGenerator::getSegmentFormation() {
  return _segment_formation;
}


/**
 * @brief Initializes Track azimuthal angles, start and end Points.
 * @details This method computes the azimuthal angles and effective track
 *          spacing to use to guarantee cyclic Track wrapping. Based on the
 *          angles and spacing, the number of Tracks per angle and the start
 *          and end Points for each Track are computed.
 */
void TrackGenerator::initializeTracks() {

  /* Make sure that the width and height of the Geometry are nonzero */
  if (_geometry->getWidthX() <= 0 || _geometry->getWidthY() <= 0)
    log_printf(ERROR, "The total height and width of the Geometry must"
               " be nonzero for Track generation. Create a CellFill which "
               "is filled by the entire geometry and bounded by XPlanes "
               "and YPlanes to enable the Geometry to determine the "
               "total width and height of the model.");

  log_printf(NORMAL, "Initializing 2D tracks...");

  /* Allocate memory for arrays */
  _tracks_2D        = new Track*[_num_azim/2];
  _num_x            = new int[_num_azim/2];
  _num_y            = new int[_num_azim/2];
  _tracks_per_azim  = new long[_num_azim/2];
  _num_2D_tracks    = 0;

  double x1, x2, y1, y2;
  double phi;
  double width  = _geometry->getWidthX();
  double height = _geometry->getWidthY();
  double dx_eff[_num_azim/2];
  double dy_eff[_num_azim/2];
  double x_min = _geometry->getMinX();
  double y_min = _geometry->getMinY();

  /* Determine angular quadrature and track spacing */
  for (int a = 0; a < _num_azim/4; a++) {

    /* Get the desired azimuthal angle */
    phi = _quadrature->getPhi(a);

    /* The number of intersections with x,y-axes */
    double module_width = width / _geometry->getNumXModules();
    double module_height = height / _geometry->getNumYModules();

    /* The number of intersections with x,y-axes */
    _num_x[a] = (int) (fabs(module_width / _azim_spacing * sin(phi))) + 1;
    _num_y[a] = (int) (fabs(module_height / _azim_spacing * cos(phi))) + 1;

    /* Extend the number of intersections with x,y axes for modules */
    _num_x[a] *= _geometry->getNumXModules();
    _num_y[a] *= _geometry->getNumYModules();

    /* Save number of intersections for supplementary angles */
    _num_x[_num_azim/2 - a - 1] = _num_x[a];
    _num_y[_num_azim/2 - a - 1] = _num_y[a];

    /* Effective/actual angle (not the angle we desire, but close) */
    double phi = atan((height * _num_x[a]) / (width * _num_y[a]));
    _quadrature->setPhi(phi, a);

    /* Effective Track spacing (not spacing we desire, but close) */
    dx_eff[a]   = (width / _num_x[a]);
    dy_eff[a]   = (height / _num_y[a]);
    double azim_spacing = dx_eff[a] * sin(phi);
    _quadrature->setAzimSpacing(azim_spacing, a);

    /* Save spacings for supplementary angles */
    dx_eff[_num_azim/2 - a - 1] = dx_eff[a];
    dy_eff[_num_azim/2 - a - 1] = dy_eff[a];
  }

  /* Generate 2D tracks */
  for (int a=0; a < _num_azim/2; a++) {

    /* Allocate memory for the 2D tracks array */
    _tracks_2D[a] = new Track[_num_x[a] + _num_y[a]];
    _num_2D_tracks += _num_x[a] + _num_y[a];

    /* Get the azimuthal angle for all tracks with this azimuthal angle */
    phi = _quadrature->getPhi(a);

    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      /* Get track and set angle and track indices */
      Track* track = &_tracks_2D[a][i];
      track->setPhi(phi);
      track->setAzimIndex(a);
      track->setXYIndex(i);

      /* Set start point */
      if (a < _num_azim/4) {
        if (i < _num_x[a])
          track->getStart()->setCoords(x_min + width - dx_eff[a] * (i + 0.5),
                                       y_min);
        else
          track->getStart()->setCoords(x_min, y_min + dy_eff[a] *
                                       (i - _num_x[a] + 0.5));
      }
      else {
        if (i < _num_x[a])
          track->getStart()->setCoords(x_min + dx_eff[a] * (i + 0.5), y_min);
        else
          track->getStart()->setCoords(x_min + width, y_min + dy_eff[a] *
                                       (i-_num_x[a] + 0.5));
      }

      /* Set end point */
      if (a < _num_azim/4) {
        if (i < _num_y[a])
          track->getEnd()->setCoords(x_min + width, y_min + dy_eff[a] *
                                     (i + 0.5));
        else
          track->getEnd()->setCoords(x_min + width - dx_eff[a] * ((i-_num_y[a])
                                     + 0.5), y_min + height);
      }
      else {
        if (i < _num_y[a])
          track->getEnd()->setCoords(x_min, y_min + dy_eff[a] * (i + 0.5));
        else
          track->getEnd()->setCoords(x_min + dx_eff[a] * (i-_num_y[a] + 0.5),
                                     y_min + height);
      }
    }
  }

  /* Set the flag indicating 2D tracks have been generated */
  _contains_2D_tracks = true;

  /* Initialize the track reflections */
  TrackGenerator::initializeTrackReflections();

  /* Initialize the 1D array of Tracks for all Tracks */
  initializeTracksArray();
}


/**
 * @brief Initializes 2D Track reflections.
 * @details This method computes the connecting Tracks for all 2D Tracks in
 *          the TrackGenerator analytically, handling both reflective and
 *          periodic boundaries.
 */
void TrackGenerator::initializeTrackReflections() {

  log_printf(NORMAL, "Initializing 2D tracks reflections...");

  /* Generate the 2D track cycles */
  for (int a=0; a < _num_azim/2; a++) {
    int ac = _num_azim/2 - a - 1;
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      /* Get current track */
      Track* track = &_tracks_2D[a][i];

      /* Set the forward boundary conditions */
      if (a < _num_azim/4) {
        if (i < _num_y[a]) {
          track->setBCFwd(_geometry->getMaxXBoundaryType());
          track->setSurfaceOut(SURFACE_X_MAX);
#ifdef MPIx
          track->setDomainFwd(_geometry->getNeighborDomain(1, 0, 0));
#endif
        }
        else {
          track->setBCFwd(_geometry->getMaxYBoundaryType());
          track->setSurfaceOut(SURFACE_Y_MAX);
#ifdef MPIx
          track->setDomainFwd(_geometry->getNeighborDomain(0, 1, 0));
#endif
        }

        if (i < _num_x[a]) {
          track->setBCBwd(_geometry->getMinYBoundaryType());
          track->setSurfaceIn(SURFACE_Y_MIN);
#ifdef MPIx
          track->setDomainBwd(_geometry->getNeighborDomain(0, -1, 0));
#endif

        }
        else {
          track->setBCBwd(_geometry->getMinXBoundaryType());
          track->setSurfaceIn(SURFACE_X_MIN);
#ifdef MPIx
          track->setDomainBwd(_geometry->getNeighborDomain(-1, 0, 0));
#endif
        }
      }

      /* Set the backward boundary conditions */
      else {
        if (i < _num_y[a]) {
          track->setBCFwd(_geometry->getMinXBoundaryType());
          track->setSurfaceOut(SURFACE_X_MIN);
#ifdef MPIx
          track->setDomainFwd(_geometry->getNeighborDomain(-1, 0, 0));
#endif
        }
        else {
          track->setBCFwd(_geometry->getMaxYBoundaryType());
          track->setSurfaceOut(SURFACE_Y_MAX);
#ifdef MPIx
          track->setDomainFwd(_geometry->getNeighborDomain(0, 1, 0));
#endif
        }

        if (i < _num_x[a]) {
          track->setBCBwd(_geometry->getMinYBoundaryType());
          track->setSurfaceIn(SURFACE_Y_MIN);
#ifdef MPIx
          track->setDomainBwd(_geometry->getNeighborDomain(0, -1, 0));
#endif
        }
        else {
          track->setBCBwd(_geometry->getMaxXBoundaryType());
          track->setSurfaceIn(SURFACE_X_MAX);
#ifdef MPIx
          track->setDomainBwd(_geometry->getNeighborDomain(1, 0, 0));
#endif
        }
      }

      /* Set connecting tracks in forward direction */
      if (i < _num_y[a]) {
        track->setNextFwdFwd(true);
        track->setTrackPrdcFwd(get2DTrackID(a, i + _num_x[a]));
        track->setTrackReflFwd(get2DTrackID(ac, i + _num_x[a]));

        if (track->getBCFwd() == PERIODIC || track->getBCFwd() == INTERFACE)
          track->setTrackNextFwd(get2DTrackID(a, i + _num_x[a]));
        else
          track->setTrackNextFwd(get2DTrackID(ac, i + _num_x[a]));
      }
      else {
        track->setTrackPrdcFwd(get2DTrackID(a, i - _num_y[a]));
        track->setTrackReflFwd
            (get2DTrackID(ac, (_num_x[a] + _num_y[a]) - (i - _num_y[a]) - 1));

        if (_geometry->getMaxYBoundaryType() == PERIODIC ||
            _geometry->getMaxYBoundaryType() == INTERFACE) {
          track->setNextFwdFwd(true);
          track->setTrackNextFwd(get2DTrackID(a, i - _num_y[a]));
        }
        else {
          track->setNextFwdFwd(false);
          track->setTrackNextFwd
              (get2DTrackID(ac, (_num_x[a] + _num_y[a]) - (i - _num_y[a]) - 1));
        }
      }

      /* Set connecting tracks in backward direction */
      if (i < _num_x[a]) {
        track->setTrackPrdcBwd(get2DTrackID(a, i + _num_y[a]));
        track->setTrackReflBwd(get2DTrackID(ac, _num_x[a] - i - 1));

        if (_geometry->getMinYBoundaryType() == PERIODIC ||
            _geometry->getMinYBoundaryType() == INTERFACE) {
          track->setNextBwdFwd(false);
          track->setTrackNextBwd(get2DTrackID(a, i + _num_y[a]));
        }
        else {
          track->setTrackNextBwd(get2DTrackID(ac, _num_x[a] - i - 1));
          track->setNextBwdFwd(true);
        }
      }
      else {
        track->setNextBwdFwd(false);
        track->setTrackPrdcBwd(get2DTrackID(a, i - _num_x[a]));
        track->setTrackReflBwd(get2DTrackID(ac, i - _num_x[a]));

        if (track->getBCBwd() == PERIODIC || track->getBCBwd() == INTERFACE)
          track->setTrackNextBwd(get2DTrackID(a, i - _num_x[a]));
        else
          track->setTrackNextBwd(get2DTrackID(ac, i - _num_x[a]));
      }
    }
  }
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize() {

  log_printf(NORMAL, "Ray tracing for 2D track segmentation...");

  /* Check to ensure the Geometry is infinite in axial direction */
  double max_z = _geometry->getRootUniverse()->getMaxZ();
  double min_z = _geometry->getRootUniverse()->getMinZ();
  if (max_z - min_z < FLT_INFINITY) {
    log_printf(WARNING_ONCE, "The Geometry was set with non-infinite "
               "z-boundaries and supplied to a 2D TrackGenerator. The min-z "
               "boundary was set to %5.2f and the max-z boundary was set to "
               "%5.2f. Z-boundaries are assumed to be infinite in 2D "
               "TrackGenerators.", min_z, max_z);

#ifdef MPIx
    /* Check that the geometry is not domain decomposed in Z */
    int domains_xyz[3];
    _geometry->getDomainStructure(domains_xyz);
    if (domains_xyz[2] > 1)
      log_printf(ERROR, "A geometry with an axial domain domain decomposition "
                 "has been supplied to a 2D ray tracer.");
#endif

    /* Check that the track generator is ray tracing within bounds */
    if (_z_coord < min_z || _z_coord > max_z) {
      log_printf(WARNING_ONCE, "The track generator z-coordinate (%.2e cm) "
                 "lies outside the geometry bounds [%.2e %.2e] cm. Z"
                 "-coordinate moved to %.2e cm", _z_coord, min_z, max_z,
                 (min_z + max_z) / 2);
      _z_coord = (min_z + max_z) / 2;
    }
  }

#ifdef MPIx
  /* If domain decomposed, add artificial infinite Z bounds */
  if (_geometry->isDomainDecomposed()) {
    ZPlane* min_z_plane = new ZPlane(-FLT_INFINITY);
    ZPlane* max_z_plane = new ZPlane(+FLT_INFINITY);
    min_z_plane->setBoundaryType(REFLECTIVE);
    max_z_plane->setBoundaryType(REFLECTIVE);

    std::map<int, Cell*>::iterator cell;
    std::map<int, Cell*> cells = _geometry->getRootUniverse()->getCells();
    for (cell = cells.begin(); cell != cells.end(); ++cell) {
      (cell->second)->addSurface(+1, min_z_plane);
      (cell->second)->addSurface(-1, max_z_plane);
    }
    _geometry->getRootUniverse()->calculateBoundaries();
  }
#endif

  /* Make sure CMFD lattice is initialized and has right offset for 2D */
  Cmfd* cmfd = _geometry->getCmfd();
  if (cmfd != NULL) {

    /* Check that CMFD has been initialized */
    if (cmfd->getLattice() == NULL)
      log_printf(ERROR, "CMFD has not been initialized before generating "
                 "tracks. A call to geometry.initializeFlatSourceRegions() "
                 "may be missing.");

    /* Re-initialize CMFD lattice with 2D dimensions */
    Point offset;
    offset.setX(cmfd->getLattice()->getOffset()->getX());
    offset.setY(cmfd->getLattice()->getOffset()->getY());
    offset.setZ(_z_coord);
    cmfd->initializeLattice(&offset, true);
  }

  /* FSR numbering can change between two ray tracing */
  _geometry->resetContainsFSRCentroids();

  std::string msg = "Segmenting 2D tracks";
  Progress progress(_num_2D_tracks, msg, 0.1, _geometry, true);

  /* Loop over all Tracks */
#pragma omp parallel for schedule(dynamic)
  for (int t=0; t < _num_2D_tracks; t++) {
    _geometry->segmentize2D(_tracks_2D_array[t], _z_coord);
    progress.incrementCounter();
  }

  _geometry->initializeFSRVectors();
  _contains_2D_segments = true;

  /* Output memory consumption of explicit ray tracing */
  printMemoryReport();
}


/**
 * @brief This method creates a directory to store Track files, and reads
 *        in ray tracing data for Tracks and segments from a Track file
 *        if one exists.
 * @details This method is called by the TrackGenerator::generateTracks()
 *          class method. If a Track file exists for this Geometry, number
 *          of azimuthal angles, and track spacing, then this method will
 *          import the ray tracing Track and segment data to fill the
 *          appropriate data structures.
 */
void TrackGenerator::initializeTrackFileDirectory() {

  struct stat buffer;
  std::stringstream directory;

  /** Create directory to store Track files with pre-generated ray tracing data
   *  if the directory does not yet exist */
  directory << get_output_directory() << "/tracks";
  struct stat st;
  if (!(stat(directory.str().c_str(), &st) == 0))
    mkdir(directory.str().c_str(), S_IRWXU);

  /* Check to see if a Track file exists for this geometry, number of azimuthal
   * angles, and track spacing, and if so, import the ray tracing data */
  _tracks_filename = getTestFilename(directory.str());
  if (!stat(_tracks_filename.c_str(), &buffer)) {
    if (readSegmentsFromFile()) {
      _use_input_file = true;
      setContainsSegments(true);
    }
  }
}


/**
 * @brief Returns the filename for writing tracking data.
 */
std::string TrackGenerator::getTestFilename(std::string directory) {

  std::stringstream test_filename;

  if (_geometry->getCmfd() != NULL)
    test_filename << directory << "/2D_"
                  << _num_azim << "_azim_"
                  << _azim_spacing << "_cm_spacing_cmfd_"
                  << _geometry->getCmfd()->getNumX()
                  << "x" << _geometry->getCmfd()->getNumY();
  else
    test_filename << directory << "/2D_"
                  << _num_azim << "_angles_"
                  << _azim_spacing << "_cm_spacing";

  test_filename << ".data";

  return test_filename.str();
}


/**
 * @brief Updates whether the TrackGenerator contains segments.
 * @param contains_segments whether the TrackGenerator contains segments
 */
void TrackGenerator::setContainsSegments(bool contains_segments) {
  _contains_2D_segments = contains_segments;
}


/**
 * @brief Writes all Track and segment data to a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 */
void TrackGenerator::dumpSegmentsToFile() {

  /* Check whether the segments should be dumped */
  if (!_dump_segments)
    return;

  log_printf(NORMAL, "Dumping segments to file...");

  if (!containsSegments())
    log_printf(ERROR, "Unable to dump Segments to a file since no Segments "
               "have been generated for %d azimuthal angles and %f track "
               "spacing", _num_azim, _azim_spacing);

  FILE* out;
  out = fopen(_tracks_filename.c_str(), "w");
  if (out == NULL)
    log_printf(ERROR, "Segments file %s could not be written.",
               &_tracks_filename[0]);

  /* Get a string representation of the Geometry's attributes. This is used to
   * check whether or not ray tracing has been performed for this Geometry */
  std::string geometry_to_string = _geometry->toString();
  int string_length = geometry_to_string.length() + 1;

  /* Write geometry metadata to the Track file */
  fwrite(&string_length, sizeof(int), 1, out);
  fwrite(geometry_to_string.c_str(), sizeof(char)*string_length, 1, out);

  /* Write segment data to Track file */
  DumpSegments dump_segments(this);
  dump_segments.setOutputFile(out);
  if (_segment_formation == EXPLICIT_2D || _segment_formation == EXPLICIT_3D)
    dump_segments.execute();

  /* Get FSR vector maps */
  ParallelHashMap<std::string, fsr_data*>& FSR_keys_map =
      _geometry->getFSRKeysMap();
  std::vector<std::string>& FSRs_to_keys = _geometry->getFSRsToKeys();
  std::vector<int>& FSRs_to_material_IDs = _geometry->getFSRsToMaterialIDs();

  std::string fsr_key;
  long fsr_id;
  double x, y, z;

  /* Write number of FSRs */
  long num_FSRs = _geometry->getNumFSRs();
  fwrite(&num_FSRs, sizeof(long), 1, out);

  /* Write FSR vector maps to file */
  std::string* fsr_key_list = FSR_keys_map.keys();
  fsr_data** fsr_data_list = FSR_keys_map.values();
  Cmfd* cmfd = _geometry->getCmfd();
  for (long i=0; i < num_FSRs; i++) {

    /* Write data to file from FSR_keys_map */
    fsr_key = fsr_key_list[i];
    string_length = fsr_key.length() + 1;
    fwrite(&string_length, sizeof(int), 1, out);
    fwrite(fsr_key.c_str(), sizeof(char)*string_length, 1, out);

    fsr_id = fsr_data_list[i]->_fsr_id;
    x = fsr_data_list[i]->_point->getX();
    y = fsr_data_list[i]->_point->getY();
    z = fsr_data_list[i]->_point->getZ();
    fwrite(&fsr_id, sizeof(long), 1, out);
    fwrite(&x, sizeof(double), 1, out);
    fwrite(&y, sizeof(double), 1, out);
    fwrite(&z, sizeof(double), 1, out);

    /* Write data to file from FSRs_to_material_IDs */
    fwrite(&(FSRs_to_material_IDs.at(i)), sizeof(int), 1, out);

    /* Write data to file from FSRs_to_keys */
    fsr_key = FSRs_to_keys.at(i);
    string_length = fsr_key.length() + 1;
    fwrite(&string_length, sizeof(int), 1, out);
    fwrite(fsr_key.c_str(), sizeof(char)*string_length, 1, out);
  }

  /* Write cmfd_fsrs vector of vectors to file */
  if (cmfd != NULL) {
    std::vector< std::vector<long> >* cell_fsrs = cmfd->getCellFSRs();
    std::vector<long>::iterator iter;
    int num_cells = cmfd->getNumCells();
    fwrite(&num_cells, sizeof(int), 1, out);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      int num_cell_FSRs = cell_fsrs->at(cell).size();
      fwrite(&num_cell_FSRs, sizeof(int), 1, out);

      /* Loop over FSRs within cell */
      for (iter = cell_fsrs->at(cell).begin();
           iter != cell_fsrs->at(cell).end(); ++iter)
        fwrite(&(*iter), sizeof(long), 1, out);
    }
  }

  /* Delete key and value lists */
  delete [] fsr_key_list;
  delete [] fsr_data_list;

  /* Write 2D basis information for 3D solvers */
  if (_segment_formation != EXPLICIT_2D && _segment_formation != EXPLICIT_3D)
    writeExtrudedFSRInfo(out);

  /* Close the Track file */
  fclose(out);

  /* Inform other the TrackGenerator::generateTracks() method that it may
   * import ray tracing data from this file if it is called and the ray
   * tracing parameters have not changed */
  _use_input_file = true;
}


/**
 * @brief Write information of all Extruded FSRs to a file
 //TODO Use implementation in 3D track generator
 * @param out file to write to
 */
void TrackGenerator::writeExtrudedFSRInfo(FILE* out) {}


/**
 * @brief Reads Tracks in from a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 * @return true if able to read Tracks in from a file; false otherwise
 */
bool TrackGenerator::readSegmentsFromFile() {

  int ret;
  FILE* in = fopen(_tracks_filename.c_str(), "r");
  if (in == NULL)
    log_printf(ERROR, "Segments file %s could not be found.",
               &_tracks_filename[0]);

  int string_length;

  /* Import Geometry metadata from the Track file */
  ret = _geometry->twiddleRead(&string_length, sizeof(int), 1, in);
  char* geometry_to_string = new char[string_length];
  ret = _geometry->twiddleRead(geometry_to_string, sizeof(char)*string_length,
                               1, in);

  /* Check if our Geometry is exactly the same as the Geometry in the
   * Track file for this number of azimuthal angles and track spacing */
  if (_geometry->toString().compare(std::string(geometry_to_string)) != 0)
    return false;

  delete [] geometry_to_string;

  log_printf(NORMAL, "Importing ray tracing data from file...");

  /* Load all segment data into Tracks */
  ReadSegments read_segments(this);
  read_segments.setInputFile(in);
  if (_segment_formation == EXPLICIT_2D || _segment_formation == EXPLICIT_3D)
    read_segments.execute();

  /* Create FSR vector maps */
  ParallelHashMap<std::string, fsr_data*>& FSR_keys_map =
    _geometry->getFSRKeysMap();
  std::vector<int>& FSRs_to_material_IDs = _geometry->getFSRsToMaterialIDs();
  std::vector<std::string>& FSRs_to_keys = _geometry->getFSRsToKeys();
  std::vector<Point*>& FSRs_to_centroids = _geometry->getFSRsToCentroids();
  std::vector<int>& FSRs_to_CMFD_cells = _geometry->getFSRsToCMFDCells();

  long num_FSRs;
  std::string fsr_key;
  long fsr_key_id;
  double x, y, z;

  /* Get number of FSRs */
  ret = _geometry->twiddleRead(&num_FSRs, sizeof(long), 1, in);

  /* Resize vectors */
  FSRs_to_centroids.resize(num_FSRs);
  FSRs_to_CMFD_cells.resize(num_FSRs);

  /* Read FSR vector maps from file */
  for (long fsr_id=0; fsr_id < num_FSRs; fsr_id++) {

    /* Read key for FSR_keys_map */
    ret = _geometry->twiddleRead(&string_length, sizeof(int), 1, in);
    char* char_buffer1 = new char[string_length];
    ret = _geometry->twiddleRead(char_buffer1, sizeof(char)*string_length, 1, in);
    fsr_key = std::string(char_buffer1);

    /* Read data from file for FSR_keys_map */
    ret = _geometry->twiddleRead(&fsr_key_id, sizeof(long), 1, in);
    ret = _geometry->twiddleRead(&x, sizeof(double), 1, in);
    ret = _geometry->twiddleRead(&y, sizeof(double), 1, in);
    ret = _geometry->twiddleRead(&z, sizeof(double), 1, in);
    fsr_data* fsr = new fsr_data;
    fsr->_fsr_id = fsr_key_id;
    Point* point = new Point();
    point->setCoords(x,y,z);
    fsr->_point = point;
    FSR_keys_map.insert(fsr_key, fsr);

    /* Read data from file for FSR_to_materials_IDs */
    int material_id;
    ret = _geometry->twiddleRead(&material_id, sizeof(int), 1, in);
    FSRs_to_material_IDs.push_back(material_id);

    /* Read data from file for FSR_to_keys */
    ret = _geometry->twiddleRead(&string_length, sizeof(int), 1, in);
    char* char_buffer2 = new char[string_length];
    ret = _geometry->twiddleRead(char_buffer2, sizeof(char)*string_length, 1, in);
    fsr_key = std::string(char_buffer2);
    FSRs_to_keys.push_back(fsr_key);
  }

  /* Read cmfd cell_fsrs vector of vectors from file */
  Cmfd* cmfd = _geometry->getCmfd();
  if (cmfd != NULL) {
    std::vector< std::vector<long> > cell_fsrs;
    int num_cells;
    ret = _geometry->twiddleRead(&num_cells, sizeof(int), 1, in);
    long fsr_id;

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      std::vector<long>* fsrs = new std::vector<long>;
      cell_fsrs.push_back(*fsrs);
      int num_cell_FSRs;
      ret = _geometry->twiddleRead(&num_cell_FSRs, sizeof(int), 1, in);

      /* Loop over FSRs within cell */
      for (int fsr = 0; fsr < num_cell_FSRs; fsr++) {
        ret = _geometry->twiddleRead(&fsr_id, sizeof(long), 1, in);
        cell_fsrs.at(cell).push_back(fsr_id);
        FSRs_to_CMFD_cells.at(fsr_id) = cell;
      }
    }

    /* Set CMFD cell_fsrs vector of vectors */
    cmfd->setCellFSRs(&cell_fsrs);
  }

  /* Read 2D basis information for OTF 3D solvers */
  if (_segment_formation != EXPLICIT_2D && _segment_formation != EXPLICIT_3D)
    readExtrudedFSRInfo(in);

  /* Close the Track file */
  fclose(in);

  return true;
}


/**
 * @brief Read information of all Extruded FSRs from a file.
 * @param in file to read from
 */
void TrackGenerator::readExtrudedFSRInfo(FILE* in) {}


/**
 * @brief Splits Track segments into sub-segments for a user-defined
 *        maximum optical length for the problem.
 * @details This routine is needed so that all segment lengths fit
 *          within the exponential interpolation table used in the MOC
 *          transport sweep.
 * @param max_optical_length the maximum optical length
 */
void TrackGenerator::splitSegments(FP_PRECISION max_optical_length) {

  if (!containsSegments())
    log_printf(ERROR, "Unable to split segments since segments have not yet "
                      "been generated");

  if (_segment_formation != EXPLICIT_3D && _segment_formation != EXPLICIT_2D)
    log_printf(ERROR, "Segments cannot be split for on-the-fly ray tracing");

  /* Split all segments along all Tracks */
  _max_optical_length = max_optical_length;
  SegmentSplitter segment_splitter(this);
  segment_splitter.execute();
}


/**
 * @brief Generates the numerical centroids of the FSRs.
 * @details This routine generates the numerical centroids of the FSRs
 *          by weighting the average x and y values of each segment in the
 *          FSR by the segment's length and azimuthal weight. The numerical
 *          centroid fomula can be found in R. Ferrer et. al. "Linear Source
 *          Approximation in CASMO 5", PHYSOR 2012.
 * @param FSR_volumes An array of FSR volumes.
 */
void TrackGenerator::generateFSRCentroids(FP_PRECISION* FSR_volumes) {

  long num_FSRs = _geometry->getNumFSRs();

  /* Create temporary array of centroids and initialize to origin */
  Point** centroids = new Point*[num_FSRs];
  for (long r=0; r < num_FSRs; r++) {
    centroids[r] = new Point();
    centroids[r]->setCoords(0.0, 0.0, 0.0);
  }

  /* Generate FSR centroids by looping over all Tracks */
  CentroidGenerator centroid_generator(this);
  centroid_generator.setCentroids(centroids);
  centroid_generator.execute();

  /* Set the centroid for the FSR */
  for (long r=0; r < num_FSRs; r++)
    _geometry->setFSRCentroid(r, centroids[r]);

  /* Recenter the segments around FSR centroid */
  if ((_segment_formation == EXPLICIT_2D || _segment_formation == EXPLICIT_3D)
      && _segments_centered == false) {
    log_printf(NORMAL, "Centering segments around FSR centroid...");
    RecenterSegments rs(this);
    rs.execute();
    _segments_centered = true;
  }

  /* Print FSR volumes, centroids and volume moments for debugging purposes */
  double total_volume[4];
  memset(total_volume, 0, 4 * sizeof(double));
  FP_PRECISION min_volume = 1e10;
  FP_PRECISION max_volume = 0.;

  for (long r=0; r < num_FSRs; r++) {
    total_volume[0] += _FSR_volumes[r];
    total_volume[1] += _FSR_volumes[r] * centroids[r]->getX();
    total_volume[2] += _FSR_volumes[r] * centroids[r]->getY();
    total_volume[3] += _FSR_volumes[r] * centroids[r]->getZ();

    min_volume = std::min(_FSR_volumes[r], min_volume);
    max_volume = std::max(_FSR_volumes[r], max_volume);

    log_printf(DEBUG, "FSR ID = %d has volume = %.6f, centroid"
               " (%.3f %.3f %.3f)", r, _FSR_volumes[r], centroids[r]->getX(),
               centroids[r]->getY(), centroids[r]->getZ());
  }

  log_printf(DEBUG, "Total volume %.6f cm3, moments of volume "
             "(%.4e %.4e %.4e).", total_volume[0], total_volume[1],
             total_volume[2], total_volume[3]);
  log_printf(DEBUG, "Average / min / max volumes of FSRs : %.2e / %.2e / %.2e",
             total_volume[0] / num_FSRs, min_volume, max_volume);
  delete [] centroids;
}


/**
 * @brief Sets the max optical path length of 3D segments for use in
 *        on-the-fly computation.
 * @param tau maximum optical path length
 */
void TrackGenerator::setMaxOpticalLength(FP_PRECISION tau) {
  _max_optical_length = tau;
}


/**
 * @brief Sets the maximum number of segments per Track.
 * @param max_num_segments the maximum number of segments per Track
 */
void TrackGenerator::setMaxNumSegments(int max_num_segments) {
  _max_num_segments = max_num_segments;
}


/**
 * @brief Retrieves the max optical path length of 3D segments for use in
 *        on-the-fly computation.
 * @return maximum optical path length
 */
FP_PRECISION TrackGenerator::retrieveMaxOpticalLength() {
  return _max_optical_length;
}


/**
 * @brief Counts the number of segments for each Track in the Geometry.
 * @details All segments are subject to the max optical path length to
 *          determine the number of segments for each track as well as the
 *          maximum number of segments per Track in the Geometry. For
 *          on-the-fly calculations, the temporary segment buffer is expanded
 *          to fit the calculated maximum number of segments per Track.
 */
void TrackGenerator::countSegments() {

  /* Count the number of segments on each track and update the maximum */
  SegmentCounter counter(this);
  counter.execute();

  /* Allocate new temporary segments if necessary */
  if (_segment_formation != EXPLICIT_3D && _segment_formation != EXPLICIT_2D)
    allocateTemporarySegments();
}


/**
 * @brief Creates a Track array by increasing uid.
 * @details An array is created which indexes Tracks by increasing uid.
 *          Parallel groups are also initialized -- groups of Tracks that can
 *          be computed in parallel without the potential of overwriting
 *          angular fluxes of connecting tracks prematurely.
 */
void TrackGenerator::initializeTracksArray() {

  log_printf(NORMAL, "Initializing 2D tracks array...");

  /* Allocate memory for tracks array */
  if (_tracks_2D_array != NULL)
    delete [] _tracks_2D_array;
  long num_2D_tracks = getNum2DTracks();
  _tracks_2D_array = new Track*[num_2D_tracks];

  /* Loop over all 2D tracks */
  long uid = 0;
  for (int a = 0; a < _num_azim / 2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      /* Get current track and azim group ids */
      Track* track = &_tracks_2D[a][i];

      track->setUid(uid);
      _tracks_2D_array[uid] = track;
      uid++;
    }
  }
}


/**
 * @brief Returns whether periodic boundaries are present in Track generation.
 * @return a boolean value - true if periodic; false otherwise
 */
bool TrackGenerator::getPeriodic() {
  return _periodic;
}


/**
 * @brief Sets a flag to record all segment information in the tracking file.
 * @param dump_segments whether or not to record segment information in the
 *        tracking file: true to record, false not to record
 */
void TrackGenerator::setDumpSegments(bool dump_segments) {
  _dump_segments = dump_segments;
}


/**
 * @brief Resets the TrackGenerator to not contain tracks or segments.
 */
void TrackGenerator::resetStatus() {
  _contains_2D_tracks = false;
  _contains_2D_segments = false;
  _use_input_file = false;
  _tracks_filename = "";
}


/**
 * @brief Allocates memory for temporary segment storage if necessary.
 * @details Temporary segments are not allocated for 2D calculations
 */
void TrackGenerator::allocateTemporarySegments() {}


/**
 * @brief Get the id of a 2D Track based on its azimuthal angle and index in the
 *        azimuthal stack.
 * @param a azimuthal angle of the Track
 * @param x index in azimuthal stack
 * @return Track unique id
 */
int TrackGenerator::get2DTrackID(int a, int x) {

  int uid = 0;

  for (int ai = 0; ai < a; ai++)
    uid += getNumX(ai) + getNumY(ai);

  uid += x;
  return uid;
}


/**
 * @brief Print the track generation timer report.
 * @param mpi_reduce whether to reduce timer across MPI processes (only once!)
 */
void TrackGenerator::printTimerReport(bool mpi_reduce) {

#ifdef MPIx
  if (_geometry->isDomainDecomposed() && mpi_reduce)
    _timer->reduceTimer(_geometry->getMPICart());
#endif

  double gen_time = _timer->getSplit("Track Generation Time");
  std::string msg_string = "Total Track Generation & Segmentation Time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), gen_time);
}


/**
 * @brief Print the track generation memory report.
 */
void TrackGenerator::printMemoryReport() {

  long track_storage = _num_2D_tracks * sizeof(Track);
  long segment_storage = getNum2DSegments() * sizeof(segment);
  long max_segment_storage = segment_storage;

#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Allreduce(&segment_storage, &max_segment_storage, 1, MPI_LONG, MPI_MAX,
                  _geometry->getMPICart());
    /* NOTE: Same tracks on every domain */
#endif

  log_printf(INFO_ONCE, "Max 2D explicit track storage per domain %.2f MB",
             track_storage / float(1e6));
  log_printf(INFO_ONCE, "Max 2D explicit segment storage per domain %.2f MB",
             max_segment_storage / float(1e6));
}
