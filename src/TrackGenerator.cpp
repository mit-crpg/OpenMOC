#include "TrackGenerator.h"
#include "TrackTraversingAlgorithms.h"

/**
 * @brief Constructor for the TrackGenerator assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param azim_spacing azimuthal track spacing (cm)
 */
TrackGenerator::TrackGenerator(Geometry* geometry, int num_azim,
                               double azim_spacing) {

  setNumThreads(1);
  _geometry = geometry;
  setNumAzim(num_azim);
  setDesiredAzimSpacing(azim_spacing);
  _quadrature = NULL;
  _user_quadrature = false;
  _contains_tracks = false;
  _contains_segments = false;
  _use_input_file = false;
  _dump_segments = true;
  _tracks_filename = "";
  _z_coord = 0.0;
  _segment_formation = EXPLICIT_2D;
  _max_optical_length = std::numeric_limits<FP_PRECISION>::max();
  _FSR_volumes = NULL;
  _FSR_locks = NULL;
  _timer = new Timer();
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator::~TrackGenerator() {

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_tracks) {
    delete [] _num_tracks;
    delete [] _num_x;
    delete [] _num_y;
    delete [] _num_tracks_by_parallel_group;

    for (int i = 0; i < _num_azim_2; i++)
      delete [] _tracks[i];

    delete [] _tracks;
    delete [] _tracks_by_parallel_group;
  }

  if (_FSR_locks != NULL)
    delete [] _FSR_locks;

  if (_FSR_volumes != NULL)
    delete [] _FSR_volumes;

  if (_quadrature != NULL && !_user_quadrature)
    delete _quadrature;

  if (_timer != NULL)
    delete _timer;
}

/**
 * @brief Return the number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @return the number of azimuthal angles in \f$ 2\pi \f$
 */
int TrackGenerator::getNumAzim() {
  return 2*_num_azim_2;
}


/**
 * @brief Return the azimuthal track spacing (cm).
 * @details This will return the user-specified azimuthal track spacing and NOT
 *          the effective track spacing which is computed and used to generate
 *          cyclic tracks.
 * @return the azimuthal track spacing (cm)
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
 * @brief Return the total number of Tracks generated.
 * @return The number of Tracks generated
 */
int TrackGenerator::getNumTracks() {
  if (!_contains_tracks)
    log_printf(ERROR, "Unable to return the total number of Tracks since "
               "Tracks have not yet been generated.");

  int num_tracks = 0;

  for (int i=0; i < _num_azim_2; i++) {
    num_tracks += _num_tracks[i];
  }

  return num_tracks;
}


/**
 * @brief Return the number of parallel track groups.
 * @return The number of parallel track groups
 */
int TrackGenerator::getNumParallelTrackGroups() {

  if (!_contains_tracks)
    log_printf(ERROR, "Unable to return the number of parallel track groups "
               "since Tracks have not yet been generated.");

  return _num_parallel_track_groups;
}


/**
 * @brief Return the number of tracks on the x-axis for a given azimuthal angle.
 * @param azim An azimuthal angle index
 * @return The number of Tracks on the x-axis
 */
int TrackGenerator::getNumX(int azim) {
  if (!_contains_tracks)
    log_printf(ERROR, "Unable to return the number of Tracks on the x-axis"
               " for azimuthal angle %d since Tracks have not yet been "
               "generated.", azim);

  return _num_x[azim];
}


/**
 * @brief Return the number of tracks on the y-axis for a given azimuthal angle.
 * @param azim An azimuthal angle index
 * @return The number of Tracks on the y-axis
 */
int TrackGenerator::getNumY(int azim) {
  if (!_contains_tracks)
    log_printf(ERROR, "Unable to return the number of Tracks on the y-axis"
               " for azimuthal angle %d since Tracks have not yet been "
               "generated.", azim);

  return _num_y[azim];
}


/**
 * @brief Return the number of tracks in a given parallel track group.
 * @return the number of tracks in a given parallel track group.
 */
int TrackGenerator::getNumTracksByParallelGroup(int group) {

  if (group < 0 || group >= _num_parallel_track_groups)
    log_printf(ERROR, "Unable to return the number of tracks in parallel track"
               " group %d which is not between 0 and the number of parallel "
               "groups, %d", group, _num_parallel_track_groups);

  if (!_contains_tracks)
    log_printf(ERROR, "Unable to return the number of tracks in the parallel "
               "track group %d since Tracks have not yet been generated."
               , group);

  return _num_tracks_by_parallel_group[group];
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator::getNumSegments() {
  if (!_contains_segments)
    log_printf(ERROR, "Unable to return the total number of segments since "
               "Tracks have not yet been generated.");

  int num_segments = 0;

  for (int i=0; i < _num_azim_2; i++) {
#pragma omp parallel for reduction(+:num_segments)
    for (int j=0; j < _num_tracks[i]; j++)
      num_segments += _tracks[i][j].getNumSegments();
  }

  return num_segments;
}


/**
 * @brief Returns a 2D jagged array of the Tracks.
 * @details The first index into the array is the azimuthal angle and the
 *          second index is the Track number for a given azimuthal angle.
 * @return the 2D jagged array of Tracks
 */
Track **TrackGenerator::getTracks() {
  if (!_contains_tracks)
    log_printf(ERROR, "Unable to return the 2D ragged array of the Tracks "
               "since Tracks have not yet been generated.");

  return _tracks;
}


/**
 * @brief Returns a 1D array of Track pointers.
 * @details The tracks in the _tracks_by_parallel_group array are organized
 *          by parallel track group. The index into the array is also the
 *          corresponding Track's UID.
 * @return The 1D array of Track pointers
 */
Track **TrackGenerator::getTracksByParallelGroup() {
  if (!_contains_tracks)
    log_printf(ERROR, "Unable to return the 1D array of the Track pointers "
               "arranged by parallel group since Tracks have not yet been "
               "generated.");

  return _tracks_by_parallel_group;
}


/**
 * @brief Returns a pointer to the Quadrature.
 * @return a pointer to the Quadrature
 */
Quadrature* TrackGenerator::getQuadrature() {

  if (_quadrature == NULL)
    log_printf(ERROR, "Unable to return the TrackGenerator's Quadrature "
               "since it has not yet been set");

  return _quadrature;
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int TrackGenerator::getNumThreads() {
  return _num_threads;
}


/**
 * @brief Returns the z-coord where the 2D Tracks should be created.
 * @return the z-coord where the 2D Tracks should be created.
 */
double TrackGenerator::getZCoord() {
  return _z_coord;
}


/**
 * @brief Returns an array of volumes indexed by FSR.
 * @return a pointer to the array of FSR volumes
 */
FP_PRECISION* TrackGenerator::getFSRVolumes() {

  if (!containsTracks())
    log_printf(ERROR, "Unable to get the FSR volumes since tracks "
               "have not yet been generated");

 if (_FSR_volumes != NULL)
   return _FSR_volumes;

#pragma omp critical
  {
    if (_FSR_volumes == NULL) {
      int num_FSRs = _geometry->getNumFSRs();
      _FSR_volumes = new FP_PRECISION[num_FSRs];
      calculateFSRVolumes();
    }
  }

  return _FSR_volumes;
}


/**
 * @brief Calculates the volume of each FSR
 * @details The _FSR_volumes array is reset to zero and then FSR
 *          volumes are re-computed.
 */
void TrackGenerator::calculateFSRVolumes() {


  if (_FSR_volumes == NULL)
    log_printf(ERROR, "Unable to calculate the FSR volumes since the FSR "
               "volumes array has not yet been allocated");

  /* Reset FSR volumes to zero */
  int num_FSRs = _geometry->getNumFSRs();
  if (_FSR_volumes != NULL)
    memset(_FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));

  /* Create volume calculator and calculate new FSR volumes */
  VolumeCalculator volume_calculator(this);
  volume_calculator.execute();
}


/**
 * @brief Deletes the memory associated with the FSR volumes and resets it NULL
 */
void TrackGenerator::resetFSRVolumes() {

#pragma omp critical
  {
    if (_FSR_volumes != NULL) {
      delete [] _FSR_volumes;
      _FSR_volumes = NULL;
    }
  }
}


/**
 * @brief Computes and returns the volume of an FSR.
 * @param fsr_id the ID for the FSR of interest
 * @return the FSR volume
 */
FP_PRECISION TrackGenerator::getFSRVolume(int fsr_id) {

  if (!containsTracks())
    log_printf(ERROR, "Unable to get the FSR %d volume since tracks "
               "have not yet been generated");

  else if (fsr_id < 0 || fsr_id > _geometry->getNumFSRs())
    log_printf(ERROR, "Unable to get the volume for FSR %d since the FSR IDs "
               "lie in the range (0, %d)", fsr_id, _geometry->getNumFSRs());

  segment* curr_segment;
  FP_PRECISION volume = 0;

  /* Calculate the FSR's "volume" by accumulating the total length of *
   * all Track segments multipled by the Track "widths" for the FSR.  */
  for (int i=0; i < _num_azim_2; i++) {
#pragma omp parallel for reduction(+:volume) private(curr_segment)
    for (int j=0; j < _num_tracks[i]; j++) {
      for (int s=0; s < _tracks[i][j].getNumSegments(); s++) {
        curr_segment = _tracks[i][j].getSegment(s);
        if (curr_segment->_region_id == fsr_id)
          volume += curr_segment->_length * _quadrature->getAzimWeight(i)
            * _quadrature->getAzimSpacing(i);
      }
    }
  }

  return volume;
}


/**
 * @breif Calculates and returns the maximum optcial length for any segment
 *        in the Geomtry.
 * @details The maximum optical length is recomputed, updated, and returned.
 *          This value determines the when segments must be split during ray
 *          tracing.
 * @return The maximum optical length of any segment in the Geometry
 */
FP_PRECISION TrackGenerator::getMaxOpticalLength() {
  MaxOpticalLength update_max_optical_length(this);
  update_max_optical_length.execute();
  return _max_optical_length;
}


/**
 * @brief Retrieves the max optical path length of 3D segments for use in
 *        on-the-fly computation
 * @return maximum optical path length
 */
FP_PRECISION TrackGenerator::retrieveMaxOpticalLength() {
  return _max_optical_length;
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

  _num_threads = num_threads;

  /* Set the number of threads for OpenMP */
  omp_set_num_threads(_num_threads);
}


/**
 * @brief Sets the z-coord where the 2D Tracks should be created.
 * @param z_coord the z-coord where the 2D Tracks should be created.
 */
void TrackGenerator::setZCoord(double z_coord) {
  _z_coord = z_coord;
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

  _num_azim_2 = num_azim/2;
  resetStatus();
}


/**
 * @brief Set the suggested azimuthal track spacing (cm).
 * @param azim_spacing the suggested azimuthal track spacing
 */
void TrackGenerator::setDesiredAzimSpacing(double azim_spacing) {
  if (azim_spacing < 0)
    log_printf(ERROR, "Unable to set a negative azimuthal track spacing %f for"
               " the TrackGenerator.", azim_spacing);

  _azim_spacing = azim_spacing;
  resetStatus();
}


/**
 * @brief Set a pointer to the Geometry to use for track generation.
 * @param geometry a pointer to the Geometry
 */
void TrackGenerator::setGeometry(Geometry* geometry) {
  _geometry = geometry;
  resetStatus();
}


/**
 * @brief Assign a Quadrature object to the Solver.
 * @details This routine allows use of a Quadrature with any polar angle
 *          quadrature. Alternatively, this routine may take in any subclass
 *          of the Quadrature parent class, including TYPolarQuad (default),
 *          LeonardPolarQuad, GLPolarQuad, etc.
 *
 *          Users may assign a Quadrature object to the Solver from
 *          Python script as follows:
 *
 * @code
 *          quadrature = openmoc.LeonardPolarQuad()
 *          quadrature.setNumPolarAngles(2)
 *          track_generator.setQuadrature(quadrature)
 * @endcode
 *
 * @param quadrature a pointer to a Quadrature object
 */
void TrackGenerator::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
  _user_quadrature = true;
}


/**
 * @brief Returns whether or not the TrackGenerator contains Track that are
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains Tracks; false otherwise
 */
bool TrackGenerator::containsTracks() {
  return _contains_tracks;
}


/**
 * @brief Returns whether or not the TrackGenerator contains segments for
 *        the current Track laydown.
 * @return true if the TrackGenerator conatains segments; false otherwise
 */
bool TrackGenerator::containsSegments() {
  return _contains_segments;
}



/**
 * @brief Fills an array with the x,y,z coordinates for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires on due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords
 *                     (num_tracks*NUM_VALUES_PER_RETRIEVED_TRACK)
 * @endcode
 *
 * @param coords an array of coords of length NUM_VALUES_PER_RETRIEVED_TRACK
 *        times the number of Tracks
 * @param length_coords the total number of Tracks times
 *                      NUM_VALUES_PER_RETRIEVED_TRACK
 */
void TrackGenerator::retrieveTrackCoords(double* coords, int length_coords) {

  if (length_coords != NUM_VALUES_PER_RETRIEVED_TRACK*getNumTracks())
    log_printf(ERROR, "Unable to retrieve the Track coordinates since the "
               "TrackGenerator contains %d Tracks with %d coordinates per track"
               " but an array of length %d was input", getNumTracks(),
               NUM_VALUES_PER_RETRIEVED_TRACK, length_coords);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  for (int i=0; i < _num_azim_2; i++) {
    for (int j=0; j < _num_tracks[i]; j++) {
      coords[counter] = _tracks[i][j].getStart()->getX();
      coords[counter+1] = _tracks[i][j].getStart()->getY();
      coords[counter+2] = _tracks[i][j].getStart()->getZ();
      coords[counter+3] = _tracks[i][j].getEnd()->getX();
      coords[counter+4] = _tracks[i][j].getEnd()->getY();
      coords[counter+5] = _tracks[i][j].getEnd()->getZ();
      counter += NUM_VALUES_PER_RETRIEVED_TRACK;
    }
  }

  return;
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
 *          coords = track_generator.retrieveSegmentCoords
 *                     (num_segments*NUM_VALUES_PER_RETRIEVED_SEGMENT)
 * @endcode
 *
 * @param coords an array of coords of length NUM_VALUES_PER_RETRIEVED_SEGMENT
 *               times the number of segments
 * @param length_coords the total number of Track segments times
 *                      NUM_VALUES_PER_RETRIEVED_SEGMENT
 */
void TrackGenerator::retrieveSegmentCoords(double* coords, int length_coords) {

  if (length_coords != NUM_VALUES_PER_RETRIEVED_SEGMENT*getNumSegments())
    log_printf(ERROR, "Unable to retrieve the Track segment coordinates since "
               "the TrackGenerator contains %d segments with %d values per "
               "segment but an array of length %d was input", getNumSegments(),
               NUM_VALUES_PER_RETRIEVED_SEGMENT, length_coords);

  segment* curr_segment = NULL;
  double x0, x1, y0, y1, z;
  double phi;
  segment* segments;

  int counter = 0;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  for (int i=0; i < _num_azim_2; i++) {
    for (int j=0; j < _num_tracks[i]; j++) {

      x0 = _tracks[i][j].getStart()->getX();
      y0 = _tracks[i][j].getStart()->getY();
      z = _tracks[i][j].getStart()->getZ();
      phi = _tracks[i][j].getPhi();

      segments = _tracks[i][j].getSegments();

      for (int s=0; s < _tracks[i][j].getNumSegments(); s++) {
        curr_segment = &segments[s];

        coords[counter] = curr_segment->_region_id;

        coords[counter+1] = x0;
        coords[counter+2] = y0;
        coords[counter+3] = z;

        x1 = x0 + cos(phi) * curr_segment->_length;
        y1 = y0 + sin(phi) * curr_segment->_length;

        coords[counter+4] = x1;
        coords[counter+5] = y1;
        coords[counter+6] = z;

        x0 = x1;
        y0 = y1;

        counter += NUM_VALUES_PER_RETRIEVED_SEGMENT;
      }
    }
  }

    return;
}


/**
 * @brief Generates tracks for some number of azimuthal angles and track spacing
 * @details Computes the effective angles and track spacing. Computes the
 *          number of Tracks for each azimuthal angle, allocates memory for
 *          all Tracks at each angle and sets each Track's starting and ending
 *          Points, azimuthal angle, and azimuthal angle quadrature weight.
 * @brief neighbor_cells whether to use neighbor cell optimizations
 * @brief store whether to store the tracks to a file for reuse
 */
void TrackGenerator::generateTracks(bool store, bool neighbor_cells) {

  if (_geometry == NULL)
    log_printf(ERROR, "Unable to generate Tracks since no Geometry "
               "has been set for the TrackGenerator");

  /* Check for valid quadrature */
  if (_quadrature != NULL) {
    if (_quadrature->getNumAzimAngles() != 2*_num_azim_2) {
      if (_user_quadrature) {
        log_printf(ERROR, "User defined quadrature does not match the "
                          " number of azimuthal angles in the TrackGenerator");
      }
      else {
        delete _quadrature;
        _quadrature = NULL;
      }
    }
  }

  /* Initialize quadrature */
  if (_quadrature == NULL) {
    _quadrature = new TYPolarQuad();
    _quadrature->setNumAzimAngles(2*_num_azim_2);
    _quadrature->setNumPolarAngles(6);
  }
  _quadrature->initialize();

  /* Clear all timing data from previous track generation */
  clearTimerSplits();

  /* Start the timer to record the total time to generate tracks */
  _timer->startTimer();

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_tracks) {
    delete [] _num_tracks;
    delete [] _num_tracks_by_parallel_group;
    delete [] _num_x;
    delete [] _num_y;
    delete [] _tracks_by_parallel_group;

    for (int i = 0; i < _num_azim_2; i++)
      delete [] _tracks[i];

    delete [] _tracks;
  }

  /* Initialize the CMFD object */
  if (_geometry->getCmfd() != NULL)
    _geometry->initializeCmfd();

  /* Initialize FSRs with pin cell discretization and neighbor cell lists */
  _geometry->initializeFSRs(neighbor_cells);

  initializeTrackFileDirectory();

  /* If not Tracks input file exists, generate Tracks */
  if (_use_input_file == false) {

    /* Allocate memory for the Tracks */
    try {
      _num_tracks = new int[_num_azim_2];
      _num_x = new int[_num_azim_2];
      _num_y = new int[_num_azim_2];
      _tracks = new Track*[_num_azim_2];
    }
    catch (std::exception &e) {
      log_printf(ERROR, "Unable to allocate memory for TrackGenerator");
    }

    /* Check to make sure that width of the Geometry in x and y are nonzero */
    if (_geometry->getWidthX() <= 0 || _geometry->getWidthY() <= 0)
      log_printf(ERROR, "The total x-width and y-width of the Geometry must be "
                 "non-zero for Track generation. Create a Cell which "
                 "is filled by the entire geometry and bounded by XPlanes "
                 "and YPlanes to enable the Geometry to determine the total "
                 "x-width and y-width of the model.");

    /* Generate Tracks, perform ray tracing across the geometry, and store
     * the data to a Track file */
    try {
      initializeTracks();
      recalibrateTracksToOrigin();
      segmentize();
      if (store)
	      dumpTracksToFile();
    }
    catch (std::exception &e) {
      log_printf(ERROR, "Unable to allocate memory for Tracks");
    }
  }
  else {

    /* Determine azimuthal spacings */
    for (int i = 0; i < _num_azim_2/2; i++) {

      /* Get the azimuthal angle */
      double phi = _quadrature->getPhi(i);

      /* Recompute azimuthal track spacing */
      double dx = _geometry->getWidthX() / _num_x[i];
      _quadrature->setAzimSpacing(dx * sin(phi), i);
    }
  }

  /* Precompute quadrature weights */
  _quadrature->precomputeWeights(false);

  /* Initialize the track boundary conditions and set the track UIDs */
  initializeBoundaryConditions();
  initializeTrackCycleIndices(PERIODIC);
  initializeTrackUids();
  initializeFSRLocks();
  initializeVolumes();

  _timer->stopTimer();
  _timer->recordSplit("Total time");

  return;
}


/**
 * @brief Create an array of OpenMP mutual exclusion locks for each FSR.
 * @details This method allocates and initializes an array of OpenMP
 *          mutual exclusion locks for each FSR for use in loops that
 *          set or update values by FSR.
 */
void TrackGenerator::initializeFSRLocks() {

  /* Delete old FSR locks, if they exist */
  if (_FSR_locks != NULL)
    delete [] _FSR_locks;

  /* Allocate array of mutex locks for each FSR */
  int num_FSRs = _geometry->getNumFSRs();
  _FSR_locks = new omp_lock_t[num_FSRs];

  /* Loop over all FSRs to initialize OpenMP locks */
#pragma omp parallel for schedule(guided)
  for (int r=0; r < num_FSRs; r++)
    omp_init_lock(&_FSR_locks[r]);
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

  std::stringstream directory;
  struct stat buffer;
  std::stringstream test_filename;

  /** Create directory to store Track files with pre-generated ray tracing data
   *  if the directory does not yet exist */

  directory << get_output_directory() << "/tracks";
  struct stat st;
  if ((!stat(directory.str().c_str(), &st)) == 0)
    mkdir(directory.str().c_str(), S_IRWXU);

  if (_geometry->getCmfd() != NULL) {
    test_filename << directory.str() << "/"
                  << 2*_num_azim_2 << "_angles_"
                  << _azim_spacing << "_cm_spacing_z_"
                  << _z_coord << "_("
                  << _geometry->getCmfd()->getNumX()
                  << "x" << _geometry->getCmfd()->getNumY()
                  << ")_cmfd.data";
    }
  else{
    test_filename << directory.str() << "/"
                  << 2*_num_azim_2 << "_angles_"
                  << _azim_spacing << "_cm_spacing_z_"
                  << _z_coord << ".data";
  }

  _tracks_filename = test_filename.str();

  /* Check to see if a Track file exists for this geometry, number of azimuthal
   * angles, and track spacing, and if so, import the ray tracing data */
  if ((!stat(_tracks_filename.c_str(), &buffer))) {
    if (readTracksFromFile()) {
      _use_input_file = true;
      _contains_tracks = true;
      _contains_segments = true;
    }
  }
}


/**
 * @brief Returns the type of ray tracing used for segment formation
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

  log_printf(INFO, "Computing azimuthal angles and track spacing...");

  /* Each element in arrays corresponds to an angle in phi_eff */
  /* Track spacing along x,y-axes, and perpendicular to each Track */
  double* dx_eff = new double[_num_azim_2];
  double* dy_eff = new double[_num_azim_2];
  double* d_eff = new double[_num_azim_2];

  double x1, x2;
  double width_x = _geometry->getWidthX();
  double width_y = _geometry->getWidthY();

  /* Determine azimuthal angles and track spacing */
  for (int i = 0; i < _num_azim_2/2; i++) {

    /* A desired azimuthal angle for the user-specified number of
     * azimuthal angles */
    double phi = M_PI / _num_azim_2 * (0.5 + i);

    /* The number of intersections with x,y-axes */
    _num_x[i] = (int) (fabs(width_x / _azim_spacing * sin(phi))) + 1;
    _num_y[i] = (int) (fabs(width_y / _azim_spacing * cos(phi))) + 1;

    /* Total number of Tracks */
    _num_tracks[i] = _num_x[i] + _num_y[i];

    /* Effective/actual angle (not the angle we desire, but close) */
    phi = atan((width_y * _num_x[i]) / (width_x * _num_y[i]));
    _quadrature->setPhi(phi, i);

    /* Effective Track spacing (not spacing we desire, but close) */
    dx_eff[i] = width_x / _num_x[i];
    dy_eff[i] = width_y / _num_y[i];
    d_eff[i] = dx_eff[i] * sin(phi);
    _quadrature->setAzimSpacing(d_eff[i], i);

    /* Set attributes for complimentary angles */
    _num_x[_num_azim_2-i-1] = _num_x[i];
    _num_y[_num_azim_2-i-1] = _num_y[i];
    _num_tracks[_num_azim_2-i-1] = _num_tracks[i];
    dx_eff[_num_azim_2-i-1] = dx_eff[i];
    dy_eff[_num_azim_2-i-1] = dy_eff[i];
    d_eff[_num_azim_2-i-1] = d_eff[i];
  }

  log_printf(INFO, "Generating Track start and end points...");

  /* Compute Track starting and end points */
  for (int i = 0; i < _num_azim_2; i++) {

    /* Extract the azimuthal angle */
    double phi = _quadrature->getPhi(i);

    /* Tracks for azimuthal angle i */
    _tracks[i] = new Track[_num_tracks[i]];

    /* Compute start points for Tracks starting on x-axis */
    for (int j = 0; j < _num_x[i]; j++) {
      if (i < _num_azim_2 / 2)
        _tracks[i][j].getStart()->setCoords(
            dx_eff[i] * (_num_x[i] - j - 0.5), 0, _z_coord);
      else
        _tracks[i][j].getStart()->setCoords(dx_eff[i] * (0.5 + j), 0, _z_coord);
    }

    /* Compute start points for Tracks starting on y-axis */
    for (int j = 0; j < _num_y[i]; j++) {

      /* If Track points to the upper right */
      if (i < _num_azim_2 / 2)
        _tracks[i][_num_x[i]+j].getStart()->setCoords(
            0, dy_eff[i] * (0.5 + j), _z_coord);

      /* If Track points to the upper left */
      else
        _tracks[i][_num_x[i]+j].getStart()->setCoords(
            width_x, dy_eff[i] * (0.5 + j), _z_coord);
    }

    /* Compute the exit points for each Track */
    for (int j = 0; j < _num_tracks[i]; j++) {

      /* Set the Track's end point */
      Point* start = _tracks[i][j].getStart();
      Point* end = _tracks[i][j].getEnd();
      computeEndPoint(start, end, phi, width_x, width_y);

      /* Set the Track's azimuthal angle */
      _tracks[i][j].setPhi(phi);
    }
  }

  _contains_tracks = true;
  _contains_segments = true;
  delete [] dx_eff;
  delete [] dy_eff;
  delete [] d_eff;
}


/**
 * @brief Recalibrates Track start and end points to the origin of the Geometry.
 * @details The origin of the Geometry is designated at its center by
 *          convention, but for track initialization the origin is assumed to be
 *          at the bottom right corner for simplicity. This method corrects
 *          for this by re-assigning the start and end Point coordinates.
 */
void TrackGenerator::recalibrateTracksToOrigin() {

  /* Recalibrate the tracks to the origin. */
  for (int a=0; a < _num_azim_2; a++) {
    for (int i=0; i < _num_tracks[a]; i++) {

      double x0 = _tracks[a][i].getStart()->getX();
      double y0 = _tracks[a][i].getStart()->getY();
      double x1 = _tracks[a][i].getEnd()->getX();
      double y1 = _tracks[a][i].getEnd()->getY();
      double new_x0 = x0 + _geometry->getMinX();
      double new_y0 = y0 + _geometry->getMinY();
      double new_x1 = x1 + _geometry->getMinX();
      double new_y1 = y1 + _geometry->getMinY();
      double phi = _tracks[a][i].getPhi();

      /* The z-coordinates don't need to be recalibrated to the origin. */
      double z0 = _tracks[a][i].getStart()->getZ();
      double z1 = _tracks[a][i].getEnd()->getZ();

      _tracks[a][i].setValues(new_x0, new_y0, z0, new_x1, new_y1, z1, phi);
      _tracks[a][i].setAzimAngleIndex(a);
    }
  }
}


/**
 * @brief Set the cycle index for each track in a PERIODIC or REFLECTIVE cycle.
 * @details The tracks can be separated into cycles as they traverse the
 *          geometry. It is important to set the periodic cycle indices for
 *          problems with periodic boundary conditions so transport sweeps can
 *          be performed in parallel over groups of tracks that do not transfer
 *          their flux into each other. In this method, all tracks are looped
 *          over; if the track's cycle index has not been set, the cycle index
 *          for the input boundaryType is incremented and the cycle index is set
 *          for that track and all tracks in its cycle.
 */
void TrackGenerator::initializeTrackCycleIndices(boundaryType bc) {

  if (bc != REFLECTIVE && bc != PERIODIC)
    log_printf(ERROR, "Cannot initialize Track cycle indices for boundaryType"
               " other than REFLECTIVE or PERIODIC");

  Track* track;
  int track_index;
  int next_i, next_a;
  bool fwd;

  /* Loop over all tracks */
  for (int a=0; a < _num_azim_2; a++) {
    for (int i=0; i < _num_tracks[a]; i++) {

      /* Get the current track */
      track = &_tracks[a][i];

      /* Set the track indices for PERIODIC boundaryType */
      if (bc == PERIODIC) {

        /* Check if track index has been set */
        if (track->getPeriodicTrackIndex() == -1) {

          /* Initialize the track index counter */
          track_index = 0;
          next_i = i;

          /* Set the periodic track indexes for all tracks in periodic cycle */
          while (track->getPeriodicTrackIndex() == -1) {

            /* Set the track periodic cycle */
            track->setPeriodicTrackIndex(track_index);

            /* Get the xy index of the next track in cycle */
            if (next_i < _num_y[a])
              next_i +=_num_x[a];
            else
              next_i -=_num_y[a];

            /* Set the next track in cycle */
            track = &_tracks[a][next_i];

            /* Increment index counter */
            track_index++;
          }
        }
      }

      /* Set the track indices for REFLECTIVE boundaryType */
      else {

        /* Check if track index has been set */
        if (track->getReflectiveTrackIndex() == -1) {

          /* Initialize the track index counter */
          track_index = 0;
          next_i = i;
          next_a = a;
          fwd = true;

          /* Set the reflective track indexes for all tracks in reflective
           * cycle */
          while (track->getReflectiveTrackIndex() == -1) {

            /* Set the track reflective cycle */
            track->setReflectiveTrackIndex(track_index);

            /* Set the azimuthal angle of the next track in the cycle */
            next_a = _num_azim_2 - next_a - 1;

            /* Set the xy index and direction of the next track in the cycle */
            if (fwd) {
              if (next_i < _num_y[a]) {
                next_i = next_i + _num_x[a];
                fwd = true;
              }
              else {
                next_i = _num_x[a] + 2*_num_y[a] - next_i - 1;
                fwd = false;
              }
            }

            else {
              if (next_i < _num_x[a]) {
                next_i = _num_x[a] - next_i - 1;
                fwd = true;
              }
              else{
                next_i = next_i - _num_x[a];
                fwd = false;
              }
            }

            /* Set the next track in cycle */
            track = &_tracks[next_a][next_i];

            /* Increment index counter */
            track_index++;
          }
        }
      }
    }
  }
}


/**
 * @brief Set the Track UIDs for all tracks and generate the 1D array of
 *        track pointers that separates the groups of tracks by parallel group.
 * @details The Solver requires the tracks to be separated into groups of tracks
 *          that can be looped over in parallel without data races. This method
 *          creates a 1D array of Track pointers where the tracks are arranged
 *          by parallel group. If the geometry has a periodic BC, 6 periodic
 *          groups are created; otherwise, 2 periodic groups are created. The
 *          Track UIDs are also set to their index in the tracks by periodic
 *          group array.
 */
void TrackGenerator::initializeTrackUids() {

  Track* track;
  bool periodic;
  int uid = 0;
  int azim_group_id, periodic_group_id;
  int track_azim_group_id, track_periodic_group_id;
  int track_periodic_index;
  int num_tracks;

  /* If periodic boundary conditions are present, set the number of parallel
   * track groups to 6; else, set the number of parallel track groups to 2. Also
   * set the periodic boolean indicating whether periodic bcs are present. */
  if (_geometry->getMinXBoundaryType() == PERIODIC ||
      _geometry->getMinYBoundaryType() == PERIODIC) {
    _num_parallel_track_groups = 6;
    periodic = true;
  }
  else {
    _num_parallel_track_groups = 2;
    periodic = false;
  }

  /* Allocate memory for the num tracks by parallel group array */
  _num_tracks_by_parallel_group = new int[_num_parallel_track_groups];

  /* Allocate memory for the tracks by parallel group array */
  _tracks_by_parallel_group = new Track*[getNumTracks()];

  /* Loop over the parallel track groups */
  for (int g = 0; g < _num_parallel_track_groups; g++) {

    /* Initialize the number of tracks counter */
    num_tracks = 0;

    /* Set the azimuthal and periodic group ids */
    azim_group_id = g % 2;
    periodic_group_id = g / 2;

    /* Loop over all tracks and add all tracks belonging to the current
     * parallel group to the tracks by parallel groups array */
    for (int a=0; a < _num_azim_2; a++) {
      for (int i=0; i < _num_tracks[a]; i++) {

        /* Get current track */
        track = &_tracks[a][i];

        /* Get the track azim group id */
        track_azim_group_id = a / (_num_azim_2 / 2);

        /* Get the track periodic group id */
        if (periodic) {
          track_periodic_index = track->getPeriodicTrackIndex();

          /* If the track is the first track in periodic cycle, assign it a
           * periodic group id of 0 */
          if (track_periodic_index == 0)
            track_periodic_group_id = 0;

          /* If the track has an odd periodic cycle index, assign it a periodic
           * group id of 1 */
          else if (track_periodic_index % 2 == 1)
            track_periodic_group_id = 1;

          /* Else the track must have an even periodic cycle index and not be
           * the first track in the periodic cycle; assign it a periodic group
           * id of 2 */
          else
            track_periodic_group_id = 2;
        }

        /* If the geometry does not have any periodic BCs, assign the track a
         * periodic group id of 0 */
        else
          track_periodic_group_id = 0;

        /* Check if track has current azim_group_id and periodic_group_id */
        if (azim_group_id == track_azim_group_id &&
            periodic_group_id == track_periodic_group_id) {
          track->setUid(uid);
          _tracks_by_parallel_group[uid] = track;
          uid++;
          num_tracks++;
        }
      }
    }

    /* Set the number of tracks in this parallel group */
    _num_tracks_by_parallel_group[g] = num_tracks;
  }
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
  resetFSRVolumes();
  FP_PRECISION* fsr_volumes = getFSRVolumes();

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
 * @brief This helper method for TrackGenerator::generateTracks() finds the end
 *        Point of a Track with a defined start Point and an angle from x-axis.
 * @details This function does not return a value but instead saves the x/y
 *          coordinates of the end Point directly within the Track's end Point.
 * @param start pointer to the Track start Point
 * @param end pointer to a Point to store the end Point coordinates
 * @param phi the azimuthal angle
 * @param width_x the x-width of the Geometry (cm)
 * @param width_y the y-width of the Geometry (cm)
 */
void TrackGenerator::computeEndPoint(Point* start, Point* end,
                                     const double phi, const double width_x,
                                     const double width_y) {

  double m = sin(phi) / cos(phi);             /* slope */
  double xin = start->getX();                 /* x-coord */
  double yin = start->getY();                 /* y-coord */
  double zin = start->getZ();                 /* z-coord */

  /* Allocate memory for the possible intersection points */
  Point *points = new Point[4];

  /* Determine all possible Points */
  points[0].setCoords(0, yin - m * xin, zin);
  points[1].setCoords(width_x, yin + m * (width_x - xin), zin);
  points[2].setCoords(xin - yin / m, 0, zin);
  points[3].setCoords(xin - (yin - width_y) / m, width_y, zin);

  /* For each of the possible intersection Points */
  for (int i = 0; i < 4; i++) {
    /* neglect the trivial Point (xin, yin) */
    if (points[i].getX() == xin && points[i].getY() == yin) { }

    /* The Point to return will be within the bounds of the cell */
    else if (points[i].getX() >= 0 && points[i].getX() <= width_x
             && points[i].getY() >= 0 && points[i].getY() <= width_y) {
      end->setCoords(points[i].getX(), points[i].getY(), zin);
    }
  }

    delete [] points;

    return;
}


/**
 * @brief Initializes boundary conditions for each Track.
 * @details Sets boundary conditions by setting the incoming and outgoing Tracks
 *          for each Track using a special indexing scheme into the 2D jagged
 *          array of Tracks.
 */
void TrackGenerator::initializeBoundaryConditions() {

  log_printf(INFO, "Initializing Track boundary conditions...");

  /* Check for symmetry of periodic boundary conditions */
  if ((_geometry->getMinXBoundaryType() == PERIODIC &&
       _geometry->getMaxXBoundaryType() != PERIODIC) ||
      (_geometry->getMinXBoundaryType() != PERIODIC &&
       _geometry->getMaxXBoundaryType() == PERIODIC))
    log_printf(ERROR, "Cannot create tracks with only one x boundary"
               " set to PERIODIC");
  else if ((_geometry->getMinYBoundaryType() == PERIODIC &&
            _geometry->getMaxYBoundaryType() != PERIODIC) ||
           (_geometry->getMinYBoundaryType() != PERIODIC &&
            _geometry->getMaxYBoundaryType() == PERIODIC))
    log_printf(ERROR, "Cannot create tracks with only one y boundary"
               " set to PERIODIC");
  Track* track;
  int ic;

  /* Loop over the all the tracks and set the incoming and outgoing tracks
   * and incoming and outgoing boundary conditions. */
  for (int i=0; i < _num_azim_2; i++) {
    ic = _num_azim_2 - i - 1;
    for (int j=0; j < _num_tracks[i]; j++) {

      /* Get current track */
      track = &_tracks[i][j];

      /* Set boundary conditions for tracks in [0, PI/2] */
      if (i < _num_azim_2/2) {
        if (j < _num_y[i])
          track->setBCOut(_geometry->getMaxXBoundaryType());
        else
          track->setBCOut(_geometry->getMaxYBoundaryType());

        if (j < _num_x[i])
          track->setBCIn(_geometry->getMinYBoundaryType());
        else
          track->setBCIn(_geometry->getMinXBoundaryType());
      }

      /* Set boundary conditions for tracks in [PI/2, PI] */
      else {
        if (j < _num_y[i])
          track->setBCOut(_geometry->getMinXBoundaryType());
        else
          track->setBCOut(_geometry->getMaxYBoundaryType());

        if (j < _num_x[i])
          track->setBCIn(_geometry->getMinYBoundaryType());
        else
          track->setBCIn(_geometry->getMaxXBoundaryType());
      }

      /* Set connecting tracks in forward direction */
      if (j < _num_y[i]) {
        track->setNextOut(false);
        if (track->getBCOut() == PERIODIC)
          track->setTrackOut(&_tracks[i][j + _num_x[i]]);
        else
          track->setTrackOut(&_tracks[ic][j + _num_x[i]]);
      }
      else {
        if (track->getBCOut() == PERIODIC) {
          track->setNextOut(false);
          track->setTrackOut(&_tracks[i][j - _num_y[i]]);
        }
        else {
          track->setNextOut(true);
          track->setTrackOut(&_tracks[ic][_num_x[i] + 2*_num_y[i] - j - 1]);
        }
      }

      /* Set connecting tracks in backward direction */
      if (j < _num_x[i]) {
        if (track->getBCIn() == PERIODIC) {
          track->setNextIn(true);
          track->setTrackIn(&_tracks[i][j + _num_y[i]]);
        }
        else {
          track->setNextIn(false);
          track->setTrackIn(&_tracks[ic][_num_x[i] - j - 1]);
        }
      }
      else {
        track->setNextIn(true);
        if (track->getBCIn() == PERIODIC)
          track->setTrackIn(&_tracks[i][j - _num_x[i]]);
        else
          track->setTrackIn(&_tracks[ic][j - _num_x[i]]);
      }
    }
  }

  return;
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize() {

  log_printf(NORMAL, "Ray tracing for track segmentation...");

  /* This section loops over all Track and segmentizes each one if the
   * Tracks were not read in from an input file */
  if (!_use_input_file) {

    /* Loop over all Tracks */
    for (int i=0; i < _num_azim_2; i++) {
#pragma omp parallel for
      for (int j=0; j < _num_tracks[i]; j++) {
        _geometry->segmentize(&_tracks[i][j]);
      }
    }
  }

  _geometry->initializeFSRVectors();

  _contains_segments = true;

  return;
}


/**
 * @brief Writes all Track and segment data to a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 */
void TrackGenerator::dumpTracksToFile() {

  /* Check whether the segments should be dumped */
  if (!_dump_segments)
    return;

  log_printf(NORMAL, "Dumping segments to file...");

  if (!containsSegments())
    log_printf(ERROR, "Unable to dump Segments to a file since no Segments "
               "have been generated for %d azimuthal angles and %f track "
               "spacing", 2*_num_azim_2, _azim_spacing);

  FILE* out;
  out = fopen(_tracks_filename.c_str(), "w");

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
  dump_segments.execute();

  /* Get FSR vector maps */
  Cmfd* cmfd = _geometry->getCmfd();
  ParallelHashMap<std::string, fsr_data*>& FSR_keys_map =
      _geometry->getFSRKeysMap();
  std::vector<std::string>& FSRs_to_keys = _geometry->getFSRsToKeys();
  std::string fsr_key;
  int fsr_id;
  double x, y, z;

  /* Write number of FSRs */
  int num_FSRs = _geometry->getNumFSRs();
  fwrite(&num_FSRs, sizeof(int), 1, out);

  /* Write FSR vector maps to file */
  std::string* fsr_key_list = FSR_keys_map.keys();
  fsr_data** fsr_data_list = FSR_keys_map.values();
  for (int i=0; i < num_FSRs; i++) {

    /* Write key to file from FSR_keys_map */
    fsr_key = fsr_key_list[i];
    string_length = fsr_key.length() + 1;
    fwrite(&string_length, sizeof(int), 1, out);
    fwrite(fsr_key.c_str(), sizeof(char)*string_length, 1, out);

    /* Write data to file from FSR_keys_map */
    fsr_id = fsr_data_list[i]->_fsr_id;
    x = fsr_data_list[i]->_point->getX();
    y = fsr_data_list[i]->_point->getY();
    z = fsr_data_list[i]->_point->getZ();
    fwrite(&fsr_id, sizeof(int), 1, out);
    fwrite(&x, sizeof(double), 1, out);
    fwrite(&y, sizeof(double), 1, out);
    fwrite(&z, sizeof(double), 1, out);

    /* Write data to file from FSRs_to_keys */
    fsr_key = FSRs_to_keys.at(i);
    string_length = fsr_key.length() + 1;
    fwrite(&string_length, sizeof(int), 1, out);
    fwrite(fsr_key.c_str(), sizeof(char)*string_length, 1, out);
  }

  /* Write cmfd_fsrs vector of vectors to file */
  if (cmfd != NULL) {
    std::vector< std::vector<int> >* cell_fsrs = cmfd->getCellFSRs();
    std::vector<int>::iterator iter;
    int num_cells = cmfd->getNumCells();
    fwrite(&num_cells, sizeof(int), 1, out);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      num_FSRs = cell_fsrs->at(cell).size();
      fwrite(&num_FSRs, sizeof(int), 1, out);

      /* Loop over FSRs within cell */
      for (iter = cell_fsrs->at(cell).begin();
           iter != cell_fsrs->at(cell).end(); ++iter)
        fwrite(&(*iter), sizeof(int), 1, out);
    }
  }

  /* Delete key and value lists */
  delete [] fsr_key_list;
  delete [] fsr_data_list;

  /* Close the Track file */
  fclose(out);

  /* Inform other the TrackGenerator::generateTracks() method that it may
   * import ray tracing data from this file if it is called and the ray
   * tracing parameters have not changed */
  _use_input_file = true;
}


/**
 * @brief Reads Tracks in from a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 * @return true if able to read Tracks in from a file; false otherwise
 */
bool TrackGenerator::readTracksFromFile() {

  int ret;
  FILE* in;
  in = fopen(_tracks_filename.c_str(), "r");

  int string_length;

  /* Import Geometry metadata from the Track file */
  ret = fread(&string_length, sizeof(int), 1, in);
  char* geometry_to_string = new char[string_length];
  ret = fread(geometry_to_string, sizeof(char)*string_length, 1, in);

  /* Check if our Geometry is exactly the same as the Geometry in the
   * Track file for this number of azimuthal angles and track spacing */
  if (_geometry->toString().compare(std::string(geometry_to_string)) != 0)
    return false;

  delete [] geometry_to_string;

  log_printf(NORMAL, "Importing ray tracing data from file...");

  /* Load all segment data into Tracks */
  ReadSegments read_segments(this);
  read_segments.setInputFile(in);
  read_segments.execute();

  /* Create FSR vector maps */
  ParallelHashMap<std::string, fsr_data*>& FSR_keys_map =
      _geometry->getFSRKeysMap();
  std::vector<std::string>& FSRs_to_keys =
      _geometry->getFSRsToKeys();
  Cmfd* cmfd = _geometry->getCmfd();
  FSR_keys_map.clear();
  FSRs_to_keys.clear();
  int num_FSRs;
  std::string fsr_key;
  int fsr_key_id;
  double x, y, z;

  /* Get number of FSRs */
  ret = fread(&num_FSRs, sizeof(int), 1, in);

  /* Read FSR vector maps from file */
  for (int fsr_id=0; fsr_id < num_FSRs; fsr_id++) {

    /* Read key for FSR_keys_map */
    ret = fread(&string_length, sizeof(int), 1, in);
    char* char_buffer1 = new char[string_length];
    ret = fread(char_buffer1, sizeof(char)*string_length, 1, in);
    fsr_key = std::string(char_buffer1);

    /* Read data from file for FSR_keys_map */
    ret = fread(&fsr_key_id, sizeof(int), 1, in);
    ret = fread(&x, sizeof(double), 1, in);
    ret = fread(&y, sizeof(double), 1, in);
    ret = fread(&z, sizeof(double), 1, in);
    fsr_data* fsr = new fsr_data;
    fsr->_fsr_id = fsr_key_id;
    Point* point = new Point();
    point->setCoords(x,y,z);
    fsr->_point = point;
    FSR_keys_map.insert(fsr_key, fsr);

    /* Read data from file for FSR_to_keys */
    ret = fread(&string_length, sizeof(int), 1, in);
    char* char_buffer2 = new char[string_length];
    ret = fread(char_buffer2, sizeof(char)*string_length, 1, in);
    fsr_key = std::string(char_buffer2);
    FSRs_to_keys.push_back(fsr_key);
  }

  /* Read cmfd cell_fsrs vector of vectors from file */
  if (cmfd != NULL) {
    std::vector< std::vector<int> > cell_fsrs;
    int num_cells, fsr_id;
    ret = fread(&num_cells, sizeof(int), 1, in);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      cell_fsrs.push_back(std::vector<int>());
      ret = fread(&num_FSRs, sizeof(int), 1, in);

      /* Loop over FRSs within cell */
      for (int fsr = 0; fsr < num_FSRs; fsr++) {
        ret = fread(&fsr_id, sizeof(int), 1, in);
        cell_fsrs.at(cell).push_back(fsr_id);
      }
    }

    /* Set CMFD cell_fsrs vector of vectors */
    cmfd->setCellFSRs(&cell_fsrs);
  }

  /* Inform the rest of the class methods that Tracks have been initialized */
  if (ret)
    _contains_tracks = true;

  /* Close the Track file */
  fclose(in);

  return true;
}


/**
 * @brief Assign a correct volume for some FSR.
 * @details This routine adjusts the length of each track segment crossing
 *          a FSR such that the integrated volume is identical to the true
 *          volume assigned by the user.
 * @param fsr_id the ID of the FSR of interest
 * @param fsr_volume the correct FSR volume to use
 */
void TrackGenerator::correctFSRVolume(int fsr_id, FP_PRECISION fsr_volume) {

  if (!_contains_tracks)
    log_printf(ERROR, "Unable to correct FSR volume since "
	       "tracks have not yet been generated");

  /* Compute the current volume approximation for the flat source region */
  FP_PRECISION curr_volume = getFSRVolume(fsr_id);

  log_printf(INFO, "Correcting FSR %d volume from %f to %f",
             fsr_id, curr_volume, fsr_volume);

  int num_segments, azim_index;
  double dx_eff, d_eff;
  double volume, corr_factor;
  segment* curr_segment;
  segment* segments;

  /* Correct volume separately for each azimuthal angle */
  for (int i=0; i < _num_azim_2; i++) {

    /* Initialize volume to zero for this azimuthal angle */
    volume = 0;

    /* Compute effective track spacing for this azimuthal angle */
    dx_eff = (_geometry->getWidthX() / _num_x[i]);
    d_eff = (dx_eff * sin(_tracks[i][0].getPhi()));

    /* Compute the current estimated volume of the FSR for this angle */
#pragma omp parallel for  private(num_segments, segments, curr_segment) \
  reduction(+:volume)
    for (int j=0; j < _num_tracks[i]; j++) {

      num_segments = _tracks[i][j].getNumSegments();
      segments = _tracks[i][j].getSegments();

      for (int s=0; s < num_segments; s++) {
        curr_segment = &segments[s];
        if (curr_segment->_region_id == fsr_id)
          volume += curr_segment->_length * d_eff;
      }
    }

    /* Compute correction factor to the volume */
    corr_factor = fsr_volume / volume;

    log_printf(DEBUG, "Volume correction factor for FSR %d and azim "
               "angle %d is %f", fsr_id, i, corr_factor);

    /* Correct the length of each segment which crosses the FSR */
#pragma omp parallel for  private(num_segments, segments, curr_segment)
    for (int j=0; j < _num_tracks[i]; j++) {

      num_segments = _tracks[i][j].getNumSegments();
      segments = _tracks[i][j].getSegments();

      for (int s=0; s < num_segments; s++) {
        curr_segment = &segments[s];
        if (curr_segment->_region_id == fsr_id)
          curr_segment->_length *= corr_factor;
      }
    }
  }
}


/**
 * @brief Generates the numerical centroids of the FSRs.
 * @details This routine generates the numerical centroids of the FSRs
 *          by weighting the average x and y values of each segment in the
 *          FSR by the segment's length and azimuthal weight. The numerical
 *          centroid fomula can be found in R. Ferrer et. al. "Linear Source
 *          Approximation in CASMO 5", PHYSOR 2012.
 */
void TrackGenerator::generateFSRCentroids() {

  /* Create temporary array of centroids and initialize to origin */
  int num_FSRs = _geometry->getNumFSRs();
  Point** centroids = new Point*[num_FSRs];
  for (int r=0; r < num_FSRs; r++) {
    centroids[r] = new Point();
    centroids[r]->setCoords(0.0, 0.0, 0.0);
  }

  /* Generate FSR centroids by looping over all Tracks */
  CentroidGenerator centroid_generator(this);
  centroid_generator.setCentroids(centroids);
  centroid_generator.execute();

  /* Set the centroid for the FSR */
  for (int r=0; r < num_FSRs; r++)
    _geometry->setFSRCentroid(r, centroids[r]);

  /* Delete temporary array of FSR volumes and centroids */
  delete [] centroids;
}


/**
 * @brief Sets the max optical path length of segments during computation
 * @param tau maximum optical path length
 */
void TrackGenerator::setMaxOpticalLength(FP_PRECISION tau) {
  _max_optical_length = tau;
}


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
 * @brief Initialize track segments with pointers to FSR Materials.
 * @details This is called by the Solver at simulation time. This
 *          initialization is necessary since Materials in each FSR
 *          may be interchanged by the user in between different
 *          simulations. This method links each segment and fsr_data
 *          struct with the current Material found in each FSR.
 */
void TrackGenerator::initializeSegments() {

  if (!_contains_segments)
    log_printf(ERROR, "Unable to initialize segments since "
	       "segments have not yet been generated");

  /* Get all of the Materials from the Geometry */
  std::map<int, Material*> materials = _geometry->getAllMaterials();

  /* Get the mappings of FSR to keys to fsr_data to update Materials */
  ParallelHashMap<std::string, fsr_data*>& FSR_keys_map =
      _geometry->getFSRKeysMap();
  std::vector<std::string>& FSRs_to_keys = _geometry->getFSRsToKeys();

#pragma omp parallel
  {

    int region_id, mat_id;
    segment* curr_segment;
    Material* mat;

    /* Set the Material for each FSR */
#pragma omp for
    for (int r=0; r < _geometry->getNumFSRs(); r++) {
      mat = _geometry->findFSRMaterial(r);
      FSR_keys_map.at(FSRs_to_keys.at(r))->_mat_id = mat->getId();
    }

    /* Set the Material for each segment */
    for (int i=0; i < _num_azim_2; i++) {
#pragma omp for
      for (int j=0; j < _num_tracks[i]; j++) {
        for (int s=0; s < _tracks[i][j].getNumSegments(); s++) {
          curr_segment = _tracks[i][j].getSegment(s);
          region_id = curr_segment->_region_id;
          mat_id = FSR_keys_map.at(FSRs_to_keys.at(region_id))->_mat_id;
          curr_segment->_material = materials[mat_id];
        }
      }
    }
  }
}


/**
 * @brief Returns the azimuthal angle for a given azimuthal angle index.
 * @param the azimuthal angle index.
 * @return the desired azimuthal angle.
 */
double TrackGenerator::getPhi(int azim) {

  if (!_contains_tracks)
    log_printf(ERROR, "Unable to get Phi since the TrackGenerator does not"
               " contain tracks.");

  if (azim < 0 || azim >= 2*_num_azim_2)
    log_printf(ERROR, "Unable to get phi for azimuthal angle %d since there"
               "are  %d azimuthal angles", azim, 2*_num_azim_2);

  return _quadrature->getPhi(azim);
}


/**
 * @brief Deletes the Timer's timing entries for each timed code section
 *        code in the loop over track for ray tracing.
 */
void TrackGenerator::clearTimerSplits() {
  _timer->clearSplits();
}


/**
 * @brief Prints a report of the timing statistics to the console.
 */
void TrackGenerator::printTimerReport() {

  std::string msg_string;

  log_printf(TITLE, "TRACKGENERATOR TIMING REPORT");

  /* Get the total runtime */
  double tot_time = _timer->getSplit("Total time");
  msg_string = "Total time for ray tracing";
  msg_string.resize(REPORT_WIDTH, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);

  /* Time per segment */
  double time_per_segment = (tot_time / getNumSegments());
  msg_string = "Time per segment";
  msg_string.resize(REPORT_WIDTH, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_segment);

  set_separator_character('-');
  log_printf(SEPARATOR, "-");
}


/**
 * @brief Resets the TrackGenerator to not contain tracks or segments
 */
void TrackGenerator::resetStatus() {
  _contains_tracks = false;
  _use_input_file = false;
  _tracks_filename = "";
}
