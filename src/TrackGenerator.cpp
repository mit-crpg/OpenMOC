#include "TrackGenerator.h"
#include "TrackTraversingAlgorithms.h"
#include <iomanip>

/**
 * @brief Constructor for the TrackGenerator assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param spacing track spacing (cm)
 */
TrackGenerator::TrackGenerator(Geometry* geometry, int num_azim, int num_polar,
                               double azim_spacing) {

  setNumThreads(1);
  _geometry = geometry;
  setNumAzim(num_azim);
  setNumPolar(num_polar);
  setDesiredAzimSpacing(azim_spacing);
  _contains_2D_tracks = false;
  _contains_2D_segments = false;
  _quadrature = NULL;
  _z_coord = 0.0;
  _segment_formation = EXPLICIT_2D;
  _max_optical_length = std::numeric_limits<FP_PRECISION>::max();
  _max_num_segments = 0;
  _FSR_volumes = NULL;
  _dump_segments = true;
  _FSR_locks = NULL;
  _tracks_2D_array = NULL;
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

    /* Delete track laydown information */
    delete [] _num_x;
    delete [] _num_y;
    delete [] _cycles_per_azim;
    delete [] _tracks_per_cycle;
    delete [] _cycle_length;
  }

  if (_FSR_locks != NULL)
    delete [] _FSR_locks;

  if (_FSR_volumes != NULL)
    delete [] _FSR_volumes;
}


/**
 * @brief Return the number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @return the number of azimuthal angles in \f$ 2\pi \f$
 */
int TrackGenerator::getNumAzim() {
  return _num_azim;
}


/**
 * @brief Return the number of polar angles in \f$ [0, \pi] \f$
 * @return the number of polar angles in \f$ \pi \f$
 */
int TrackGenerator::getNumPolar() {
  return _num_polar;
}


/**
 * @brief Return the track azimuthal spacing (cm).
 * @ditails This will return the user-specified track spacing and NOT the
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
 * @brief Return the array used to store the FSR volumes
 * @return _FSR_volumes the FSR volumes array indexed by FSR ID
 */
FP_PRECISION* TrackGenerator::getFSRVolumesBuffer() {
#pragma omp critical
  {
    if (_FSR_volumes == NULL) {
      int num_FSRs = _geometry->getNumFSRs();
      _FSR_volumes = new FP_PRECISION[num_FSRs];
      memset(_FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));
    }
  }

  return _FSR_volumes;
}


/**
 * @brief Return the total number of Tracks across the Geometry.
 * @return the total number of Tracks
 */
int TrackGenerator::getNumTracks() {
  return getNum2DTracks();
}


/**
 * @brief Return the total number of 2D Tracks across the Geometry.
 * @return the total number of 2D Tracks
 */
int TrackGenerator::getNum2DTracks() {

  int num_2D_tracks = 0;

  for (int a=0; a < _num_azim/2; a++)
    num_2D_tracks += _num_x[a] + _num_y[a];

  return num_2D_tracks;
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator::getNumSegments() {
  return getNum2DSegments();
}


/**
 * @brief Return the total number of 2D Track segments across the Geometry.
 * @return the total number of 2D Track segments
 */
int TrackGenerator::getNum2DSegments() {

  if (!TrackGenerator::containsSegments())
    log_printf(ERROR, "Cannot get the number of 2D segments since they "
               "have not been generated.");

  int num_2D_segments = 0;

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
 * @breif Calculates and returns the maximum optcial length for any segment
 *        in the Geomtry.
 * @details The _max_optical_length value is recomputed, updated, and returned.
 *          This value determines the when segments must be split during ray
 *          tracing.
 * @return _max_optical_length the maximum optical length of any segment in the
 *         Geometry
 */
FP_PRECISION TrackGenerator::getMaxOpticalLength() {
  MaxOpticalLength update_max_optical_length(this);
  update_max_optical_length.execute();
  return _max_optical_length;
}


/**
 * @brief Returns the maximum number of segments along a single track
 * @details The TrackGenerator::countSegments routine must be called before
 *          this function will return a correct value
 * @return the maximum number of segments
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
 * @brief Returns an array of the number of 2D Tracks in a cycle
 * @details The number of Tracks in a 2D cycle depends on the azimuthal angle
 *          index. This function returns an array of the number of 2D Tracks in
 *          each cycle, indexed by azimuthal anlge index. NOTE: all 2D cycles
 *          with the same azimuthal angle have the same number of Tracks.
 * @return the array of cycle lengths
 */
int* TrackGenerator::getTracksPerCycle() {
  return _tracks_per_cycle;
}


/**
 * @brief Returns an array describing the number of cycles per azimuthal angle
 * @details An array of the number of cycles per azimuthal angle is returned,
 *          indexed by azimuthal index.
 * @return the number of cycles per azimuthal angle
 */
int* TrackGenerator::getCyclesPerAzim() {
  return _cycles_per_azim;
}


/**
 * @brief Returns the number of 2D Tracks in a cycle for a given azimuthal
 *        angle index
 * @details The number of Tracks in a 2D cycle depends on the azimuthal angle
 *          index. This function returns the number of 2D Tracks in for a cycle
 *          with a given azimuthal angle index.
 * @param azim the azimuthal angle index in the first quadrant
 * @return the number of 2D Tracks in the cycle
 */
double TrackGenerator::getCycleLength(int azim) {
  if (azim > _num_azim/4)
    log_printf(ERROR, "Azimuthal angle index %d refers to an angle that is "
                      "not in the first quadrant", azim);

  return _cycle_length[azim];
}


/**
 * @brief Returns the number of 2D Tracks in the x-direction for a given
 *        azimuthal angle index
 * @param azim the azimuthal angle index
 * @return the number of 2D Tracks in the x-direction of the Geometry
 */
int TrackGenerator::getNumX(int azim) {
  return _num_x[azim];
}


/**
 * @brief Returns the number of 2D Tracks in the y-direction for a given
 *        azimuthal angle index
 * @param azim the azimuthal angle index
 * @return the number of 2D Tracks in the y-direction of the Geometry
 */
int TrackGenerator::getNumY(int azim) {
  return _num_y[azim];
}


/**
 * @brief FSR volumes are coppied to an array input by the user
 * @param out_volumes The array to which FSR volumes are coppied
 * @param num_fsrs The number of FSR volumes to copy. The first num_fsrs
 *        volumes stored in the FSR volumes array are coppied.
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
  int num_FSRs = _geometry->getNumFSRs();
  if (_FSR_volumes != NULL)
    memset(_FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));

  /* Create volume calculator and calculate new FSR volumes */
  VolumeCalculator volume_calculator(this);
  volume_calculator.execute();

  /* Check to ensure all FSRs are crossed by at least one track */
  for (int i=0; i < num_FSRs; i++) {
    if (_FSR_volumes[i] == 0.0) {
      log_printf(NORMAL, "Zero volume calculated for FSR %d, point (%f, %f, %f)",
                 i, _geometry->getFSRPoint(i)->getX(), _geometry->getFSRPoint(i)->getY(),
                 _geometry->getFSRPoint(i)->getZ());
      log_printf(ERROR, "Zero volume calculated in an FSR region since no "
               "track traversed the FSR. Use a finer track laydown to ensure "
               "every FSR is traversed.");
    }
  }

  return _FSR_volumes;
}


/**
 * @brief Returns the volume of an FSR.
 * @param fsr_id the ID for the FSR of interest
 * @return the FSR volume
 */
FP_PRECISION TrackGenerator::getFSRVolume(int fsr_id) {

  if (_FSR_volumes == NULL)
    log_printf(ERROR, "Unable to get the FSR volume since FSR volumes "
               "have not yet been generated");

  else if (fsr_id < 0 || fsr_id > _geometry->getNumFSRs())
    log_printf(ERROR, "Unable to get the volume for FSR %d since the FSR IDs "
               "lie in the range (0, %d)", fsr_id, _geometry->getNumFSRs());
  return _FSR_volumes[fsr_id];
}


/**
 * @brief Returns the z-coord of the radial plane used in 2D calcualtions
 * @return the z-coord of the 2D calculation
 */
double TrackGenerator::getZCoord() {
  return _z_coord;
}


/**
 * @brief Returns the Quadrature object
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

  _num_threads = num_threads;

  /* Set the number of threads for OpenMP */
  omp_set_num_threads(_num_threads);
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
 * @brief Set the number of polar angles in \f$ [0, \pi] \f$.
 * @param num_polar the number of polar angles in \f$ \pi \f$
 */
void TrackGenerator::setNumPolar(int num_polar) {

  if (num_polar < 0)
    log_printf(ERROR, "Unable to set a negative number of polar angles "
               "%d for the TrackGenerator.", num_polar);

  if (num_polar % 2 != 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d for the "
               "TrackGenerator since it is not a multiple of 2", num_polar);

  _num_polar = num_polar;
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
  resetStatus();
}


/**
 * @brief Sets the z-coord of the raidal plane used in 2D calculations
 * @param z_coord the z-coord of the radial plane
 */
void TrackGenerator::setZCoord(double z_coord) {
  _z_coord = z_coord;
}


/**
 * @brief sets the Quadrature used for integrating the MOC equations
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
 *        geometry for it's current segmentation type.
 * @return true if the TrackGenerator conatains segments; false otherwise
 */
bool TrackGenerator::containsSegments() {
  return _contains_2D_segments;
}


/**
 * @brief Fills an array with the x,y coordinates for each Track.
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
 * @param coords an array of coords of length 4 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator::retrieveTrackCoords(double* coords, int num_tracks) {
  retrieve2DTrackCoords(coords, num_tracks);
}


/**
 * @brief Fills an array with the x,y coordinates for each Track.
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
 * @param coords an array of coords of length 4 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator::retrieve2DTrackCoords(double* coords, int num_tracks) {

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
      coords[counter+2] = _tracks_2D[a][i].getEnd()->getX();
      coords[counter+3] = _tracks_2D[a][i].getEnd()->getY();

      counter += 4;
    }
  }

  return;
}


/**
 * @brief Fills an array with the x,y coordinates and the periodic cycle ID
 *        for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords(num_tracks*5)
 * @endcode
 *
 * @param coords an array of coords of length 5 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator::retrieve2DPeriodicCycleCoords(double* coords,
                                                   int num_tracks) {

  if (num_tracks != NUM_VALUES_PER_RETRIEVED_TRACK * getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the 2D Track periodic cycle "
               "coordinates since the TrackGenerator contains %d Tracks with "
               "%d coordinates but an array of length %d was input",
               getNum2DTracks(), NUM_VALUES_PER_RETRIEVED_TRACK *
               getNum2DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      coords[counter]   = _tracks_2D[a][i].getStart()->getX();
      coords[counter+1] = _tracks_2D[a][i].getStart()->getY();
      coords[counter+2] = _tracks_2D[a][i].getEnd()->getX();
      coords[counter+3] = _tracks_2D[a][i].getEnd()->getY();
      coords[counter+4] = _tracks_2D[a][i].getPeriodicCycleId();

      counter += 5;
    }
  }
}


/**
 * @brief Fills an array with the x,y coordinates and the reflective cycle ID
 *        for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords(num_tracks*5)
 * @endcode
 *
 * @param coords an array of coords of length 5 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator::retrieve2DReflectiveCycleCoords(double* coords,
                                                     int num_tracks) {

  if (num_tracks != NUM_VALUES_PER_RETRIEVED_TRACK * getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the 2D Track reflective cycle "
               "coordinates since the TrackGenerator contains %d Tracks with "
               "%d coordinates but an array of length %d was input",
               getNum2DTracks(), NUM_VALUES_PER_RETRIEVED_TRACK *
               getNum2DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      coords[counter]   = _tracks_2D[a][i].getStart()->getX();
      coords[counter+1] = _tracks_2D[a][i].getStart()->getY();
      coords[counter+2] = _tracks_2D[a][i].getEnd()->getX();
      coords[counter+3] = _tracks_2D[a][i].getEnd()->getY();
      coords[counter+4] = _tracks_2D[a][i].getReflectiveCycleId();

      counter += 5;
    }
  }
}


/**
 * @brief Fills an array with the x,y coordinates for each Track segment.
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
 * @param coords an array of coords of length 5 times the number of segments
 * @param num_segments the total number of Track segments
 */
void TrackGenerator::retrieveSegmentCoords(double* coords, int num_segments) {
  retrieve2DSegmentCoords(coords, num_segments);
}


/**
 * @brief Fills an array with the x,y coordinates for each Track segment.
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
 * @param coords an array of coords of length 5 times the number of segments
 * @param num_segments the total number of Track segments
 */
void TrackGenerator::retrieve2DSegmentCoords(double* coords, int num_segments) {

  if (num_segments != NUM_VALUES_PER_RETRIEVED_SEGMENT * getNum2DSegments())
    log_printf(ERROR, "Unable to retrieve the Track segment coordinates since "
               "the TrackGenerator contains %d segments with %d coordinates "
               "but an array of length %d was input",
               getNum2DSegments(), NUM_VALUES_PER_RETRIEVED_SEGMENT *
               getNum2DSegments(), num_segments);

  segment* curr_segment = NULL;
  double x0, x1, y0, y1;
  double phi;
  segment* segments;

  int counter = 0;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      x0    = _tracks_2D[a][i].getStart()->getX();
      y0    = _tracks_2D[a][i].getStart()->getY();
      phi   = _tracks_2D[a][i].getPhi();

      segments = _tracks_2D[a][i].getSegments();

      for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++) {
        curr_segment = &segments[s];

        coords[counter] = curr_segment->_region_id;

        coords[counter+1] = x0;
        coords[counter+2] = y0;

        x1 = x0 + cos(phi) * curr_segment->_length;
        y1 = y0 + sin(phi) * curr_segment->_length;

        coords[counter+3] = x1;
        coords[counter+4] = y1;

        x0 = x1;
        y0 = y1;

        counter += 5;
      }
    }
  }

  return;
}


/**
 * @brief Checks the boundary conditions for all 2D surfaces for inconsistent
 *        periodic boundary conditions
 */
void TrackGenerator::checkBoundaryConditions() {

  /* Check X and Y boundaries for consistency */
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

  /* Check for correct track method if a PERIODIC bc is present */
  if (_geometry->getMinXBoundaryType() == PERIODIC ||
      _geometry->getMinYBoundaryType() == PERIODIC)

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

  /* Generate Tracks, perform ray tracing across the geometry, and store
   * the data to a Track file */
  try {

    /* Create default quadrature set if user one has not been set */
    if (_quadrature == NULL)
      initializeDefaultQuadrature();

    /* Initialize the quadrature set */
    _quadrature->setNumPolarAngles(_num_polar);
    _quadrature->setNumAzimAngles(_num_azim);
    _quadrature->initialize();

    /* Check periodic BCs for symmetry */
    checkBoundaryConditions();

    /* Lay down Tracks accross the Geometry */
    if (_geometry == NULL)
    log_printf(ERROR, "Unable to lay down Tracks since no Geometry "
               "has been set for the TrackGenerator");

    /* Initialize the Tracks */
    initializeTracks();

    /* Recalibrate the Tracks back to the geometry origin */
    recalibrateTracksToOrigin();

    /* Initialize the 1D array of Tracks for all Tracks */
    initializeTracksArray();

    /* Initialize the track file directory and read in tracks if they exist */
    initializeTrackFileDirectory();

    /* If track file not present, generate segments */
    if (_use_input_file == false) {

      /* Segmentize the tracks */
      segmentize();
      if (_segment_formation == EXPLICIT_2D ||
          _segment_formation == EXPLICIT_3D)
        dumpSegmentsToFile();
    }

    /* Allocate array of mutex locks for each FSR */
    int num_FSRs = _geometry->getNumFSRs();
    _FSR_locks = new omp_lock_t[num_FSRs];

    /* Loop over all FSRs to initialize OpenMP locks */
#pragma omp parallel for schedule(guided)
    for (int r=0; r < num_FSRs; r++)
      omp_init_lock(&_FSR_locks[r]);

    /* Precompute the quadrature weights */
    _quadrature->precomputeWeights(_segment_formation != EXPLICIT_2D);

  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to allocate memory needed to generate "
               "Tracks. Backtrace:\n%s", e.what());
  }

  return;
}


/**
 * @brief Allocates a new Quadrature with the default Quadrature
 * @details The defualt quadrature for 2D calculations is the TY quadrature
 */
void TrackGenerator::initializeDefaultQuadrature() {
  if (_quadrature != NULL)
    delete _quadrature;
  _quadrature = new TYPolarQuad();
}


/**
 * @brief calcualtes the least common multiple of two numbers a and b
 * @param first number a
 * @param second nuber b (order does not matter)
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

  /* Make sure that the width and height of the Geometry are nonzero */
  if (_geometry->getWidthX() <= 0 || _geometry->getWidthY() <= 0)
    log_printf(ERROR, "The total height and width of the Geometry must"
               " be nonzero for Track generation. Create a CellFill which "
               "is filled by the entire geometry and bounded by XPlanes "
               "and YPlanes to enable the Geometry to determine the "
               "total width and height of the model.");

  log_printf(NORMAL, "Initializing 2D tracks...");

  /* Allocate memory for arrays */
  _tracks_per_cycle = new int[_num_azim/4];
  _cycles_per_azim  = new int[_num_azim/4];
  _tracks_2D        = new Track*[_num_azim/2];
  _num_x            = new int[_num_azim/2];
  _num_y            = new int[_num_azim/2];
  _cycle_length     = new double[_num_azim/4];
  _num_2D_tracks    = 0;

  double x1, x2, y1, y2;
  double phi;
  double width  = _geometry->getWidthX();
  double height = _geometry->getWidthY();
  double dx_eff[_num_azim/2];
  double dy_eff[_num_azim/2];

  /* Determine angular quadrature and track spacing */
  for (int a = 0; a < _num_azim/4; a++) {

    /* Get the desired azimuthal angle */
    phi = _quadrature->getPhi(a);

    /* The number of intersections with x,y-axes */
    _num_x[a] = (int) (fabs(width / _azim_spacing * sin(phi))) + 1;
    _num_y[a] = (int) (fabs(height / _azim_spacing * cos(phi))) + 1;

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

    /* The length of all tracks in a 2D cycle */
    _cycle_length[a] = dx_eff[a] / cos(_quadrature->getPhi(a)) *
      leastCommonMultiple(2 * _num_x[a], 2 * height /
                          (tan(_quadrature->getPhi(a)) * dx_eff[a]));

    /* Get the number of tracks per cycle */
    _tracks_per_cycle[a] = (int)
      (round(_cycle_length[a] * sin(_quadrature->getPhi(a)) / width) +
       round(_cycle_length[a] * cos(_quadrature->getPhi(a)) / height));

    /* Compute the number of cycles */
    _cycles_per_azim[a] = (_num_x[a] + _num_y[a]) * 2 / _tracks_per_cycle[a];
  }

  Track* track;

  /* Generate 2D tracks */
  for (int a=0; a < _num_azim/2; a++) {

    /* Allocate memory for the 2D tracks array */
    _tracks_2D[a] = new Track[_num_x[a] + _num_y[a]];
    _num_2D_tracks += _num_x[a] + _num_y[a];

    /* Get the azimuthal angle for all tracks with this azimuthal angle */
    phi = _quadrature->getPhi(a);

    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      /* Get track and set angle and track indices */
      track = (&_tracks_2D[a][i]);
      track->setPhi(phi);
      track->setAzimIndex(a);
      track->setXYIndex(i);

      /* Set start point */
      if (a < _num_azim/4) {
        if (i < _num_x[a])
          track->getStart()->setCoords(width - dx_eff[a] * (i + 0.5), 0.0);
        else
          track->getStart()->setCoords(0.0, dy_eff[a] * (i-_num_x[a] + 0.5));
      }
      else {
        if (i < _num_x[a])
          track->getStart()->setCoords(dx_eff[a] * (i + 0.5), 0.0);
        else
          track->getStart()->setCoords(width, dy_eff[a] * (i-_num_x[a] + 0.5));
      }

      /* Set end point */
      if (a < _num_azim/4) {
        if (i < _num_y[a])
          track->getEnd()->setCoords(width, dy_eff[a] * (i + 0.5));
        else
          track->getEnd()->setCoords(width - dx_eff[a] * ((i-_num_y[a]) + 0.5),
                                     height);
      }
      else {
        if (i < _num_y[a])
          track->getEnd()->setCoords(0.0, dy_eff[a] * (i + 0.5));
        else
          track->getEnd()->setCoords(dx_eff[a] * (i-_num_y[a] + 0.5), height);
      }
    }
  }

  /* Set the flag indicating 2D tracks have been generated */
  _contains_2D_tracks = true;

  /* Initialize the track reflections and cycle ids */
  TrackGenerator::initializeTrackReflections();
  TrackGenerator::initializeTrackCycleIds();
  TrackGenerator::initializeTrackPeriodicIndices();
  initializeTrackCycles();
}


/**
 * @brief Initializes 2D Track cycles array
 * @details This method creates an array of 2D Tracks ordered by azimuthal
 *          angle, cycle index, and train index.
 */
void TrackGenerator::initializeTrackCycles() {

  _tracks_2D_cycle  = new Track***[_num_azim/4];
  for (int a=0; a < _num_azim/4; a++) {
    _tracks_2D_cycle[a] = new Track**[_cycles_per_azim[a]];
    for (int c=0; c < _cycles_per_azim[a]; c++) {
      _tracks_2D_cycle[a][c] = new Track*[_tracks_per_cycle[a]];
    }
  }

  bool fwd;
  Track* track;
  Track* track_prev;

  for (int a=0; a < _num_azim/4; a++) {
    for (int c=0; c < _cycles_per_azim[a]; c++) {
      track = &_tracks_2D[a][c];
      fwd = true;
      for (int i=0; i < _tracks_per_cycle[a]; i++) {

        track_prev = track;

        /* Add Track to 2D Track cycles array and set direction in cycle */
        _tracks_2D_cycle[a][c][i] = track;
        track->setDirectionInCycle(fwd);

        if (fwd) {
          track = static_cast<Track*>(track_prev->getTrackReflFwd());
          fwd = track_prev->getReflFwdFwd();
        }
        else {
          track = static_cast<Track*>(track_prev->getTrackReflBwd());
          fwd = track_prev->getReflBwdFwd();
        }
      }
    }
  }
}

/**
 * @brief Initializes 2D Track reflections
 * @details This method computes the connecting Tracks for all 2D Tracks in
 *          the TrackGenerator analytically, handling both reflective and
 *          periodic boundaries.
 */
void TrackGenerator::initializeTrackReflections() {

  log_printf(NORMAL, "Initializing 2D tracks reflections...");

  Track* track;
  int ac;

  /* Generate the 2D track cycles */
  for (int a=0; a < _num_azim/2; a++) {
    ac = _num_azim/2 - a - 1;
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      /* Get current track */
      track = &_tracks_2D[a][i];

      /* Set connecting tracks in forward direction */
      if (i < _num_y[a]) {
        track->setReflFwdFwd(true);
        track->setTrackReflFwd(&_tracks_2D[ac][i + _num_x[a]]);
        track->setTrackPrdcFwd(&_tracks_2D[a][i + _num_x[a]]);
      }
      else {
        track->setReflFwdFwd(false);
        track->setTrackReflFwd
          (&_tracks_2D[ac][_num_x[a] + _num_y[a] - (i - _num_y[a]) - 1]);
        track->setTrackPrdcFwd(&_tracks_2D[a][i - _num_y[a]]);
      }

      /* Set connecting tracks in backward direction */
      if (i < _num_x[a]) {
        track->setReflBwdFwd(true);
        track->setTrackReflBwd(&_tracks_2D[ac][_num_x[a] - i - 1]);
        track->setTrackPrdcBwd(&_tracks_2D[a][i + _num_y[a]]);
      }
      else {
        track->setReflBwdFwd(false);
        track->setTrackReflBwd(&_tracks_2D[ac][i - _num_x[a]]);
        track->setTrackPrdcBwd(&_tracks_2D[a][i - _num_x[a]]);
      }

      /* Set the foward boundary conditions */
      if (a < _num_azim/4) {
        if (i < _num_y[a])
          track->setBCFwd(_geometry->getMaxXBoundaryType());
        else
          track->setBCFwd(_geometry->getMaxYBoundaryType());

        if (i < _num_x[a])
          track->setBCBwd(_geometry->getMinYBoundaryType());
        else
          track->setBCBwd(_geometry->getMinXBoundaryType());
      }

      /* Set the backward boundary conditions */
      else {
        if (i < _num_y[a])
          track->setBCFwd(_geometry->getMinXBoundaryType());
        else
          track->setBCFwd(_geometry->getMaxYBoundaryType());

        if (i < _num_x[a])
          track->setBCBwd(_geometry->getMinYBoundaryType());
        else
          track->setBCBwd(_geometry->getMaxXBoundaryType());
      }
    }
  }
}


/**
 * @brief Cycle IDs are created for all 2D cycles and assigned to 2D Tracks
 * @details All tracks are traversed through connecting tracks, assigning cycle
 *          numbers until all 2D Tracks are traversed. This is done for both
 *          periodic and reflective connections.
 */
void TrackGenerator::initializeTrackCycleIds() {

  log_printf(NORMAL, "Initializing 2D track cycle ids...");

  int id = 0;
  bool fwd;
  Track* track;

  /* Set the periodic track cycle ids */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      track = &_tracks_2D[a][i];

      if (track->getPeriodicCycleId() == -1) {
        while (track->getPeriodicCycleId() == -1) {

          /* Set the periodic cycle id */
          track->setPeriodicCycleId(id);
          track = track->getTrackPrdcFwd();
        }
        id++;
      }
    }
  }

  id = 0;

  /* Set the reflective track cycle ids */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      track = &_tracks_2D[a][i];
      fwd = true;

      if (track->getReflectiveCycleId() == -1) {
        while (track->getReflectiveCycleId() == -1) {

          /* Set the reflective cycle id */
          track->setReflectiveCycleId(id);
          if (fwd) {
            fwd = track->getReflFwdFwd();
            track = track->getTrackReflFwd();
          }
          else {
            fwd = track->getReflBwdFwd();
            track = track->getTrackReflBwd();
          }
        }
        id++;
      }
    }
  }
}



/**
 * @brief Recalibrates Track start and end points to the origin of the Geometry.
 * @details The origin of the Geometry is designated at its center by
 *          convention, but for track initialization the origin is assumed to be
 *          at the bottom right corner for simplicity. This method corrects
 *          for this by re-assigning the start and end Point coordinates.
 */
void TrackGenerator::recalibrateTracksToOrigin() {

  /* Recalibrate the tracks to the origin and set the uid. Note that the
   * loop structure is unconventional in order to preserve an increasing
   * track uid value in the Solver's tracks array. The tracks array is
   * oriented with this loop structure in order to maintain reproducability
   * for parallel runs.
   */

  /* Loop over azim reflective halfspaces */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      double x0 = _tracks_2D[a][i].getStart()->getX();
      double y0 = _tracks_2D[a][i].getStart()->getY();
      double x1 = _tracks_2D[a][i].getEnd()->getX();
      double y1 = _tracks_2D[a][i].getEnd()->getY();
      double new_x0 = x0 + _geometry->getMinX();
      double new_y0 = y0 + _geometry->getMinY();
      double new_x1 = x1 + _geometry->getMinX();
      double new_y1 = y1 + _geometry->getMinY();

      _tracks_2D[a][i].setCoords(new_x0, new_y0, new_x1, new_y1);
    }
  }
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize() {

  log_printf(NORMAL, "Ray tracing for 2D track segmentation...");

  int tracks_segmented = 0;
  int num_2D_tracks = getNum2DTracks();

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    log_printf(NORMAL, "segmenting 2D tracks - Percent complete: %5.2f %%",
               double(tracks_segmented) / num_2D_tracks * 100.0);
#pragma omp parallel for
    for (int i=0; i < _num_x[a] + _num_y[a]; i++)
      _geometry->segmentize2D(&_tracks_2D[a][i], _z_coord);

    tracks_segmented += _num_x[a] + _num_y[a];
  }

  _geometry->initializeFSRVectors();
  _contains_2D_segments = true;

  return;
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
  if (!stat(directory.str().c_str(), &st) == 0)
    mkdir(directory.str().c_str(), S_IRWXU);

  /* Check to see if a Track file exists for this geometry, number of azimuthal
   * angles, and track spacing, and if so, import the ray tracing data */
  _tracks_filename = getTestFilename(directory.str());
  if (!stat(_tracks_filename.c_str(), &buffer)) {
    if (_segment_formation == EXPLICIT_3D || _segment_formation == EXPLICIT_2D) {
      if (readSegmentsFromFile()) {
        _use_input_file = true;
        setContainsSegments(true);
      }
    }
  }
}


/**
 * @brief Returns the filename for writing tracking data
 */
std::string TrackGenerator::getTestFilename(std::string directory) {

  std::stringstream test_filename;

  if (_geometry->getCmfd() != NULL)
    test_filename << directory << "/2D_"
                  << _num_azim << "_azim_"
                  << _azim_spacing << "_cm_spacing_cmfd_"
                  << _geometry->getCmfd()->getNumX()
                  << "x" << _geometry->getCmfd()->getNumY()
                  << ".data";
  else
    test_filename << directory << "/2D_"
                  << _num_azim << "_angles_"
                  << _azim_spacing << "_cm_spacing.data";

  return test_filename.str();
}


/**
 * @brief Updates whether the TrackGenerator contains segments
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
  ParallelHashMap<std::size_t, fsr_data*>* FSR_keys_map =
      _geometry->getFSRKeysMap();
  std::vector<std::size_t>* FSRs_to_keys = _geometry->getFSRsToKeys();
  std::vector<int>* FSRs_to_material_IDs = _geometry->getFSRsToMaterialIDs();
  std::size_t fsr_key;
  int fsr_id;
  int fsr_counter = 0;
  double x, y, z;

  /* Write number of FSRs */
  int num_FSRs = _geometry->getNumFSRs();
  fwrite(&num_FSRs, sizeof(int), 1, out);

  /* Write FSR vector maps to file */
  std::size_t* fsr_key_list = FSR_keys_map->keys();
  fsr_data** fsr_data_list = FSR_keys_map->values();
  Cmfd* cmfd = _geometry->getCmfd();
  for (int i=0; i < num_FSRs; i++) {

    /* Write data to file from FSR_keys_map */
    fsr_key = fsr_key_list[i];
    fsr_id = fsr_data_list[i]->_fsr_id;
    x = fsr_data_list[i]->_point->getX();
    y = fsr_data_list[i]->_point->getY();
    z = fsr_data_list[i]->_point->getZ();
    fwrite(&fsr_key, sizeof(std::size_t), 1, out);
    fwrite(&fsr_id, sizeof(int), 1, out);
    fwrite(&x, sizeof(double), 1, out);
    fwrite(&y, sizeof(double), 1, out);
    fwrite(&z, sizeof(double), 1, out);

    /* Write data to file from FSRs_to_material_IDs */
    fwrite(&(FSRs_to_material_IDs->at(fsr_counter)), sizeof(int), 1, out);

    /* Write data to file from FSRs_to_keys */
    fwrite(&(FSRs_to_keys->at(fsr_counter)), sizeof(std::size_t), 1, out);

    /* Increment FSR ID counter */
    fsr_counter++;
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

  return;
}


/**
 * @brief Reads Tracks in from a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 * @return true if able to read Tracks in from a file; false otherwise
 */
bool TrackGenerator::readSegmentsFromFile() {

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
  ParallelHashMap<std::size_t, fsr_data*>* FSR_keys_map =
      new ParallelHashMap<std::size_t, fsr_data*>;
  std::vector<int>* FSRs_to_material_IDs
    = new std::vector<int>;
  std::vector<std::size_t>* FSRs_to_keys
    = new std::vector<std::size_t>;
  int num_FSRs;
  std::size_t fsr_key;
  int fsr_key_id;
  double x, y, z;

  /* Get number of FSRs */
  ret = fread(&num_FSRs, sizeof(int), 1, in);

  /* Read FSR vector maps from file */
  for (int fsr_id=0; fsr_id < num_FSRs; fsr_id++) {

    /* Read data from file for FSR_keys_map */
    ret = fread(&fsr_key, sizeof(std::size_t), 1, in);
    ret = fread(&fsr_key_id, sizeof(int), 1, in);
    ret = fread(&x, sizeof(double), 1, in);
    ret = fread(&y, sizeof(double), 1, in);
    ret = fread(&z, sizeof(double), 1, in);
    fsr_data* fsr = new fsr_data;
    fsr->_fsr_id = fsr_key_id;
    Point* point = new Point();
    point->setCoords(x,y,z);
    fsr->_point = point;
    FSR_keys_map->insert(fsr_key, fsr);

    /* Read data from file for FSR_to_materials_IDs */
    int material_id;
    ret = fread(&material_id, sizeof(int), 1, in);
    FSRs_to_material_IDs->push_back(material_id);

    /* Read data from file for FSR_to_keys */
    ret = fread(&fsr_key, sizeof(std::size_t), 1, in);
    FSRs_to_keys->push_back(fsr_key);
  }

  /* Set FSR vector maps */
  _geometry->setFSRKeysMap(FSR_keys_map);
  _geometry->setFSRsToMaterialIDs(FSRs_to_material_IDs);
  _geometry->setFSRsToKeys(FSRs_to_keys);

  /* Read cmfd cell_fsrs vector of vectors from file */
  Cmfd* cmfd = _geometry->getCmfd();
  if (cmfd != NULL) {
    std::vector< std::vector<int> > cell_fsrs;
    int num_cells, fsr_id;
    ret = fread(&num_cells, sizeof(int), 1, in);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      std::vector<int>* fsrs = new std::vector<int>;
      cell_fsrs.push_back(*fsrs);
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

  /* Close the Track file */
  fclose(in);

  return true;
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
 * @brief Generates the numerical centroids of the FSRs.
 * @details This routine generates the numerical centroids of the FSRs
 *          by weighting the average x and y values of each segment in the
 *          FSR by the segment's length and azimuthal weight. The numerical
 *          centroid fomula can be found in R. Ferrer et. al. "Linear Source
 *          Approximation in CASMO 5", PHYSOR 2012.
 * @param FSR_volumes An array of FSR volumes.
 */
void TrackGenerator::generateFSRCentroids(FP_PRECISION* FSR_volumes) {

  int num_FSRs = _geometry->getNumFSRs();

  /* Create temporary array of centroids and initialize to origin */
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

  delete [] centroids;
}


/**
 * @brief Returns the direction in the cycle of the Track indexed by azimuthal
 *        index, cycle, and track index in the cycle
 * @details The 2D Track cycle indicated by the azimuthal angle and cycle
 *          number is traversed across track_index tracks, returning the
 *          direction of the Track at that position
 * @param azim The azimuthal index
 * @param cycle The 2D cycle number
 * @param track_index The track index into the cycle
 * @return the direction of the matching Track
 */
bool TrackGenerator::getCycleDirection(int azim, int cycle, int track_index) {

  return _tracks_2D_cycle[azim][cycle][track_index]->getDirectionInCycle();
}


/**
 * @brief Sets the max optical path length of 3D segments for use in
 *        on-the-fly computation
 * @param tau maximum optical path length
 */
void TrackGenerator::setMaxOpticalLength(FP_PRECISION tau) {
  _max_optical_length = tau;
}


/**
 * @breif Sets the maximum number of segments per Track
 * @param max_num_segments the maximum number of segments per Track
 */
void TrackGenerator::setMaxNumSegments(int max_num_segments) {
  _max_num_segments = max_num_segments;
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
 * @brief Counts the number of segments for each Track in the Geomtery
 * @details All segments are subject to the max optical path length to
 *          determine the number of segments for each track as well as the
 *          maximum number of segments per Track in the Geometry. For
 *          on-the-fly calculations, the temporary segment buffer is expanded
 *          to fit the calculated maximum number of segments per Track.
 */
void TrackGenerator::countSegments() {

  std::string msg = "Counting segments";
  Progress progress(_num_2D_tracks, msg);

  /* Count the number of segments on each track and update the maximium */
  SegmentCounter counter(this);
  counter.execute();

  /* Allocate new temporary segments if necessary */
  if (_segment_formation != EXPLICIT_3D && _segment_formation != EXPLICIT_2D)
    allocateTemporarySegments();
}


/**
 * @brief Sets the track periodic indices of all 2D Tracks
 * @details Periodic cylces are traversed until all 2D Tracks are visited and
 *          their periodic indices are set
 */
void TrackGenerator::initializeTrackPeriodicIndices() {

  log_printf(NORMAL, "Initializing track periodic indices...");

  if (!_periodic)
    return;

  Track* track;
  int track_index;

  /* Set the track periodic cycle indices for 2D tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {

      /* Get the current track */
      track = &_tracks_2D[a][i];

      /* Check if periodic track index has been set */
      if (track->getPeriodicTrackIndex() == -1) {

        /* Initialize the track index counter */
        track_index = 0;

        /* Set the periodic track indexes for all tracks in periodic cycle */
        while (track->getPeriodicTrackIndex() == -1) {

          /* Set the track periodic cycle */
          track->setPeriodicTrackIndex(track_index);

          /* Get the next track in cycle */
          track = track->getTrackPrdcFwd();

          /* Increment index counter */
          track_index++;
        }
      }
    }
  }
}


/**
 * @brief Creates a Track array by increasing uid
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
  int num_2D_tracks = getNum2DTracks();
  _tracks_2D_array = new Track*[num_2D_tracks];

  /* Loop over all 2D tracks */
  int uid = 0;
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
 * @brief returns whether periodic boundaries are present in Track generation
 * @return a boolean value - true if periodic; false otherwise
 */
bool TrackGenerator::getPeriodic() {
  return _periodic;
}


/**
 * @brief Sets a flag to record all segment information in the tracking file
 * @param A boolean value to determine whether or not to record segment
 *        information in the tracking file: true to record, false not to record
 */
void TrackGenerator::setDumpSegments(bool dump_segments) {
  _dump_segments = dump_segments;
}


/**
 * @brief Resets the TrackGenerator to not contain tracks or segments
 */
void TrackGenerator::resetStatus() {
  _contains_2D_tracks = false;
  _contains_2D_segments = false;
  _use_input_file = false;
  _tracks_filename = "";
}


/**
 * @brief Allocates memory for temporary segment storage if necessary
 * @details Temporary segments are not allocated for 2D calculations
 */
void TrackGenerator::allocateTemporarySegments() {}
