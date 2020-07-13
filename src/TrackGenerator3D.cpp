#include "TrackGenerator3D.h"
#include "TrackTraversingAlgorithms.h"


/**
 * @brief Constructor for the TrackGenerator3D assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param num_polar number of polar angles in \f$ [0, \pi] \f$
 * @param azim_spacing track azimuthal spacing (cm)
 * @param z_spacing track axial spacing (cm)
 */
TrackGenerator3D::TrackGenerator3D(Geometry* geometry, int num_azim,
                                   int num_polar, double azim_spacing,
                                   double z_spacing) :
                    TrackGenerator(geometry, num_azim, azim_spacing) {
  setNumPolar(num_polar);
  setDesiredZSpacing(z_spacing);
  _contains_3D_tracks = false;
  _contains_3D_segments = false;
  _contains_global_z_mesh = false;
  _contains_segmentation_heights = false;
  _contains_temporary_segments = false;
  _contains_temporary_tracks = false;
  _segment_formation = EXPLICIT_3D;
  _max_num_tracks_per_stack = 0;
  _num_seg_matrix_rows = 0;
  _num_seg_matrix_columns = 0;
  _tracks_3D = NULL;
  _tracks_2D_chains = NULL;

  _cum_tracks_per_stack = NULL;
  _cum_tracks_per_xy = NULL;
  _tracks_per_stack = NULL;
  _first_lz_of_stack = NULL;
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator3D::~TrackGenerator3D() {

  /* Delete 2D chains if created */
  if (_tracks_2D_chains != NULL) {
    for (int a=0; a < _num_azim/2; a++) {
      for (int x=0; x < _num_x[a]; x++)
        delete [] _tracks_2D_chains[a][x];
      delete [] _tracks_2D_chains[a];
    }
    delete [] _tracks_2D_chains;
  }

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_3D_tracks) {

    /* Delete 3D tracks */
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
        delete [] _cum_tracks_per_stack[a][i];
        delete [] _tracks_per_stack[a][i];
        delete [] _first_lz_of_stack[a][i];
      }
      delete [] _cum_tracks_per_stack[a];
      delete [] _tracks_per_stack[a];
      delete [] _cum_tracks_per_xy[a];
      delete [] _first_lz_of_stack[a];
    }
    delete [] _cum_tracks_per_stack;
    delete [] _tracks_per_stack;
    delete [] _cum_tracks_per_xy;
    delete [] _first_lz_of_stack;

    if (_tracks_3D != NULL) {
      for (int a=0; a < _num_azim/2; a++) {
        for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
          for (int p=0; p < _num_polar; p++)
            delete [] _tracks_3D[a][i][p];
          delete [] _tracks_3D[a][i];
        }
        delete [] _tracks_3D[a];
      }
      delete [] _tracks_3D;
      _tracks_3D = NULL;
    }

    /* Delete book keeping for 3D tracks */
    for (int a=0; a < _num_azim/2; a++) {
      delete [] _num_l[a];
      delete [] _num_z[a];
      delete [] _dz_eff[a];
      delete [] _dl_eff[a];
    }
    delete [] _num_l;
    delete [] _num_z;
    delete [] _dz_eff;
    delete [] _dl_eff;

    /* Delete temporary segments if they exist */
    if (_contains_temporary_segments) {
      for (int t = 0; t < _num_threads; t++) {
        delete [] _temporary_segments.at(t);
      }
    }

    /* Delete temporary Tracks if they exist */
    if (_contains_temporary_tracks) {
      for (int t = 0; t < _num_threads; t++) {
        delete [] _temporary_3D_tracks.at(t);
        delete [] _temporary_tracks_array.at(t);
      }
    }
  }
}


/**
 * @brief Return the number of polar angles in \f$ [0, \pi] \f$.
 * @return the number of polar angles in \f$ [0, \pi] \f$
 */
int TrackGenerator3D::getNumPolar() {
  return _num_polar;
}


/**
 * @brief Return the track polar spacing (cm).
 * @details This will return the user-specified track spacing and NOT the
 *          effective track spacing which is computed and used to generate
 *          cyclic tracks.
 * @return the track polar spacing (cm)
 */
double TrackGenerator3D::getDesiredZSpacing() {
  return _z_spacing;
}


/**
 * @brief Return the total number of 3D Tracks across the Geometry.
 * @return the total number of 3D Tracks
 */
long TrackGenerator3D::getNumTracks() {
  return getNum3DTracks();
}


/**
 * @brief Return the total number of 3D Tracks across the Geometry.
 * @return the total number of 3D Tracks
 */
long TrackGenerator3D::getNum3DTracks() {

  int a = _num_azim/2 - 1;
  int xy = _num_x[a] + _num_y[a] - 1;
  int p = _num_polar - 1;
  return _cum_tracks_per_stack[a][xy][p] + _tracks_per_stack[a][xy][p];
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
long TrackGenerator3D::getNumSegments() {
  return getNum3DSegments();
}


/**
 * @brief Return the total number of 3D Track segments across the Geometry.
 * @return the total number of 3D Track segments
 */
long TrackGenerator3D::getNum3DSegments() {

  if (!containsSegments())
    log_printf(ERROR, "Cannot get the number of 3D segments since they "
               "have not been generated.");


  /* Loop over all Tracks and count segments */
  long num_3D_segments = 0;
  if (_segment_formation == EXPLICIT_3D) {
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
            num_3D_segments += _tracks_3D[a][i][p][z].getNumSegments();
        }
      }
    }
  }
  else {
    SegmentCounter counter(this);
    counter.countTotalNumSegments();
    counter.execute();
    num_3D_segments = counter.getTotalNumSegments();
  }

  return num_3D_segments;
}


/**
 * @brief Returns the spacing between tracks in the axial direction for the
 *        requested azimuthal angle index and polar angle index.
 * @param azim the requested azimuthal angle index
 * @param polar the requested polar angle index
 * @return the effective axial spacing
 */
double TrackGenerator3D::getZSpacing(int azim, int polar) {
  return _dz_eff[azim][polar];
}


/**
 * @brief Returns the maximum number of tracks in a single stack.
 * @return the maximum number of tracks
 */
int TrackGenerator3D::getMaxNumTracksPerStack() {
  return _max_num_tracks_per_stack;
}


/**
 * @brief Returns the number of rows in the temporary segment storage matrix.
 * @details For on-the-fly computation, a matrix of temporary segments is
 *          allocated for each thread. This matrix is indexed by the z-stack
 *          index (row) and the segment number (column). For ray tracing by
 *          individual Tracks the number of rows is always one since temporary
 *          segments only need to be stored for one Track at a time.
 * @return the number of rows in the temporary segment storage matrix
 */
int TrackGenerator3D::getNumRows() {
  _num_seg_matrix_rows = 1;
  return _num_seg_matrix_rows;
}


/**
 * @brief Returns the number of columns in the temporary segment storage matrix.
 * @details For on-the-fly computation, a matrix of temporary segments is
 *          allocated for each thread. This matrix is indexed by the z-stack
 *          index (row) and the segment number (column). The number of columns
 *          is equal to the maximum number of segments per Track at the time of
 *          allocation.
 * @return the number of columns in the temporary segment storage matrix
 */
int TrackGenerator3D::getNumColumns() {
  return _num_seg_matrix_columns;
}


/**
 * @brief Returns an array of temporary segments for use in on-the-fly
 *        computations.
 * @details For on-the-fly computation, a matrix of temporary segments is
 *          allocated for each thread. This matrix is indexed by the z-stack
 *          index (row) and the segment number (column). The row index should
 *          be one if temporary segments are only required for one Track at a
 *          given time. If the segmentation method is not an on-the-fly method,
 *          NULL is returned.
 * @param thread_id The thread ID, as assigned by OpenMP
 * @return a pointer to the array of temporary segments
 */
segment* TrackGenerator3D::getTemporarySegments(int thread_id) {
  if (_contains_temporary_segments)
    return _temporary_segments.at(thread_id);
  else
    return NULL;
}


/**
 * @brief Returns an array of temporary 3D Tracks for use in on-the-fly
 *        computations.
 * @param thread_id The thread ID, as assigned by OpenMP
 */
Track3D* TrackGenerator3D::getTemporary3DTracks(int thread_id) {
  if (_contains_temporary_tracks)
    return _temporary_3D_tracks.at(thread_id);
  else
    return NULL;
}


/**
 * @brief Returns an array of pointers to temporary Tracks for use in
 *        on-the-fly computations.
 * @param thread_id The thread ID, as assigned by OpenMP
 */
Track** TrackGenerator3D::getTemporaryTracksArray(int thread_id) {
  if (_contains_temporary_tracks)
    return _temporary_tracks_array.at(thread_id);
  else
    return NULL;
}


/**
 * @brief Returns whether or not the TrackGenerator contains an allocation
 *        for temporary Tracks to be filled on-the-fly.
 */
bool TrackGenerator3D::containsTemporaryTracks() {
  return _contains_temporary_tracks;
}


/**
 * @brief Returns a 3D array of the number of 3D Tracks in each z-stack.
 * @details A 3D array is returned indexed first by azimuthal angle, second by
 *          2D track number, and third by polar angle. This array describes
 *          the number of tracks in each z-stack.
 * @return A 3D array of the number of tracks in each z-stack
 */
int*** TrackGenerator3D::getTracksPerStack() {
  return _tracks_per_stack;
}


/**
 * @brief Returns the number of 3D Tracks in the z-direction for a given
 *        azimuthal angle index and polar angle index.
 * @param azim the azimuthal angle index
 * @param polar the polar angle index
 * @return the number of 3D Tracks in the z-direction of the Geometry
 */
int TrackGenerator3D::getNumZ(int azim, int polar) {
  return _num_z[azim][polar];
}


/**
 * @brief Returns the number of 3D Tracks in the radial direction for a given
 *        azimuthal angle index and polar angle index.
 * @param azim the azimuthal angle index
 * @param polar the polar angle index
 * @return the number of 3D Tracks in the radial direction of the Geometry
 */
int TrackGenerator3D::getNumL(int azim, int polar) {
  return _num_l[azim][polar];
}


/**
 * @brief Set the number of polar angles in \f$ [0, \pi] \f$.
 * @param num_polar the number of polar angles in \f$ \pi \f$
 */
void TrackGenerator3D::setNumPolar(int num_polar) {

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
 * @brief Set the suggested track polar spacing (cm).
 * @param spacing the suggested track polar spacing
 */
void TrackGenerator3D::setDesiredZSpacing(double spacing) {
  if (spacing < 0)
    log_printf(ERROR, "Unable to set a negative track z-spacing "
               "%f for the TrackGenerator.", spacing);

  _z_spacing = spacing;
  resetStatus();
}


/**
 * @brief Set the type of segmentation used for segment formation.
 * @param segmentation_type a segmentationType defining the type of
 *        segmentation to be used in segment formation. Options are:
 *          - EXPLICIT_3D: explicit 2D/3D segment formation
 *          - OTF_TRACKS: axial on-the-fly ray tracing by individual tracks
 *          - OTF_STACKS: axial on-the-fly ray tracing by entire z-stacks
 */
void TrackGenerator3D::setSegmentFormation(segmentationType segmentation_type) {
  _segment_formation = segmentation_type;
}


/**
 * @brief Sets the z-planes over which 2D segmentation is performed for
 *        on-the-fly calculations
 * @param zones the z-coordinates defining the height of the radial
 *        segmentation planes
 */
void TrackGenerator3D::setSegmentationZones(std::vector<double> zones) {

  /* Clear any previous segmentation heights and note their existence */
  _contains_segmentation_heights = true;
  _segmentation_heights.clear();

  /* Check that the zones contain the Geometry minimum and maximum */
  double z_min = _geometry->getRootUniverse()->getMinZ();
  double z_max = _geometry->getRootUniverse()->getMaxZ();
  if (zones.size() >= 2) {
    if (fabs(zones.at(0) - z_min) > TINY_MOVE)
      log_printf(ERROR, "Segmentation zones must contain the Geometry minimum."
                 " The first value of the segmentation heights is %f and the "
                 "Geometry minimum is %f.", zones.at(0), z_min);
    if (fabs(zones.at(zones.size()-1) - z_max) > TINY_MOVE)
      log_printf(ERROR, "Segmentation zones must contain the Geometry maximum."
                 " The last value of the segmentation heights is %f and the "
                 "Geometry maximum is %f.", zones.at(zones.size()-1), z_max);
  }
  else {
    log_printf(ERROR, "Segmentation zones must contain the Geometry minimum "
               "and maximum. Therefore the segmentation zones must be of "
               "length greater than or equal to 2.");
  }

  /* Check that zones are monotonic and within total Geometry bounds */
  for (int i=1; i < zones.size(); i++)
    if (zones.at(i) < zones.at(i-1))
      log_printf(ERROR, "Segmentation zones must be monotonically "
                 "increasing. Axial level %f is less than axial level %f",
                 zones.at(i), zones.at(i-1));

  /* Find minimum and maximum indexes for this domain */
  int min_idx = 0;
  int max_idx = 0;
  double domain_z_min = _geometry->getMinZ();
  double domain_z_max = _geometry->getMaxZ();
  for (int i=0; i < zones.size(); i++) {
    if (zones.at(i) - TINY_MOVE < domain_z_min)
      min_idx++;
    if (zones.at(i) + TINY_MOVE < domain_z_max)
      max_idx++;
  }

  /* Form segmentation heights */
  if (min_idx == max_idx) {
    _segmentation_heights.push_back((domain_z_min + domain_z_max)/2);
  }
  else {
    _segmentation_heights.push_back((domain_z_min + zones.at(min_idx))/2);
    for (int i=min_idx+1; i < max_idx; i++)
      _segmentation_heights.push_back((zones.at(i-1) + zones.at(i))/2);
    _segmentation_heights.push_back((zones.at(max_idx-1) + domain_z_max)/2);
  }
}


/**
 * @brief Sets a global z-mesh to use during axial on-the-fly ray tracing.
 * @details In axial on-the-fly ray tracing, normally each extruded FSR
 *          contains a z-mesh. During on-the-fly segmentation when a new
 *          extruded FSR is entered, a binary search must be conducted to
 *          determine the axial cell. Alternatively, this function can be
 *          called which creates a global z-mesh from the geometry so that
 *          binary searches must only be conducted at the beginning of the
 *          track.
 */
void TrackGenerator3D::useGlobalZMesh() {
  _contains_global_z_mesh = true;
  _global_z_mesh = _geometry->getUniqueZHeights(true);
}


/**
 * @brief Provides the global z-mesh and size if available.
 * @details For some cases, a global z-mesh is generated for the Geometry. If
 *          so, a pointer to the associated mesh (array) is updated as well as
 *          the number of FSRs in the mesh. If no global z-mesh has been
 *          generated, a null pointer is given to z_mesh and the number of FSRs
 *          is assigned to be zero.
 * @param z_mesh The global z-mesh to be updated
 * @param num_fsrs The number of FSRs in the z-mesh
 */
void TrackGenerator3D::retrieveGlobalZMesh(double*& z_mesh, int& num_fsrs) {
  if (_contains_global_z_mesh) {
    z_mesh = &_global_z_mesh[0];
    num_fsrs = _global_z_mesh.size() - 1;
  }
  else {
    z_mesh = NULL;
    num_fsrs = 0;
  }
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
 *          coords = track_generator.retrieveTrackCoords(num_tracks*6)
 * @endcode
 *
 * @param coords an array of coords of length 6 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator3D::retrieveTrackCoords(double* coords, long num_tracks) {
  retrieve3DTrackCoords(coords, num_tracks);
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
 *          num_tracks = track_generator.getNum3DTracks()
 *          coords = track_generator.retrieve3DTrackCoords(num_tracks*6)
 * @endcode
 *
 * @param coords an array of coords of length 6 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator3D::retrieve3DTrackCoords(double* coords, long num_tracks) {

  if (num_tracks != NUM_VALUES_PER_RETRIEVED_TRACK * getNum3DTracks())
    log_printf(ERROR, "Unable to retrieve the Track coordinates since the "
               "TrackGenerator contains %d Tracks with %d coordinates but an "
               "array of length %d was input",
               getNum3DTracks(), NUM_VALUES_PER_RETRIEVED_TRACK *
               getNum3DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          TrackStackIndexes tsi;
          tsi._azim = a;
          tsi._xy = i;
          tsi._polar = p;
          tsi._z = z;

          Track3D track;
          getTrackOTF(&track, &tsi);
          coords[counter]   = track.getStart()->getX();
          coords[counter+1] = track.getStart()->getY();
          coords[counter+2] = track.getStart()->getZ();
          coords[counter+3] = track.getEnd()->getX();
          coords[counter+4] = track.getEnd()->getY();
          coords[counter+5] = track.getEnd()->getZ();
          counter += NUM_VALUES_PER_RETRIEVED_TRACK;
        }
      }
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
 *          coords = track_generator.retrieveSegmentCoords(num_segments*7)
 * @endcode
 *
 * @param coords an array of coords of length 7 times the number of segments
 * @param num_segments the total number of Track segments
 */
void TrackGenerator3D::retrieveSegmentCoords(double* coords,
                                             long num_segments) {
  retrieve3DSegmentCoords(coords, num_segments);
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
 *          num_segments = track_generator.getNum3DSegments()
 *          coords = track_generator.retrieve3DSegmentCoords(num_segments*7)
 * @endcode
 *
 * @param coords an array of coords of length 7 times the number of segments
 * @param num_segments the total number of Track segments
 */
void TrackGenerator3D::retrieve3DSegmentCoords(double* coords,
                                               long num_segments) {

  if (num_segments != NUM_VALUES_PER_RETRIEVED_SEGMENT * getNum3DSegments())
    log_printf(ERROR, "Unable to retrieve the Track segment coordinates since "
               "the TrackGenerator contains %d segments with %d coordinates "
               "but an array of length %d was input",
               getNum3DSegments(), NUM_VALUES_PER_RETRIEVED_SEGMENT *
               getNum3DSegments(), num_segments);

  segment* curr_segment = NULL;
  double x0, x1, y0, y1, z0, z1;
  double phi, theta;
  segment* segments;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  int counter = 0;
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          TrackStackIndexes tsi;
          tsi._azim = a;
          tsi._xy = i;
          tsi._polar = p;
          tsi._z = z;

          Track3D track;
          getTrackOTF(&track, &tsi);

          x0    = track.getStart()->getX();
          y0    = track.getStart()->getY();
          z0    = track.getStart()->getZ();
          phi   = track.getPhi();
          theta = track.getTheta();

          _geometry->segmentize3D(&track);
          segments = track.getSegments();

          for (int s=0; s < track.getNumSegments(); s++) {

            curr_segment = &segments[s];

            coords[counter] = curr_segment->_region_id;

            coords[counter+1] = x0;
            coords[counter+2] = y0;
            coords[counter+3] = z0;

            x1 = x0 + cos(phi) * sin(theta) * curr_segment->_length;
            y1 = y0 + sin(phi) * sin(theta) * curr_segment->_length;
            z1 = z0 + cos(theta) * curr_segment->_length;

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
  }
}


/**
 * @brief Initializes Track azimuthal angles, start and end Points.
 * @details This method computes the azimuthal angles and effective track
 *          spacing to use to guarantee cyclic Track wrapping. Based on the
 *          angles and spacing, the number of Tracks per angle and the start
 *          and end Points for each Track are computed.
 */
void TrackGenerator3D::initializeTracks() {

  /* Initialize the 2D Tracks */
  TrackGenerator::initializeTracks();

  /* Initialize the 2D Track chains */
  initialize2DTrackChains();

  /* Make sure that the depth of the Geometry is nonzero */
  if (_geometry->getWidthZ() <= 0)
    log_printf(ERROR, "The total depth of the Geometry must"
               " be nonzero for 3D Track generation. Create a CellFill which "
               "is filled by the entire geometry and bounded by XPlanes, "
               "YPlanes, and ZPlanes to enable the Geometry to determine the "
               "total width, height, and depth of the model.");

  if (!_contains_2D_tracks)
    log_printf(ERROR, "Cannot initialize 3D tracks since the 2D tracks "
               "have not been created");

  log_printf(NORMAL, "Initializing 3D tracks...");

  /* Allocate arrays for 3D Track data */
  _dz_eff = new double*[_num_azim/2];
  _dl_eff = new double*[_num_azim/2];
  _num_z            = new int*[_num_azim/2];
  _num_l            = new int*[_num_azim/2];
  _num_3D_tracks    = 0;

  for (int i=0; i < _num_azim/2; i++) {
    _dz_eff[i]         = new double[_num_polar];
    _dl_eff[i]         = new double[_num_polar];
    _num_z[i]          = new int[_num_polar];
    _num_l[i]          = new int[_num_polar];
  }

  /* Allocate memory for tracks per stack */
  _cum_tracks_per_stack = new long**[_num_azim/2];
  _cum_tracks_per_xy = new long*[_num_azim/2];
  _tracks_per_stack = new int**[_num_azim/2];
  _first_lz_of_stack = new int**[_num_azim/2];
  for (int a=0; a < _num_azim/2; a++) {
    _cum_tracks_per_stack[a] = new long*[_num_x[a] + _num_y[a]];
    _cum_tracks_per_xy[a] = new long[_num_x[a] + _num_y[a]];
    _tracks_per_stack[a] = new int*[_num_x[a] + _num_y[a]];
    _first_lz_of_stack[a] = new int*[_num_x[a] + _num_y[a]];
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      _cum_tracks_per_stack[a][i] = new long[_num_polar];
      _cum_tracks_per_xy[a][i] = 0;
      _tracks_per_stack[a][i] = new int[_num_polar];
      _first_lz_of_stack[a][i] = new int[_num_polar];
      for (int p=0; p < _num_polar; p++) {
        _cum_tracks_per_stack[a][i][p] = 0;
        _tracks_per_stack[a][i][p] = 0;
        _first_lz_of_stack[a][i][p] = -1;
      }
    }
  }

  double x1, x2, y1, y2, z1, z2;
  double width_x = _geometry->getWidthX();
  double width_y = _geometry->getWidthY();
  double width_z = _geometry->getWidthZ();
  double avg_polar_spacing = 0.0;
  double sum_correction = 0.0;
  double max_correction = 0.0;

  /* Determine angular quadrature and track spacing */
  for (int i = 0; i < _num_azim/4; i++) {

    double phi = _quadrature->getPhi(i);

    /* Determine the polar angles and spacing for this azimuthal angle */
    for (int j=0; j < _num_polar/2; j++) {

      /* Compute the cosine weighted average angle */
      double theta = _quadrature->getTheta(i, j);
      double theta_new;

      /* Compute the length to traverse one domain y-width */
      double module_width_y = width_y / _geometry->getNumYModules();
      double length = module_width_y / sin(phi);

      /* The number of intersections with xy (denoted "l") plane */
      _num_l[i][j] = int(ceil(length * tan(M_PI_2 - theta) / _z_spacing));

      /* Number of crossings along the z axis */
      double module_width_z = width_z / _geometry->getNumZModules();
      _num_z[i][j] = (int) ceil(module_width_z * _num_l[i][j] * tan(theta)
                                / length);
      _num_l[i][j] *= _geometry->getNumYModules();
      _num_z[i][j] *= _geometry->getNumZModules();

      /* Effective track spacing */
      _dl_eff[i][j] = width_y / (sin(phi) * _num_l[i][j]);
      _dz_eff[i][j] = width_z / _num_z[i][j];

      /* Evaluate the polar correction */
      theta_new = atan(_dl_eff[i][j] / _dz_eff[i][j]);
      sum_correction += pow(((theta_new - theta)/theta), 2);
      max_correction = std::max(max_correction, std::abs((theta_new - theta)/theta));

      /* Set the corrected polar angle */
      theta = theta_new;
      _quadrature->setTheta(theta, i, j);
      _quadrature->setPolarSpacing(_dz_eff[i][j] * sin(theta), i, j);

      /* Save information for supplementary angles */
      int supp_azim = _num_azim/2 - i - 1;
      _num_z[supp_azim][j] = _num_z[i][j];
      _dz_eff[supp_azim][j] = _dz_eff[i][j];
      _num_l[supp_azim][j] = _num_l[i][j];
      _dl_eff[supp_azim][j] = _dl_eff[i][j];

      /* Save information for complementary angles */
      int comp_polar = _num_polar - j - 1;
      _num_z[i][comp_polar] = _num_z[i][j];
      _dz_eff[i][comp_polar] = _dz_eff[i][j];
      _num_l[i][comp_polar] = _num_l[i][j];
      _dl_eff[i][comp_polar] = _dl_eff[i][j];

      _num_z[supp_azim][comp_polar] = _num_z[i][j];
      _dz_eff[supp_azim][comp_polar] = _dz_eff[i][j];
      _num_l[supp_azim][comp_polar] = _num_l[i][j];
      _dl_eff[supp_azim][comp_polar] = _dl_eff[i][j];
    }
  }

  /* Report polar correction information */
  log_printf(NORMAL, "RMS polar angle correction %f %%", sqrt(sum_correction)
                      / _num_azim / _num_polar * 8 * 100);
  log_printf(NORMAL, "Max polar angle correction %f %%", max_correction * 100);

  /* Compute the total number of chains */
  int num_chains = 0;
  for (int a = 0; a < _num_azim/2; a++) {
    for (int x = 0; x < _num_x[a]; x++)
      num_chains += _num_polar;
  }

  /* Create array of track chain indices */
  TrackChainIndexes* tcis = new TrackChainIndexes[num_chains];

  num_chains = 0;
  for (int a = 0; a < _num_azim/2; a++) {
    for (int x = 0; x < _num_x[a]; x++) {
      for (int p = 0; p < _num_polar; p++) {
        tcis[num_chains]._azim  = a;
        tcis[num_chains]._x = x;
        tcis[num_chains]._polar = p;
        num_chains++;
      }
    }
  }

  /* Create the 3D tracks for each lz plane */
  /*
   *                A sample 2D track layout
   *       _____________________________________________
   *      |        !    !         !       !             |
   *      |        !    !         !       !             |
   * ^    |        !    !         !       !             |
   * |    |        !    !         !       !             |
   * z+   |        !    !         !       !             |
   *       ------->!--->!-------->!------>!------------>|
   * l+ ->
   */
  getCycleTrackData(tcis, num_chains, false);

  /* Allocate memory for 3D track stacks */
  create3DTracksArrays();

  /* Initialize tracks in _tracks_3D array, save tracks */
  if (_segment_formation == EXPLICIT_3D)
    getCycleTrackData(tcis, num_chains, true);

  /* Delete the array of chain track indexes */
  delete [] tcis;

  /* Record the number of tracks, maximum number of tracks in a single stack */
  for (int a=0; a < _num_azim/2; a++)
    for (int i=0; i < _num_x[a] + _num_y[a]; i++)
      for (int p=0; p < _num_polar; p++) {
        _num_3D_tracks += _tracks_per_stack[a][i][p];
        if (_tracks_per_stack[a][i][p] > _max_num_tracks_per_stack)
          _max_num_tracks_per_stack = _tracks_per_stack[a][i][p];
      }

  /* Allocate temporary Tracks if necessary */
  if (_segment_formation == OTF_STACKS)
    allocateTemporaryTracks();

  log_printf(NORMAL, "3D Tracks in domain: %ld", getNum3DTracks());

  _contains_3D_tracks = true;
}


/**
 * @brief Initializes 2D Track chains array.
 * @details This method creates an array of 2D Tracks ordered by azimuthal
 *          angle, x index, and link index.
 */
void TrackGenerator3D::initialize2DTrackChains() {

  Track* track;
  int link_index;

  _tracks_2D_chains = new Track***[_num_azim/2];
  for (int a=0; a < _num_azim/2; a++) {
    _tracks_2D_chains[a] = new Track**[_num_x[a]];
    for (int x=0; x < _num_x[a]; x++) {

      /* Get the first track in the 2D chain */
      link_index = 0;
      track = &_tracks_2D[a][x];
      track->setLinkIndex(link_index);

      /* Cycle through 2D chain's tracks, set their index, get chain length */
      while (track->getXYIndex() < _num_y[a]) {
        link_index++;
        track = _tracks_2D_array[track->getTrackPrdcFwd()];
        track->setLinkIndex(link_index);
      }

      /* Allocate memory for the track chains */
      _tracks_2D_chains[a][x] = new Track*[link_index + 1];

      /* Assign tracks to the chains array */
      link_index = 0;
      track = &_tracks_2D[a][x];
      _tracks_2D_chains[a][x][link_index] = track;

      while (track->getXYIndex() < _num_y[a]) {
        link_index++;
        track = _tracks_2D_array[track->getTrackPrdcFwd()];
        _tracks_2D_chains[a][x][link_index] = track;
      }
    }
  }
}


/**
 * @brief Save or initialize 3D tracks by traversing all the 3D chain tracks,
 *        and split each 3D chain track into 3D tracks belonging to
 *        the z-stacks.
 * @param tcis array of Track chain indexes, sorted by chain
 * @param num_chains number of 3D Track chains to examine
 * @param save_tracks whether saving tracks or initializing them
 */
void TrackGenerator3D::getCycleTrackData(TrackChainIndexes* tcis,
                                         int num_chains, bool save_tracks) {

  /* Determine if initializing or allocating tracks */
  std::string msg;
  if (save_tracks)
    msg = "Saving 3D Tracks";
  else
    msg = "Initializing 3D Tracks";
  Progress progress(num_chains, msg);

  /* Loop over 3D track chains */
#pragma omp parallel
  {
    Track3D track_3D;
    TrackChainIndexes* tci;
    int nl, nz;
#pragma omp for
    for (int chain=0; chain < num_chains; chain++) {
      tci = &tcis[chain];
      nl = _num_l[tci->_azim][tci->_polar];
      nz = _num_z[tci->_azim][tci->_polar];
      for (int lz=0; lz < nl + nz; lz++) {
        tci->_lz = lz;
        tci->_link = -1;
        set3DTrackData(tci, &track_3D, true, save_tracks);
      }
      progress.incrementCounter();
    }
  }
}


/**
 * @brief Compute the starting coordinates of a 3D chain track, and the link
          number of the 2D track where its start point reside on.
 * @param tci the track chain index of a chain track
 * @param track_3D the chain track
 * @return the link number of the 2D track where the track_3D start point
 *         reside on
 */
int TrackGenerator3D::getFirst2DTrackLinkIndex(TrackChainIndexes* tci,
                                               Track3D* track_3D) {

  /* Initialize variables */
  Track* track_2D = _tracks_2D_chains[tci->_azim][tci->_x][0];

  double phi      = track_2D->getPhi();
  double cos_phi  = cos(phi);
  double sin_phi  = sin(phi);

  Point* start    = track_2D->getStart();
  double x_start  = start->getX();
  double y_start  = start->getY();

  int lz          = tci->_lz;
  int nl          = _num_l[tci->_azim][tci->_polar];
  int nz          = _num_z[tci->_azim][tci->_polar];
  double dz       = _dz_eff[tci->_azim][tci->_polar];
  double dl       = _dl_eff[tci->_azim][tci->_polar];

  /* Get Geometry information */
  double width_x = _x_max - _x_min;
  double width_y = _y_max - _y_min;

  /* Get the length from the beginning of the chain to the Track start point */
  double l_start  = 0.0;
  if (tci->_polar < _num_polar / 2 && lz < nl)
    l_start = width_y / sin_phi - (lz + 0.5) * dl;
  else if (tci->_polar >= _num_polar / 2 && lz >= nz)
    l_start = dl * (lz - nz + 0.5);

  double x_ext = x_start - _x_min + l_start * cos_phi;
  double y_ext = y_start - _y_min + l_start * sin_phi;

  /* If the start point is on a double reflection, nudge it to avoid track
   * with the same start and end points */
  bool nudged = false;
  if (fabs(l_start) > FLT_EPSILON) {
    /* x_ext and y_ext are both 0 or n*width_x/y at a double reflection */
    //FIXME Understand why sides instead of corners
    if (fabs(round(x_ext / width_x) * width_x - x_ext) < TINY_MOVE ||
        fabs(round(y_ext / width_y) * width_y - y_ext) < TINY_MOVE) {
      l_start += 10 * TINY_MOVE;
      x_ext = x_start - _x_min + l_start * cos_phi;
      y_ext = y_start - _y_min + l_start * sin_phi;
      nudged = true;
    }
  }

  /* Get the index of the first 2D Track link */
  int link_index = abs(int(floor(x_ext / width_x)));

  /* Get the starting x-coord */
  double x1;
  if (x_ext < 0.0)
    x1 = fmod(x_ext, width_x) + _x_max;
  else
    x1 = fmod(x_ext, width_x) + _x_min;

  /* Get the starting y-coord */
  double y1 = fmod(y_ext, width_y) + _y_min;

  /* Get the track's starting and ending z-coords */
  /* lz index of track in [0, nz+nl]
   *             nl tracks
   *             cutting top plane
   *           _____________
   *          |/ / / / / / /|
   *  nz      | / / / / / / | nz
   *  tracks  |/ / / / / / /| tracks cutting Z axis
   *  cutting | / / / / / / | on other side
   *  Z axis  |/ / / / / / /|
   *          |_/_/_/_/_/_/_|
   *           / / / / / /
   *          / / / / / /
   *           / / / / /
   *  nl      / / / / /
   *  tracks   / / / /
   *  cutting / / / /
   *  bottom   / / /
   *  horiz.  / / /
   *  plane    / /
   *          / /
   *           /
   *          /
   */
  double z1;
  double z2;
  if (tci->_polar < _num_polar / 2) {
    z1 = _z_min + std::max(0., (lz - nl + 0.5)) * dz;
    z2 = _z_max + std::min(0., (lz - nz + 0.5)) * dz;
  }
  else {
    z1 = _z_max + std::min(0., (lz - nz + 0.5)) * dz;
    z2 = _z_min + std::max(0., (lz - nl + 0.5)) * dz;
  }

  /* Reverse nudging of point */
  if (nudged) {
    x1 -= 10 * TINY_MOVE * cos_phi;
    y1 -= 10 * TINY_MOVE * sin_phi;
  }

  /* Set the Track start point and ending z-coord */
  track_3D->getStart()->setCoords(x1, y1, z1);
  track_3D->getEnd()->setZ(z2);

  return link_index;
}


/**
 * @brief A 3D Track is decomposed by intersections with x and y boundaries of
 *        z-stacks.
 * @details A 3D Track which initially neglects intersections with x and y
 *          boundaries is split into multiple tracks by finding those
 *          x and y intersections. If create_arrays is True, the number of
 *          tracks in the associated z-stacks is incremented.
 * @param tci pointer to a track chain index representing a chain track
 * @param track The 3D track to be decomposed
 * @param create_arrays whether creating the stack array of Tracks
 * @param save_tracks whether to save the Track3D or just initialize it. It may
 *        only be true for EXPLICIT_3D ray tracing.
 */
void TrackGenerator3D::set3DTrackData(TrackChainIndexes* tci,
                                      Track3D* track,
                                      bool create_arrays,
                                      bool save_tracks) {

  /* Initialize variables */
  double theta = _quadrature->getTheta(tci->_azim, tci->_polar);
  bool end_of_chain = false;

  /* Get the index of the first 2D Track link */
  int link = getFirst2DTrackLinkIndex(tci, track);
  int first_link = link;

  /* Record start and end points */
  double x1;
  double y1;
  double z1;

  double x2 = track->getStart()->getX();
  double y2 = track->getStart()->getY();
  double z2 = track->getStart()->getZ();

  /* Set the start and end point for each 3D track */
  while (!end_of_chain) {

    /* Get the 2D track associated with this 3D track */
    Track* track_2D = _tracks_2D_chains[tci->_azim][tci->_x][link];
    double phi = track_2D->getPhi();
    int azim = track_2D->getAzimIndex();
    int xy = track_2D->getXYIndex();

    /* Set the start and end points */
    x1 = x2;
    y1 = y2;
    z1 = z2;

    /* Get the distance to the nearest X or Y boundary */
    double dx;
    double dy;
    double dl_xy;
    /* For first link, track starts at x1,y1,z1, in middle of Z-boundary */
    if (link == first_link) {
      dx = track_2D->getEnd()->getX() - x1;
      dy = track_2D->getEnd()->getY() - y1;
      dl_xy = sqrt(dx*dx + dy*dy);
    }
    /* For other links, track starts on the X or Y boundaries */
    else
      dl_xy = track_2D->getLength();

    /* Get the distance to the nearest Z boundary */
    double dl_z;
    if (tci->_polar < _num_polar / 2)
      dl_z = (track->getEnd()->getZ() - z1) * tan(theta);
    else
      dl_z = (z1 - track->getEnd()->getZ()) / tan(theta - M_PI_2);

    double dl = std::min(dl_z, dl_xy);

    /* Set the end point coords */
    int polar;
    x2 = x1 + dl * cos(phi);
    y2 = y1 + dl * sin(phi);
    polar = tci->_polar;

    if (tci->_polar < _num_polar / 2)
      z2 = z1 + dl / tan(theta);
    else
      z2 = z1 - dl * tan(theta - M_PI_2);

    /* If a track has the same start and end point, ignore it */
    if (fabs(x2 - x1) < TINY_MOVE ||
        fabs(y2 - y1) < TINY_MOVE ||
        fabs(z2 - z1) < TINY_MOVE)
      break;

    /* If the Z boundary or last link or the desired link has been reached,
     * exit. tci->_link appoints the desired track */
    if (dl_z < dl_xy || track_2D->getXYIndex() >= _num_y[tci->_azim] ||
        tci->_link == link - first_link)
      end_of_chain = true;

    if (create_arrays || save_tracks) {

      /* Get the index in the z-stack */
      int z = _tracks_per_stack[azim][xy][polar];
      _tracks_per_stack[azim][xy][polar]++;

      if (z == 0)
        _first_lz_of_stack[azim][xy][polar] = tci->_lz;

      if (_segment_formation == EXPLICIT_3D && save_tracks) {

        /* Get a pointer to this 3D track in the global tracks array */
        Track3D* track_3D = &_tracks_3D[azim][xy][polar][z];

        /* Ensure the track does not lie outside the geometry bounds,
           due to floating point errors */
        x1 = std::max(_x_min, std::min(_x_max, x1));
        y1 = std::max(_y_min, std::min(_y_max, y1));
        z1 = std::max(_z_min, std::min(_z_max, z1));
        x2 = std::max(_x_min, std::min(_x_max, x2));
        y2 = std::max(_y_min, std::min(_y_max, y2));
        z2 = std::max(_z_min, std::min(_z_max, z2));

        /* Set the start and end points */
        track_3D->getStart()->setCoords(x1, y1, z1);
        track_3D->getEnd()->setCoords(x2, y2, z2);

        /* Set track attributes */
        track_3D->setAzimIndex(azim);
        track_3D->setXYIndex(xy);
        track_3D->setPolarIndex(polar);
        track_3D->setTheta(theta);
        track_3D->setPhi(phi);
      }
    }

    /* Increment the link index */
    link++;

    /* Move the starting x-coord to account for periodic BCs for chains */
    if (!end_of_chain) {
      if (tci->_azim < _num_azim / 4)
        x2 = _x_min;
      else
        x2 = _x_max;
    }
  }

  /* Ensure the track does not lie outside the geometry bounds,
     due to floating point errors */
  x1 = std::max(_x_min, std::min(_x_max, x1));
  y1 = std::max(_y_min, std::min(_y_max, y1));
  z1 = std::max(_z_min, std::min(_z_max, z1));
  x2 = std::max(_x_min, std::min(_x_max, x2));
  y2 = std::max(_y_min, std::min(_y_max, y2));
  z2 = std::max(_z_min, std::min(_z_max, z2));

  /* Set the start and end points */
  track->getStart()->setCoords(x1, y1, z1);
  track->getEnd()->setCoords(x2, y2, z2);

  /* Set the link index, if it has not been set */
  if (tci->_link == -1)
    tci->_link = link - first_link;
}


/**
 * @brief Generates 2D segments for each extruded track across the Geometry,
 *        initializing axially extruded regions as well as 3D FSRs.
 */
void TrackGenerator3D::segmentizeExtruded() {

  log_printf(NORMAL, "Ray tracing for axially extruded track segmentation...");

  /* Get all unique z-coords at which 2D radial segmentation is performed */
  std::vector<double> z_coords;
  if (_contains_segmentation_heights)
    z_coords = _segmentation_heights;
  else
    z_coords = _geometry->getUniqueZPlanes();

  /* Loop over all extruded Tracks */
  Progress progress(_num_2D_tracks, "Segmenting 2D Tracks", 0.1, _geometry,
                    true);
#pragma omp parallel for schedule(dynamic)
  for (int index=0; index < _num_2D_tracks; index++) {
    progress.incrementCounter();
    _geometry->segmentizeExtruded(_tracks_2D_array[index], z_coords);
  }

  /* Output memory consumption of 2D explicit ray tracing */
  _contains_2D_segments = true;
  printMemoryReport();

  /* Initialize 3D FSRs and their associated vectors */
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Barrier(_geometry->getMPICart());
#endif
  log_printf(NORMAL, "Initializing FSRs axially...");
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Barrier(_geometry->getMPICart());
#endif
  _geometry->initializeAxialFSRs(_global_z_mesh);
  _geometry->initializeFSRVectors();

  /* Count the number of segments in each track */
  countSegments();
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator3D::segmentize() {

  /* Check for on-the-fly methods */
  if (_segment_formation != EXPLICIT_3D) {
    segmentizeExtruded();
    return;
  }

  log_printf(NORMAL, "Ray tracing for 3D track segmentation...");

  long num_segments = 0;
  Progress progress(_num_3D_tracks, "Segmenting 3D Tracks", 0.1, _geometry,
                    true);

  /* Loop over all Tracks */  //FIXME Move openmp section over all tracks
  for (int a=0; a < _num_azim/2; a++) {

#pragma omp parallel for schedule(dynamic)
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
          progress.incrementCounter();
          _geometry->segmentize3D(&_tracks_3D[a][i][p][z]);
          TrackStackIndexes tsi;
          tsi._azim = a;
          tsi._xy = i;
          tsi._polar = p;
          tsi._z = z;
          _tracks_3D[a][i][p][z].setUid(get3DTrackID(&tsi));

          /* Set boundary conditions and linking Track data */
          //FIXME Move to track initialization rather than segmentation
          TrackChainIndexes tci;
          convertTSItoTCI(&tsi, &tci);
          setLinkingTracks(&tsi, &tci, true, &_tracks_3D[a][i][p][z]);
          setLinkingTracks(&tsi, &tci, false, &_tracks_3D[a][i][p][z]);

#pragma omp atomic update
          num_segments += _tracks_3D[a][i][p][z].getNumSegments();
        }
      }
    }
  }
  _geometry->initializeFSRVectors();
  _contains_3D_segments = true;

  log_printf(INFO, "Explicit 3D segments storage = %.2f MB", num_segments *
             sizeof(segment) / 1e6);
}


/**
 * @brief Fills an array with the x,y,z coordinates for a given track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          coords = track_generator.retrieveTrackCoords(track_id)
 * @endcode
 *
 * @param coords an array of coords of length 6
 * @param track_id the ID of the requested track
 */
void TrackGenerator3D::retrieveSingle3DTrackCoords(double coords[6],
                                                   long track_id) {

  /* Find 3D track associated with track_id */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          TrackStackIndexes tsi;
          tsi._azim = a;
          tsi._xy = i;
          tsi._polar = p;
          tsi._z = z;

          Track3D track;
          getTrackOTF(&track, &tsi);
          if (track.getUid() == track_id) {
            coords[0] = track.getStart()->getX();
            coords[1] = track.getStart()->getY();
            coords[2] = track.getStart()->getZ();
            coords[3] = track.getEnd()->getX();
            coords[4] = track.getEnd()->getY();
            coords[5] = track.getEnd()->getZ();
            return;
          }
        }
      }
    }
  }
  log_printf(ERROR, "Unable to find a 3D track associated with the given track"
                    "ID %ld during coordinate retrieval", track_id);
  return;
}


/**
 * @brief Allocates memory for 3D Tracks.
 * @details Before calling this function, the number of tracks per z-stack
 *          should be known and initialized in the _tracks_per_stack 3D array
 */
void TrackGenerator3D::create3DTracksArrays() {

  long num_tracks = 0;
  long tracks_per_azim = 0;

  /* Count number of tracks : total, per azimuthal and polar angle */
  for (int a=0; a < _num_azim/2; a++) {
    tracks_per_azim = 0;
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++) {
        _cum_tracks_per_stack[a][i][p] = num_tracks;
        num_tracks += _tracks_per_stack[a][i][p];
        tracks_per_azim += _tracks_per_stack[a][i][p];
      }
      _cum_tracks_per_xy[a][i] = num_tracks;
    }
    _tracks_per_azim[a] = tracks_per_azim;
  }

  /* Allocate tracks arrays if using explicit ray tracing */
  if (_segment_formation == EXPLICIT_3D) {

    log_printf(NORMAL, "Explicit 3D Track storage = %.2f MB", num_tracks *
               sizeof(Track3D) / 1e6);

    _tracks_3D = new Track3D***[_num_azim/2];
    for (int a=0; a < _num_azim/2; a++) {
      _tracks_3D[a] = new Track3D**[_num_x[a] + _num_y[a]];
      for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
        _tracks_3D[a][i] = new Track3D*[_num_polar];
        for (int p=0; p < _num_polar; p++) {
          _tracks_3D[a][i][p] = new Track3D[_tracks_per_stack[a][i][p]];
        }
      }
    }

    /* Reset the number of tracks per stack */
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
        for (int p=0; p < _num_polar; p++)
          _tracks_per_stack[a][i][p] = 0;
      }
    }
  }
}


/**
 * @brief Returns the 3D Track ID based on its indexes.
 * @param tsi The indexes of the Track of interest (azimuthal, xy track, polar,
 *        and z stack)
 */
long TrackGenerator3D::get3DTrackID(TrackStackIndexes* tsi) {
  return _cum_tracks_per_stack[tsi->_azim][tsi->_xy][tsi->_polar] + tsi->_z;
}


/**
 * @brief Allocates a new Quadrature with the default Quadrature.
 * @details The default quadrature for 3D calculations is equal weight
 */
void TrackGenerator3D::initializeDefaultQuadrature() {

  if (_quadrature != NULL)
    delete _quadrature;

  log_printf(NORMAL, "Initializing a default angular quadrature...");
  _quadrature = new GLPolarQuad();
  _quadrature->setNumAzimAngles(_num_azim);
  _quadrature->setNumPolarAngles(_num_polar);
}


/**
 * @brief Returns the filename for writing tracking data.
 */
std::string TrackGenerator3D::getTestFilename(std::string directory) {

  std::stringstream test_filename;

  /* Get the quadrature and track method types */
  std::string quad_type;
  std::string track_method;

  if (_quadrature->getQuadratureType() == EQUAL_WEIGHT)
    quad_type = "EQ_WGT";
  else if (_quadrature->getQuadratureType() == EQUAL_ANGLE)
    quad_type = "EQ_ANG";
  if (_quadrature->getQuadratureType() == TABUCHI_YAMAMOTO)
    quad_type = "TABUCHI_YAMAMOTO";
  else if (_quadrature->getQuadratureType() == LEONARD)
    quad_type = "LEONARD";
  else if (_quadrature->getQuadratureType() == GAUSS_LEGENDRE)
    quad_type = "GAUSS_LEGENDRE";

  if (_geometry->getCmfd() != NULL)
    test_filename << directory << "/3D_"
                  << _num_azim << "_azim_"
                  << _num_polar << "_polar_"
                  << _azim_spacing << "x" << _z_spacing
                  << "_cm_spacing_cmfd_"
                  << _geometry->getCmfd()->getNumX()
                  << "x" << _geometry->getCmfd()->getNumY()
                  << "x" << _geometry->getCmfd()->getNumZ()
                  << "_quad_" << quad_type << "_track_" << track_method;
  else
    test_filename << directory << "/3D_"
                  << _num_azim << "_azim_"
                  << _num_polar << "_polar_"
                  << _azim_spacing << "x" << _z_spacing
                  << "_cm_spacing_quad_"
                  << quad_type << "_track_" << track_method;

  if (_geometry->isDomainDecomposed()) {
    test_filename << "_dd";
    struct stat st;
    if (!(stat(test_filename.str().c_str(), &st) == 0))
      mkdir(test_filename.str().c_str(), S_IRWXU);
    int indexes[3];
    _geometry->getDomainIndexes(indexes);
    test_filename << "/domain";
    for (int i=0; i<3; i++)
      test_filename << "_" << indexes[i];
  }

  test_filename << ".data";

  return test_filename.str();
}


/**
 * @brief Updates whether the TrackGenerator contains segments.
 * @param contains_segments whether the TrackGenerator contains segments
 */
void TrackGenerator3D::setContainsSegments(bool contains_segments) {
  if (_segment_formation == EXPLICIT_3D)
    _contains_3D_segments = contains_segments;
  else
    _contains_2D_segments = contains_segments;
}


/**
 * @brief Checks the boundary conditions for all 3D surfaces for inconsistent
 *        periodic boundary conditions.
 */
void TrackGenerator3D::checkBoundaryConditions() {

  /* Check X and Y boundaries */
  TrackGenerator::checkBoundaryConditions();

  /* Check Z boundaries for consistency */
  if ((_geometry->getMinZBoundaryType() == PERIODIC &&
        _geometry->getMaxZBoundaryType() != PERIODIC) ||
       (_geometry->getMinZBoundaryType() != PERIODIC &&
        _geometry->getMaxZBoundaryType() == PERIODIC))
    log_printf(ERROR, "Cannot create tracks with only one z boundary"
               " set to PERIODIC");

  /* Check that there are no periodic boundaries if domain decomposed */
  if (_geometry->isDomainDecomposed())
    if (_geometry->getMinZBoundaryType() == PERIODIC)
      log_printf(ERROR, "Periodic boundaries are not supported for domain "
                 "decomposition");

  /* Check for correct track method if a PERIODIC bc is present */
  if (_geometry->getMinZBoundaryType() == PERIODIC ||
      _geometry->getMinZBoundaryType() == INTERFACE ||
      _geometry->getMaxZBoundaryType() == INTERFACE)
    _periodic = true;
}


/**
 * @brief Allocates memory for temporary segment storage if necessary.
 * @details New memory is only allocated if _max_num_segments exceeds
 *          _num_seg_matrix_columns (the maximum when the segments were allocated)
 */
void TrackGenerator3D::allocateTemporarySegments() {

  /* Check if a resize is unnecessary */
  if (_max_num_segments <= _num_seg_matrix_columns)
    return;

  /* Widen the number of columns in the temporary segments matrix */
  _num_seg_matrix_columns = _max_num_segments;

  /* Delete temporary segments if already allocated */
  if (_contains_temporary_segments) {
    for (int t = 0; t < _num_threads; t++) {
      delete [] _temporary_segments.at(t);
    }
  }
  else {
    _temporary_segments.resize(_num_threads);
    _contains_temporary_segments = true;
  }

  /* Determine storage size of temporary segments */
  long max_size = _num_seg_matrix_columns;
#ifdef MPIX
  if (_geometry->isDomainDecomposed())
    MPI_Allreduce(&_num_seg_matrix_columns, &max_size, 1, MPI_LONG, MPI_MAX,
                  _geometry->getMPICart());
#endif
  double max_size_mb = (double) (max_size * _num_threads * sizeof(segment))
      / (double) (1e6);

  log_printf(INFO_ONCE, "Max temporary segment storage per domain = %6.2f MB",
             max_size_mb);

  /* Allocate new temporary segments */
  for (int t = 0; t < _num_threads; t++)
    _temporary_segments.at(t) = new segment[_num_seg_matrix_columns];
}


/**
 * @brief Allocates memory for temporary Track storage if necessary.
 */
void TrackGenerator3D::allocateTemporaryTracks() {

  /* Delete temporary tracks if already allocated */
  if (_contains_temporary_tracks) {
    for (int t = 0; t < _num_threads; t++) {
      delete [] _temporary_3D_tracks.at(t);
      delete [] _temporary_tracks_array.at(t);
    }
  }
  else {
    _temporary_3D_tracks.resize(_num_threads);
    _temporary_tracks_array.resize(_num_threads);
    _contains_temporary_tracks = true;
  }

  /* Report memory usage */
  double size_mb = (double) (_num_threads * _max_num_tracks_per_stack
        * sizeof(Track3D)) / (double) 1e6;
  log_printf(INFO_ONCE, "Temporary Track storage per domain = %6.2f MB",
             size_mb);

  /* Allocate new temporary tracks arrays */
  for (int t = 0; t < _num_threads; t++) {
    _temporary_3D_tracks.at(t) = new Track3D[_max_num_tracks_per_stack];
    _temporary_tracks_array.at(t) = new Track*[_max_num_tracks_per_stack];
    for (int i=0; i < _max_num_tracks_per_stack; i++)
      _temporary_tracks_array.at(t)[i] = &_temporary_3D_tracks.at(t)[i];
  }
}


/**
 * @brief Resets the TrackGenerator to not contain tracks or segments.
 */
void TrackGenerator3D::resetStatus() {
  TrackGenerator::resetStatus();
  _contains_3D_tracks = false;
  _contains_3D_segments = false;
}


/**
 * @brief Returns whether or not the TrackGenerator contains Tracks
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains Tracks; false otherwise
 */
bool TrackGenerator3D::containsTracks() {
  return _contains_3D_tracks;
}


/**
 * @brief Returns whether or not the TrackGenerator contains 3D segments
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains 3D segments; false otherwise
 */
bool TrackGenerator3D::containsSegments() {
  if (_segment_formation == EXPLICIT_3D)
    return _contains_3D_segments;
  else
    return _contains_2D_segments;
}


/**
 * @brief Updates the provided Track with the correct information based on the
 *        Track indexes.
 * @param track The 3D Track whose information is updated
 * @param tsi The indexes of the 3D Track
 */
void TrackGenerator3D::getTrackOTF(Track3D* track, TrackStackIndexes* tsi) {

  try {
    Track* track_2D = &_tracks_2D[tsi->_azim][tsi->_xy];
    TrackChainIndexes tci;
    convertTSItoTCI(tsi, &tci);

    /* Set the start and end points */
    if (_segment_formation == EXPLICIT_3D) {
      Track3D* track_3D =
        &_tracks_3D[tsi->_azim][tsi->_xy][tsi->_polar][tsi->_z];
      Point* start_3d_1 = track_3D->getStart();
      Point* end_3d_1 = track_3D->getEnd();
      Point* start_3d_2 = track->getStart();
      Point* end_3d_2 = track->getEnd();

      double x1 = start_3d_1->getX();
      double y1 = start_3d_1->getY();
      double z1 = start_3d_1->getZ();
      double x2 = end_3d_1->getX();
      double y2 = end_3d_1->getY();
      double z2 = end_3d_1->getZ();

      start_3d_2->setCoords(x1, y1, z1);
      end_3d_2->setCoords(x2, y2, z2);
    }
    else
      set3DTrackData(&tci, track, false, false);

    /* Set the angles */
    track->setPhi(track_2D->getPhi());
    track->setTheta(_quadrature->getTheta(tsi->_azim, tsi->_polar));
    track->setAzimIndex(tsi->_azim);
    track->setXYIndex(tsi->_xy);
    track->setPolarIndex(tsi->_polar);

    /* Set the linking Track data */
    setLinkingTracks(tsi, &tci, true, track);
    setLinkingTracks(tsi, &tci, false, track);

    /* Set the UID for the current Track */
    track->setUid(get3DTrackID(tsi));
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to get track otf (%d, %d, %d, %d)."
               " Backtrace:\n%s", tsi->_azim, tsi->_xy, tsi->_polar,
               tsi->_z, e.what());
  }
}


/**
 * @brief Get the array of 3D Tracks, only used in EXPLICT-3D ray tracing.
 * @return a pointer to the array of 3D Tracks
 */
Track3D**** TrackGenerator3D::get3DTracks() {
  return _tracks_3D;
}


/**
 * @brief Converts TrackStackIndexes to TrackChainIndexes.
 * @param tsi The TrackStackIndexes of the 3D Track
 * @param tci The TrackChainIndexes of the 3D Track to be updated
 */
void TrackGenerator3D::convertTSItoTCI(TrackStackIndexes* tsi,
                                       TrackChainIndexes* tci) {

  tci->_azim = tsi->_azim;
  tci->_x = tsi->_xy % _num_x[tsi->_azim];
  tci->_polar = tsi->_polar;
  tci->_lz = _first_lz_of_stack[tsi->_azim][tsi->_xy][tsi->_polar] + tsi->_z;

  setLinkIndex(tci, tsi);
}


/**
 * @brief Converts TrackChainIndexes to TrackStackIndexes.
 * @param tci The TrackChainIndexes of the 3D Track
 * @param tsi The TrackStackIndexes of the 3D Track to be updated
 */
void TrackGenerator3D::convertTCItoTSI(TrackChainIndexes* tci,
                                       TrackStackIndexes* tsi) {

  int link = getLinkIndex(tci);
  Track* track_2D = _tracks_2D_chains[tci->_azim][tci->_x][link];
  tsi->_azim = tci->_azim;
  tsi->_xy = track_2D->getXYIndex();
  tsi->_polar = tci->_polar;
  tsi->_z = tci->_lz - _first_lz_of_stack[tsi->_azim][tsi->_xy][tsi->_polar];
}


/**
 * @brief Returns the index into the chain based on chain indexes.
 * @param tci The chain indexes
 * @return the index into the 3D z-stack
 */
int TrackGenerator3D::getLinkIndex(TrackChainIndexes* tci) {
  Track3D track;
  return getFirst2DTrackLinkIndex(tci, &track) + tci->_link;
}


/**
 * @brief Updates the index into the 3D chain based on chain and stack indexes.
 * @param tci The chain indexes to be updated
 * @param tsi The stack indexes
 */
void TrackGenerator3D::setLinkIndex(TrackChainIndexes* tci,
                                    TrackStackIndexes* tsi) {
  Track3D track;
  int first_link = getFirst2DTrackLinkIndex(tci, &track);
  Track* track_2D = &_tracks_2D[tsi->_azim][tsi->_xy];
  tci->_link = track_2D->getLinkIndex() - first_link;
}


/**
 * @brief Calculates the number of Tracks in the l-z traversal.
 * @param tci The chain indexes
 * @return the number of Tracks in the l-z traversal
 */
int TrackGenerator3D::getNum3DTrackChainLinks(TrackChainIndexes* tci) {

  Track3D track_3D;
  int first_link = getFirst2DTrackLinkIndex(tci, &track_3D);
  int link = first_link;
  Track* track_2D = _tracks_2D_chains[tci->_azim][tci->_x][link];
  int azim = track_2D->getAzimIndex();
  int xy = track_2D->getXYIndex();
  int polar = tci->_polar;
  int min_lz, max_lz;
  int lz = tci->_lz;

  while (true) {
    track_2D = _tracks_2D_chains[tci->_azim][tci->_x][link];
    azim = track_2D->getAzimIndex();
    xy = track_2D->getXYIndex();

    min_lz = _first_lz_of_stack[azim][xy][polar];
    max_lz = _tracks_per_stack[azim][xy][polar] + min_lz - 1;

    if (tci->_polar < _num_polar / 2 && tci->_lz > max_lz)
      break;
    else if (tci->_polar >= _num_polar / 2 && tci->_lz < min_lz)
      break;

    link++;

    if (track_2D->getXYIndex() >= _num_y[azim])
      break;
  }

  return link - first_link;
}


/**
 * @brief Fills the provided 3D Track with its linking information.
 * @param tsi The stack indexes
 * @param tci The chain indexes
 * @param outgoing A boolean indicating the direction of the Track
 *        (True = forward, False = Backward)
 * @param track The 3D Track whose information is updated
 */
void TrackGenerator3D::setLinkingTracks(TrackStackIndexes* tsi,
                                        TrackChainIndexes* tci,
                                        bool outgoing,
                                        Track3D* track) {

  Track* track_2D = &_tracks_2D[tsi->_azim][tsi->_xy];
  TrackChainIndexes tci_next;
  TrackChainIndexes tci_prdc;
  TrackChainIndexes tci_refl;
  TrackStackIndexes tsi_next;
  TrackStackIndexes tsi_prdc;
  TrackStackIndexes tsi_refl;

  /* Set the next TCI to the current TCI */
  tci_next._azim  = tci->_azim;
  tci_next._x = tci->_x;
  tci_next._polar = tci->_polar;
  tci_next._lz = tci->_lz;
  tci_next._link = 0;

  /* Set the periodic TCI to the current TCI */
  tci_prdc._azim  = tci->_azim;
  tci_prdc._x = tci->_x;
  tci_prdc._polar = tci->_polar;
  tci_prdc._lz = tci->_lz;
  tci_prdc._link = 0;

  /* Set the reflective TCI to the current TCI */
  tci_refl._azim  = tci->_azim;
  tci_refl._x     = tci->_x;
  tci_refl._polar = tci->_polar;
  tci_refl._lz    = tci->_lz;
  tci_refl._link = 0;

  int nz = _num_z[tci->_azim][tci->_polar];
  int nl = _num_l[tci->_azim][tci->_polar];
  int lz = tci->_lz;
  int ac = _num_azim / 2 - tci->_azim - 1;
  int pc = _num_polar - tci->_polar - 1;
  int xc = _num_x[tci->_azim] - tci->_x - 1;

  int num_links = getNum3DTrackChainLinks(tci);
  bool next_fwd = outgoing;
  bool refl_fwd = outgoing;
  boundaryType bc;

  int surface_2D;
  if (outgoing) {
    bc = track_2D->getBCFwd();
    surface_2D = track_2D->getSurfaceOut();
  }
  else {
    bc = track_2D->getBCBwd();
    surface_2D = track_2D->getSurfaceIn();
  }

  /* Initialize the Track's connecting domains with their 2D connections */
#ifdef MPIx
  int domain_delta_x = (surface_2D % 3 == 0) * (2 * (surface_2D/3) - 1);
  int domain_delta_y = (surface_2D % 3 == 1) * (2 * (surface_2D/3) - 1);
  int domain_delta_z = 0;
#endif

  /* Tracks pointing in the positive z direction in the lz plane */
  if (tci->_polar < _num_polar / 2) {

    /* SURFACE_Z_MAX */
    if (tci->_link == num_links - 1 && lz >= nz && outgoing) {

      bc = _geometry->getMaxZBoundaryType();

      tci_prdc._lz    = lz - nz;
      tci_refl._polar = pc;
      tci_refl._lz    = nl + 2 * nz - lz - 1;

      /* PERIODIC or INTERFACE BC */
      if (_geometry->getMaxZBoundaryType() == PERIODIC ||
          _geometry->getMaxZBoundaryType() == INTERFACE)
        tci_next._lz    = lz - nz;

      /* REFLECTIVE OR VACUUM BC */
      else {
        tci_next._polar = pc;
        tci_next._lz    = nl + 2 * nz - lz - 1;
      }
#ifdef MPIx
      domain_delta_x = 0;
      domain_delta_y = 0;
      domain_delta_z = +1;
#endif

      /* Check for a double reflection */
      convertTCItoTSI(&tci_prdc, &tsi_prdc);

      if (tsi_prdc._xy != tsi->_xy) {
#ifdef MPIx
        if (track_2D->getBCFwd() == INTERFACE) {
          if (tsi->_azim < _num_azim / 4)
            domain_delta_x = +1;
          else
            domain_delta_x = -1;
        }
#endif
        tci_refl._azim = ac;

        /* Set the next Track */
        boundaryType bc_xy = track_2D->getBCFwd();
        if (bc_xy != PERIODIC && bc_xy != INTERFACE) {
          tci_next._azim = ac;
        }
        if (bc_xy == INTERFACE && bc != VACUUM) {
#ifdef MPIx
          if (bc != INTERFACE)
            domain_delta_z = 0;
#endif
          bc = INTERFACE;
        }
        else if (bc_xy == VACUUM)
          bc = VACUUM;
      }
    }

    /* SURFACE_Z_MIN */
    else if (tci->_link == 0 && lz < nl && !outgoing) {

      bc = _geometry->getMinZBoundaryType();

      tci_prdc._lz    = lz + nz;
      tci_prdc._link = getNum3DTrackChainLinks(&tci_prdc) - 1;
      tci_refl._polar = pc;
      tci_refl._lz    = nl - lz - 1;
      tci_refl._link = getNum3DTrackChainLinks(&tci_refl) - 1;

      /* PERIODIC or INTERFACE BC */
      if (_geometry->getMinZBoundaryType() == PERIODIC ||
          _geometry->getMinZBoundaryType() == INTERFACE) {
        tci_next._lz    = lz + nz;
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
      }

      /* REFLECTIVE OR VACUUM BC */
      else {
        tci_next._polar = pc;
        tci_next._lz    = nl - lz - 1;
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
      }
#ifdef MPIx
      domain_delta_x = 0;
      domain_delta_y = 0;
      domain_delta_z = -1;
#endif

      /* Check for a double reflection */
      convertTCItoTSI(&tci_prdc, &tsi_prdc);

      if (tsi_prdc._xy != tsi->_xy) {
#ifdef MPIx
        if (track_2D->getBCBwd() == INTERFACE) {
          if (tsi->_azim < _num_azim / 4)
            domain_delta_x = -1;
          else
            domain_delta_x = +1;
        }
#endif
        tci_refl._azim = ac;

        /* Set the next Track */
        boundaryType bc_xy = track_2D->getBCBwd();
        if (bc_xy != PERIODIC && bc_xy != INTERFACE) {
          tci_next._azim = ac;
        }
        if (bc_xy == INTERFACE && bc != VACUUM) {
#ifdef MPIx
          if (bc != INTERFACE)
            domain_delta_z = 0;
#endif
          bc = INTERFACE;
        }
        else if (bc_xy == VACUUM)
          bc = VACUUM;
      }
    }

    /* SURFACE_Y_MIN */
    else if (tci->_link == 0 && lz >= nl && !outgoing) {

      tci_prdc._lz    = lz - nl;
      tci_prdc._x     = _tracks_2D_array[track_2D->getTrackPrdcBwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_prdc._link = getNum3DTrackChainLinks(&tci_prdc) - 1;
      tci_refl._azim  = ac;
      tci_refl._x     = _tracks_2D_array[track_2D->getTrackReflBwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_refl._polar = pc;
      refl_fwd = true;
      tci_refl._lz = lz - nl;
      tci_next._lz = lz - nl;

      if (_geometry->getMinYBoundaryType() == PERIODIC ||
          _geometry->getMinYBoundaryType() == INTERFACE) {
        tci_next._x     = _tracks_2D_array[track_2D->getTrackPrdcBwd()]
            ->getXYIndex() % _num_x[tci->_azim];
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
      }
      else {
        tci_next._azim  = ac;
        tci_next._x     = _tracks_2D_array[track_2D->getTrackReflBwd()]
            ->getXYIndex() % _num_x[tci->_azim];
        tci_next._polar = pc;
        next_fwd = true;
      }
    }

    /* SURFACE_Y_MAX */
    else if (tci->_link == num_links - 1 && lz < nz && outgoing) {

      tci_prdc._lz    = nl + lz;
      tci_prdc._x     = _tracks_2D_array[track_2D->getTrackPrdcFwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_refl._azim  = ac;
      tci_refl._x     = _tracks_2D_array[track_2D->getTrackReflFwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_refl._polar = pc;
      tci_refl._lz = nl + lz;
      tci_refl._link = getNum3DTrackChainLinks(&tci_refl) - 1;
      tci_next._lz = nl + lz;
      refl_fwd = false;

      if (_geometry->getMaxYBoundaryType() == PERIODIC ||
          _geometry->getMaxYBoundaryType() == INTERFACE)
        tci_next._x     = _tracks_2D_array[track_2D->getTrackPrdcFwd()]
            ->getXYIndex() % _num_x[tci->_azim];
      else {
        tci_next._azim  = ac;
        tci_next._x     = _tracks_2D_array[track_2D->getTrackReflFwd()]
            ->getXYIndex() % _num_x[tci->_azim];
        tci_next._polar = pc;
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
        next_fwd = false;
      }
    }

    /* SURFACE_X_MIN or SURFACE_X_MAX */
    else if (outgoing) {

      /* Set the link index */
      tci_prdc._link = tci->_link + 1;
      tci_next._link = tci->_link + 1;
      tci_refl._link = tci->_link + 1;
      tci_refl._azim  = ac;

      /* Set the next track */
      if (track_2D->getBCFwd() != PERIODIC &&
          track_2D->getBCFwd() != INTERFACE)
        tci_next._azim  = ac;
    }

    /* Tracks hitting any of the four x or y surfaces */
    else {

      /* Set the link index */
      tci_prdc._link = tci->_link - 1;
      tci_next._link = tci->_link - 1;
      tci_refl._link = tci->_link - 1;
      tci_refl._azim  = ac;

      /* Set the next track */
      if (track_2D->getBCBwd() != PERIODIC &&
          track_2D->getBCBwd() != INTERFACE)
        tci_next._azim  = ac;
    }
  }
  else {

    /* SURFACE_Z_MAX */
    if (tci->_link == 0 && lz >= nz && !outgoing) {

      bc = _geometry->getMaxZBoundaryType();

      tci_prdc._lz    = lz - nz;
      tci_prdc._link = getNum3DTrackChainLinks(&tci_prdc) - 1;
      tci_refl._polar = pc;
      tci_refl._lz    = nl + 2 * nz - lz - 1;
      tci_refl._link = getNum3DTrackChainLinks(&tci_refl) - 1;

      /* PERIODIC or INTERFACE BC */
      if (_geometry->getMaxZBoundaryType() == PERIODIC ||
          _geometry->getMaxZBoundaryType() == INTERFACE) {
        tci_next._lz    = lz - nz;
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
      }
      /* REFLECTIVE OR VACUUM BC */
      else {
        tci_next._polar = pc;
        tci_next._lz    = nl + 2 * nz - lz - 1;
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
      }
#ifdef MPIx
      domain_delta_x = 0;
      domain_delta_y = 0;
      domain_delta_z = +1;
#endif

      /* Check for a double reflection */
      convertTCItoTSI(&tci_prdc, &tsi_prdc);

      if (tsi_prdc._xy != tsi->_xy) {
#ifdef MPIx
        if (track_2D->getBCBwd() == INTERFACE) {
          if (tsi->_azim < _num_azim / 4)
            domain_delta_x = -1;
          else
            domain_delta_x = +1;
        }
#endif
        tci_refl._azim = ac;

        /* Set the next Track */
        boundaryType bc_xy = track_2D->getBCBwd();
        if (bc_xy != PERIODIC && bc_xy != INTERFACE) {
          tci_next._azim = ac;
        }
        if (bc_xy == INTERFACE && bc != VACUUM) {
#ifdef MPIx
          if (bc != INTERFACE)
            domain_delta_z = 0;
#endif
          bc = INTERFACE;
        }
        else if (bc_xy == VACUUM)
          bc = VACUUM;
      }
    }

    /* SURFACE_Z_MIN */
    else if (tci->_link == num_links - 1 && lz < nl && outgoing) {

      bc = _geometry->getMinZBoundaryType();

      tci_prdc._lz    = lz + nz;
      tci_refl._polar = pc;
      tci_refl._lz    = nl - lz - 1;

      /* PERIODIC or INTERFACE BC */
      if (_geometry->getMinZBoundaryType() == PERIODIC ||
          _geometry->getMinZBoundaryType() == INTERFACE)
        tci_next._lz    = lz + nz;

      /* REFLECTIVE OR VACUUM BC */
      else {
        tci_next._polar = pc;
        tci_next._lz    = nl - lz - 1;
      }
#ifdef MPIx
      domain_delta_x = 0;
      domain_delta_y = 0;
      domain_delta_z = -1;
#endif

      /* Check for a double reflection */
      convertTCItoTSI(&tci_prdc, &tsi_prdc);

      if (tsi_prdc._xy != tsi->_xy) {
#ifdef MPIx
        if (track_2D->getBCFwd() == INTERFACE) {
          if (tsi->_azim < _num_azim / 4)
            domain_delta_x = +1;
          else
            domain_delta_x = -1;
        }
#endif
        tci_refl._azim = ac;

        /* Set the next Track */
        boundaryType bc_xy = track_2D->getBCFwd();
        if (bc_xy != PERIODIC && bc_xy != INTERFACE) {
          tci_next._azim = ac;
        }
        if (bc_xy == INTERFACE && bc != VACUUM) {
#ifdef MPIx
          if (bc != INTERFACE)
            domain_delta_z = 0;
#endif
          bc = INTERFACE;
        }
        else if (bc_xy == VACUUM)
          bc = VACUUM;
      }
    }

    /* SURFACE_Y_MIN */
    else if (tci->_link == 0 && lz < nz && !outgoing) {

      tci_prdc._lz    = lz + nl;
      tci_prdc._x     = _tracks_2D_array[track_2D->getTrackPrdcBwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_prdc._link = getNum3DTrackChainLinks(&tci_prdc) - 1;
      tci_refl._azim  = ac;
      tci_refl._x     = _tracks_2D_array[track_2D->getTrackReflBwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_refl._polar = pc;
      tci_next._lz    = lz + nl;
      tci_refl._lz    = lz + nl;
      refl_fwd = true;

      if (_geometry->getMinYBoundaryType() == PERIODIC ||
          _geometry->getMinYBoundaryType() == INTERFACE) {
        tci_next._x     = _tracks_2D_array[track_2D->getTrackPrdcBwd()]
            ->getXYIndex() % _num_x[tci->_azim];
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
      }
      else {
        tci_next._azim  = ac;
        tci_next._x     = _tracks_2D_array[track_2D->getTrackReflBwd()]
            ->getXYIndex() % _num_x[tci->_azim];
        tci_next._polar = pc;
        next_fwd = true;
      }
    }

    /* SURFACE_Y_MAX */
    else if (tci->_link == num_links - 1 && lz >= nl && outgoing) {

      tci_prdc._lz    = lz - nl;
      tci_prdc._x     = _tracks_2D_array[track_2D->getTrackPrdcFwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_refl._azim  = ac;
      tci_refl._x     = _tracks_2D_array[track_2D->getTrackReflFwd()]
          ->getXYIndex() % _num_x[tci->_azim];
      tci_refl._polar = pc;
      tci_refl._lz    = lz - nl;
      tci_refl._link = getNum3DTrackChainLinks(&tci_refl) - 1;
      tci_next._lz    = lz - nl;
      refl_fwd = false;

      if (_geometry->getMaxYBoundaryType() == PERIODIC ||
          _geometry->getMaxYBoundaryType() == INTERFACE)
        tci_next._x     = _tracks_2D_array[track_2D->getTrackPrdcFwd()]
            ->getXYIndex() % _num_x[tci->_azim];
      else {
        tci_next._azim  = ac;
        tci_next._x     = _tracks_2D_array[track_2D->getTrackReflFwd()]
            ->getXYIndex() % _num_x[tci->_azim];
        tci_next._polar = pc;
        tci_next._link = getNum3DTrackChainLinks(&tci_next) - 1;
        next_fwd = false;
      }
    }

    /* SURFACE_X_MIN or SURFACE_X_MAX */
    else if (outgoing) {

      /* Set the link index */
      tci_prdc._link = tci->_link + 1;
      tci_next._link = tci->_link + 1;
      tci_refl._link = tci->_link + 1;
      tci_refl._azim  = ac;

      /* Set the next track */
      if (track_2D->getBCFwd() != PERIODIC &&
          track_2D->getBCFwd() != INTERFACE)
        tci_next._azim  = ac;
    }

    /* Tracks hitting any of the four x or y surfaces */
    else {

      /* Set the link index */
      tci_prdc._link = tci->_link - 1;
      tci_next._link = tci->_link - 1;
      tci_refl._link = tci->_link - 1;
      tci_refl._azim  = ac;

      /* Set the next track */
      if (track_2D->getBCBwd() != PERIODIC &&
          track_2D->getBCBwd() != INTERFACE)
        tci_next._azim  = ac;
    }
  }

  convertTCItoTSI(&tci_next, &tsi_next);
  convertTCItoTSI(&tci_prdc, &tsi_prdc);
  convertTCItoTSI(&tci_refl, &tsi_refl);

  if (outgoing) {
    track->setTrackNextFwd(get3DTrackID(&tsi_next));
    track->setTrackPrdcFwd(get3DTrackID(&tsi_prdc));
    track->setTrackReflFwd(get3DTrackID(&tsi_refl));
    track->setNextFwdFwd(next_fwd);
    track->setBCFwd(bc);
#ifdef MPIx
    int domain_fwd = _geometry->getNeighborDomain(domain_delta_x,
                                                  domain_delta_y,
                                                  domain_delta_z);
    track->setDomainFwd(domain_fwd);
#endif
  }
  else {
    track->setTrackNextBwd(get3DTrackID(&tsi_next));
    track->setTrackPrdcBwd(get3DTrackID(&tsi_prdc));
    track->setTrackReflBwd(get3DTrackID(&tsi_refl));
    track->setNextBwdFwd(next_fwd);
    track->setBCBwd(bc);
#ifdef MPIx
    int domain_bwd = _geometry->getNeighborDomain(domain_delta_x,
                                                  domain_delta_y,
                                                  domain_delta_z);
    track->setDomainBwd(domain_bwd);
#endif
  }
}


/**
 * @brief Updates stack indexes to reflect those from the given Track ID.
 * @param id The 3D Track ID
 * @param tsi The stack indexes to be updated
 */
void TrackGenerator3D::getTSIByIndex(long id, TrackStackIndexes* tsi) {

  long cum_track_index = 0;
  int a;
  for (a=0; a < _num_azim/2; a++) {
    cum_track_index += _tracks_per_azim[a];
    if (id < cum_track_index)
      break;
  }
  tsi->_azim = a;
  tsi->_xy = -1;

  for (int xy=0; xy < _num_x[tsi->_azim] + _num_y[tsi->_azim]; xy++) {
    if (id < _cum_tracks_per_xy[tsi->_azim][xy]) {
      tsi->_xy = xy;
      break;
    }
  }

  if (tsi->_xy == -1)
    log_printf(ERROR, "could not generate TSI xy from track ID: %ld", id);

  for (int p=0; p < _num_polar; p++) {
    if (id < _cum_tracks_per_stack[tsi->_azim][tsi->_xy][p] +
        _tracks_per_stack[tsi->_azim][tsi->_xy][p]) {
      tsi->_polar = p;
      tsi->_z = id - _cum_tracks_per_stack[tsi->_azim][tsi->_xy][p];
      return;
    }
  }

  log_printf(ERROR, "could not generate TSI from track ID: %ld", id);
}


/**
 * @brief Write information of all Extruded FSRs to a file.
 * @param out file to write to
 */
void TrackGenerator3D::writeExtrudedFSRInfo(FILE* out) {

  /* Module to write track info */
  DumpSegments dump_segments(this);
  dump_segments.setOutputFile(out);

  /* Write extruded FSR data */
    ParallelHashMap<std::string, ExtrudedFSR*>& extruded_FSR_keys_map =
        _geometry->getExtrudedFSRKeysMap();
    std::string* extruded_fsr_key_list = extruded_FSR_keys_map.keys();
    ExtrudedFSR** extruded_fsr_list = extruded_FSR_keys_map.values();

    /* Write number of extruded FSRs */
    int num_extruded_FSRs = extruded_FSR_keys_map.size();
    fwrite(&num_extruded_FSRs, sizeof(int), 1, out);

    /* Write extruded FSR data */
    for (int i=0; i < num_extruded_FSRs; i++) {

      std::string key = extruded_fsr_key_list[i];
      int string_length = key.length() + 1;
      fwrite(&string_length, sizeof(int), 1, out);
      fwrite(key.c_str(), sizeof(char)*string_length, 1, out);

      ExtrudedFSR* extruded_fsr = extruded_fsr_list[i];

      int extruded_fsr_id = extruded_fsr->_fsr_id;
      fwrite(&extruded_fsr_id, sizeof(int), 1, out);

      double x = extruded_fsr->_coords->getX();
      double y = extruded_fsr->_coords->getY();
      double z = extruded_fsr->_coords->getZ();
      fwrite(&x, sizeof(double), 1, out);
      fwrite(&y, sizeof(double), 1, out);
      fwrite(&z, sizeof(double), 1, out);

      int num_fsrs = extruded_fsr->_num_fsrs;
      fwrite(&num_fsrs, sizeof(int), 1, out);

      double init_mesh_val = extruded_fsr->_mesh[0];
      fwrite(&init_mesh_val, sizeof(double), 1, out);

      for (int j=0; j < num_fsrs; j++) {
        long fsr_id = extruded_fsr->_fsr_ids[j];
        fwrite(&fsr_id, sizeof(long), 1, out);
        double mesh_val = extruded_fsr->_mesh[j+1];
        fwrite(&mesh_val, sizeof(double), 1, out);
      }
    }

    /* Delete extruded FSR key and value lists */
    delete [] extruded_fsr_key_list;
    delete [] extruded_fsr_list;

    /* Record 2D track info */
    int num_2D_tracks = getNum2DTracks();
    fwrite(&num_2D_tracks, sizeof(int), 1, out);
    for (int i=0; i < num_2D_tracks; i++) {
      Track* track = _tracks_2D_array[i];
      segment* segments = track->getSegments();
      dump_segments.onTrack(track, segments);
    }

    /* Record maximum number of segments */
    fwrite(&_max_num_segments, sizeof(int), 1, out);
}


/**
 * @brief Read information of all Extruded FSRs from a file.
 * @param in file to read from
 */
void TrackGenerator3D::readExtrudedFSRInfo(FILE* in) {

  /* Module to read track info */
  ReadSegments read_segments(this);
  read_segments.setInputFile(in);

    /* Read number of extruded FSRs */
    ParallelHashMap<std::string, ExtrudedFSR*>& extruded_FSR_keys_map =
        _geometry->getExtrudedFSRKeysMap();
    int num_extruded_FSRs;
    int ret = _geometry->twiddleRead(&num_extruded_FSRs, sizeof(int), 1, in);

    /* Resize for the number of extruded FSRs */
    std::vector<ExtrudedFSR*>& extruded_FSR_lookup =
        _geometry->getExtrudedFSRLookup();
    extruded_FSR_lookup.resize(num_extruded_FSRs);
    extruded_FSR_keys_map.realloc(2*num_extruded_FSRs);

    /* Write extruded FSR data */
    for (int i=0; i < num_extruded_FSRs; i++) {

      /* Read the extruded FSR key */
      int string_length;
      ret = _geometry->twiddleRead(&string_length, sizeof(int), 1, in);
      char* char_buffer2 = new char[string_length];
      ret = _geometry->twiddleRead(char_buffer2, sizeof(char)*string_length, 1,
                                   in);
      std::string key = std::string(char_buffer2);

      /* Create new extruded FSR and add to map */
      ExtrudedFSR* extruded_fsr = new ExtrudedFSR;
      extruded_FSR_keys_map.insert(key, extruded_fsr);

      /* Read ID */
      int extruded_fsr_id;
      ret = _geometry->twiddleRead(&extruded_fsr_id, sizeof(int), 1, in);
      extruded_fsr->_fsr_id = extruded_fsr_id;

      /* Read coordinates */
      double x, y, z;
      ret = _geometry->twiddleRead(&x, sizeof(double), 1, in);
      ret = _geometry->twiddleRead(&y, sizeof(double), 1, in);
      ret = _geometry->twiddleRead(&z, sizeof(double), 1, in);
      LocalCoords* coords = new LocalCoords(x,y,z);
      coords->setUniverse(_geometry->getRootUniverse());
      extruded_fsr->_coords = coords;

      /* Read the number of FSRs and allocate memory */
      int num_fsrs;
      ret = _geometry->twiddleRead(&num_fsrs, sizeof(int), 1, in);
      extruded_fsr->_num_fsrs = num_fsrs;
      extruded_fsr->_materials = new Material*[num_fsrs];
      extruded_fsr->_fsr_ids = new long[num_fsrs];
      extruded_fsr->_mesh = new double[num_fsrs+1];

      /* Read the mesh values and FSR IDs */
      double init_mesh_val;
      ret = _geometry->twiddleRead(&init_mesh_val, sizeof(double), 1, in);
      extruded_fsr->_mesh[0] = init_mesh_val;

      for (int j=0; j < num_fsrs; j++) {
        long fsr_id;
        ret = _geometry->twiddleRead(&fsr_id, sizeof(long), 1, in);
        double mesh_val;
        ret = _geometry->twiddleRead(&mesh_val, sizeof(double), 1, in);
        extruded_fsr->_fsr_ids[j] = fsr_id;
        extruded_fsr->_materials[j] = _geometry->findFSRMaterial(fsr_id);
        extruded_fsr->_mesh[j+1] = mesh_val;
      }

      /* Setup reverse lookup */
      extruded_FSR_lookup[extruded_fsr_id] = extruded_fsr;
    }

    /* Record 2D track info */
    int num_2D_tracks;
    ret = _geometry->twiddleRead(&num_2D_tracks, sizeof(int), 1, in);
    for (int i=0; i < num_2D_tracks; i++) {
      Track* track = _tracks_2D_array[i];
      read_segments.onTrack(track, NULL);
    }
    _contains_2D_segments = true;

    /* Record maximum number of segments */
    ret = _geometry->twiddleRead(&_max_num_segments, sizeof(int), 1, in);

    /* Allocate temporary segments */
    allocateTemporarySegments();
}
