#include "TrackGenerator3D.h"
#include "TrackTraversingAlgorithms.h"


/**
 * @brief Constructor for the TrackGenerator assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param spacing track spacing (cm)
 */
TrackGenerator3D::TrackGenerator3D(Geometry* geometry, int num_azim,
                                   int num_polar, double azim_spacing,
                                   double polar_spacing) :
                    TrackGenerator(geometry, num_azim, num_polar,
                                   azim_spacing) {
  setDesiredPolarSpacing(polar_spacing);
  _contains_3D_tracks = false;
  _contains_3D_segments = false;
  _contains_global_z_mesh = false;
  _contains_segmentation_heights = false;
  _contains_temporary_segments = false;
  _equal_z_spacing = false;
  _segment_formation = EXPLICIT_3D;
  _max_num_tracks_per_stack = 0;
  _num_seg_matrix_rows = 0;
  _num_seg_matrix_columns = 0;
  _track_generation_method = GLOBAL_TRACKING;
  _tracks_3D = NULL;

  _cum_tracks_per_stack = NULL;
  _cum_tracks_per_xy = NULL;
  _tracks_per_stack = NULL;
  _first_lz_of_stack = NULL;
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator3D::~TrackGenerator3D() {

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
  }
}


/**
 * @brief Return the track polar spacing (cm).
 * @details This will return the user-specified track spacing and NOT the
 *          effective track spacing which is computed and used to generate
 *          cyclic tracks.
 * @return the track polar spacing (cm)
 */
double TrackGenerator3D::getDesiredPolarSpacing() {
  return _polar_spacing;
}


/**
 * @brief Return the total number of 3D Tracks across the Geometry.
 * @return the total number of 3D Tracks
 */
int TrackGenerator3D::getNumTracks() {
  return getNum3DTracks();
}


/**
 * @brief Return the total number of 3D Tracks across the Geometry.
 * @return the total number of 3D Tracks
 */
int TrackGenerator3D::getNum3DTracks() {

  int a = _num_azim/2 - 1;
  int xy = _num_x[a] + _num_y[a] - 1;
  int p = _num_polar - 1;
  return _cum_tracks_per_stack[a][xy][p] + _tracks_per_stack[a][xy][p];
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator3D::getNumSegments() {
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
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
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
 *        requested azimuthal angle index and polar angle index
 * @param azim the requested azimuthal angle index
 * @param polar the requested polar angle index
 * @return the requested axial spacing
 */
double TrackGenerator3D::getZSpacing(int azim, int polar) {
  return _dz_eff[azim][polar];
}


/**
 * @brief Returns the maximum number of tracks in a single stack
 * @return the maximum number of tracks
 */
int TrackGenerator3D::getMaxNumTracksPerStack() {
  return _max_num_tracks_per_stack;
}


/**
 * @brief Returns the number of rows in the temporary segment storage matrix
 * @details For on-the-fly computation, a matrix of temporary segments is
 *          allocated for each thread. This matrix is indexed by the z-stack
 *          index (row) and the segment number (column). For ray tracing by
 *          individual Tracks the number of rows is always one since temporary
 *          segments only need to be stored for one Track at a time.
 * @return _num_seg_matrix_rows the number of rows in the temporary segment storage matrix
 */
int TrackGenerator3D::getNumRows() {
  _num_seg_matrix_rows = 1;
  return _num_seg_matrix_rows;
}


/**
 * @brief Returns the number of columns in the temporary segment storage matrix
 * @details For on-the-fly computation, a matrix of temporary segments is
 *          allocated for each thread. This matrix is indexed by the z-stack
 *          index (row) and the segment number (column). The number of columns
 *          is equal to the maximum number of segments per Track at the time of
 *          allocation.
 * @return _num_seg_matrix_columns the number of columns in the temporary segment storage
 *         matrix
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
 * @brief Returns a 3D array of the number of 3D Tracks in each z-stack
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
 *        azimuthal angle index and polar angle index
 * @param azim the azimuthal angle index
 * @param polar the polar angle index
 * @return the number of 3D Tracks in the z-direction of the Geometry
 */
int TrackGenerator3D::getNumZ(int azim, int polar) {
  return _num_z[azim][polar];
}


/**
 * @brief Returns the number of 3D Tracks in the radial direction for a given
 *        azimuthal angle index and polar angle index
 * @param azim the azimuthal angle index
 * @param polar the polar angle index
 * @return the number of 3D Tracks in the radial direction of the Geometry
 */
int TrackGenerator3D::getNumL(int azim, int polar) {
  return _num_l[azim][polar];
}


/**
 * @brief Set the suggested track polar spacing (cm).
 * @param spacing the suggested track polar spacing
 */
void TrackGenerator3D::setDesiredPolarSpacing(double spacing) {
  if (spacing < 0)
    log_printf(ERROR, "Unable to set a negative track polar spacing "
               "%f for the TrackGenerator.", spacing);

  _polar_spacing = spacing;
  resetStatus();
}


/**
 * @brief sets the type of segmentation used for segment formation
 * @param segmentation_type a segmentationType defining the type of
 *        segmentation to be used in segment formation. Options are:
 *          - EXPLICIT_3D: explicit 2D/3D segment formation
 *          - OTF_TRACKS: axial on-the-fly ray tracing by individaul tracks
 *          - OTF_STACKS: axial on-the-fly ray tracing by entire z-stacks
 */
void TrackGenerator3D::setSegmentFormation(segmentationType segmentation_type) {
  _segment_formation = segmentation_type;
}


/**
 * @brief Sets the z-planes over which 2D segmentation is performed for
 *        on-the-fly calculations
 * @param z_mesh the z-coordinates defining the height of the radial
 *        segmentation planes
 */
void TrackGenerator3D::setSegmentationHeights(std::vector<FP_PRECISION>
                                              z_mesh) {
  _contains_segmentation_heights = true;
  _segmentation_heights = z_mesh;
}


/**
 * @brief Sets a global z-mesh to use during axial on-the-fly ray tracing
 * @details In axial on-the-fly ray tracing, normally each extruded FSR
 *          contians a z-mesh. During on-the-fly segmentation when a new
 *          extruded FSR is entered, a binary search must be conducted to
 *          determine the axial cell. Alternatively, this function can be
 *          called which creates a global z-mesh from the geometry so that
 *          binary searches must only be conducted at the beginning of the
 *          track.
 */
void TrackGenerator3D::useGlobalZMesh() {
  _contains_global_z_mesh = true;
  _global_z_mesh = _geometry->getUniqueZHeights();
}


/**
 * @brief Provides the global z-mesh and size if available
 * @details For some cases, a global z-mesh is generated for the Geometry. If
 *          so, a pointer to the assocaited mesh (array) is updated as well as
 *          the number of FSRs in the mesh. If no global z-mesh has been
 *          generated, a null pointer is given to z_mesh and the number of FSRs
 *          is assigned to be zero.
 * @param z_mesh The global z-mesh to be updated
 * @param num_fsrs The number of FSRs in the z-mesh
 */
void TrackGenerator3D::retrieveGlobalZMesh(FP_PRECISION*& z_mesh,
                                           int& num_fsrs) {
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
 * @brief Fills an array with the x,y coordinates for each Track.
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
void TrackGenerator3D::retrieveTrackCoords(double* coords, int num_tracks) {
  retrieve3DTrackCoords(coords, num_tracks);
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
 *          num_tracks = track_generator.getNum3DTracks()
 *          coords = track_generator.retrieve3DTrackCoords(num_tracks*6)
 * @endcode
 *
 * @param coords an array of coords of length 6 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator3D::retrieve3DTrackCoords(double* coords, int num_tracks) {

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

          StackTrackIndexes sti;
          sti._azim = a;
          sti._xy = i;
          sti._polar = p;
          sti._z = z;

          Track3D track;
          getTrackOTF(&track, &sti);
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
 * @brief Fills an array with the x,y coordinates for each Track segment.
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
void TrackGenerator3D::retrieveSegmentCoords(double* coords, int num_segments) {
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
void TrackGenerator3D::retrieve3DSegmentCoords(double* coords, int num_segments) {

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

          StackTrackIndexes sti;
          sti._azim = a;
          sti._xy = i;
          sti._polar = p;
          sti._z = z;

          Track3D track;
          getTrackOTF(&track, &sti);

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

  /* Initialize the 2D tracks */
  TrackGenerator::initializeTracks();

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

  /* Allocate memory for arrays */
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
  double theta;
  double width_x  = _geometry->getWidthX();
  double width_y = _geometry->getWidthY();
  double width_z  = _geometry->getWidthZ();
  double avg_polar_spacing = 0.0;

  /* Determine angular quadrature and track spacing */
  for (int i = 0; i < _num_azim/4; i++) {

    /* Determine the polar angles and spacing for this azimuthal angle */
    for (int j=0; j < _num_polar/2; j++) {

      /* Compute the cosine weighted average angle */
      theta = _quadrature->getTheta(i, j);

      /* The number of intersections with xy (denoted "l") plane */
      if (_track_generation_method == GLOBAL_TRACKING) {
        if (_equal_z_spacing)
          _num_l[i][j] = (int) (fabs(_cycle_length[i] * tan(M_PI_2 - theta)
                                   / _polar_spacing)) + 1;
        else
          _num_l[i][j] = (int) (fabs(_cycle_length[i] * tan(M_PI_2 - theta)
                                   * sin(theta) / _polar_spacing)) + 1;
      }
      else if (_track_generation_method == MODULAR_RAY_TRACING) {
        if (_equal_z_spacing)
          _num_l[i][j] = 2 * (int) (fabs(_cycle_length[i] * 0.5 *
                tan(M_PI_2 - theta) / _polar_spacing) + 1);
        else
          _num_l[i][j] = 2 * (int) (fabs(_cycle_length[i] * 0.5 *
                tan(M_PI_2 - theta) * sin(theta) / _polar_spacing) + 1);
      }
      else if (_track_generation_method == SIMPLIFIED_MODULAR_RAY_TRACING) {
        double dx = width_x / _num_x[i];
        double dy = width_y / _num_y[i];
        double dl = sqrt(dx*dx + dy*dy);
        if (_equal_z_spacing)
          _num_l[i][j] =  (int) round(_cycle_length[i] / dl * ceil(fabs(dl *
                          tan(M_PI_2 - theta) / _polar_spacing)));
        else
          _num_l[i][j] =  (int) round(_cycle_length[i] / dl * ceil(fabs(dl *
                          tan(M_PI_2 - theta) * sin(theta) / _polar_spacing)));
      }

      /* Number of crossings along the z axis */
      _num_z[i][j] = (int) (fabs(width_z * _num_l[i][j] * tan(theta)
                                 / _cycle_length[i])) + 1;

      /* Effective track spacing */
      _dl_eff[i][j]          = _cycle_length[i] / _num_l[i][j];
      _dz_eff[i][j]         = width_z / _num_z[i][j];
      FP_PRECISION theta = atan(_dl_eff[i][j] / _dz_eff[i][j]);
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


  /* Compute the total number of cycles */
  int num_cycles = 0;
  for (int a = 0; a < _num_azim/4; a++) {
    for (int c = 0; c < _cycles_per_azim[a]; c++)
      num_cycles += _num_polar;
  }

  /* Create array of cycle track indices */
  CycleTrackIndexes* ctis = new CycleTrackIndexes[num_cycles];

  num_cycles = 0;
  for (int a = 0; a < _num_azim/4; a++) {
    for (int c = 0; c < _cycles_per_azim[a]; c++) {
      for (int p = 0; p < _num_polar; p++) {
        ctis[num_cycles]._azim  = a;
        ctis[num_cycles]._cycle = c;
        ctis[num_cycles]._polar = p;
        num_cycles++;
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
  getCycleTrackData(ctis, num_cycles, false);

  /* Allocate memory for 3D track stacks */
  create3DTracksArrays();

  /* Save explicit Track data if necessary */
  if (_segment_formation == EXPLICIT_3D)
    getCycleTrackData(ctis, num_cycles, true);

  /* Delete the array of cycle track indexes */
  delete [] ctis;

  /* Record the maximum number of tracks in a single stack */
  for (int a=0; a < _num_azim/2; a++)
    for (int i=0; i < _num_x[a] + _num_y[a]; i++)
      for (int p=0; p < _num_polar; p++)
        if (_tracks_per_stack[a][i][p] > _max_num_tracks_per_stack)
          _max_num_tracks_per_stack = _tracks_per_stack[a][i][p];

  _contains_3D_tracks = true;
}


//FIXME: description
void TrackGenerator3D::getCycleTrackData(CycleTrackIndexes* ctis,
                                         int num_cycles, bool save_tracks) {

  /* Determine if initializing or allocating tracks */
  std::string msg;
  if (save_tracks)
    msg = "Saving 3D Tracks";
  else
    msg = "Initializing 3D Tracks";
  Progress progress(num_cycles, msg);

  /* Loop over 3D track cycles */
#pragma omp parallel
  {
    Track3D track_3D;
    CycleTrackIndexes* cti;
    int nlz;
#pragma omp for
    for (int cycle=0; cycle < num_cycles; cycle++) {
      cti = &ctis[cycle];
      nlz = _num_l[cti->_azim][cti->_polar] + _num_z[cti->_azim][cti->_polar];
      std::cout << "Cycle " << cycle << " NLZ = " << nlz << std::endl;
      for (int lz=0; lz < nlz; lz++) {
        cti->_lz = lz;
        /*
        if (lz == 94) {
          std::cout << "SPOT 1" << std::endl;
          exit(1);
        }
        */
        cti->_train = -1;
        get3DTrack(cti, &track_3D, true, save_tracks);
      }
      progress.incrementCounter();
    }
  }
}


/* Create tracks starting on Z_MIN and L_MIN surfaces */
/*
 *             The track layout in the lz plane
 *       _____________________________________________
 *      | /    /    /    /    /    /    /    /    /   |
 *      |/    /    /    /    /    /    /    /    /    |
 * ^  9 |    /    /    /    /    /    /    /    /    /|
 * |    |   /    /    /    /    /    /    /    /    / |
 * z+   |__/____/____/____/____/____/____/____/____/__|
 *         8    7    6    5    4    3    2    1    0
 * l+ ->
 */

/* Create tracks starting on L_MIN and Z_MAX surfaces */
/*          1    2    3    4     5    6    7    8   9
 *       ______________________________________________
 *      |   \    \    \     \    \    \    \    \    \ |
 *      |    \    \    \     \    \    \    \    \    \|
 * ^  0 |\    \    \    \     \    \    \    \    \    |
 * |    | \    \    \    \     \    \    \    \    \   |
 * z+   |__\____\____\____\ ____\____\____\____\____\__|
 *
 * l+ ->
 */
//FIXME
double TrackGenerator3D::getLStart(CycleTrackIndexes* cti) {

  double l_start;
  int lz = cti->_lz;
  int nl = _num_l[cti->_azim][cti->_polar];
  int nz = _num_z[cti->_azim][cti->_polar];
  double dl = _dl_eff[cti->_azim][cti->_polar];

  if (cti->_polar < _num_polar / 2) {
    if (lz < nl)
      l_start = _cycle_length[cti->_azim] - (lz + 0.5) * dl;
    else
      l_start = 0.0;
  }
  else {
    if (lz < nz)
      l_start = 0.0;
    else
      l_start = dl * (lz - nz + 0.5);
  }

  /* If the track starts on an edge nudge it forward */
  if (l_start != 0.0) {
    Track* track_2D = _tracks_2D_cycle[cti->_azim][cti->_cycle][0];
    double phi      = track_2D->getPhi();
    double cos_phi  = cos(phi);
    double sin_phi  = sin(phi);
    double x_start  = track_2D->getStart()->getX();
    double y_start  = track_2D->getStart()->getY();
    double x_ext = x_start - _geometry->getMinX() + l_start * cos_phi;
    double y_ext = y_start - _geometry->getMinY() + l_start * sin_phi;
    int x_cell = int(round(x_ext / _geometry->getWidthX()));
    int y_cell = int(round(y_ext / _geometry->getWidthY()));

    /* If the start point is on a double reflection, nudge it to avoid track
     * with the same start and end points */
    if (fabs(x_cell * _geometry->getWidthX() - x_ext) < TINY_MOVE ||
        fabs(y_cell * _geometry->getWidthY() - y_ext) < TINY_MOVE )
      l_start += TINY_MOVE * 10;
  }

  return l_start;
}


//FIXME: description
int TrackGenerator3D::getFirstStack(CycleTrackIndexes* cti, Track3D* track_3D) {

  /* Initialize variables */
  Track* track_2D = _tracks_2D_cycle[cti->_azim][cti->_cycle][0];

  double phi      = track_2D->getPhi();
  double cos_phi  = cos(phi);
  double sin_phi  = sin(phi);

  double x_start  = track_2D->getStart()->getX();
  double y_start  = track_2D->getStart()->getY();
  double l_start  = getLStart(cti);

  int lz          = cti->_lz;
  int nl          = getNumL(cti->_azim, cti->_polar);
  int nz          = getNumZ(cti->_azim, cti->_polar);
  double dz       = _dz_eff[cti->_azim][cti->_polar];

  double x_min = _geometry->getMinX();
  double y_min = _geometry->getMinY();
  double z_min = _geometry->getMinZ();

  double x_max = _geometry->getMaxX();
  double y_max = _geometry->getMaxY();
  double z_max = _geometry->getMaxZ();

  double width_x = _geometry->getWidthX();
  double width_y = _geometry->getWidthY();

  double x_ext = x_start - x_min + l_start * cos_phi;
  double y_ext = y_start - y_min + l_start * sin_phi;

  /* Get the index of the first stack */
  int x_cell = int(floor(x_ext / width_x));
  int y_cell = int(floor(y_ext / width_y));

  /* Get the starting x-coord */
  double x1 = fmod(x_ext, width_x) + x_min;

  /* If the x_cell is odd, invert the x1 value in the geometry */
  if (x_cell % 2 == 1)
    x1 = x_max - (x1 - x_min);

  /* Get the starting y-coord */
  double y1 = fmod(y_ext, width_y) + y_min;

  /* If the y_cell is odd, invert the y1 value in the geometry */
  if (y_cell % 2 == 1)
    y1 = y_max - (y1 - y_min);

  /* Get the z-coords */
  double z1;
  double z2;
  if (cti->_polar < _num_polar / 2) {
    z1 = z_min + std::max(0., (lz - nl + 0.5)) * dz;
    z2 = z_max + std::min(0., (lz - nz + 0.5)) * dz;
  }
  else {
    z1 = z_max + std::min(0., (lz - nz + 0.5)) * dz;
    z2 = z_min + std::max(0., (lz - nl + 0.5)) * dz;
  }

  track_3D->getStart()->setCoords(x1, y1, z1);
  track_3D->getEnd()->setZ(z2);
  return x_cell + y_cell;
}


/**
 * @brief A 3D Track is decomposed by intersections with x and y boundaries
 * @details A 3D Track which initially neglects intersections with x and y
 *          boundaries is split into multiple tracks by finding those
 *          x and y intersections. If create_tracks is set to true, the Tracks
 *          will be altered to represent the correct 3D Tracks. If not, the
 *          number of tracks in the train will simply be counted. Whenever this
 *          function is called, the number of tracks in the associated z-stacks
 *          are incremented.
 * @param track The 3D track to be decomposed
 * @param l_start The distance accross the radial or "l" direction to where the
 *        track starts
 * @param l_end The distance accross the radial or "l" direction to where the
 *        track end
 * @param azim The azimuthal index of the track
 * @param cycle The cycle index of the track
 * @param polar The polar index of the track
 * @param lz_index The lz index into the Track3D cycles array
 * @param create_tracks Boolean value determining whether to create the
 *        associated decomposed tracks (true) or simply count them (false)
 */
void TrackGenerator3D::get3DTrack(CycleTrackIndexes* cti,
                                  Track3D* track,
                                  bool create_arrays,
                                  bool save_tracks) {

  /* Initialize variables */
  double theta = _quadrature->getTheta(cti->_azim, cti->_polar);
  bool end_of_train = false;

  /* Get the index of the first stack */
  int stack = getFirstStack(cti, track);
  int first_stack = stack;
  int train = 0;

  /* Record start and end points */
  double x1;
  double y1;
  double z1;

  double x2 = track->getStart()->getX();
  double y2 = track->getStart()->getY();
  double z2 = track->getStart()->getZ();

  bool fwd;

  /* Set the start and end point for each 3D track */
  while (!end_of_train) {

    /* Get the 2D track associated with this 3D track */
    Track* track_2D = _tracks_2D_cycle[cti->_azim][cti->_cycle][stack];
    double phi = track_2D->getPhi();
    fwd = track_2D->getDirectionInCycle();
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
    if (stack == first_stack && fwd) {
      dx = track_2D->getEnd()->getX() - x1;
      dy = track_2D->getEnd()->getY() - y1;
      dl_xy = sqrt(dx*dx + dy*dy);
    }
    else if (stack == first_stack && !fwd) {
      dx = track_2D->getStart()->getX() - x1;
      dy = track_2D->getStart()->getY() - y1;
      dl_xy = sqrt(dx*dx + dy*dy);
    }
    else
      dl_xy = track_2D->getLength();

    /* Get the distance to the nearest Z boundary */
    double dl_z;
    if (cti->_polar < _num_polar / 2)
      dl_z = (track->getEnd()->getZ() - z1) * tan(theta);
    else
      dl_z = (z1 - track->getEnd()->getZ()) / tan(theta - M_PI_2);

    double dl = std::min(dl_z, dl_xy);

    /* Set the end point coords */
    int polar;
    double track_theta;
    if (fwd) {
      x2 = x1 + dl * cos(phi);
      y2 = y1 + dl * sin(phi);
      polar = cti->_polar;
      track_theta = theta;
    }
    else {
      x2 = x1 - dl * cos(phi);
      y2 = y1 - dl * sin(phi);
      polar = _num_polar - cti->_polar - 1;
      track_theta = M_PI - theta;
    }

    if (cti->_polar < _num_polar / 2)
      z2 = z1 + dl / tan(theta);
    else
      z2 = z1 - dl * tan(theta - M_PI_2);

    /* If a track has the same start and end point, ignore it */
    if (fabs(x2 - x1) < TINY_MOVE ||
        fabs(y2 - y1) < TINY_MOVE ||
        fabs(z2 - z1) < TINY_MOVE)
      break;

    /* If the Z boundary or last stack or the desired
     * train has been reached, exit */
    if (dl_z < dl_xy || stack == _tracks_per_cycle[cti->_azim] - 1 ||
        cti->_train == stack - first_stack)
      end_of_train = true;

    if (create_arrays || save_tracks) {

      /* Get the index in the z-stack */
      int z = _tracks_per_stack[azim][xy][polar];
      _tracks_per_stack[azim][xy][polar]++;

      if (z == 0)
        _first_lz_of_stack[azim][xy][polar] = cti->_lz;

      if (save_tracks && _segment_formation == EXPLICIT_3D) {

        /* Get this 3D track */
        Track3D* track_3D = &_tracks_3D[azim][xy][polar][z];

        /* Set the start and end points */
        if (fwd) {
          track_3D->getStart()->setCoords(x1, y1, z1);
          track_3D->getEnd()->setCoords(x2, y2, z2);
        }
        else {
          track_3D->getEnd()->setCoords(x1, y1, z1);
          track_3D->getStart()->setCoords(x2, y2, z2);
        }

        /* Set track attributes */
        track_3D->setAzimIndex(azim);
        track_3D->setXYIndex(xy);
        track_3D->setPolarIndex(polar);
        track_3D->setTheta(track_theta);
        track_3D->setPhi(phi);
        track_3D->setDirectionInCycle(fwd);
      }
    }

    /* Increment the stack index */
    stack++;
  }

  /* Extract Geometry information */
  double x_min = _geometry->getMinX();
  double x_max = _geometry->getMaxX();
  double y_min = _geometry->getMinY();
  double y_max = _geometry->getMaxY();
  double z_min = _geometry->getMinZ();
  double z_max = _geometry->getMaxZ();

  /* Ensure the track does not lie outside the geometry bounds */
  x1 = std::max(x_min, std::min(x_max, x1));
  y1 = std::max(y_min, std::min(y_max, y1));
  z1 = std::max(z_min, std::min(z_max, z1));
  x2 = std::max(x_min, std::min(x_max, x2));
  y2 = std::max(y_min, std::min(y_max, y2));
  z2 = std::max(z_min, std::min(z_max, z2));

  /* Set the start and end points */
  if (fwd) {
    track->getStart()->setCoords(x1, y1, z1);
    track->getEnd()->setCoords(x2, y2, z2);
  }
  else {
    track->getEnd()->setCoords(x1, y1, z1);
    track->getStart()->setCoords(x2, y2, z2);
  }

  /* Set the train index, if it has not been set */
  if (cti->_train == -1)
    cti->_train = stack - first_stack;
}


/**
 * @brief Generates 2D segments for each extruded track across the Geometry,
 *        initializing axially extruded regions as well as 3D FSRs.
 */
void TrackGenerator3D::segmentizeExtruded() {

  log_printf(NORMAL, "Ray tracing for axially extruded track segmentation...");

  /* Get all unique z-coords at which 2D radial segementation is performed */
  std::vector<FP_PRECISION> z_coords;
  if (_contains_segmentation_heights)
    z_coords = _segmentation_heights;
  else
    z_coords = _geometry->getUniqueZPlanes();

  /* Loop over all extruded Tracks */
#pragma omp parallel for
  for (int index=0; index < _num_2D_tracks; index++)
    _geometry->segmentizeExtruded(_tracks_2D_array[index], z_coords);

  /* Initialize 3D FSRs and their associated vectors*/
  log_printf(NORMAL, "Initializing FSRs axially...");
  _geometry->initializeAxialFSRs(_global_z_mesh);
  _geometry->initializeFSRVectors();
  _contains_2D_segments = true;

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

  int tracks_segmented = 0;
  int num_3D_tracks = getNum3DTracks();

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++) {

    log_printf(NORMAL, "segmenting 3D tracks - Percent complete: %5.2f %%",
               double(tracks_segmented) / num_3D_tracks * 100.0);

#pragma omp parallel for
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          _geometry->segmentize3D(&_tracks_3D[a][i][p][z]);
      }
    }

    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++)
        tracks_segmented += _tracks_per_stack[a][i][p];
    }
  }

  _geometry->initializeFSRVectors();
  _contains_3D_segments = true;
}


/**
 * @brief Sets the track laydown method for generation of 3D Tracks
 * @details Options for the track laydown are GLOBAL_TRACKING,
 *          MODULAR_RAY_TRACING, and SIMPLIFIED_MODULAR_RAY_TRACING
 * @param method The track laydown method
 */
void TrackGenerator3D::setTrackGenerationMethod(int method) {

  if (method != GLOBAL_TRACKING && method != MODULAR_RAY_TRACING &&
      method != SIMPLIFIED_MODULAR_RAY_TRACING)
    log_printf(ERROR, "Unable to set Track Generation Method to %i. Valid"
               " methods include GLOBAL_TRACKING, MODULAR_RAY_TRACING, "
               "and SIMPLIFIED_MODULAR_RAY_TRACING", method);

  _track_generation_method = method;
}


/**
 * @brief Returns the track laydown method used for generating 3D Tracks
 * @return The track generation method: GLOBAL_TRACKING, MODULAR_RAY_TRACING,
 *         or SIMPLIFIED_MODULAR_RAY_TRACING
 */
int TrackGenerator3D::getTrackGenerationMethod() {
  return _track_generation_method;
}


/**
 * @brief Fills an array with the x,y coordinates for a given track.
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
                                                 int track_id) {

  /* Find 3D track associated with track_id */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          StackTrackIndexes sti;
          sti._azim = a;
          sti._xy = i;
          sti._polar = p;
          sti._z = z;

          Track3D track;
          getTrackOTF(&track, &sti);
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
                    "ID during coordinate retrieval");
  return;
}


/**
 * @brief allocates memory for 3D Tracks
 * @details Before calling this function, the number of tracks per z-stack
 *          should be known and initialized in the _tracks_per_stack 3D array
 */
void TrackGenerator3D::create3DTracksArrays() {

  long num_tracks = 0;
  long tracks_per_azim = 0;

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

  log_printf(NORMAL, "Total number of Tracks = %d", num_tracks);

  if (_segment_formation == EXPLICIT_3D) {
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

    /* Reset the number of tracks per stack....TODO: why? */
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < _num_x[a] + _num_y[a]; i++) {
        for (int p=0; p < _num_polar; p++)
          _tracks_per_stack[a][i][p] = 0;
      }
    }
  }
}


//FIXME
long TrackGenerator3D::get3DTrackID(StackTrackIndexes* sti) {
  return _cum_tracks_per_stack[sti->_azim][sti->_xy][sti->_polar] + sti->_z;
}


/**
 * @brief Allocates a new Quadrature with the default Quadrature
 * @details The defualt quadrature for 3D calculations is equal weight
 */
void TrackGenerator3D::initializeDefaultQuadrature() {
  if (_quadrature != NULL)
    delete _quadrature;
  _quadrature = new EqualWeightPolarQuad();
}


/**
 * @brief Returns the filename for writing tracking data
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

  if (_track_generation_method == GLOBAL_TRACKING)
    track_method = "GT";
  else if (_track_generation_method == MODULAR_RAY_TRACING)
    track_method = "MRT";
  else if (_track_generation_method == SIMPLIFIED_MODULAR_RAY_TRACING)
    track_method = "sMRT";

  if (_geometry->getCmfd() != NULL)
    test_filename << directory << "/3D_"
                  << _num_azim << "_azim_"
                  << _num_polar << "_polar_"
                  << _azim_spacing << "x" << _polar_spacing
                  << "_cm_spacing_cmfd_"
                  << _geometry->getCmfd()->getNumX()
                  << "x" << _geometry->getCmfd()->getNumY()
                  << "x" << _geometry->getCmfd()->getNumZ()
                  << "_quad_" << quad_type << "_track_" << track_method
                  << ".data";
  else
    test_filename << directory << "/3D_"
                  << _num_azim << "_azim_"
                  << _num_polar << "_polar_"
                  << _azim_spacing << "x" << _polar_spacing
                  << "_cm_spacing_quad_"
                  << quad_type << "_track_" << track_method << ".data";

  return test_filename.str();
}


/**
 * @brief Updates whether the TrackGenerator contains segments
 * @param contains_segments whether the TrackGenerator contains segments
 */
void TrackGenerator3D::setContainsSegments(bool contains_segments) {
  _contains_3D_segments = contains_segments;
}


/**
 * @brief Checks the boundary conditions for all 3D surfaces for inconsistent
 *        periodic boundary conditions
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

  /* Check for correct track method if a PERIODIC bc is present */
  if (_geometry->getMinZBoundaryType() == PERIODIC)
    _periodic = true;

  /* Check if 3D track generation method allows periodic boundaries */
  if (_periodic && _track_generation_method != MODULAR_RAY_TRACING &&
        _track_generation_method != SIMPLIFIED_MODULAR_RAY_TRACING)
      log_printf(ERROR, "Cannot create tracks for a geometry containing a"
                 " periodic BC with a track generation method that is not"
                 " modular");
}


/**
 * @brief Allocates memory for temporary segment storage if necessary
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

  /* Allocate new temporary segments */
  for (int t = 0; t < _num_threads; t++)
    _temporary_segments.at(t) = new segment[_num_seg_matrix_columns];
}


/**
 * @brief Resets the TrackGenerator to not contain tracks or segments
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


//FIXME
void TrackGenerator3D::deleteTemporarySegments() {
  /* Delete temporary segments if they exist */
  if (_contains_temporary_segments) {
    for (int t = 0; t < _num_threads; t++) {
      delete [] _temporary_segments.at(t);
    }
  }
  _contains_temporary_segments = false;
}


//FIXME
void TrackGenerator3D::getTrackOTF(Track3D* track, StackTrackIndexes* sti) {

  try {
    Track* track_2D = &_tracks_2D[sti->_azim][sti->_xy];
    double x1, x2, y1, y2, z1, z2;
    CycleTrackIndexes cti;
    convertSTItoCTI(sti, &cti);

    /* Set the start and end points */
    if (_segment_formation == EXPLICIT_3D) {
      Track3D* track_3D =
        &_tracks_3D[sti->_azim][sti->_xy][sti->_polar][sti->_z];
      Point* start_3d_1 = track_3D->getStart();
      Point* end_3d_1 = track_3D->getEnd();
      Point* start_3d_2 = track->getStart();
      Point* end_3d_2 = track->getEnd();

      x1 = start_3d_1->getX();
      y1 = start_3d_1->getY();
      z1 = start_3d_1->getZ();
      x2 = end_3d_1->getX();
      y2 = end_3d_1->getY();
      z2 = end_3d_1->getZ();

      start_3d_2->setCoords(x1, y1, z1);
      end_3d_2->setCoords(x2, y2, z2);
    }
    else
      get3DTrack(&cti, track, false, false);

    /* Set the angles */
    track->setDirectionInCycle(track_2D->getDirectionInCycle());
    track->setPhi(track_2D->getPhi());
    track->setTheta(_quadrature->getTheta(sti->_azim, sti->_polar));
    track->setAzimIndex(sti->_azim);
    track->setXYIndex(sti->_xy);
    track->setPolarIndex(sti->_polar);
    get3DTrackData(sti, &cti, true, track);
    get3DTrackData(sti, &cti, false, track);
    track->setUid(get3DTrackID(sti));
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to get track otf (%d, %d, %d, %d)."
               " Backtrace:\n%s", sti->_azim, sti->_xy, sti->_polar,
               sti->_z, e.what());
  }
}


//FIXME
Track3D**** TrackGenerator3D::get3DTracks() {
  return _tracks_3D;
}


//FIXME
void TrackGenerator3D::convertSTItoCTI(StackTrackIndexes* sti,
                                       CycleTrackIndexes* cti) {

  std::cout << "STI: " << sti->_azim << ", " << sti->_xy << ", " << sti->_polar
    << ", " << sti->_z << std::endl;
  Track* track_2D = &_tracks_2D[sti->_azim][sti->_xy];
  cti->_azim = sti->_azim;
  cti->_cycle = track_2D->getCycleIndex();

  if (track_2D->getDirectionInCycle())
    cti->_polar = sti->_polar;
  else
    cti->_polar = _num_polar - sti->_polar - 1;

  cti->_lz = _first_lz_of_stack[sti->_azim][sti->_xy][sti->_polar] + sti->_z;
/*
  if (cti->_lz == 94) {
      std::cout << "AZIM = " << sti->_azim << std::endl;
      std::cout << "XY = " << sti->_xy << std::endl;
      std::cout << "POLAR = " << sti->_polar << std::endl;
      std::cout << "Z = " << sti->_z << std::endl;
      std::cout << "FIRST LZ OF STACK = " <<
        _first_lz_of_stack[sti->_azim][sti->_xy][sti->_polar] << std::endl;
      exit(1);
  }
  */

  getTrainIndex(cti, sti);
}


//FIXME
void TrackGenerator3D::convertCTItoSTI(CycleTrackIndexes* cti,
                                       StackTrackIndexes* sti) {

  std::cout << "CTI: " << cti->_azim << ", " << cti->_cycle << ", " << cti->_polar
    << ", " << cti->_lz << ", " << cti->_train << std::endl;
  int stack = getStackIndex(cti);
  Track* track_2D = _tracks_2D_cycle[cti->_azim][cti->_cycle][stack];
  //FIXME HERE
/*
  std::cout << "c2s AZIM = " << cti->_azim << std::endl;
  std::cout << "c2s CYCLE = " << cti->_cycle << std::endl;
  std::cout << "c2s stack = " << stack << std::endl;
  std::cout << "Max stack should be " << _tracks_per_cycle[cti->_azim]
    << std::endl;
*/
  sti->_azim = track_2D->getAzimIndex();
//  std::cout << "c2s COMPLETE " << std::endl;
  sti->_xy = track_2D->getXYIndex();

  if (track_2D->getDirectionInCycle())
    sti->_polar = cti->_polar;
  else
    sti->_polar = _num_polar - cti->_polar - 1;

  sti->_z = cti->_lz - _first_lz_of_stack[sti->_azim][sti->_xy][sti->_polar];
}


//FIXME
int TrackGenerator3D::getStackIndex(CycleTrackIndexes* cti) {

  Track3D track;
  int first_stack = getFirstStack(cti, &track);
/*
  if (first_stack + cti->_train > 4) {
    std::cout << "First stack = " << first_stack << std::endl;
    std::cout << "Train = " << cti->_train << std::endl;
  }
*/
  return first_stack + cti->_train;
}


//FIXME
void TrackGenerator3D::getTrainIndex(CycleTrackIndexes* cti,
                                     StackTrackIndexes* sti) {

  Track3D track;
  int first_stack = getFirstStack(cti, &track);
  Track* track_2D = &_tracks_2D[sti->_azim][sti->_xy];
  cti->_train = track_2D->getStackIndex() - first_stack;
}


//FIXME
int TrackGenerator3D::getNumTracksPerLZ(CycleTrackIndexes* cti) {

  Track3D track_3D;
  int first_stack = getFirstStack(cti, &track_3D);
  int stack = first_stack;
  Track* track_2D = _tracks_2D_cycle[cti->_azim][cti->_cycle][stack];
  int azim = track_2D->getAzimIndex();
  int xy = track_2D->getXYIndex();
  int polar, min_lz, max_lz;
  int lz = cti->_lz;

  if (track_2D->getDirectionInCycle())
    polar = cti->_polar;
  else
    polar = _num_polar - cti->_polar - 1;

  while (stack < _tracks_per_cycle[cti->_azim]) {

    track_2D = _tracks_2D_cycle[cti->_azim][cti->_cycle][stack];
    azim = track_2D->getAzimIndex();
    xy = track_2D->getXYIndex();

    if (track_2D->getDirectionInCycle())
      polar = cti->_polar;
    else
      polar = _num_polar - cti->_polar - 1;

    min_lz = _first_lz_of_stack[azim][xy][polar];
    max_lz = _tracks_per_stack[azim][xy][polar] + min_lz - 1;

    if (cti->_polar < _num_polar / 2 && cti->_lz > max_lz)
      break;
    else if (cti->_polar >= _num_polar / 2 && cti->_lz < min_lz)
      break;

    stack++;
  }

  return stack - first_stack;
}


//FIXME description
void TrackGenerator3D::get3DTrackData(StackTrackIndexes* sti,
                                      CycleTrackIndexes* cti,
                                      bool outgoing,
                                      Track3D* track) {

  Track* track_2D = &_tracks_2D[sti->_azim][sti->_xy];
  bool cycle_fwd = _tracks_2D[sti->_azim][sti->_xy].getDirectionInCycle();
  CycleTrackIndexes cti_next;
  cti_next._azim  = cti->_azim;
  cti_next._cycle = cti->_cycle;
  int nz = _num_z[cti->_azim][cti->_polar];
  int nl = _num_l[cti->_azim][cti->_polar];
  int lz = cti->_lz;
  int pc = _num_polar - cti->_polar - 1;
  int tracks_per_lz = getNumTracksPerLZ(cti);
  bool check_double_reflection = false;
  bool next_fwd;
  boundaryType bc;

  if (outgoing) {
    next_fwd = _tracks_2D[sti->_azim][sti->_xy].getNextFwdFwd();
    bc = track_2D->getBCFwd();
  }
  else {
    next_fwd = _tracks_2D[sti->_azim][sti->_xy].getNextBwdFwd();
    bc = track_2D->getBCBwd();
  }

  /* Tracks pointing in the positive z direction in the lz plane */
  if (cti->_polar < _num_polar/2) {

    /* SURFACE_Z_MAX */
    if (cti->_train == tracks_per_lz - 1 && lz >= nz &&
        outgoing == cycle_fwd) {

      bc = _geometry->getMaxZBoundaryType();
      check_double_reflection = true;

      if (cycle_fwd)
        next_fwd = true;
      else
        next_fwd = false;

      /* PERIODIC BC */
      if (_geometry->getMaxZBoundaryType() == PERIODIC) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz - nz;
        cti_next._train = 0;
      }

      /* REFLECTIVE OR VACUUM BC */
      else {
        cti_next._polar = pc;
        cti_next._lz    = nl + 2 * nz - lz - 1;
        cti_next._train = 0;
      }
    }

    /* SURFACE_Z_MIN */
    else if (cti->_train == 0 && lz < nl && outgoing != cycle_fwd) {

      bc = _geometry->getMinZBoundaryType();
      check_double_reflection = true;

      if (cycle_fwd)
        next_fwd = false;
     else
        next_fwd = true;

      /* PERIODIC BC */
      if (_geometry->getMinZBoundaryType() == PERIODIC) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz + nz;
        cti_next._train = getNumTracksPerLZ(&cti_next) - 1;
      }

      /* REFLECTIVE OR VACUUM BC */
      else {
        cti_next._polar = pc;
        cti_next._lz    = nl - lz - 1;
        cti_next._train = getNumTracksPerLZ(&cti_next) - 1;
      }
    }

    /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURAFCE_Y_MAX */
    else {

      if (cti->_train == 0 && lz >= nl && outgoing != cycle_fwd) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz - nl;
        cti_next._train = getNumTracksPerLZ(&cti_next) - 1;
      }

      else if (cti->_train == tracks_per_lz - 1 && lz < nz &&
               outgoing == cycle_fwd) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = nl + lz;
        cti_next._train = 0;
      }

      else {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz;

        if (outgoing == cycle_fwd) {
          cti_next._train = cti->_train + 1; //TODO
          //std::cout << "Spot 4" << std::endl;
          //std::cout << "Computed " << cti_next._train << std::endl;
        }
        else
          cti_next._train = cti->_train - 1;
      }
    }
  }

  /* Tracks pointing in the negative z direction in the lz plane */
  else if (cti->_polar >= _num_polar/2) {

    /* SURFACE_Z_MAX */
    if (cti->_train == 0 && lz >= nz && outgoing != cycle_fwd) {

      bc = _geometry->getMaxZBoundaryType();
      check_double_reflection = true;

      if (cycle_fwd)
        next_fwd = false;
      else
        next_fwd = true;

      /* PERIODIC BC */
      if (_geometry->getMaxZBoundaryType() == PERIODIC) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz - nz;
        cti_next._train = getNumTracksPerLZ(&cti_next) - 1;
      }

      /* REFLECTIVE OR VACUUM BC */
      else {
        cti_next._polar = pc;
        cti_next._lz    = nl + 2 * nz - lz - 1;
        cti_next._train = getNumTracksPerLZ(&cti_next) - 1;
      }
    }

    /* SURFACE_Z_MIN */
    else if (cti->_train == tracks_per_lz - 1 && lz < nl &&
             outgoing == cycle_fwd) {

      bc = _geometry->getMinZBoundaryType();
      check_double_reflection = true;

      if (cycle_fwd)
        next_fwd = true;
      else
        next_fwd = false;

      /* PERIODIC BC */
      if (_geometry->getMinZBoundaryType() == PERIODIC) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz + nz;
        cti_next._train = 0;
      }

      /* REFLECTIVE OR VACUUM BC */
      else {
        cti_next._polar = pc;
        cti_next._lz    = nl - lz - 1;
        cti_next._train = 0;
      }
    }

    /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURAFCE_Y_MAX */
    else {

      if (cti->_train == tracks_per_lz - 1 && lz >= nl &&
          outgoing == cycle_fwd) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz - nl;
        cti_next._train = 0;
      }

      else if (cti->_train == 0 && lz < nz && outgoing != cycle_fwd) {
        cti_next._polar = cti->_polar;
        cti_next._lz    = lz + nl;
        cti_next._train = getNumTracksPerLZ(&cti_next) - 1;
      }

      else {
        cti_next._polar = cti->_polar;
        cti_next._lz = lz;

        if (outgoing == cycle_fwd)
          cti_next._train = cti->_train + 1;
        else
          cti_next._train = cti->_train - 1;
      }
    }
  }

  StackTrackIndexes sti_next;
  convertCTItoSTI(&cti_next, &sti_next);
  bool next_cycle_fwd =
    _tracks_2D[sti_next._azim][sti_next._xy].getDirectionInCycle();

  if (check_double_reflection && cycle_fwd != next_cycle_fwd)
    next_fwd = !next_fwd;

  if (outgoing) {
    track->setTrackNextFwd(get3DTrackID(&sti_next));
    track->setNextFwdFwd(next_fwd);
    track->setBCFwd(bc);
  }
  else {
    track->setTrackNextBwd(get3DTrackID(&sti_next));
    track->setNextBwdFwd(next_fwd);
    track->setBCBwd(bc);
  }
}


//FIXME
void TrackGenerator3D::getSTIByIndex(long id, StackTrackIndexes* sti) {

  long int cum_track_index = 0;
  int a;
  for (a=0; a < _num_azim/2; a++) {
    cum_track_index += _tracks_per_azim[a];
    if (id < cum_track_index)
      break;
  }
  sti->_azim = a;
  sti->_xy = -1;

  for (int xy=0; xy < _num_x[sti->_azim] + _num_y[sti->_azim]; xy++) {
    if (id < _cum_tracks_per_xy[sti->_azim][xy]) {
      sti->_xy = xy;
      break;
    }
  }

  if (sti->_xy == -1)
    log_printf(ERROR, "could not generate STI xy from track ID: %ld", id);

  for (int p=0; p < _num_polar; p++) {
    if (id < _cum_tracks_per_stack[sti->_azim][sti->_xy][p] +
        _tracks_per_stack[sti->_azim][sti->_xy][p]) {
      sti->_polar = p;
      sti->_z = id - _cum_tracks_per_stack[sti->_azim][sti->_xy][p];
      return;
    }
  }

  log_printf(ERROR, "could not generate STI from track ID: %ld", id);
}


//FIXME
void TrackGenerator3D::useEqualZSpacing() {
  _equal_z_spacing = true;
}
