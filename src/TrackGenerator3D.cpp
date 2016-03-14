#include "TrackGenerator3D.h"


/**
 * @brief Constructor for the TrackGenerator assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param spacing track spacing (cm)
 */
TrackGenerator3D::TrackGenerator3D(Geometry* geometry, const int num_azim,
                                   const int num_polar,
                                   const double azim_spacing,
                                   const double polar_spacing) :
                    TrackGenerator(geometry, num_azim, num_polar,
                                   azim_spacing) {
  setDesiredPolarSpacing(polar_spacing);
  _contains_3D_tracks = false;
  _contains_3D_segments = false;
  _contains_global_z_mesh = false;
  _contains_segmentation_heights = false;
  _contains_temporary_segments = false;
  _segment_formation = EXPLICIT_3D;
  _max_num_tracks_per_stack = 0;
  _num_seg_matrix_rows = 0;
  _num_seg_matrix_columns = 0;
  _track_generation_method = GLOBAL_TRACKING;
  _tracks_3D_array = NULL;
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator3D::~TrackGenerator3D() {

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_3D_tracks) {

    /* Delete 3D tracks */
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          delete [] _tracks_3D[a][i][p];
        }
        delete [] _tracks_3D[a][i];
        delete [] _tracks_per_stack[a][i];
      }
      delete [] _tracks_3D[a];
      delete [] _tracks_per_stack[a];
    }
    delete [] _tracks_3D;
    delete [] _tracks_per_stack;

    /* Delete book keeping for 3D tracks */
    for (int a = 0; a < _num_azim/4; a++) {
      for (int c = 0; c < _cycles_per_azim[a]; c++) {
        for (int p=0; p < _num_polar; p++)
          delete [] _tracks_per_train[a][c][p];
        delete [] _tracks_per_train[a][c];
      }
      delete [] _tracks_per_train[a];
    }
    delete [] _tracks_per_train;

    /* Delete book keeping for 3D tracks */
    if (_tracks_3D_cycle != NULL) {
      for (int a = 0; a < _num_azim/4; a++) {
        for (int c = 0; c < _cycles_per_azim[a]; c++) {
          for (int p=0; p < _num_polar; p++) {
            for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++)
              delete [] _tracks_3D_cycle[a][c][p][i];
            delete [] _tracks_3D_cycle[a][c][p];
          }
          delete [] _tracks_3D_cycle[a][c];
        }
        delete [] _tracks_3D_cycle[a];
      }
      delete [] _tracks_3D_cycle;
    }

    /* Delete book keeping for 3D tracks */
    for (int a=0; a < _num_azim/4; a++) {
      delete [] _num_l[a];
      delete [] _num_z[a];
      delete [] _dz_eff[a];
      delete [] _polar_spacings[a];
    }
    delete [] _num_l;
    delete [] _num_z;
    delete [] _dz_eff;
    delete [] _polar_spacings;

    /* Delete temporary segments if they exist */
    if (_contains_temporary_segments) {
      for (int t = 0; t < _num_threads; t++) {
        for (int z = 0; z < _num_seg_matrix_rows; z++)
          delete [] _temporary_segments.at(t)[z];
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
int TrackGenerator3D::getNum3DTracks() {

  int num_3D_tracks = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++)
        num_3D_tracks += _tracks_per_stack[a][i][p];
    }
  }

  return num_3D_tracks;
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator3D::getNum3DSegments() {

  if (!containsSegments())
    log_printf(ERROR, "Cannot get the number of 3D segments since they "
               "have not been generated.");

  int num_3D_segments = 0;

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          num_3D_segments += _tracks_3D[a][i][p][z].getNumSegments();
      }
    }
  }

  return num_3D_segments;
}


/**
 * @brief Returns an array of the Track pointers by increasing UID.
 * @details An array of pointers to all Track objects in the Geometry is
 *          returned, arranged by increasing unique identifier (UID).
 * @return the array of Track pointers
 */
Track** TrackGenerator3D::getTracksArray() {

  if (!containsTracks())
    log_printf(ERROR, "Unable to return the 1D array of Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_3D_array;
}


/**
 * @brief Returns a 4D jagged array of the 3D Tracks.
 * @details The first index into the array is the azimuthal angle, the second
 *          index is the 2D Track number, the third index is the polar angle,
 *          and the fourth index is the z-stack number.
 * @return the 4D jagged array of 3D Tracks
 */
Track3D**** TrackGenerator3D::get3DTracks() {

  if (!containsTracks())
    log_printf(ERROR, "Unable to return the 3D ragged array of the 3D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_3D;
}


/**
 * @brief Returns a 2D array of adjusted polar spacings
 * @details An array of polar spacings after adjustment is returned,
 *          indexed first by azimuthal angle and then by polar angle
 * @return the 2D array of polar spacings
 */
double** TrackGenerator3D::getPolarSpacings() {
  return _polar_spacings;
}


/**
 * @brief Returns the adjusted polar spacing at the requested azimuthal
 *        angle index and polar angle index
 * @details The polar spacing depends on the azimuthal angle and the polar
 *          angle. This function returns the azimuthal spacing used at the
 *          desired azimuthal angle and polar angle indexes.
 * @param azim the requested azimuthal angle index
 * @param polar the requested polar angle index
 * @return the requested polar spacing
 */
double TrackGenerator3D::getPolarSpacing(int azim, int polar) {
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
  return _polar_spacings[azim][polar];
}


/**
 * @brief Returns the spacing between tracks in the axial direction for the
 *        requested azimuthal angle index and polar angle index
 * @param azim the requested azimuthal angle index
 * @param polar the requested polar angle index
 * @return the requested axial spacing
 */
double TrackGenerator3D::getZSpacing(int azim, int polar) {
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
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
  if (_segment_formation == OTF_STACKS)
    _num_seg_matrix_rows = _max_num_tracks_per_stack;
  else
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
 * @param num_row The requested row number associated with the z-stack index
 * @return a pointer to the array of temporary segments
 */
segment* TrackGenerator3D::getTemporarySegments(int thread_id, int row_num) {
  if (_contains_temporary_segments)
    return _temporary_segments.at(thread_id)[row_num];
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
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
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
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
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
  _contains_2D_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _use_input_file = false;
  _tracks_filename = "";
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
void TrackGenerator3D::setSegmentationHeights(std::vector<double> z_mesh) {
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
void TrackGenerator3D::setGlobalZMesh() {
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
 * @brief Fills an array with the x,y,z coordinates and the periodic cycle ID
 *        for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords(num_tracks*7)
 * @endcode
 *
 * @param coords an array of coords of length 7 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator3D::retrieve3DPeriodicCycleCoords(double* coords,
                                                   int num_tracks) {

  if (num_tracks != 7*getNum3DTracks())
    log_printf(ERROR, "Unable to retrieve the 3D Track periodic cycle "
               "coordinates since the TrackGenerator contains %d Tracks with "
               "%d coordinates but an array of length %d was input",
               getNum3DTracks(), 7*getNum3DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          coords[counter]   = _tracks_3D[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D[a][i][p][z].getEnd()->getZ();
          coords[counter+6] = _tracks_3D[a][i][p][z].getPeriodicCycleId();
          counter += 7;
        }
      }
    }
  }

  return;
}


/**
 * @brief Fills an array with the x,y,z coordinates and the reflective cycle ID
 *        for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_tracks = track_generator.getNumTracks()
 *          coords = track_generator.retrieveTrackCoords(num_tracks*7)
 * @endcode
 *
 * @param coords an array of coords of length 7 times the number of Tracks
 * @param num_tracks the total number of Tracks
 */
void TrackGenerator3D::retrieve3DReflectiveCycleCoords(double* coords,
                                                     int num_tracks) {

  if (num_tracks != 7*getNum3DTracks())
    log_printf(ERROR, "Unable to retrieve the 3D Track reflective cycle "
               "coordinates since the TrackGenerator contains %d Tracks with "
               "%d coordinates but an array of length %d was input",
               getNum3DTracks(), 7*getNum3DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          coords[counter]   = _tracks_3D[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D[a][i][p][z].getEnd()->getZ();
          coords[counter+6] = _tracks_3D[a][i][p][z].getReflectiveCycleId();
          counter += 7;
        }
      }
    }
  }

  return;
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
void TrackGenerator3D::retrieve3DTrackCoords(double* coords, int num_tracks) {

  if (num_tracks != 6*getNum3DTracks())
    log_printf(ERROR, "Unable to retrieve the Track coordinates since the "
               "TrackGenerator contains %d Tracks with %d coordinates but an "
               "array of length %d was input",
               getNum3DTracks(), 6*getNum3DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          coords[counter]   = _tracks_3D[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D[a][i][p][z].getEnd()->getZ();
          counter += 6;
        }
      }
    }
  }

  return;
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
void TrackGenerator3D::retrieve3DSegmentCoords(double* coords, int num_segments) {

  if (num_segments != 7*getNum3DSegments())
    log_printf(ERROR, "Unable to retrieve the Track segment coordinates since "
               "the TrackGenerator contains %d segments with %d coordinates "
               "but an array of length %d was input",
               getNum3DSegments(), 7*getNum3DSegments(), num_segments);

  segment* curr_segment = NULL;
  double x0, x1, y0, y1, z0, z1;
  double phi, theta;
  segment* segments;

  int counter = 0;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          x0    = _tracks_3D[a][i][p][z].getStart()->getX();
          y0    = _tracks_3D[a][i][p][z].getStart()->getY();
          z0    = _tracks_3D[a][i][p][z].getStart()->getZ();
          phi   = _tracks_3D[a][i][p][z].getPhi();
          theta = _tracks_3D[a][i][p][z].getTheta();

          segments = _tracks_3D[a][i][p][z].getSegments();

          for (int s=0; s < _tracks_3D[a][i][p][z].getNumSegments(); s++) {

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

            counter += 7;
          }
        }
      }
    }
  }

  return;
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
  _dz_eff    = new double*[_num_azim/4];
  _num_z            = new int*[_num_azim/4];
  _num_l            = new int*[_num_azim/4];
  _polar_spacings   = new double*[_num_azim/4];
  _tracks_3D_cycle  = new Track3D*****[_num_azim/4];
  _tracks_per_train = new int***[_num_azim/4];
  _num_3D_tracks    = 0;

  for (int i=0; i < _num_azim/4; i++) {
    _dz_eff[i]         = new double[_num_polar/2];
    _num_z[i]          = new int[_num_polar/2];
    _num_l[i]          = new int[_num_polar/2];
    _polar_spacings[i] = new double[_num_polar/2];
  }

  /* Allocate memory for tracks per stack */
  _tracks_per_stack = new int**[_num_azim/2];
  for (int a=0; a < _num_azim/2; a++) {
    _tracks_per_stack[a] = new int*[getNumX(a) + getNumY(a)];
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      _tracks_per_stack[a][i] = new int[_num_polar];
      for (int p=0; p < _num_polar; p++)
        _tracks_per_stack[a][i][p] = 0;
    }
  }

  double x1, x2, y1, y2, z1, z2;
  double theta;
  double width  = _geometry->getWidthX();
  double height = _geometry->getWidthY();
  double depth  = _geometry->getWidthZ();
  double dl_eff[_num_azim/4][_num_polar/2];

  /* Determine angular quadrature and track spacing */
  for (int i = 0; i < _num_azim/4; i++) {

    /* Determine the polar angles and spacing for this azimuthal angle */
    for (int j=0; j < _num_polar/2; j++) {

      /* Compute the cosine weighted average angle */
      theta = _quadrature->getTheta(i, j);

      /* The number of intersections with xy (denoted "l") plane */
      if (_track_generation_method == GLOBAL_TRACKING) {
        _num_l[i][j] = (int) (fabs(_cycle_length[i] * tan(M_PI_2 - theta)
                                   * sin(theta) / _polar_spacing)) + 1;
      }
      else if (_track_generation_method == MODULAR_RAY_TRACING) {
        _num_l[i][j] = 2 * (int)
          (fabs(_cycle_length[i] * 0.5 * tan(M_PI_2 - theta)
                * sin(theta) / _polar_spacing) + 1);
      }
      else if (_track_generation_method == SIMPLIFIED_MODULAR_RAY_TRACING) {
        double dx = width / _num_x[i];
        double dy = height / _num_y[i];
        double dl = sqrt(dx*dx + dy*dy);
        _num_l[i][j] =  (int)
          round(_cycle_length[i] / dl *
                ceil(fabs(dl * tan(M_PI_2 - theta)
                          * sin(theta) / _polar_spacing)));
      }

      /* Number of crossings along the z axis */
      _num_z[i][j] = (int) (fabs(depth * _num_l[i][j] * tan(theta)
                                 / _cycle_length[i])) + 1;

      /* Effective track spacing */
      dl_eff[i][j]          = _cycle_length[i] / _num_l[i][j];
      _dz_eff[i][j]          = depth            / _num_z[i][j];
      _quadrature->setTheta(atan(dl_eff[i][j] / _dz_eff[i][j]), i, j);
      _polar_spacings[i][j] = _dz_eff[i][j] * sin(_quadrature->getTheta(i, j));
    }
  }

  /* Allocate memory for tracks per train */
  for (int a = 0; a < _num_azim/4; a++) {
    _tracks_per_train[a] = new int**[_cycles_per_azim[a]];
    for (int c = 0; c < _cycles_per_azim[a]; c++) {
      _tracks_per_train[a][c] = new int*[_num_polar];
      for (int p=0; p < _num_polar; p++)
        _tracks_per_train[a][c][p] = new int[getNumL(a,p) + getNumZ(a,p)];
    }
  }

  Track3D track_3D;
  Track* track_2D;
  int a, c, p, d, pc;
  double l_start, l_end;

  std::vector<std::tuple<int, int, int, int> > cycle_tuples;
  int tot_num_cycles = 0;
  for (int direction = 0; direction < 2; direction++) {
    for (a = 0; a < _num_azim/4; a++) {
      for (c = 0; c < _cycles_per_azim[a]; c++) {
        for (p = 0; p < _num_polar/2; p++) {
          cycle_tuples.push_back(std::make_tuple(direction, a, c, p));
          tot_num_cycles++;
        }
      }
    }
  }

  std::string msg = "initializing 3D Tracks";
  Progress progress(tot_num_cycles, msg);

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
  for (int create_tracks = 0; create_tracks < 2; create_tracks++) {

    progress.reset();

    /* Allocate memory for 3D track stacks */
    if (create_tracks)
      create3DTracksArrays();

    /* Loop over 3D track cycles */
    #pragma omp parallel for private(l_start, x1, y1, z1, track_2D, x2, y2, \
      z2, l_end, track_3D, pc, a, c, p, d)
    for (int ac = 0; ac < tot_num_cycles; ac++) {

      d = std::get<0>(cycle_tuples[ac]);
      a = std::get<1>(cycle_tuples[ac]);
      c = std::get<2>(cycle_tuples[ac]);
      p = std::get<3>(cycle_tuples[ac]);

      if (d == 0) {

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
        for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++) {

          /* Get the starting point */
          if (i < _num_l[a][p]) {
            l_start = _cycle_length[a] - (i + 0.5) * dl_eff[a][p];
            x1 = convertLtoX(l_start, a, c);
            y1 = convertLtoY(l_start, a, c);
            z1 = 0.0;
          }
          else{
            l_start = 0.0;
            track_2D = _tracks_2D_cycle[a][c][0];
            x1 = track_2D->getStart()->getX();
            y1 = track_2D->getStart()->getY();
            z1 = _dz_eff[a][p] * (i - _num_l[a][p] + 0.5);
          }

          /* Get the end point */
          if (i < _num_z[a][p]) {
            l_end = _cycle_length[a];
            track_2D = _tracks_2D_cycle[a][c][0];
            x2 = track_2D->getStart()->getX();
            y2 = track_2D->getStart()->getY();
            z2 = _dz_eff[a][p] * (i + 0.5);
          }
          else{
            l_end = _cycle_length[a] - dl_eff[a][p] *
                (i - _num_z[a][p] + 0.5);
            x2 = convertLtoX(l_end, a, c);
            y2 = convertLtoY(l_end, a, c);
            z2 = depth;
          }

          /* Set start and end points and save polar angle */
          track_3D.getStart()->setCoords(x1, y1, z1);
          track_3D.getEnd()->setCoords(x2, y2, z2);
          track_3D.setTheta(_quadrature->getTheta(a,p));

          /* Decompose the track in the LZ plane by splitting it
           * based on the x and y geometry boundaries */
          decomposeLZTrack(&track_3D, l_start, l_end, a, c, p, i,
                           create_tracks);
        }
      }
      else {

        pc = _num_polar-p-1;

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
        for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++) {

          /* Get the starting point */
          if (i < _num_z[a][p]) {
            l_start = 0.0;
            track_2D = _tracks_2D_cycle[a][c][0];
            x1 = track_2D->getStart()->getX();
            y1 = track_2D->getStart()->getY();
            z1 = _dz_eff[a][p] * (i + 0.5);
          }
          else{
            l_start = dl_eff[a][p] * (i - _num_z[a][p] + 0.5);
            x1 = convertLtoX(l_start, a, c);
            y1 = convertLtoY(l_start, a, c);
            z1 = depth;
          }

          /* Get the end point */
          if (i < _num_l[a][p]) {
            l_end = dl_eff[a][p] * (i + 0.5);
            x2 = convertLtoX(l_end, a, c);
            y2 = convertLtoY(l_end, a, c);
            z2 = 0.0;
          }
          else{
            l_end = _cycle_length[a];
            track_2D = _tracks_2D_cycle[a][c][0];
            x2 = track_2D->getStart()->getX();
            y2 = track_2D->getStart()->getY();
            z2 = _dz_eff[a][p] * (i - _num_l[a][p] + 0.5);
          }

          /* Set start and end points and save polar angle */
          track_3D.getStart()->setCoords(x1, y1, z1);
          track_3D.getEnd()->setCoords(x2, y2, z2);
          track_3D.setTheta(_quadrature->getTheta(a,pc));

          /* Decompose the track in the LZ plane by splitting it
           * based on the x and y geometry boundaries */
          decomposeLZTrack(&track_3D, l_start, l_end, a, c, pc, i,
                           create_tracks);
        }
      }

      progress.incrementCounter();

    }
  }

  /* Record the maximum number of tracks in a single stack */
  for (int a=0; a < _num_azim/2; a++)
    for (int i=0; i < getNumX(a) + getNumY(a); i++)
      for (int p=0; p < _num_polar; p++)
        if (_tracks_per_stack[a][i][p] > _max_num_tracks_per_stack)
          _max_num_tracks_per_stack = _tracks_per_stack[a][i][p];

  _contains_3D_tracks = true;

  /* Initialize the 3D track reflections and cycle ids */
  initializeTrackReflections();
  initializeTrackCycleIds();
  initializeTrackPeriodicIndices();
}


/**
 * @brief Initializes 3D Track reflections
 * @details This method computes the connecting Tracks for all 3D Tracks in
 *          the TrackGenerator analytically, handling both reflective and
 *          periodic boundaries.
 */
void TrackGenerator3D::initializeTrackReflections() {

  int pc;
  int ai, xi, pi, zi, xp, zp, pp, ci;

  /* Set reflective tracks and periodic top and bottom indices */
  for (int a = 0; a < _num_azim/4; a++) {
    for (int c = 0; c < _cycles_per_azim[a]; c++) {

      /* Loop over polar angles < PI/2 */
      for (int p=0; p < _num_polar/2; p++) {

        /* Set the complementary polar angle */
        pc = _num_polar-p-1;
        Track3D *** polar_group = _tracks_3D_cycle[a][c][p];

        for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++) {
          for (int t=0; t < _tracks_per_train[a][c][p][i]; t++) {

            if (t == _tracks_per_train[a][c][p][i]-1) {

              /* SURFACE_Y_MIN */
              if (i < _num_z[a][p]) {
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCFwd(_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflFwdFwd
                    (polar_group[_num_l[a][p] + i][0]->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (polar_group[_num_l[a][p] + i][0]);

                  /* PERIODIC */
                  if (_periodic) {
                    ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    zi = polar_group[i][t]->getZIndex();
                    ci = polar_group[i][t]->getCycleTrackIndex();
                    xp = _tracks_2D_cycle[a][c][ci]
                        ->getTrackPrdcFwd()->getXYIndex();
                    zp = polar_group[_num_l[a][p] + i][0]->getZIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                        (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
                else {

                  polar_group[i][t]->setBCBwd(_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflBwdFwd
                    (polar_group[_num_l[a][p] + i][0]->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (polar_group[_num_l[a][p] + i][0]);

                  /* PERIODIC */
                  if (_periodic) {
                    ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    zi = polar_group[i][t]->getZIndex();
                    ci = polar_group[i][t]->getCycleTrackIndex();
                    xp = _tracks_2D_cycle[a][c][ci]
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[_num_l[a][p] + i][0]->getZIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
              }

              /* SURFACE_Z_MAX */
              else{
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCFwd
                    (_geometry->getMaxZBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]);

                  /* PERIODIC */
                  if (_periodic) {
                    polar_group[i][t]->setTrackPrdcFwd
                      (polar_group[i - _num_z[a][p]][0]);
                  }
                }
                else{

                  polar_group[i][t]->setBCBwd
                    (_geometry->getMaxZBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflBwdFwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]);

                  /* PERIODIC */
                  if (_periodic) {
                    polar_group[i][t]->setTrackPrdcBwd
                      (polar_group[i - _num_z[a][p]][0]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURAFCE_Y_MAX */
            else{
              if (polar_group[i][t]->getCycleFwd()) {

                polar_group[i][t]->setBCFwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCFwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflFwdFwd
                  (polar_group[i][t+1]->getCycleFwd());
                polar_group[i][t]->setTrackReflFwd
                  (polar_group[i][t+1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]
                    ->getTrackPrdcFwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCBwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCBwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflBwdFwd
                  (polar_group[i][t+1]->getCycleFwd());
                polar_group[i][t]->setTrackReflBwd
                  (polar_group[i][t+1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]
                    ->getTrackPrdcBwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
            }

            if (t == 0) {

              /* SURFACE_Z_MIN */
              if (i < _num_l[a][p]) {
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCBwd
                    (_geometry->getMinZBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflBwdFwd
                    (!_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]
                     ->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]);

                  /* PERIODIC */
                  if (_periodic) {
                    polar_group[i][t]->setTrackPrdcBwd
                      (polar_group[i + _num_z[a][p]]
                       [_tracks_per_train[a][c][p][i + _num_z[a][p]]-1]);
                  }
                }
                else {

                  polar_group[i][t]->setBCFwd
                    (_geometry->getMinZBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflFwdFwd
                    (!_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]
                     ->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]);

                  /* PERIODIC */
                  if (_periodic) {
                    polar_group[i][t]->setTrackPrdcFwd
                      (polar_group[i + _num_z[a][p]]
                       [_tracks_per_train[a][c][p][i + _num_z[a][p]]-1]);
                  }
                }
              }

              /* SURFACE_Y_MIN */
              else{
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCBwd
                    (_geometry->getMinYBoundaryType());

                  /* REFFLECTIVE */
                  polar_group[i][t]->setReflBwdFwd
                    (!polar_group[i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (polar_group[i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]);

                  /* PERIODIC */
                  if (_periodic) {
                    ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    zi = polar_group[i][t]->getZIndex();
                    ci = polar_group[i][t]->getCycleTrackIndex();
                    xp = _tracks_2D_cycle[a][c][ci]
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][p]]
                      [_tracks_per_train[a][c][p]
                       [i - _num_l[a][p]] - 1]->getZIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
                else {

                  polar_group[i][t]->setBCFwd
                    (_geometry->getMinYBoundaryType());

                  /* REFFLECTIVE */
                  polar_group[i][t]->setReflFwdFwd
                    (polar_group[i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (polar_group[i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]);

                  /* PERIODIC */
                  if (_periodic) {
                    ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    zi = polar_group[i][t]->getZIndex();
                    xp = _tracks_2D_cycle[a][c][xi]
                      ->getTrackPrdcFwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][p]]
                      [_tracks_per_train[a][c][p]
                       [i - _num_l[a][p]] - 1]->getZIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURAFCE_Y_MAX */
            else{
              if (polar_group[i][t]->getCycleFwd()) {

                polar_group[i][t]->setBCBwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCBwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflBwdFwd
                  (!polar_group[i][t-1]->getCycleFwd());
                polar_group[i][t]->setTrackReflBwd
                  (polar_group[i][t-1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]
                    ->getTrackPrdcBwd()->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCFwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCFwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflFwdFwd
                  (!polar_group[i][t-1]->getCycleFwd());
                polar_group[i][t]->setTrackReflFwd
                  (polar_group[i][t-1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]
                    ->getTrackPrdcFwd()->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
            }
          }
        }
      }

      /* Loop over polar angles > PI/2 */
      for (int p = _num_polar/2; p < _num_polar; p++) {

        pc = _num_polar-p-1;
        Track3D *** polar_group = _tracks_3D_cycle[a][c][p];

        for (int i=0; i < _num_l[a][pc] + _num_z[a][pc]; i++) {
          for (int t=0; t < _tracks_per_train[a][c][p][i]; t++) {

            if (t == _tracks_per_train[a][c][p][i]-1) {

              /* SURFACE_Z_MIN */
              if (i < _num_l[a][pc]) {
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCFwd
                    (_geometry->getMinZBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]
                     ->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]);

                  /* PERIODIC */
                 if (_periodic) {
                   polar_group[i][t]->setTrackPrdcFwd
                     (polar_group[i + _num_z[a][pc]][0]);
                 }
                }
                else {

                  polar_group[i][t]->setBCBwd
                    (_geometry->getMinZBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflBwdFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]
                     ->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]);

                  /* PERIODIC */
                  if (_periodic) {
                    polar_group[i][t]->setTrackPrdcBwd
                      (polar_group[i + _num_z[a][pc]][0]);
                  }
                }
              }

              /* SURFACE_Y_MIN */
              else{
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCFwd
                    (_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflFwdFwd
                    (polar_group[i - _num_l[a][pc]][0]->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (polar_group[i - _num_l[a][pc]][0]);

                  /* PERIODIC */
                  if (_periodic) {
                    ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    zi = polar_group[i][t]->getZIndex();
                    ci = polar_group[i][t]->getCycleTrackIndex();
                    xp = _tracks_2D_cycle[a][c][ci]
                      ->getTrackPrdcFwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][pc]][0]->getZIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
                else {

                  polar_group[i][t]->setBCBwd
                    (_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflBwdFwd
                    (polar_group[i - _num_l[a][pc]][0]->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (polar_group[i - _num_l[a][pc]][0]);

                  /* PERIODIC */
                  if (_periodic) {
                    ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    zi = polar_group[i][t]->getZIndex();
                    ci = polar_group[i][t]->getCycleTrackIndex();
                    xp = _tracks_2D_cycle[a][c][ci]
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][pc]][0]->getZIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, or SURFACE_Y_MAX */
            else{
              if (polar_group[i][t]->getCycleFwd()) {

                polar_group[i][t]->setBCFwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCFwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflFwdFwd
                  (polar_group[i][t+1]->getCycleFwd());
                polar_group[i][t]->setTrackReflFwd
                  (polar_group[i][t+1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]
                    ->getTrackPrdcFwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCBwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCBwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflBwdFwd
                  (polar_group[i][t+1]->getCycleFwd());
                polar_group[i][t]->setTrackReflBwd
                  (polar_group[i][t+1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]
                    ->getTrackPrdcBwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
            }

            if (t == 0) {

              /* SURFACE_Y_MIN */
              if (i < _num_z[a][pc]) {
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCBwd
                    (_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflBwdFwd
                    (!polar_group[_num_l[a][pc] + i]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] + i] - 1]->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (polar_group[_num_l[a][pc] + i]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] + i] - 1]);

                  /* PERIODIC */
                  if (_periodic) {
                    ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    zi = polar_group[i][t]->getZIndex();
                    ci = polar_group[i][t]->getCycleTrackIndex();
                    xp = _tracks_2D_cycle[a][c][ci]
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[_num_l[a][pc] + i]
                      [_tracks_per_train[a][c][p]
                       [_num_l[a][pc] + i] - 1]->getZIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
                else {

                  polar_group[i][t]->setBCFwd
                    (_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  polar_group[i][t]->setReflFwdFwd
                    (polar_group[_num_l[a][pc] + _num_z[a][pc] - i - 1]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]
                      ->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (polar_group[_num_l[a][pc] + _num_z[a][pc] - i - 1]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]);

                  /* PERIODIC */
                  if (_periodic) {ai = polar_group[i][t]->getAzimIndex();
                    xi = polar_group[i][t]->getXYIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    zi = polar_group[i][t]->getZIndex();
                    ci = polar_group[i][t]->getCycleTrackIndex();
                    xp = _tracks_2D_cycle[a][c][ci]->getTrackPrdcFwd()
                      ->getXYIndex();
                    zp = polar_group[_num_l[a][pc] + _num_z[a][pc] - i - 1]
                      [_tracks_per_train[a][c][p]
                       [_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]->getZIndex();
                    _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D[ai][xp][pi][zp]);
                  }
                }
              }

              /* SURFACE_Z_MAX */
              else{
                if (polar_group[i][t]->getCycleFwd()) {

                  polar_group[i][t]->setBCBwd
                    (_geometry->getMaxZBoundaryType());

                  /* REFLECTIVE */
                  int idx = 2 * _num_z[a][pc] + _num_l[a][pc] - i - 1;
                  polar_group[i][t]->setReflBwdFwd
                    (!_tracks_3D_cycle[a][c][pc][idx]
                     [_tracks_per_train[a][c][pc][idx] - 1]->getCycleFwd());
                  polar_group[i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc][idx]
                     [_tracks_per_train[a][c][pc][idx] - 1]);

                  /* PERIODIC */
                  if (_periodic) {
                    polar_group[i][t]->setTrackPrdcBwd
                      (polar_group[i - _num_z[a][pc]]
                       [_tracks_per_train[a][c][p][i - _num_z[a][pc]]-1]);
                  }
                }
                else {

                  polar_group[i][t]->setBCFwd
                    (_geometry->getMaxZBoundaryType());

                  /* REFLECTIVE */
                  int idx = 2 * _num_z[a][pc] + _num_l[a][pc] - i - 1;
                  polar_group[i][t]->setReflFwdFwd
                    (!_tracks_3D_cycle[a][c][pc][idx]
                     [_tracks_per_train[a][c][pc][idx] - 1]->getCycleFwd());
                  polar_group[i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc][idx]
                     [_tracks_per_train[a][c][pc][idx] - 1]);

                  /* PERIODIC */
                  if (_periodic) {
                    polar_group[i][t]->setTrackPrdcFwd
                      (polar_group[i - _num_z[a][pc]]
                       [_tracks_per_train[a][c][p][i - _num_z[a][pc]]-1]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, or SURFACE_Y_MAX */
            else{
              if (polar_group[i][t]->getCycleFwd()) {

                polar_group[i][t]->setBCBwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCBwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflBwdFwd
                  (!polar_group[i][t-1]->getCycleFwd());
                polar_group[i][t]->setTrackReflBwd
                  (polar_group[i][t-1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]->getTrackPrdcBwd()
                    ->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCFwd
                    (_tracks_2D_cycle[a][c]
                     [polar_group[i][t]->getCycleTrackIndex()]->getBCFwd());

                /* REFLECTIVE */
                polar_group[i][t]->setReflFwdFwd
                  (!polar_group[i][t-1]->getCycleFwd());
                polar_group[i][t]->setTrackReflFwd
                  (polar_group[i][t-1]);

                /* PERIODIC */
                if (_periodic) {
                  ai = polar_group[i][t]->getAzimIndex();
                  xi = polar_group[i][t]->getXYIndex();
                  pi = polar_group[i][t]->getPolarIndex();
                  zi = polar_group[i][t]->getZIndex();
                  ci = polar_group[i][t]->getCycleTrackIndex();
                  xp = _tracks_2D_cycle[a][c][ci]->getTrackPrdcFwd()
                    ->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D[ai][xp][pi][zp]);
                }
              }
            }
          }
        }
      }
    }
  }
}


/**
 * @brief Cycle IDs are created for all 3D cycles and assigned to 3D Tracks
 * @details All tracks are traversed through connecting tracks, assigning cycle
 *          numbers until all 3D Tracks are traversed. This is done for both
 *          periodic and reflective connections.
 */
void TrackGenerator3D::initializeTrackCycleIds() {

  int id = 0;
  Track* track;
  bool fwd;

  /* Set the periodic track cycle ids */
  if (_periodic) {
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

            track = &_tracks_3D[a][i][p][z];

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
      }
    }
  }

  id = 0;

  /* Set the reflective track cycle ids */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          track = &_tracks_3D[a][i][p][z];
          fwd = true;

          if (track->getReflectiveCycleId() == -1) {
            while (track->getReflectiveCycleId() == -1) {

              /* Set the periodic cycle id */
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
  }
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
void TrackGenerator3D::decomposeLZTrack(Track3D* track, double l_start,
                                      double l_end, int azim, int cycle,
                                      int polar, int lz_index,
                                      bool create_tracks) {

  if (l_start > _cycle_length[azim] || l_start < 0.0)
    log_printf(ERROR, "Unable to tally tracks per stack for track with "
               "starting length not within "
               "(0, %f), l_start: %f", _cycle_length[azim], l_start);
  else if  (l_start > _cycle_length[azim] || l_start < 0.0)
    log_printf(ERROR, "Unable to tally tracks per stack for track with "
               "ending length not within "
               "(0, %f), l_end: %f", _cycle_length[azim], l_end);

  /* Initialize variables */
  double length_sum = 0.0;
  int first_stack = 0;
  int last_stack = _tracks_per_cycle[azim] - 1;
  Track* track_2d;
  double theta = track->getTheta();
  double phi = track->getPhi();
  double dx = cos(phi) * sin(theta) * TINY_MOVE;
  double dy = sin(phi) * sin(theta) * TINY_MOVE;
  double nudge = sqrt(dx*dx + dy*dy);
  double track_theta;

  /* Find the last cycle index */
  for (int i=0; i < _tracks_per_cycle[azim]; i++) {

    track_2d = _tracks_2D_cycle[azim][cycle][i];

    if (l_end < length_sum + track_2d->getLength() + 1.05 * nudge) {
      last_stack = i;
      break;
    }
    else{
      length_sum += track_2d->getLength();
    }
  }

  length_sum = 0.0;

  /* Find the first cycle index */
  for (int i=0; i < _tracks_per_cycle[azim]; i++) {

    track_2d = _tracks_2D_cycle[azim][cycle][i];

    if (l_start < length_sum + track_2d->getLength() - 1.05 * nudge) {
      first_stack = i;
      break;
    }
    else{
      length_sum += track_2d->getLength();
    }
  }

  /* Check to make sure at least one track is created */
  if (last_stack < first_stack)
    log_printf(ERROR, "Unable to tally tracks per stack for track with "
               "last stack less than first stack. "
               "first stack: %i, last stack: %i", first_stack, last_stack);

  int ti, zi, pi, ai;
  double l_start_2 = l_start;
  double length_sum_2 = length_sum;
  bool fwd;
  double x1, y1, z1;
  double x2, y2, z2;

  /* If tracks are to be created during this call to decomposeLZTrack,
   * decompose the track, compute endpoints, and set necessary attributes of
   * the 3D tracks
   */
  if (create_tracks) {

    Track3D* track_3d;
    int t=0;

    /* Set the start and end point for each 3D track */
    for (int i=first_stack; i <= last_stack; i++) {

      /* Get the 2D track associated with this 3D track */
      track_2d = _tracks_2D_cycle[azim][cycle][i];
      fwd = getCycleDirection(azim, cycle, i);
      ai = track_2d->getAzimIndex();
      ti = track_2d->getXYIndex();

      /* Set the start and end points */
      if (i == first_stack) {
        x1 = track->getStart()->getX();
        y1 = track->getStart()->getY();
        z1 = track->getStart()->getZ();
      }
      else{
        x1 = x2;
        y1 = y2;
        z1 = z2;
      }

      if (i == last_stack) {
        x2 = track->getEnd()->getX();
        y2 = track->getEnd()->getY();
        z2 = track->getEnd()->getZ();
      }
      else{
        length_sum += track_2d->getLength();
        if (fwd) {
          x2 = track_2d->getEnd()->getX();
          y2 = track_2d->getEnd()->getY();
        }
        else {
          x2 = track_2d->getStart()->getX();
          y2 = track_2d->getStart()->getY();
        }
        z2 = z1 + (length_sum - l_start) * tan(M_PI_2 - theta);
        l_start = length_sum;
      }

      /* Get the polar angle index. Note that the polar index in the
       * _tracks_3D_cycle array is not always the same as the polar
       * index in the _tracks_3D array since we want the tracks
       * in the _tracks_3D array to all be pointing in the positive-y
       * direction whereas tracks in the _tracks_3D array will point in
       * all directions.
       */
      if (fwd) {
        pi = polar;
        track_theta = theta;
      }
      else {
        pi = _num_polar - polar - 1;
        track_theta = M_PI - theta;
      }

      /* Get the index in the z-stack */
      zi = _tracks_per_stack[ai][ti][pi];

      /* Get this 3D track */
      track_3d = &_tracks_3D[ai][ti][pi][zi];

      /* Set pointer to track in 3D track cycles array */
      _tracks_3D_cycle[azim][cycle][polar][lz_index][t] = track_3d;

      /* Set the start and end points */
      if (fwd) {
        track_3d->getStart()->setCoords(x1, y1, z1);
        track_3d->getEnd()->setCoords(x2, y2, z2);
      }
      else {
        track_3d->getEnd()->setCoords(x1, y1, z1);
        track_3d->getStart()->setCoords(x2, y2, z2);
      }

      /* Set track attributes */
      track_3d->setAzimIndex(ai);
      track_3d->setXYIndex(ti);
      track_3d->setZIndex(zi);
      track_3d->setPolarIndex(pi);
      track_3d->setLZIndex(lz_index);
      track_3d->setCycleIndex(cycle);
      track_3d->setTrainIndex(t);
      track_3d->setTheta(track_theta);
      track_3d->setPhi(track_2d->getPhi());
      track_3d->setCycleFwd(fwd);
      track_3d->setCycleTrackIndex(i);

      /* Increment track train index */
      t++;
    }
  }

  /* If tracks are not being created, just increment the _tracks_per_train
   * array
   */
  else {

    /* Set the number of tracks in the lz track */
    _tracks_per_train[azim][cycle][polar][lz_index] =
      last_stack - first_stack + 1;
  }

  /* Increment the tracks per stack */
  for (int i=first_stack; i <= last_stack; i++) {
    track_2d = _tracks_2D_cycle[azim][cycle][i];
    fwd = getCycleDirection(azim, cycle, i);
    ti = track_2d->getXYIndex();
    ai = track_2d->getAzimIndex();

    /* Computing the endpoint z value */
    if (i == first_stack)
      z1 = track->getStart()->getZ();
    else
      z1 = z2;

    if (i == last_stack) {
      z2 = track->getEnd()->getZ();
    }
    else{
      length_sum_2 += track_2d->getLength();
      z2 = z1 + (length_sum_2 - l_start_2) * tan(M_PI_2 - theta);
      l_start_2 = length_sum_2;
    }

    /* Get the polar angle index */
    if (fwd)
      pi = polar;
    else
      pi = _num_polar - polar - 1;

    /* Tally another track in the _tracks_per_stack array */
    _tracks_per_stack[ai][ti][pi]++;
  }
}


/**
 * @brief Recalibrates Track start and end points to the origin of the Geometry.
 * @details The origin of the Geometry is designated at its center by
 *          convention, but for track initialization the origin is assumed to be
 *          at the bottom right corner for simplicity. This method corrects
 *          for this by re-assigning the start and end Point coordinates.
 */
void TrackGenerator3D::recalibrateTracksToOrigin() {

  /* Recalibrate the tracks to the origin and set the uid. Note that the
   * loop structure is unconventional in order to preserve an increasing
   * track uid value in the Solver's tracks array. The tracks array is
   * oriented with this loop structure in order to maintain reproducability
   * for parallel runs.
   */

  /* Recalibrate the 2D tracks back to the geometry origin */
  TrackGenerator::recalibrateTracksToOrigin();

  /* Loop over azim reflective halfspaces */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          double x0 = _tracks_3D[a][i][p][z].getStart()->getX();
          double y0 = _tracks_3D[a][i][p][z].getStart()->getY();
          double z0 = _tracks_3D[a][i][p][z].getStart()->getZ();
          double x1 = _tracks_3D[a][i][p][z].getEnd()->getX();
          double y1 = _tracks_3D[a][i][p][z].getEnd()->getY();
          double z1 = _tracks_3D[a][i][p][z].getEnd()->getZ();
          double new_x0 = x0 + _geometry->getMinX();
          double new_y0 = y0 + _geometry->getMinY();
          double new_z0 = z0 + _geometry->getMinZ();
          double new_x1 = x1 + _geometry->getMinX();
          double new_y1 = y1 + _geometry->getMinY();
          double new_z1 = z1 + _geometry->getMinZ();

          _tracks_3D[a][i][p][z]
            .setCoords(new_x0, new_y0, new_z0, new_x1, new_y1, new_z1);
        }
      }
    }
  }

  /* Enusre that all tracks reside within the geometry */
  double max_z = _geometry->getMaxZ();
  double min_z = _geometry->getMinZ();
  for (int a=0; a < _num_azim/2; a++) {
    #pragma omp parallel for
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          double start_z =
            _tracks_3D[a][i][p][z].getStart()->getZ();
          if (start_z > max_z)
            _tracks_3D[a][i][p][z].getStart()->setZ(max_z);
          else if (start_z < min_z)
            _tracks_3D[a][i][p][z].getStart()->setZ(min_z);
        }
      }
    }
  }
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

  return;
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
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          _geometry->segmentize3D(&_tracks_3D[a][i][p][z]);
      }
    }

    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          tracks_segmented++;
      }
    }
  }

  _geometry->initializeFSRVectors();
  _contains_3D_segments = true;

  return;
}


/**
 * @brief Converts a length traveled along a 2D Track cycle to a displacement
 *        in the x-direction
 * @param l The distance traveled in the 2D Track cycle
 * @param azim The azimuthal angle index of the cycle
 * @param cycle The ID of the cycle
 * @return The displacement in the x-direction
 */
double TrackGenerator3D::convertLtoX(double l, int azim, int cycle) {

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to X since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);

  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++) {
    if (l <= length_sum + _tracks_2D_cycle[azim][cycle][i]->getLength()) {
      track_index = i;
      break;
    }
    else{
      length_sum += _tracks_2D_cycle[azim][cycle][i]->getLength();
    }
  }

  if (l - length_sum < 0.0)
    log_printf(ERROR, "found negative length residual in converting l to x");

  double x1 = _tracks_2D_cycle[azim][cycle][track_index]->getStart()->getX();
  double x2 = _tracks_2D_cycle[azim][cycle][track_index]->getEnd()->getX();
  double l_rel = (l - length_sum) / _tracks_2D_cycle[azim][cycle][track_index]
    ->getLength();

  double x;

  if (getCycleDirection(azim, cycle, track_index))
    x = l_rel * x2 + (1.0 - l_rel) * x1;
  else
    x = l_rel * x1 + (1.0 - l_rel) * x2;

  return x;
}


/**
 * @brief Converts a length traveled along a 2D Track cycle to a displacement
 *        in the y-direction
 * @param l The distance traveled in the 2D Track cycle
 * @param azim The azimuthal angle index of the cycle
 * @param cycle The ID of the cycle
 * @return The displacement in the y-direction
 */
double TrackGenerator3D::convertLtoY(double l, int azim, int cycle) {

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to Y since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);

  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++) {
    if (l <= length_sum + _tracks_2D_cycle[azim][cycle][i]->getLength()) {
      track_index = i;
      break;
    }
    else{
      length_sum += _tracks_2D_cycle[azim][cycle][i]->getLength();
    }
  }

  if (l - length_sum < 0.0)
    log_printf(ERROR, "found negative length residual in converting l to y");

  double y1 = _tracks_2D_cycle[azim][cycle][track_index]->getStart()->getY();
  double y2 = _tracks_2D_cycle[azim][cycle][track_index]->getEnd()->getY();
  double l_rel = (l - length_sum) / _tracks_2D_cycle[azim][cycle][track_index]
    ->getLength();

  double y;

  if (getCycleDirection(azim, cycle, track_index))
    y = l_rel * y2 + (1.0 - l_rel) * y1;
  else
    y = l_rel * y1 + (1.0 - l_rel) * y2;

  return y;
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
  for (int a=0; a < _num_azim/2; a++)
    for (int i=0; i < getNumX(a) + getNumY(a); i++)
      for (int p=0; p < _num_polar; p++)
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          if (_tracks_3D[a][i][p][z].getUid() == track_id) {

            coords[0] = _tracks_3D[a][i][p][z].getStart()->getX();
            coords[1] = _tracks_3D[a][i][p][z].getStart()->getY();
            coords[2] = _tracks_3D[a][i][p][z].getStart()->getZ();
            coords[3] = _tracks_3D[a][i][p][z].getEnd()->getX();
            coords[4] = _tracks_3D[a][i][p][z].getEnd()->getY();
            coords[5] = _tracks_3D[a][i][p][z].getEnd()->getZ();
            return;
          }
  log_printf(ERROR, "Unable to find a 3D track associated with the given track"
                    "ID during coordinate retrieval");
  return;
}


/**
 * @brief Sets the track periodic indices of all 3D Tracks
 * @details Periodic cylces are traversed until all 3D Tracks are visited and
 *          their periodic indices are set
 */
void TrackGenerator3D::initializeTrackPeriodicIndices() {

  if (!_periodic)
    return;

  Track* track;
  int track_index;

  /* Set the track periodic cycle indices for 3D tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          track = &_tracks_3D[a][i][p][z];

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
  }
}


/**
 * @brief allocates memory for 3D Tracks
 * @details Before calling this function, the number of tracks per z-stack
 *          should be known and initialized in the _tracks_per_stack 3D array
 */
void TrackGenerator3D::create3DTracksArrays() {

  _tracks_3D = new Track3D***[_num_azim/2];
  for (int a=0; a < _num_azim/2; a++) {
    _tracks_3D[a] = new Track3D**[getNumX(a) + getNumY(a)];
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      _tracks_3D[a][i] = new Track3D*[_num_polar];
      for (int p=0; p < _num_polar; p++) {
        _tracks_3D[a][i][p] = new Track3D[_tracks_per_stack[a][i][p]];
        _num_3D_tracks += _tracks_per_stack[a][i][p];
        _tracks_per_stack[a][i][p] = 0;
      }
    }
  }

  /* Allocate memory for 3D tracks cycles */
  for (int a = 0; a < _num_azim/4; a++) {
    _tracks_3D_cycle[a] = new Track3D****[_cycles_per_azim[a]];
    for (int c = 0; c < _cycles_per_azim[a]; c++) {
      _tracks_3D_cycle[a][c] = new Track3D***[_num_polar];
      for (int p=0; p < _num_polar; p++) {
        _tracks_3D_cycle[a][c][p] =
          new Track3D**[getNumZ(a,p) + getNumL(a,p)];
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++)
          _tracks_3D_cycle[a][c][p][i] =
            new Track3D*[_tracks_per_train[a][c][p][i]];
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
void TrackGenerator3D::initializeTracksArray() {

  /* First initialize 2D Tracks */
  TrackGenerator::initializeTracksArray();

  log_printf(NORMAL, "Initializing 3D tracks array...");

  Track* track;
  int uid = 0;
  int azim_group_id, periodic_group_id, polar_group_id;
  int track_azim_group_id, track_periodic_group_id, track_polar_group_id;
  int track_periodic_index;

  /* Set the number of parallel track groups */
  if (_periodic)
    _num_parallel_track_groups = 12;
  else
    _num_parallel_track_groups = 4;

  /* Create the array of track ids separating the parallel groups */
  if (_num_tracks_by_parallel_group != NULL)
    delete [] _num_tracks_by_parallel_group;
  _num_tracks_by_parallel_group = new int[_num_parallel_track_groups + 1];

  /* Set the first index in the num tracks by parallel group array to 0 */
  _num_tracks_by_parallel_group[0] = 0;

  /* Reset UID in case 2D tracks were intiialized */
  uid = 0;

  /* Allocate memory for tracks array */
  if (_tracks_3D_array != NULL)
    delete [] _tracks_3D_array;
  int num_3D_tracks = getNum3DTracks();
  _tracks_3D_array = new Track*[num_3D_tracks];

  for (int g = 0; g < _num_parallel_track_groups; g++) {

    /* Set the azimuthal, polar, and periodic group ids */
    azim_group_id = g % 2;
    polar_group_id = g % 4 / 2;
    periodic_group_id = g / 4;

    for (int a = 0; a < _num_azim / 2; a++) {
      for (int i = 0; i < getNumX(a) + getNumY(a); i++) {
        for (int p = 0; p < _num_polar; p++) {
          for (int z = 0; z < _tracks_per_stack[a][i][p]; z++) {

            /* Get current track and azim group ids */
            track = &_tracks_3D[a][i][p][z];

            /* Get the track azim group id */
            track_azim_group_id = a / (_num_azim / 4);

            /* Get the track polar group id */
            track_polar_group_id = p / (_num_polar / 2);

            /* Get the track periodic group id */
            if (_periodic) {
              track_periodic_index = track->getPeriodicTrackIndex();

              if (track_periodic_index == 0)
                track_periodic_group_id = 0;
              else if (track_periodic_index % 2 == 1)
                track_periodic_group_id = 1;
              else
                track_periodic_group_id = 2;
            }
            else
              track_periodic_group_id = 0;

            /* Check if track has current azim_group_id
               and periodic_group_id */
            if (azim_group_id == track_azim_group_id &&
                polar_group_id == track_polar_group_id &&
                periodic_group_id == track_periodic_group_id) {
              int azim_index = _quadrature->getFirstOctantAzim(a);
              int polar_index = _quadrature->getFirstOctantPolar(p);
              track->setUid(uid);
              _tracks_3D_array[uid] = track;
              uid++;
            }
          }
        }
      }
    }

    /* Set the track index boundary for this parallel group */
    _num_tracks_by_parallel_group[g + 1] = uid;
  }
}


/**
 * @brief Allocates a new Quadrature with the default Quadrature
 * @details The defualt quadrature for 3D calculations is equal weight
 */
void TrackGenerator3D::initializeDefaultQuadrature() {
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
 * @brief Calculates and assigns the weight for every Track
 */
void TrackGenerator3D::setTotalWeights() {
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          int azim_index = _quadrature->getFirstOctantAzim(a);
          int polar_index = _quadrature->getFirstOctantPolar(p);
          FP_PRECISION weight =
                  _quadrature->getPolarWeight(azim_index, polar_index)
                  * getPolarSpacing(azim_index, polar_index)
                  * _quadrature->getAzimWeight(azim_index)
                  * getAzimSpacing(azim_index);
          _tracks_3D[a][i][p][z].setWeight(weight);
        }
      }
    }
  }
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
    for (int t = 0; t < _num_threads; t++)
      for (int z = 0; z < _num_seg_matrix_rows; z++)
        delete [] _temporary_segments.at(t)[z];
  }
  else {
    _temporary_segments.resize(_num_threads);
    for (int t = 0; t < _num_threads; t++)
      _temporary_segments.at(t) = new segment*[_num_seg_matrix_rows];
    _contains_temporary_segments = true;
  }

  /* Allocate new temporary segments */
  for (int t = 0; t < _num_threads; t++)
    for (int z = 0; z < _num_seg_matrix_rows; z++)
      _temporary_segments.at(t)[z] = new segment[_num_seg_matrix_columns];
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
