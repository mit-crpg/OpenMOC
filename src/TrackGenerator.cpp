#include "TrackGenerator.h"
#include <iomanip>

/**
 * @brief Constructor for the TrackGenerator assigns default values.
 * @param geometry a pointer to a Geometry object
 * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
 * @param spacing track spacing (cm)
 */
TrackGenerator::TrackGenerator(Geometry* geometry, const int num_azim,
                               const int num_polar,
                               const double azim_spacing,
                               const double polar_spacing) {

  setNumThreads(1);

  _geometry = geometry;
  setNumAzim(num_azim);
  setNumPolar(num_polar);
  setDesiredAzimSpacing(azim_spacing);
  setDesiredPolarSpacing(polar_spacing);
  _contains_2D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_tracks = false;
  _contains_3D_segments = false;
  _contains_global_z_mesh = false;
  _contains_segmentation_heights = false;
  _quadrature = NULL;
  _z_coord = 0.0;
  _solve_3D = true;
  _OTF = false;
  _max_optical_length = std::numeric_limits<FP_PRECISION>::max();
  _max_num_segments = 0;
  _track_generation_method = GLOBAL_TRACKING;
  _dump_segments = true;
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator::~TrackGenerator() {

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_3D_tracks) {
    
    /* Delete 3D tracks */
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          delete [] _tracks_3D_stack[a][i][p];
        }
        delete [] _tracks_3D_stack[a][i];
        delete [] _tracks_per_stack[a][i];
      }
      delete [] _tracks_3D_stack[a];
      delete [] _tracks_per_stack[a];
    }
    delete [] _tracks_3D_stack;
    delete [] _tracks_per_stack;

    /* Delete book keeping for 3D tracks */
    for (int a = 0; a < _num_azim/4; a++) {
      for (int c = 0; c < _cycles_per_azim[a]; c++) {
        for (int p=0; p < _num_polar; p++) {
          for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++)
            delete [] _tracks_3D_cycle[a][c][p][i];
          delete [] _tracks_3D_cycle[a][c][p];
          delete [] _tracks_per_train[a][c][p];
        }
        delete [] _tracks_3D_cycle[a][c];
        delete [] _tracks_per_train[a][c];
      }
      delete [] _tracks_3D_cycle[a];
      delete [] _tracks_per_train[a];
    }
    delete [] _tracks_3D_cycle;
    delete [] _tracks_per_train;

    /* Delete book keeping for 3D tracks */
    for (int a=0; a < _num_azim/4; a++) {
      delete [] _num_l[a];
      delete [] _num_z[a];
      delete [] _dl_eff[a];
      delete [] _dz_eff[a];
      delete [] _polar_spacings[a];
    }
    delete [] _num_l;
    delete [] _num_z;
    delete [] _dl_eff;
    delete [] _dz_eff;
    delete [] _polar_spacings;
  }

  if (_contains_2D_tracks) {

    /* Delete 2D tracks */
    for (int a=0; a < _num_azim/2; a++)
      delete [] _tracks_2D[a];
    delete [] _tracks_2D;

    /* Delete track laydown information */
    delete [] _num_x;
    delete [] _num_y;
    delete [] _dx_eff;
    delete [] _dy_eff;
    delete [] _azim_spacings;
    delete [] _cycles_per_azim;
    delete [] _tracks_per_cycle;
    delete [] _cycle_length;
  }

  /* Delete flattened tracks used in OTF calculations */
  if (_contains_flattened_tracks)
    delete [] _flattened_tracks;
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
 * @brief Return the track polar spacing (cm).
 * @details This will return the user-specified track spacing and NOT the
 *          effective track spacing which is computed and used to generate
 *          cyclic tracks.
 * @return the track polar spacing (cm)
 */
double TrackGenerator::getDesiredPolarSpacing() {
  return _polar_spacing;
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
 * @brief Return the total number of 2D Tracks across the Geometry.
 * @return the total number of 2D Tracks
 */
int TrackGenerator::getNum2DTracks() {

  int num_2D_tracks = 0;

  for (int a=0; a < _num_azim/2; a++)
    num_2D_tracks += getNumX(a) + getNumY(a);

  return num_2D_tracks;
}


/**
 * @brief Return the total number of 3D Tracks across the Geometry.
 * @return the total number of 3D Tracks
 */
int TrackGenerator::getNum3DTracks() {

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
int TrackGenerator::getNum2DSegments() {

  if (!contains2DSegments())
    log_printf(ERROR, "Cannot get the number of 2D segments since they "
               "have not been generated.");

  int num_2D_segments = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      num_2D_segments += _tracks_2D[a][i].getNumSegments();
    }
  }

  return num_2D_segments;
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator::getNum3DSegments() {

  if ((!_OTF && !contains3DSegments()) ||
      (_OTF && !contains2DSegments()))
    log_printf(ERROR, "Cannot get the number of 3D segments since they "
               "have not been generated.");

  int num_3D_segments = 0;

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          num_3D_segments += _tracks_3D_stack[a][i][p][z].getNumSegments();
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
Track** TrackGenerator::getTracksArray() {

  if (!contains2DTracks() && !contains3DTracks())
    log_printf(ERROR, "Unable to return the 1D array of Tracks "
               "since Tracks have not yet been generated.");

  return _tracks;
}


/**
 * @brief Returns an array of the flattend Track pointers by increasing UID.
 * @details An array of pointers to all 2D Track objects for 3D on-the-fly
 *          calculation in the Geometry is returned, arranged by increasing 
 *          2D Track unique identifier (UID).
 * @return the array of flattened Track pointers
 */
Track** TrackGenerator::getFlattenedTracksArray() {

  if (!containsFlattenedTracks())
    log_printf(ERROR, "Unable to return the 1D array of Tracks "
               "since Tracks have not yet been generated.");

  return _flattened_tracks;
}


/**
 * @brief Returns a 2D jagged array of the 2D Tracks.
 * @details The first index into the array is the azimuthal angle and the
 *          second index is the Track number.
 * @return the 2D jagged array of 2D Tracks
 */
Track2D** TrackGenerator::get2DTracks() {

  if (!contains2DTracks())
    log_printf(ERROR, "Unable to return the 3D ragged array of the 2D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_2D;
}


/**
 * @brief Returns a 4D jagged array of the 3D Tracks.
 * @details The first index into the array is the azimuthal angle, the second
 *          index is the 2D Track number, the third index is the polar angle,
 *          and the fourth index is the z-stack number.
 * @return the 4D jagged array of 3D Tracks
 */
Track3D**** TrackGenerator::get3DTracks() {

  if (!contains3DTracks())
    log_printf(ERROR, "Unable to return the 3D ragged array of the 3D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_3D_stack;
}


/**
 * @brief Returns an array of adjusted azimuthal spacings
 * @details An array of azimuthal spacings after adjustment is returned,
 *          indexed by azimuthal angle
 * @return the array of azimuthal spacings
 */
double* TrackGenerator::getAzimSpacings() {
  return _azim_spacings;
}


/**
 * @brief Returns the adjusted azimuthal spacing at the requested azimuthal
 *        angle index
 * @details The aziumthal spacing depends on the azimuthal angle. This function
 *          returns the azimuthal spacing used at the desired azimuthal angle
 *          index.
 * @param azim the requested azimuthal angle index
 * @return the requested azimuthal spacing
 */
double TrackGenerator::getAzimSpacing(int azim) {
  azim = _quadrature->getFirstOctantAzim(azim);
  return _azim_spacings[azim];
}


/**
 * @brief Returns a 2D array of adjusted polar spacings
 * @details An array of polar spacings after adjustment is returned,
 *          indexed first by azimuthal angle and then by polar angle
 * @return the 2D array of polar spacings
 */
double** TrackGenerator::getPolarSpacings() {
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
double TrackGenerator::getPolarSpacing(int azim, int polar) {
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
  return _polar_spacings[azim][polar];
}


/**
 * @brief Get the maximum allowable optical length for a track segment
 * @return The max optical length
 */
FP_PRECISION TrackGenerator::getMaxOpticalLength() {

  segment* curr_segment;
  FP_PRECISION length;
  Material* material;
  FP_PRECISION* sigma_t;
  FP_PRECISION max_optical_length = 0.;

  if (_solve_3D) {

    /* Allocate array for 3D segments for OTF computation */
    segment* segments_3D;
    if (_OTF)
      segments_3D = new segment[_max_num_segments];

    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

            /* Extract 3D track and initialize segments pointer */
            Track* track_3D = &_tracks_3D_stack[a][i][p][z];
            int num_segments = track_3D->getNumSegments();

            /* Get the segments corresponding to the 3D track */
            if (_OTF) {

              Point* start = _tracks_3D_stack[a][i][p][z].getStart();
              double theta = _tracks_3D_stack[a][i][p][z].getTheta();
              Track2D* flattened_track = &_tracks_2D[a][i];

              SegmentationKernel kernel;
              kernel.setSegments(segments_3D);
              traceSegmentsOTF(flattened_track, start, theta, &kernel);
            }
            else
              segments_3D = track_3D->getSegments();

            /* Look through all segments for max optical path length */
            for (int s=0; s < num_segments; s++) {
              curr_segment = &segments_3D[s];
              length = curr_segment->_length;
              material = curr_segment->_material;
              sigma_t = material->getSigmaT();

              for (int e=0; e < material->getNumEnergyGroups(); e++)
                max_optical_length = std::max(max_optical_length,
                                              length*sigma_t[e]);
            }
          }
        }
      }
    }
    if (_OTF)
      delete[] segments_3D;
  }
  else{
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++) {
          curr_segment = _tracks_2D[a][i].getSegment(s);
          length = curr_segment->_length;
          material = curr_segment->_material;
          sigma_t = material->getSigmaT();

          for (int e=0; e < material->getNumEnergyGroups(); e++)
            max_optical_length = std::max(max_optical_length,
                                          length*sigma_t[e]);
        }
      }
    }
  }
  return max_optical_length;
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
 * @brief Returns a 3D array of the number of 3D Tracks in each z-stack
 * @details A 3D array is returned indexed first by azimuthal angle, second by
 *          2D track number, and third by polar angle. This array describes
 *          the number of tracks in each z-stack.
 * @return A 3D array of the number of tracks in each z-stack
 */
int*** TrackGenerator::getTracksPerStack() {
  return _tracks_per_stack;
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
 * @param azim the azimuthal angle index
 * @return the number of 2D Tracks in the cycle
 */
double TrackGenerator::getCycleLength(int azim) {
  azim = _quadrature->getFirstOctantAzim(azim);
  return _cycle_length[azim];
}


/**
 * @brief Returns the number of 2D Tracks in the x-direction for a given
 *        azimuthal angle index
 * @param azim the azimuthal angle index
 * @return the number of 2D Tracks in the x-direction of the Geometry
 */
int TrackGenerator::getNumX(int azim) {
  azim = _quadrature->getFirstOctantAzim(azim);
  return _num_x[azim];
}


/**
 * @brief Returns the number of 2D Tracks in the y-direction for a given
 *        azimuthal angle index
 * @param azim the azimuthal angle index
 * @return the number of 2D Tracks in the y-direction of the Geometry
 */
int TrackGenerator::getNumY(int azim) {
  azim = _quadrature->getFirstOctantAzim(azim);
  return _num_y[azim];
}


/**
 * @brief Returns the number of 3D Tracks in the z-direction for a given
 *        azimuthal angle index and polar angle index
 * @param azim the azimuthal angle index
 * @param polar the polar angle index
 * @return the number of 3D Tracks in the z-direction of the Geometry
 */
int TrackGenerator::getNumZ(int azim, int polar) {
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
int TrackGenerator::getNumL(int azim, int polar) {
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
  return _num_l[azim][polar];
}


/**
 * @brief Returns the spacing between Tracks in the x-direction for a given
 *        azimuthal angle index
 * @param azim the azimuthal angle index
 * @return the spacing between Tracks in the x-direction
 */
double TrackGenerator::getDxEff(int azim) {
  azim = _quadrature->getFirstOctantAzim(azim);
  return _dx_eff[azim];
}


/**
 * @brief Returns the spacing between Tracks in the y-direction for a given
 *        azimuthal angle index
 * @param azim the azimuthal angle index
 * @return the spacing between Tracks in the y-direction
 */
double TrackGenerator::getDyEff(int azim) {
  azim = _quadrature->getFirstOctantAzim(azim);
  return _dy_eff[azim];
}


/**
 * @brief Computes and returns an array of volumes indexed by FSR.
 * @details Note: It is the function caller's responsibility to deallocate
 *          the memory reserved for the FSR volume array.
 * @return a pointer to the array of FSR volumes
 */
FP_PRECISION* TrackGenerator::get2DFSRVolumes() {

  if (!contains2DSegments())
    log_printf(ERROR, "Unable to get the FSR volumes since 2D tracks "
               "have not yet been generated");

  int num_FSRs = _geometry->getNumFSRs();
  FP_PRECISION *FSR_volumes = new FP_PRECISION[num_FSRs];
  memset(FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));

  segment* segment;
  FP_PRECISION volume;

  /* Calculate each FSR's "volume" by accumulating the total length of *
   * all Track segments multiplied by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++) {
        segment = _tracks_2D[a][i].getSegment(s);
        volume = segment->_length * _quadrature->getAzimWeight(a)
          * getAzimSpacing(a);
        FSR_volumes[segment->_region_id] += volume;
      }
    }
  }

  return FSR_volumes;
}


/**
 * @brief Computes and returns an array of volumes indexed by FSR.
 * @details Note: It is the function caller's responsibility to deallocate
 *          the memory reserved for the FSR volume array.
 * @return a pointer to the array of FSR volumes
 */
FP_PRECISION* TrackGenerator::get3DFSRVolumes() {

  /* Determine whether to calculate on-the-fly */
  if (_OTF)
    return get3DFSRVolumesOTF();

  if (!contains3DSegments())
    log_printf(ERROR, "Unable to get the FSR volumes since 3D tracks "
               "have not yet been generated");

  int num_FSRs = _geometry->getNumFSRs();
  FP_PRECISION *FSR_volumes = new FP_PRECISION[num_FSRs];
  memset(FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));

  segment* segment;
  FP_PRECISION volume;

  /* Calculate each FSR's "volume" by accumulating the total length of *
   * all Track segments multiplied by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments();
              s++) {

            segment = _tracks_3D_stack[a][i][p][z].getSegment(s);
            volume = segment->_length * _quadrature->getAzimWeight(a)
              * _quadrature->getPolarWeight(a, p) * getAzimSpacing(a)
              * getPolarSpacing(a,p);
            FSR_volumes[segment->_region_id] += volume;

          }
        }
      }
    }
  }

  return FSR_volumes;
}


/**
 * @brief Computes and returns the volume of an FSR.
 * @param fsr_id the ID for the FSR of interest
 * @return the FSR volume
 */
FP_PRECISION TrackGenerator::get2DFSRVolume(int fsr_id) {

  if (!contains2DSegments())
    log_printf(ERROR, "Unable to get the FSR volume since 2D tracks "
               "have not yet been generated");

  else if (fsr_id < 0 || fsr_id > _geometry->getNumFSRs())
    log_printf(ERROR, "Unable to get the volume for FSR %d since the FSR IDs "
               "lie in the range (0, %d)", fsr_id, _geometry->getNumFSRs());

  segment* segment;
  FP_PRECISION volume = 0.0;

  /* Calculate each FSR's "volume" by accumulating the total length of *
   * all Track segments multiplied by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++) {
        segment = _tracks_2D[a][i].getSegment(s);
        if (segment->_region_id == fsr_id)
          volume += segment->_length * _quadrature->getAzimWeight(a)
            * getAzimSpacing(a);
      }
    }
  }

  return volume;
}


/**
 * @brief Computes and returns the volume of an FSR.
 * @param fsr_id the ID for the FSR of interest
 * @return the FSR volume
 */
FP_PRECISION TrackGenerator::get3DFSRVolume(int fsr_id) {

  if (!contains3DSegments())
    log_printf(ERROR, "Unable to get the FSR volume since 3D tracks "
               "have not yet been generated");

  else if (fsr_id < 0 || fsr_id > _geometry->getNumFSRs())
    log_printf(ERROR, "Unable to get the volume for FSR %d since the FSR IDs "
               "lie in the range (0, %d)", fsr_id, _geometry->getNumFSRs());

  segment* segment;
  FP_PRECISION volume = 0.0;

  /* Calculate each FSR's "volume" by accumulating the total length of *
   * all Track segments multiplied by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments();
              s++) {
            segment = _tracks_3D_stack[a][i][p][z].getSegment(s);
            if (segment->_region_id == fsr_id)
              volume += segment->_length * _quadrature->getAzimWeight(a)
                * _quadrature->getPolarWeight(a, p) * getAzimSpacing(a)
                * getPolarSpacing(a,p);
          }
        }
      }
    }
  }

  return volume;
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
  _contains_2D_tracks = false;
  _contains_flattened_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _use_input_file = false;
  _tracks_filename = "";
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
  _contains_2D_tracks = false;
  _contains_flattened_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _use_input_file = false;
  _tracks_filename = "";
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
  _contains_2D_tracks = false;
  _contains_flattened_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _use_input_file = false;
  _tracks_filename = "";
}


/**
 * @brief Set the suggested track polar spacing (cm).
 * @param spacing the suggested track polar spacing
 */
void TrackGenerator::setDesiredPolarSpacing(double spacing) {
  if (spacing < 0)
    log_printf(ERROR, "Unable to set a negative track polar spacing "
               "%f for the TrackGenerator.", spacing);

  _polar_spacing = spacing;
  _contains_2D_tracks = false;
  _contains_flattened_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _use_input_file = false;
  _tracks_filename = "";
}


/**
 * @brief Set a pointer to the Geometry to use for track generation.
 * @param geometry a pointer to the Geometry
 */
void TrackGenerator::setGeometry(Geometry* geometry) {
  _geometry = geometry;
  _contains_2D_tracks = false;
  _contains_flattened_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _contains_global_z_mesh = false;
  _use_input_file = false;
  _tracks_filename = "";
}


/**
 * @brief sets a flag to solve the problem in two dimensions
 */
void TrackGenerator::setSolve2D() {
  _solve_3D = false;
}


/**
 * @brief sets a flag to solve the problem in three dimensions
 */
void TrackGenerator::setSolve3D() {
  _solve_3D = true;
}


/**
 * @brief sets a flag to calculate 3D segments on-the-fly in the axial
 *        direction
 */
void TrackGenerator::setOTF() {
  _OTF = true;
}


/**
 * @brief Sets the z-coord of the raidal plane used in 2D calculations
 * @param z_coord the z-coord of the radial plane
 */
void TrackGenerator::setZCoord(double z_coord) {
  _z_coord = z_coord;
}

/**
 * @brief Sets the z-planes over which 2D segmentation is performed for
 *        on-the-fly calculations
 * @param z_mesh the z-coordinates defining the height of the radial
 *        segmentation planes
 */
void TrackGenerator::setSegmentationHeights(std::vector<double> z_mesh) {
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
void TrackGenerator::setGlobalZMesh() {
  _contains_global_z_mesh = true;
  _global_z_mesh = _geometry->getUniqueZHeights();
}


/**
 * @brief sets the Quadrature used for integrating the MOC equations
 * @param quadrature a pointer to the Quadrature object used in calculation
 */
void TrackGenerator::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
}


/**
 * @brief Returns whether or not the TrackGenerator contains 2D Tracks
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains 2D Tracks; false otherwise
 */
bool TrackGenerator::contains2DTracks() {
  return _contains_2D_tracks;
}


/**
 * @brief Returns whether or not the TrackGenerator contains flattened Tracks
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains flattened Tracks
 *         false otherwise
 */
bool TrackGenerator::containsFlattenedTracks() {
  return _contains_flattened_tracks;
}


/**
 * @brief Returns whether or not the TrackGenerator contains 3D Tracks
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains 3D Tracks; false otherwise
 */
bool TrackGenerator::contains3DTracks() {
  return _contains_3D_tracks;
}


/**
 * @brief Returns whether or not the TrackGenerator contains 2D segments
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains 2D segments; false otherwise
 */
bool TrackGenerator::contains2DSegments() {
  return _contains_2D_segments;
}


/**
 * @brief Returns whether or not the TrackGenerator contains 3D segments
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains 3D segments; false otherwise
 */
bool TrackGenerator::contains3DSegments() {
  return _contains_3D_segments;
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
void TrackGenerator::retrieve2DTrackCoords(double* coords, int num_tracks) {

  if (num_tracks != 4*getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the Track coordinates since the "
               "TrackGenerator contains %d Tracks with %d coordinates but an "
               "array of length %d was input",
               getNum2DTracks(), 4*getNum2DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
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

  if (num_tracks != 5*getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the 2D Track periodic cycle "
               "coordinates since the TrackGenerator contains %d Tracks with "
               "%d coordinates but an array of length %d was input",
               getNum2DTracks(), 5*getNum2DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      coords[counter]   = _tracks_2D[a][i].getStart()->getX();
      coords[counter+1] = _tracks_2D[a][i].getStart()->getY();
      coords[counter+2] = _tracks_2D[a][i].getEnd()->getX();
      coords[counter+3] = _tracks_2D[a][i].getEnd()->getY();
      coords[counter+4] = _tracks_2D[a][i].getPeriodicCycleId();

      counter += 5;
    }
  }

  return;
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

  if (num_tracks != 5*getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the 2D Track reflective cycle "
               "coordinates since the TrackGenerator contains %d Tracks with "
               "%d coordinates but an array of length %d was input",
               getNum2DTracks(), 5*getNum2DTracks(), num_tracks);

  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;

  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      coords[counter]   = _tracks_2D[a][i].getStart()->getX();
      coords[counter+1] = _tracks_2D[a][i].getStart()->getY();
      coords[counter+2] = _tracks_2D[a][i].getEnd()->getX();
      coords[counter+3] = _tracks_2D[a][i].getEnd()->getY();
      coords[counter+4] = _tracks_2D[a][i].getReflectiveCycleId();

      counter += 5;
    }
  }

  return;
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
void TrackGenerator::retrieve3DPeriodicCycleCoords(double* coords,
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
          coords[counter]   = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D_stack[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D_stack[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D_stack[a][i][p][z].getEnd()->getZ();
          coords[counter+6] =
            _tracks_3D_stack[a][i][p][z].getPeriodicCycleId();
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
void TrackGenerator::retrieve3DReflectiveCycleCoords(double* coords,
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
          coords[counter]   = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D_stack[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D_stack[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D_stack[a][i][p][z].getEnd()->getZ();
          coords[counter+6] =
            _tracks_3D_stack[a][i][p][z].getReflectiveCycleId();
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
void TrackGenerator::retrieve3DTrackCoords(double* coords, int num_tracks) {

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
          coords[counter]   = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D_stack[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D_stack[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D_stack[a][i][p][z].getEnd()->getZ();
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
void TrackGenerator::retrieve3DSegmentCoords(double* coords, int num_segments) {

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

          x0    = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          y0    = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          z0    = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          phi   = _tracks_3D_stack[a][i][p][z].getPhi();
          theta = _tracks_3D_stack[a][i][p][z].getTheta();

          segments = _tracks_3D_stack[a][i][p][z].getSegments();

          for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments();
              s++) {

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
void TrackGenerator::retrieve2DSegmentCoords(double* coords, int num_segments) {

  if (num_segments != 5*getNum2DSegments())
    log_printf(ERROR, "Unable to retrieve the Track segment coordinates since "
               "the TrackGenerator contains %d segments with %d coordinates "
               "but an array of length %d was input",
               getNum2DSegments(), 5*getNum2DSegments(), num_segments);

  segment* curr_segment = NULL;
  double x0, x1, y0, y1;
  double phi;
  segment* segments;

  int counter = 0;

  /* Loop over Track segments and populate array with their FSR ID and *
   * start/end points */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

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
 * @brief Checks boundary conditions for inconsistent periodic boundary
 *        conditions
 */
void TrackGenerator::checkBoundaryConditions() {

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
  else if ((_geometry->getMinZBoundaryType() == PERIODIC &&
            _geometry->getMaxZBoundaryType() != PERIODIC) ||
           (_geometry->getMinZBoundaryType() != PERIODIC &&
            _geometry->getMaxZBoundaryType() == PERIODIC))
    log_printf(ERROR, "Cannot create tracks with only one z boundary"
               " set to PERIODIC");

  /* Check for correct track method if a PERIODIC bc is present */
  if (_geometry->getMinXBoundaryType() == PERIODIC ||
      _geometry->getMinYBoundaryType() == PERIODIC ||
      _geometry->getMinZBoundaryType() == PERIODIC) {

    _periodic = true;

    if (_track_generation_method != MODULAR_RAY_TRACING &&
        _track_generation_method != SIMPLIFIED_MODULAR_RAY_TRACING &&
        _solve_3D)
      log_printf(ERROR, "Cannot create tracks for a geometry containing a"
                 " periodic BC with a track generation method that is not"
                 " modular");
  }
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

  if (_geometry == NULL)
    log_printf(ERROR, "Unable to generate Tracks since no Geometry "
               "has been set for the TrackGenerator");

  /* Check to make sure that height, width of the Geometry are nonzero */
  if (_geometry->getHeight() <= 0 || _geometry->getHeight() <= 0 ||
      _geometry->getDepth() <= 0)
    log_printf(ERROR, "The total height, width, and depth of the Geometry must"
               " be nonzero for Track generation. Create a CellFill which "
               "is filled by the entire geometry and bounded by XPlanes, "
               "YPlanes, and ZPlanes to enable the Geometry to determine the "
               "total width, height, and depth of the model.");

  /* Generate Tracks, perform ray tracing across the geometry, and store
   * the data to a Track file */
  try {

    /* Create default quadrature set if user one has not been set */
    if (_quadrature == NULL) {
      if (_solve_3D)
        _quadrature = new EqualWeightPolarQuad();
      else
        _quadrature = new TYPolarQuad();
    }

    /* Initialize the quadrature set */
    _quadrature->setNumPolarAngles(_num_polar);
    _quadrature->setNumAzimAngles(_num_azim);
    _quadrature->initialize();

    /* Check periodic BCs for symmetry */
    checkBoundaryConditions();

    /* Initialize the 2D tracks */
    initialize2DTracks();

    /* If 3D problem, initialize the 3D tracks */
    if (_solve_3D)
      initialize3DTracks();

    /* Recalibrate the 2D tracks back to the geometry origin */
    recalibrate2DTracksToOrigin();

    /* If 3D problem, recalibrate the 3D tracks back to the geometry origin */
    if (_solve_3D)
      recalibrate3DTracksToOrigin();

    /* Initialize the track file directory and read in tracks if they exist */
    initializeTrackFileDirectory();

    /* Initialize the 1D array of Tracks */
    initializeTracksArray();

    /* If track file not present, generater segments */
    if (_use_input_file == false) {

      /* Segmentize the tracks */
      if (_solve_3D) {
        if (_OTF)
          segmentizeExtruded();
        else {
          segmentize3D();
          dump3DSegmentsToFile();
        }
      }
      else{
        segmentize2D();
        dump2DSegmentsToFile();
      }
    }

    /* Precompute the quadrature weights */
    _quadrature->precomputeWeights(_solve_3D);
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to allocate memory needed to generate "
               "Tracks. Backtrace:\n%s", e.what());
  }

  return;
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
 * @brief Returns whether or not the solver is set to 2D
 * @return true if the solver is set to 2D; false otherwise
 */
bool TrackGenerator::isSolve2D() {
  return !_solve_3D;
}


/**
 * @brief Returns whether or not the solver is set to 3D
 * @return true if the solver is set to 3D; false otherwise
 */
bool TrackGenerator::isSolve3D() {
  return _solve_3D;
}


/**
 * @brief Returns whether or not the solver is set to forming 3D segments
 *        on-the-fly
 * @return true if the solver is set to axial on-the-fly segmentation; 
 *         false otherwise
 */
bool TrackGenerator::isOTF() {
  return _OTF;
}


/**
 * @brief Initializes Track azimuthal angles, start and end Points.
 * @details This method computes the azimuthal angles and effective track
 *          spacing to use to guarantee cyclic Track wrapping. Based on the
 *          angles and spacing, the number of Tracks per angle and the start
 *          and end Points for each Track are computed.
 */
void TrackGenerator::initialize2DTracks() {

  log_printf(NORMAL, "Initializing 2D tracks...");

  /* Allocate memory for arrays */
  _dx_eff    = new double[_num_azim/4];
  _dy_eff    = new double[_num_azim/4];
  _tracks_per_cycle = new int[_num_azim/4];
  _cycles_per_azim  = new int[_num_azim/4];
  _tracks_2D        = new Track2D*[_num_azim/2];
  _num_x            = new int[_num_azim/4];
  _num_y            = new int[_num_azim/4];
  _cycle_length     = new double[_num_azim/4];
  _azim_spacings    = new double[_num_azim/4];
  _num_2D_tracks    = 0;

  double x1, x2, y1, y2;
  double phi;
  double width  = _geometry->getWidth();
  double height = _geometry->getHeight();

  /* Determine angular quadrature and track spacing */
  for (int a = 0; a < _num_azim/4; a++) {

    /* Get the desired azimuthal angle */
    phi = _quadrature->getPhi(a);

    /* The number of intersections with x,y-axes */
    _num_x[a] = (int) (fabs(width / _azim_spacing * sin(phi))) + 1;
    _num_y[a] = (int) (fabs(height / _azim_spacing * cos(phi))) + 1;

    /* Effective/actual angle (not the angle we desire, but close) */
    _quadrature->setPhi(atan((height * _num_x[a]) / (width * _num_y[a])), a);

    /* Effective Track spacing (not spacing we desire, but close) */
    _dx_eff[a]   = (width / _num_x[a]);
    _dy_eff[a]   = (height / _num_y[a]);
    _azim_spacings[a] = (_dx_eff[a] * sin(_quadrature->getPhi(a)));

    /* The length of all tracks in a 2D cycle */
    _cycle_length[a] = _dx_eff[a] / cos(_quadrature->getPhi(a)) *
      leastCommonMultiple(2 * _num_x[a], 2 * height /
                          (tan(_quadrature->getPhi(a)) * _dx_eff[a]));

    /* Get the number of tracks per cycle */
    _tracks_per_cycle[a] = (int)
      (round(_cycle_length[a] * sin(_quadrature->getPhi(a)) / width) +
       round(_cycle_length[a] * cos(_quadrature->getPhi(a)) / height));

    /* Compute the number of cycles */
    _cycles_per_azim[a] = (_num_x[a] + _num_y[a]) * 2 / _tracks_per_cycle[a];
  }

  /* Generate 2D tracks */
  for (int a=0; a < _num_azim/2; a++) {

    /* Allocate memory for the 2D tracks array */
    _tracks_2D[a] = new Track2D[getNumX(a) + getNumY(a)];
    _num_2D_tracks += getNumX(a) + getNumY(a);

    /* Get the azimuthal angle for all tracks with this azimuthal angle */
    phi = _quadrature->getPhi(a);

    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

      /* Get track and set angle and track indices */
      Track2D* track = (&_tracks_2D[a][i]);
      track->setPhi(phi);
      track->setAzimIndex(a);
      track->setXYIndex(i);

      /* Set start point */
      if (a < _num_azim/4) {
        if (i < getNumX(a))
          track->getStart()->setCoords(width - getDxEff(a) * (i + 0.5), 0.0);
        else
          track->getStart()->setCoords(0.0, getDyEff(a) *
                                        ((i-getNumX(a)) + 0.5));
      }
      else{
        if (i < getNumX(a))
          track->getStart()->setCoords(getDxEff(a) * (i + 0.5), 0.0);
        else
          track->getStart()->setCoords(width, getDyEff(a) *
                                       ((i-getNumX(a)) + 0.5));
      }

      /* Set end point */
      if (a < _num_azim/4) {
        if (i < getNumY(a))
          track->getEnd()->setCoords(width, getDyEff(a) * (i + 0.5));
        else
          track->getEnd()->setCoords(width - getDxEff(a) *
                                      ((i-getNumY(a)) + 0.5), height);
      }
      else{
        if (i < getNumY(a))
          track->getEnd()->setCoords(0.0, getDyEff(a) * (i + 0.5));
        else
          track->getEnd()->setCoords(getDxEff(a) * ((i-getNumY(a)) + 0.5),
                                     height);
      }
    }
  }

  /* Set the flag indicating 2D tracks have been generated */
  _contains_2D_tracks = true;

  /* Initialize the track reflections and cycle ids */
  initialize2DTrackReflections();
  initialize2DTrackCycleIds();
  initialize2DTrackPeriodicIndices();
}


/**
 * @brief Initializes 2D Track reflections
 * @details This method computes the connecting Tracks for all 2D Tracks in
 *          the TrackGenerator analytically, handling both reflective and
 *          periodic boundaries.
 */
void TrackGenerator::initialize2DTrackReflections() {

  log_printf(NORMAL, "Initializing 2D tracks reflections...");

  Track* track;
  int ac;

  /* Generate the 2D track cycles */
  for (int a=0; a < _num_azim/2; a++) {
    ac = _num_azim/2 - a - 1;
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

      /* Get current track */
      track = &_tracks_2D[a][i];

      /* Set connecting tracks in forward direction */
      if (i < getNumY(a)) {
        track->setReflFwdFwd(true);
        track->setTrackReflFwd(&_tracks_2D[ac][i + getNumX(a)]);
        track->setTrackPrdcFwd(&_tracks_2D[a][i + getNumX(a)]);
      }
      else{
        track->setReflFwdFwd(false);
        track->setTrackReflFwd
          (&_tracks_2D[ac][(getNumX(a) + getNumY(a)) - (i - getNumY(a)) - 1]);
        track->setTrackPrdcFwd(&_tracks_2D[a][i - getNumY(a)]);
      }

      /* Set connecting tracks in backward direction */
      if (i < getNumX(a)) {
        track->setReflBwdFwd(true);
        track->setTrackReflBwd(&_tracks_2D[ac][getNumX(a) - i - 1]);
        track->setTrackPrdcBwd(&_tracks_2D[a][i + getNumY(a)]);
      }
      else{
        track->setReflBwdFwd(false);
        track->setTrackReflBwd(&_tracks_2D[ac][i - getNumX(a)]);
        track->setTrackPrdcBwd(&_tracks_2D[a][i - getNumX(a)]);
      }

      /* Set the foward boundary conditions */
      if (a < _num_azim/4) {
        if (i < getNumY(a))
          track->setBCFwd(_geometry->getMaxXBoundaryType());
        else
          track->setBCFwd(_geometry->getMaxYBoundaryType());

        if (i < getNumX(a))
          track->setBCBwd(_geometry->getMinYBoundaryType());
        else
          track->setBCBwd(_geometry->getMinXBoundaryType());
      }

      /* Set the backward boundary conditions */
      else{
        if (i < getNumY(a))
          track->setBCFwd(_geometry->getMinXBoundaryType());
        else
          track->setBCFwd(_geometry->getMaxYBoundaryType());

        if (i < getNumX(a))
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
void TrackGenerator::initialize2DTrackCycleIds() {

  log_printf(NORMAL, "Initializing 2D track cycle ids...");

  int id = 0;
  bool fwd;
  Track* track;

  /* Set the periodic track cycle ids */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

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
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

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
 * @brief Initializes Track azimuthal angles, start and end Points.
 * @details This method computes the azimuthal angles and effective track
 *          spacing to use to guarantee cyclic Track wrapping. Based on the
 *          angles and spacing, the number of Tracks per angle and the start
 *          and end Points for each Track are computed.
 */
void TrackGenerator::initialize3DTracks() {

  if (!_contains_2D_tracks)
    log_printf(ERROR, "Cannot initialize 3D tracks since the 2D tracks "
               "have not been created");

  else if (_quadrature->getQuadratureType() == TABUCHI_YAMAMOTO ||
           _quadrature->getQuadratureType() == LEONARD ||
           _quadrature->getQuadratureType() == GAUSS_LEGENDRE)
    log_printf(ERROR, "Cannot initialize 3D tracks with a quadrature "
               "type of TABUCHI_YAMAMOTO, LEONARD, or GAUSS_LEGENDRE");

  log_printf(NORMAL, "Initializing 3D tracks...");

  /* Allocate memory for arrays */
  _dz_eff    = new double*[_num_azim/4];
  _dl_eff    = new double*[_num_azim/4];
  _num_z            = new int*[_num_azim/4];
  _num_l            = new int*[_num_azim/4];
  _polar_spacings   = new double*[_num_azim/4];
  _tracks_3D_cycle  = new Track3D*****[_num_azim/4];
  _tracks_per_train = new int***[_num_azim/4];
  _num_3D_tracks    = 0;

  for (int i=0; i < _num_azim/4; i++) {
    _dz_eff[i]         = new double[_num_polar/2];
    _dl_eff[i]         = new double[_num_polar/2];
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
  double depth  = _geometry->getDepth();

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
        double cycle_length = sqrt(_dx_eff[i]*_dx_eff[i]
                                   + _dy_eff[i]*_dy_eff[i]);
        _num_l[i][j] =  (int)
          round(_cycle_length[i] / cycle_length *
                ceil(fabs(cycle_length * tan(M_PI_2 - theta)
                          * sin(theta) / _polar_spacing)));
      }

      /* Number of crossings along the z axis */
      _num_z[i][j] = (int) (fabs(depth * _num_l[i][j] * tan(theta)
                                 / _cycle_length[i])) + 1;

      /* Effective track spacing */
      _dl_eff[i][j]          = _cycle_length[i] / _num_l[i][j];
      _dz_eff[i][j]          = depth            / _num_z[i][j];
      _quadrature->setTheta(atan(_dl_eff[i][j] / _dz_eff[i][j]), i, j);
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
  int pc;
  double l_start, l_end;

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

    /* Allocate memory for 3D track stacks */
    if (create_tracks)
      create3DTracksArrays();

    /* Loop over 3D track cycles */
    #pragma omp parallel for private(track_2D, track_3D, pc, l_start, l_end, \
                                     x1, y1, z1, x2, y2, z2)
    for (int a = 0; a < _num_azim/4; a++) {
      for (int c = 0; c < _cycles_per_azim[a]; c++) {

        /* Loop over polar angles < PI/2 */
        for (int p=0; p < _num_polar/2; p++) {

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
              l_start = _cycle_length[a] - (i + 0.5) * _dl_eff[a][p];
              x1 = convertLtoX(l_start, a, c);
              y1 = convertLtoY(l_start, a, c);
              z1 = 0.0;
            }
            else{
              l_start = 0.0;
              track_2D = getTrack2DByCycle(a, c, 0);
              x1 = track_2D->getStart()->getX();
              y1 = track_2D->getStart()->getY();
              z1 = _dz_eff[a][p] * (i - _num_l[a][p] + 0.5);
            }

            /* Get the end point */
            if (i < _num_z[a][p]) {
              l_end = _cycle_length[a];
              track_2D = getTrack2DByCycle(a, c, 0);
              x2 = track_2D->getStart()->getX();
              y2 = track_2D->getStart()->getY();
              z2 = _dz_eff[a][p] * (i + 0.5);
            }
            else{
              l_end = _cycle_length[a] - _dl_eff[a][p] *
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

        /* Create tracks for polar angles [PI/2,PI] */
        for (int p=0; p < _num_polar/2; p++) {

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
              track_2D = getTrack2DByCycle(a, c, 0);
              x1 = track_2D->getStart()->getX();
              y1 = track_2D->getStart()->getY();
              z1 = _dz_eff[a][p] * (i + 0.5);
            }
            else{
              l_start = _dl_eff[a][p] * (i - _num_z[a][p] + 0.5);
              x1 = convertLtoX(l_start, a, c);
              y1 = convertLtoY(l_start, a, c);
              z1 = depth;
            }

            /* Get the end point */
            if (i < _num_l[a][p]) {
              l_end = _dl_eff[a][p] * (i + 0.5);
              x2 = convertLtoX(l_end, a, c);
              y2 = convertLtoY(l_end, a, c);
              z2 = 0.0;
            }
            else{
              l_end = _cycle_length[a];
              track_2D = getTrack2DByCycle(a, c, 0);
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
      }
    }
  }

  _contains_3D_tracks = true;

  /* Initialize the 3D track reflections and cycle ids */
  initialize3DTrackReflections();
  initialize3DTrackCycleIds();
  initialize3DTrackPeriodicIndices();
}


/**
 * @brief Initializes 3D Track reflections
 * @details This method computes the connecting Tracks for all 3D Tracks in
 *          the TrackGenerator analytically, handling both reflective and
 *          periodic boundaries.
 */
void TrackGenerator::initialize3DTrackReflections() {

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
                    xp = getTrack2DByCycle(a, c, ci)
                      ->getTrackPrdcFwd()->getXYIndex();
                    zp = polar_group[_num_l[a][p] + i][0]->getZIndex();
                    pi = polar_group[i][t]->getPolarIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                    xp = getTrack2DByCycle(a, c, ci)
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[_num_l[a][p] + i][0]->getZIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCFwd());

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
                  xp = getTrack2DByCycle(a, c, ci)
                    ->getTrackPrdcFwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCBwd
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCBwd());

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
                  xp = getTrack2DByCycle(a, c, ci)
                    ->getTrackPrdcBwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                    xp = getTrack2DByCycle(a, c, ci)
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][p]]
                      [_tracks_per_train[a][c][p]
                       [i - _num_l[a][p]] - 1]->getZIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                    xp = getTrack2DByCycle(a, c, xi)
                      ->getTrackPrdcFwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][p]]
                      [_tracks_per_train[a][c][p]
                       [i - _num_l[a][p]] - 1]->getZIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURAFCE_Y_MAX */
            else{
              if (polar_group[i][t]->getCycleFwd()) {

                polar_group[i][t]->setBCBwd
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCBwd());

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
                  xp = getTrack2DByCycle(a, c, ci)
                    ->getTrackPrdcBwd()->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCFwd
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCFwd());

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
                  xp = getTrack2DByCycle(a, c, ci)
                    ->getTrackPrdcFwd()->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                    xp = getTrack2DByCycle(a, c, ci)
                      ->getTrackPrdcFwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][pc]][0]->getZIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                    xp = getTrack2DByCycle(a, c, ci)
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[i - _num_l[a][pc]][0]->getZIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, or SURFACE_Y_MAX */
            else{
              if (polar_group[i][t]->getCycleFwd()) {
                
                polar_group[i][t]->setBCFwd
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCFwd());

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
                  xp = getTrack2DByCycle(a, c, ci)
                    ->getTrackPrdcFwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCBwd
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCBwd());

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
                  xp = getTrack2DByCycle(a, c, ci)
                    ->getTrackPrdcBwd()->getXYIndex();
                  zp = polar_group[i][t+1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                    xp = getTrack2DByCycle(a, c, ci)
                      ->getTrackPrdcBwd()->getXYIndex();
                    zp = polar_group[_num_l[a][pc] + i]
                      [_tracks_per_train[a][c][p]
                       [_num_l[a][pc] + i] - 1]->getZIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()
                      ->getXYIndex();
                    zp = polar_group[_num_l[a][pc] + _num_z[a][pc] - i - 1]
                      [_tracks_per_train[a][c][p]
                       [_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]->getZIndex();
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
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
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCBwd());

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
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()
                    ->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                polar_group[i][t]->setBCFwd
                  (getTrack2DByCycle(a, c, polar_group[i][t]
                                     ->getCycleTrackIndex())->getBCFwd());

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
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()
                    ->getXYIndex();
                  zp = polar_group[i][t-1]->getZIndex();
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
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
void TrackGenerator::initialize3DTrackCycleIds() {

  int id = 0;
  Track* track;
  bool fwd;

  /* Set the periodic track cycle ids */
  if (_periodic) {
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

            track = &_tracks_3D_stack[a][i][p][z];

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

          track = &_tracks_3D_stack[a][i][p][z];
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
void TrackGenerator::decomposeLZTrack(Track3D* track, double l_start,
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

    track_2d = getTrack2DByCycle(azim, cycle, i);

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

    track_2d = getTrack2DByCycle(azim, cycle, i);

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
      track_2d = getTrack2DByCycle(azim, cycle, i);
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
       * index in the _tracks_3D_stack array since we want the tracks
       * in the _tracks_3D_stack array to all be pointing in the positive-y
       * direction whereas tracks in the _tracks_3D_stack array will point in
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
      track_3d = &_tracks_3D_stack[ai][ti][pi][zi];

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
    track_2d = getTrack2DByCycle(azim, cycle, i);
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
    _tracks_per_stack[ai][ti][pi] += 1;
  }
}


/**
 * @brief Recalibrates Track start and end points to the origin of the Geometry.
 * @details The origin of the Geometry is designated at its center by
 *          convention, but for track initialization the origin is assumed to be
 *          at the bottom right corner for simplicity. This method corrects
 *          for this by re-assigning the start and end Point coordinates.
 */
void TrackGenerator::recalibrate2DTracksToOrigin() {

  /* Recalibrate the tracks to the origin and set the uid. Note that the
   * loop structure is unconventional in order to preserve an increasing
   * track uid value in the Solver's tracks array. The tracks array is
   * oriented with this loop structure in order to maintain reproducability
   * for parallel runs.
   */

  /* Loop over azim reflective halfspaces */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

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
 * @brief Recalibrates Track start and end points to the origin of the Geometry.
 * @details The origin of the Geometry is designated at its center by
 *          convention, but for track initialization the origin is assumed to be
 *          at the bottom right corner for simplicity. This method corrects
 *          for this by re-assigning the start and end Point coordinates.
 */
void TrackGenerator::recalibrate3DTracksToOrigin() {

  /* Recalibrate the tracks to the origin and set the uid. Note that the
   * loop structure is unconventional in order to preserve an increasing
   * track uid value in the Solver's tracks array. The tracks array is
   * oriented with this loop structure in order to maintain reproducability
   * for parallel runs.
   */

  /* Loop over azim reflective halfspaces */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          double x0 = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          double y0 = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          double z0 = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          double x1 = _tracks_3D_stack[a][i][p][z].getEnd()->getX();
          double y1 = _tracks_3D_stack[a][i][p][z].getEnd()->getY();
          double z1 = _tracks_3D_stack[a][i][p][z].getEnd()->getZ();
          double new_x0 = x0 + _geometry->getMinX();
          double new_y0 = y0 + _geometry->getMinY();
          double new_z0 = z0 + _geometry->getMinZ();
          double new_x1 = x1 + _geometry->getMinX();
          double new_y1 = y1 + _geometry->getMinY();
          double new_z1 = z1 + _geometry->getMinZ();

          _tracks_3D_stack[a][i][p][z]
            .setCoords(new_x0, new_y0, new_z0, new_x1, new_y1, new_z1);
        }
      }
    }
  }

  /* Enusre that all tracks reside within the geometry */
  FP_PRECISION max_z = _geometry->getMaxZ();
  FP_PRECISION min_z = _geometry->getMinZ();
  for (int a=0; a < _num_azim/2; a++) {
    #pragma omp parallel for
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          FP_PRECISION start_z =
            _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          if (start_z > max_z)
            _tracks_3D_stack[a][i][p][z].getStart()->setZ(max_z);
          else if (start_z < min_z)
            _tracks_3D_stack[a][i][p][z].getStart()->setZ(min_z);
        }
      }
    }
  }
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize2D() {

  log_printf(NORMAL, "Ray tracing for 2D track segmentation...");

  int tracks_segmented = 0;
  int num_2D_tracks = getNum2DTracks();

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    log_printf(NORMAL, "segmenting 2D tracks - Percent complete: %5.2f %%",
               double(tracks_segmented) / num_2D_tracks * 100.0);
    #pragma omp parallel for
    for (int i=0; i < getNumX(a) + getNumY(a); i++)
      _geometry->segmentize2D(&_tracks_2D[a][i], _z_coord);

    tracks_segmented += getNumX(a) + getNumY(a);
  }

  _geometry->initializeFSRVectors();
  _contains_2D_segments = true;

  return;
}

/**
 * @brief Generates 2D segments for each extruded track across the Geometry,
 *        initializing axially extruded regions as well as 3D FSRs.
 */
void TrackGenerator::segmentizeExtruded() {

  log_printf(NORMAL, "Ray tracing for axially extruded track segmentation...");

  /* Allocate and initialize the flattened tracks array */
  _flattened_tracks = new Track*[_num_2D_tracks];
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      int id = _tracks_2D[a][i].getUid();
      _flattened_tracks[id] = &_tracks_2D[a][i];
    }
  }

  /* Get all unique z-coords at which 2D radial segementation is performed */
  std::vector<FP_PRECISION> z_coords;
  if (_contains_segmentation_heights)
    z_coords = _segmentation_heights;
  else
    z_coords = _geometry->getUniqueZPlanes();

  /* Loop over all extruded Tracks */
  #pragma omp parallel for
  for (int index=0; index < _num_2D_tracks; index++)
    _geometry->segmentizeExtruded(_flattened_tracks[index], z_coords);

  /* Initialize 3D FSRs and their associated vectors*/
  log_printf(NORMAL, "Initializing FSRs axially...");
  _geometry->initializeAxialFSRs(_global_z_mesh);
  _geometry->initializeFSRVectors();

  /* Count the number of segments in each track */
  countSegments();
  _contains_flattened_tracks = true;
  _contains_2D_segments = true;
  
  return;
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize3D() {

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
          _geometry->segmentize3D(&_tracks_3D_stack[a][i][p][z]);
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
double TrackGenerator::convertLtoX(double l, int azim, int cycle) {

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to X since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);

  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++) {
    if (l <= length_sum + getTrack2DByCycle(azim, cycle, i)->getLength()) {
      track_index = i;
      break;
    }
    else{
      length_sum += getTrack2DByCycle(azim, cycle, i)->getLength();
    }
  }

  if (l - length_sum < 0.0)
    log_printf(ERROR, "found negative length residual in converting l to x");

  double x1 = getTrack2DByCycle(azim, cycle, track_index)->getStart()->getX();
  double x2 = getTrack2DByCycle(azim, cycle, track_index)->getEnd()->getX();
  double l_rel = (l - length_sum) / getTrack2DByCycle(azim, cycle, track_index)
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
double TrackGenerator::convertLtoY(double l, int azim, int cycle) {

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to Y since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);

  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++) {
    if (l <= length_sum + getTrack2DByCycle(azim, cycle, i)->getLength()) {
      track_index = i;
      break;
    }
    else{
      length_sum += getTrack2DByCycle(azim, cycle, i)->getLength();
    }
  }

  if (l - length_sum < 0.0)
    log_printf(ERROR, "found negative length residual in converting l to y");

  double y1 = getTrack2DByCycle(azim, cycle, track_index)->getStart()->getY();
  double y2 = getTrack2DByCycle(azim, cycle, track_index)->getEnd()->getY();
  double l_rel = (l - length_sum) / getTrack2DByCycle(azim, cycle, track_index)
    ->getLength();

  double y;

  if (getCycleDirection(azim, cycle, track_index))
    y = l_rel * y2 + (1.0 - l_rel) * y1;
  else
    y = l_rel * y1 + (1.0 - l_rel) * y2;

  return y;
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
  if (!stat(directory.str().c_str(), &st) == 0)
    mkdir(directory.str().c_str(), S_IRWXU);

  if (_solve_3D) {

    /* Get the quadrature and track method types */
    std::string quad_type;
    std::string track_method;

    if (_quadrature->getQuadratureType() == EQUAL_WEIGHT)
      quad_type = "EQ_WGT";
    else if (_quadrature->getQuadratureType() == EQUAL_ANGLE)
      quad_type = "EQ_ANG";
    else
      log_printf(ERROR, "Unable to solve 3D problem with quadrature type"
                 " other than EQUAL_WEIGHT or EQUAL_ANGLE");

    if (_track_generation_method == GLOBAL_TRACKING)
      track_method = "GT";
    else if (_track_generation_method == MODULAR_RAY_TRACING)
      track_method = "MRT";
    else if (_track_generation_method == SIMPLIFIED_MODULAR_RAY_TRACING)
      track_method = "sMRT";

    if (_geometry->getCmfd() != NULL) {
      test_filename << directory.str() << "/3D_"
                    << _num_azim << "_azim_"
                    << _num_polar << "_polar_"
                    << _azim_spacing << "x" << _polar_spacing
                    << "_cm_spacing_cmfd_"
                    << _geometry->getCmfd()->getNumX()
                    << "x" << _geometry->getCmfd()->getNumY()
                    << "x" << _geometry->getCmfd()->getNumZ()
                    << "_quad_" << quad_type << "_track_" << track_method
                    << ".data";
    }
    else{
      test_filename << directory.str() << "/3D_"
                    << _num_azim << "_azim_"
                    << _num_polar << "_polar_"
                    << _azim_spacing << "x" << _polar_spacing
                    << "_cm_spacing_quad_"
                    << quad_type << "_track_" << track_method << ".data";
    }
  }
  else{
    if (_geometry->getCmfd() != NULL) {
      test_filename << directory.str() << "/2D_"
                    << _num_azim << "_azim_"
                    << _azim_spacing << "_cm_spacing_cmfd_"
                    << _geometry->getCmfd()->getNumX()
                    << "x" << _geometry->getCmfd()->getNumY()
                    << ".data";
    }
    else{
      test_filename << directory.str() << "/2D_"
                    << _num_azim << "_angles_"
                    << _azim_spacing << "_cm_spacing.data";
    }
  }

  _tracks_filename = test_filename.str();

  /* Check to see if a Track file exists for this geometry, number of azimuthal
   * angles, and track spacing, and if so, import the ray tracing data */
  if (!stat(_tracks_filename.c_str(), &buffer)) {
    if (_solve_3D && !_OTF){
      if (read3DSegmentsFromFile()) {
        _use_input_file = true;
        _contains_3D_segments = true;
      }
    }
    else{
      if (read2DSegmentsFromFile()) {
        _use_input_file = true;
        _contains_2D_segments = true;
      }
    }
  }
}


/**
 * @brief Writes all Track and segment data to a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 */
void TrackGenerator::dump2DSegmentsToFile() {

  /* Check whether the segments should be dumped */
  if (!_dump_segments)
    return;

  log_printf(NORMAL, "Dumping 2D segments to file...");

  if (!_contains_2D_segments)
    log_printf(ERROR, "Unable to dump 2D Segments to a file since no Segments "
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
  fwrite(&_z_coord, sizeof(double), 1, out);

  Track2D* curr_track;
  int num_segments;
  std::vector<segment*> _segments;
  Cmfd* cmfd = _geometry->getCmfd();

  segment* curr_segment;
  double length;
  int material_id;
  int region_id;
  int cmfd_surface_fwd;
  int cmfd_surface_bwd;

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

      /* Get data for this Track */
      curr_track = &_tracks_2D[a][i];
      num_segments = curr_track->getNumSegments();

      /* Write data for this Track to the Track file */
      fwrite(&num_segments, sizeof(int), 1, out);

      /* Loop over all segments for this Track */
      for (int s=0; s < num_segments; s++) {

        /* Get data for this segment */
        curr_segment = curr_track->getSegment(s);
        length = curr_segment->_length;
        material_id = curr_segment->_material->getId();
        region_id = curr_segment->_region_id;

        /* Write data for this segment to the Track file */
        fwrite(&length, sizeof(double), 1, out);
        fwrite(&material_id, sizeof(int), 1, out);
        fwrite(&region_id, sizeof(int), 1, out);

        /* Write CMFD-related data for the Track if needed */
        if (cmfd != NULL) {
          cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
          cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;
          fwrite(&cmfd_surface_fwd, sizeof(int), 1, out);
          fwrite(&cmfd_surface_bwd, sizeof(int), 1, out);
        }
      }
    }
  }

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
    std::vector< std::vector<int> > cell_fsrs = cmfd->getCellFSRs();
    std::vector<int>::iterator iter;
    int num_cells = cmfd->getNumCells();
    fwrite(&num_cells, sizeof(int), 1, out);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      num_FSRs = cell_fsrs.at(cell).size();
      fwrite(&num_FSRs, sizeof(int), 1, out);

      /* Loop over FSRs within cell */
      for (iter = cell_fsrs.at(cell).begin(); iter != cell_fsrs.at(cell).end();
           ++iter)
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
 * @brief Writes all Track and segment data to a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 */
void TrackGenerator::dump3DSegmentsToFile() {

  /* Check whether the segments should be dumped */
  if (!_dump_segments)
    return;

  log_printf(NORMAL, "Dumping 3D segments to file...");

  if (!_contains_3D_segments)
    log_printf(ERROR, "Unable to dump 3D Segments to a file since no Segments "
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

  Track3D* curr_track;
  int num_segments;
  std::vector<segment*> _segments;
  Cmfd* cmfd = _geometry->getCmfd();

  segment* curr_segment;
  double length;
  int material_id;
  int region_id;
  int cmfd_surface_fwd;
  int cmfd_surface_bwd;

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          /* Get data for this Track */
          curr_track = &_tracks_3D_stack[a][i][p][z];
          num_segments = curr_track->getNumSegments();

          /* Write data for this Track to the Track file */
          fwrite(&num_segments, sizeof(int), 1, out);

          /* Loop over all segments for this Track */
          for (int s=0; s < num_segments; s++) {

            /* Get data for this segment */
            curr_segment = curr_track->getSegment(s);
            length = curr_segment->_length;
            material_id = curr_segment->_material->getId();
            region_id = curr_segment->_region_id;

            /* Write data for this segment to the Track file */
            fwrite(&length, sizeof(double), 1, out);
            fwrite(&material_id, sizeof(int), 1, out);
            fwrite(&region_id, sizeof(int), 1, out);

            /* Write CMFD-related data for the Track if needed */
            if (cmfd != NULL) {
              cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
              cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;
              fwrite(&cmfd_surface_fwd, sizeof(int), 1, out);
              fwrite(&cmfd_surface_bwd, sizeof(int), 1, out);
            }
          }
        }
      }
    }
  }

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
    std::vector< std::vector<int> > cell_fsrs = cmfd->getCellFSRs();
    std::vector<int>::iterator iter;
    int num_cells = cmfd->getNumCells();
    fwrite(&num_cells, sizeof(int), 1, out);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      num_FSRs = cell_fsrs.at(cell).size();
      fwrite(&num_FSRs, sizeof(int), 1, out);

      /* Loop over FSRs within cell */
      for (iter = cell_fsrs.at(cell).begin(); iter != cell_fsrs.at(cell).end();
           ++iter)
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
bool TrackGenerator::read2DSegmentsFromFile() {

  int ret;
  FILE* in;
  in = fopen(_tracks_filename.c_str(), "r");

  int string_length;
  double z_coord;

  /* Import Geometry metadata from the Track file */
  ret = fread(&string_length, sizeof(int), 1, in);
  char* geometry_to_string = new char[string_length];
  ret = fread(geometry_to_string, sizeof(char)*string_length, 1, in);
  ret = fread(&z_coord, sizeof(double), 1, in);

  /* Check if our Geometry is exactly the same as the Geometry in the
   * Track file for this number of azimuthal angles and track spacing */
  if (_geometry->toString().compare(std::string(geometry_to_string)) != 0 ||
      _z_coord != z_coord)
    return false;

  delete [] geometry_to_string;

  log_printf(NORMAL, "Importing ray tracing data from file...");

  Track2D* curr_track;
  int num_segments;
  Cmfd* cmfd = _geometry->getCmfd();

  double length;
  int material_id;
  int region_id;

  int cmfd_surface_fwd;
  int cmfd_surface_bwd;

  std::map<int, Material*> materials = _geometry->getAllMaterials();

  int uid = 0;

  /* Loop over Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

      /* Import data for this Track from Track file */
      ret = fread(&num_segments, sizeof(int), 1, in);

      /* Initialize a Track with this data */
      curr_track = &_tracks_2D[a][i];

      /* Loop over all segments in this Track */
      for (int s=0; s < num_segments; s++) {

        /* Import data for this segment from Track file */
        ret = fread(&length, sizeof(double), 1, in);
        ret = fread(&material_id, sizeof(int), 1, in);
        ret = fread(&region_id, sizeof(int), 1, in);

        /* Initialize segment with the data */
        segment* curr_segment = new segment;
        curr_segment->_length = length;
        curr_segment->_material = materials[material_id];
        curr_segment->_region_id = region_id;

        /* Import CMFD-related data if needed */
        if (cmfd != NULL) {
          ret = fread(&cmfd_surface_fwd, sizeof(int), 1, in);
          ret = fread(&cmfd_surface_bwd, sizeof(int), 1, in);
          curr_segment->_cmfd_surface_fwd = cmfd_surface_fwd;
          curr_segment->_cmfd_surface_bwd = cmfd_surface_bwd;
        }

        /* Add this segment to the Track */
        curr_track->addSegment(curr_segment);
      }

      uid++;
    }
  }

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
  if (cmfd != NULL) {
    std::vector< std::vector<int> > cell_fsrs;
    int num_cells, fsr_id;
    ret = fread(&num_cells, sizeof(int), 1, in);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++) {
      std::vector<int> *fsrs = new std::vector<int>;
      cell_fsrs.push_back(*fsrs);
      ret = fread(&num_FSRs, sizeof(int), 1, in);

      /* Loop over FRSs within cell */
      for (int fsr = 0; fsr < num_FSRs; fsr++) {
        ret = fread(&fsr_id, sizeof(int), 1, in);
        cell_fsrs.at(cell).push_back(fsr_id);
      }
    }

    /* Set CMFD cell_fsrs vector of vectors */
    cmfd->setCellFSRs(cell_fsrs);
  }

  /* Close the Track file */
  fclose(in);

  return true;
}


/**
 * @brief Reads Tracks in from a "*.tracks" binary file.
 * @details Storing Tracks in a binary file saves time by eliminating ray
 *          tracing for Track segmentation in commonly simulated geometries.
 * @return true if able to read Tracks in from a file; false otherwise
 */
bool TrackGenerator::read3DSegmentsFromFile() {

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

  Track3D* curr_track;
  int num_segments;
  Cmfd* cmfd = _geometry->getCmfd();

  double length;
  int material_id;
  int region_id;

  int cmfd_surface_fwd;
  int cmfd_surface_bwd;

  std::map<int, Material*> materials = _geometry->getAllMaterials();

  int uid = 0;

  /* Loop over Tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

          /* Import data for this Track from Track file */
          ret = fread(&num_segments, sizeof(int), 1, in);

          /* Get data for this Track */
          curr_track = &_tracks_3D_stack[a][i][p][z];

          /* Loop over all segments in this Track */
          for (int s=0; s < num_segments; s++) {

            /* Import data for this segment from Track file */
            ret = fread(&length, sizeof(double), 1, in);
            ret = fread(&material_id, sizeof(int), 1, in);
            ret = fread(&region_id, sizeof(int), 1, in);

            /* Initialize segment with the data */
            segment* curr_segment = new segment;
            curr_segment->_length = length;
            curr_segment->_material = materials[material_id];
            curr_segment->_region_id = region_id;

            /* Import CMFD-related data if needed */
            if (cmfd != NULL) {
              ret = fread(&cmfd_surface_fwd, sizeof(int), 1, in);
              ret = fread(&cmfd_surface_bwd, sizeof(int), 1, in);
              curr_segment->_cmfd_surface_fwd = cmfd_surface_fwd;
              curr_segment->_cmfd_surface_bwd = cmfd_surface_bwd;
            }

            /* Add this segment to the Track */
            curr_track->addSegment(curr_segment);
          }

          uid++;
        }
      }
    }
  }

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
    cmfd->setCellFSRs(cell_fsrs);
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

  if (!_contains_2D_segments && !_contains_3D_segments)
    log_printf(ERROR, "Unable to split segments since "
	       "segments have not yet been generated");

  if (_OTF)
    log_printf(ERROR, "Segments cannot be split for on-the-fly ray tracing");

  int num_cuts, min_num_cuts;
  segment* curr_segment;

  FP_PRECISION length, tau;
  int fsr_id;
  Material* material;
  FP_PRECISION* sigma_t;
  int num_groups;

  if (_solve_3D) {
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
            for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments();
                s++) {
                
              /* Extract data from this segment to compute its optical 
               * length */
              curr_segment = _tracks_3D_stack[a][i][p][z].getSegment(s);
              material = curr_segment->_material;
              length = curr_segment->_length;
              fsr_id = curr_segment->_region_id;

              /* Compute number of segments to split this segment into */
              min_num_cuts = 1;
              num_groups = material->getNumEnergyGroups();
              sigma_t = material->getSigmaT();

              for (int g=0; g < num_groups; g++) {
                tau = length * sigma_t[g];
                num_cuts = ceil(tau / max_optical_length);
                min_num_cuts = std::max(num_cuts, min_num_cuts);
              }

              /* If the segment does not need subdivisions, go to next
               * segment */
              if (min_num_cuts == 1)
                continue;

              /* Record the CMFD surfaces */
              int cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
              int cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;

              /* Split the segment into sub-segments */
              for (int k=0; k < min_num_cuts; k++) {

                /* Create a new Track segment */
                segment* new_segment = new segment;
                new_segment->_material = material;
                new_segment->_length = length / FP_PRECISION(min_num_cuts);
                new_segment->_region_id = fsr_id;

                /* Assign CMFD surface boundaries */
                if (k == 0)
                  new_segment->_cmfd_surface_bwd = cmfd_surface_bwd;

                if (k == min_num_cuts-1)
                  new_segment->_cmfd_surface_fwd = cmfd_surface_fwd;

                /* Insert the new segment to the Track */
                _tracks_3D_stack[a][i][p][z].insertSegment(s+k+1,
                    new_segment);
              }

              /* Remove the original segment from the Track */
              _tracks_3D_stack[a][i][p][z].removeSegment(s);
            }
          }
        }
      }
    }
  }
  else{
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++) {

          /* Extract data from this segment to compute it optical length */
          curr_segment = _tracks_2D[a][i].getSegment(s);
          material = curr_segment->_material;
          length = curr_segment->_length;
          fsr_id = curr_segment->_region_id;

          /* Compute number of segments to split this segment into */
          min_num_cuts = 1;
          num_groups = material->getNumEnergyGroups();
          sigma_t = material->getSigmaT();

          for (int g=0; g < num_groups; g++) {
            tau = length * sigma_t[g];
            num_cuts = ceil(tau / max_optical_length);
            min_num_cuts = std::max(num_cuts, min_num_cuts);
          }

          /* If the segment does not need subdivisions, go to next segment */
          if (min_num_cuts == 1)
            continue;

          /* Split the segment into sub-segments */
          for (int k=0; k < min_num_cuts; k++) {

            /* Create a new Track segment */
            segment* new_segment = new segment;
            new_segment->_material = material;
            new_segment->_length = length / FP_PRECISION(min_num_cuts);
            new_segment->_region_id = fsr_id;

            /* Assign CMFD surface boundaries */
            if (k == 0)
              new_segment->_cmfd_surface_bwd =
                curr_segment->_cmfd_surface_bwd;

            if (k == min_num_cuts-1)
              new_segment->_cmfd_surface_fwd =
                curr_segment->_cmfd_surface_fwd;

            /* Insert the new segment to the Track */
            _tracks_2D[a][i].insertSegment(s+k+1, new_segment);
          }

          /* Remove the original segment from the Track */
          _tracks_2D[a][i].removeSegment(s);
        }
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

  int num_FSRs = _geometry->getNumFSRs();
  FP_PRECISION* FSR_volumes;

  /* Create temporary array of centroids and initialize to origin */
  Point** centroids = new Point*[num_FSRs];
  for (int r=0; r < num_FSRs; r++) {
    centroids[r] = new Point();
    centroids[r]->setCoords(0.0, 0.0, 0.0);
  }

  if (_solve_3D) {

    FSR_volumes = get3DFSRVolumes();

    /* Allocate array for 3D segments for OTF computation */
    segment* segments;
    if (_OTF)
      segments = new segment[_max_num_segments];

    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

            int num_segments = _tracks_3D_stack[a][i][p][z].getNumSegments();
            double xx = _tracks_3D_stack[a][i][p][z].getStart()->getX();
            double yy = _tracks_3D_stack[a][i][p][z].getStart()->getY();
            double zz = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
            double phi = _tracks_3D_stack[a][i][p][z].getPhi();
            double theta = _quadrature->getTheta(a, p);
            double wgt = _quadrature->getAzimWeight(a) *
              _quadrature->getPolarWeight(a, p) * getAzimSpacing(a)
              * getPolarSpacing(a,p);
            
            if (_OTF) {

              Point* start = _tracks_3D_stack[a][i][p][z].getStart();
              double theta = _tracks_3D_stack[a][i][p][z].getTheta();
              int ext_id = _tracks_2D[a][i].getUid();
              Track* flattened_track = _flattened_tracks[ext_id];

              SegmentationKernel kernel;
              kernel.setSegments(segments);
              kernel.setMaxVal(_max_optical_length);
              traceSegmentsOTF(flattened_track, start, theta, &kernel);
            }
            else
              segments = _tracks_3D_stack[a][i][p][z].getSegments();

            for (int s=0; s < num_segments; s++) {
              segment* curr_segment = &segments[s];
              int fsr = curr_segment->_region_id;
              double volume = FSR_volumes[fsr];
              centroids[fsr]->
                setX(centroids[fsr]->getX() + wgt *
                     (xx + cos(phi) * sin(theta) * curr_segment->_length / 2.0)
                     * curr_segment->_length / FSR_volumes[fsr]);

              centroids[fsr]->
                setY(centroids[fsr]->getY() + wgt *
                     (yy + sin(phi) * sin(theta) * curr_segment->_length / 2.0)
                     * curr_segment->_length / FSR_volumes[fsr]);

              centroids[fsr]->
                setZ(centroids[fsr]->getZ() + wgt *
                     (zz + cos(theta) * curr_segment->_length / 2.0) *
                     curr_segment->_length / FSR_volumes[fsr]);

              xx += cos(phi) * sin(theta) * curr_segment->_length;
              yy += sin(phi) * sin(theta) * curr_segment->_length;
              zz += cos(theta) * curr_segment->_length;
            }
          }
        }
      }
    }
    if (_OTF)
      delete[] segments;
  }
  else{

    FSR_volumes = get2DFSRVolumes();

    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {

        int num_segments = _tracks_2D[a][i].getNumSegments();
        segment* segments = _tracks_2D[a][i].getSegments();
        double x = _tracks_2D[a][i].getStart()->getX();
        double y = _tracks_2D[a][i].getStart()->getY();
        double phi = _tracks_2D[a][i].getPhi();
        double wgt = _quadrature->getAzimWeight(a) *
          getAzimSpacing(a);

        for (int s=0; s < num_segments; s++) {
          segment* curr_segment = &segments[s];
          int fsr = curr_segment->_region_id;
          double volume = FSR_volumes[fsr];
          centroids[fsr]->
            setX(centroids[fsr]->getX() + wgt *
                 (x + cos(phi) * curr_segment->_length / 2.0) *
                 curr_segment->_length / FSR_volumes[fsr]);

          centroids[fsr]->
            setY(centroids[fsr]->getY() + wgt *
                 (y + sin(phi) * curr_segment->_length / 2.0) *
                 curr_segment->_length / FSR_volumes[fsr]);

          x += cos(phi) * curr_segment->_length;
          y += sin(phi) * curr_segment->_length;
        }
      }
    }
  }

  /* Set the centroid for the FSR */
  for (int r=0; r < num_FSRs; r++) {
    _geometry->setFSRCentroid(r, centroids[r]);
  }

  /* Delete temporary array of centroids and FSR volumes */
  delete [] FSR_volumes;
}


/**
 * @brief Sets the track laydown method for generation of 3D Tracks
 * @details Options for the track laydown are GLOBAL_TRACKING, 
 *          MODULAR_RAY_TRACING, and SIMPLIFIED_MODULAR_RAY_TRACING
 * @param method The track laydown method
 */
void TrackGenerator::setTrackGenerationMethod(int method) {

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
int TrackGenerator::getTrackGenerationMethod() {
  return _track_generation_method;
}


/**
 * @brief Returns a pointer to the Track indexed by azimuthal index, cycle,
 *        and track index in the cycle
 * @details The 2D Track cycle indicated by the azimuthal angle and cycle
 *          number is traversed across track_index tracks, returning the Track
 *          at that position if valid
 * @param azim The azimuthal index
 * @param cycle The 2D cycle number
 * @param track_index The track index into the cycle
 * @return the matching Track, if found
 */
Track* TrackGenerator::getTrack2DByCycle(int azim, int cycle, int track_index) {

  azim = _quadrature->getFirstOctantAzim(azim);
  Track* track = &_tracks_2D[azim][cycle];
  Track* track_prev;
  bool fwd = true;

  for (int i=0; i < track_index; i++) {
    track_prev = track;

    if (fwd) {
      track = track_prev->getTrackReflFwd();
      fwd = track_prev->getReflFwdFwd();
    }
    else {
      track = track_prev->getTrackReflBwd();
      fwd = track_prev->getReflBwdFwd();
    }
  }

  if (track == NULL)
    log_printf(ERROR, "Could not find track 2d by cycle for %i, %i, %i",
               azim, cycle, track_index);

  return track;
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

  azim = _quadrature->getFirstOctantAzim(azim);
  Track* track = &_tracks_2D[azim][cycle];
  Track* track_prev;
  bool fwd = true;

  for (int i=0; i < track_index; i++) {
    track_prev = track;

    if (fwd) {
      track = track_prev->getTrackReflFwd();
      fwd = track_prev->getReflFwdFwd();
    }
    else {
      track = track_prev->getTrackReflBwd();
      fwd = track_prev->getReflBwdFwd();
    }
  }

  return fwd;
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
 * @brief Retrieves the max optical path length of 3D segments for use in
 *        on-the-fly computation
 * @return maximum optical path length
 */
FP_PRECISION TrackGenerator::retrieveMaxOpticalLength() {
  return _max_optical_length;
}


/**
 * @brief A function that searches for the index into a values mesh using a
 *        binary search.
 * @details A binary search is used to calculate the index into a mesh of where
 *          the value val resides. If a mesh boundary is hit, the upper region
 *          is selected for positive-z traversing rays and the lower region is
 *          selected for negative-z traversing rays.
 * @param values an array of monotonically increasing values
 * @param size the size of the values array
 * @param val the level to be searched for in the mesh
 * @param sign the direction of the ray in the z-direction
 */
int TrackGenerator::binarySearch(FP_PRECISION* values, int size,
                                 FP_PRECISION val, int sign) {

  /* Initialize indexes into the values array */
  int imin = 0;
  int imax = size-1;

  /* Check if val is outside the range */
  if (val < values[imin] or val > values[imax]) {
    log_printf(ERROR, "Value out of the mesh range in binary search");
    return -1;
  }

  /* Search for interval containing val */
  while (imax - imin > 1) {

    int imid = (imin + imax) / 2;

    if (val > values[imid])
      imin = imid;
    else if (val < values[imid])
      imax = imid;
    else {
      if (sign > 0)
        return imid;
      else
        return imid-1;
    }
  }
  return imin;
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
void TrackGenerator::retrieveSingle3DTrackCoords(double coords[6],
                                                 int track_id) {

  /* Find 3D track associated with track_id */
  for (int a=0; a < _num_azim/2; a++)
    for (int i=0; i < getNumX(a) + getNumY(a); i++)
      for (int p=0; p < _num_polar; p++)
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          if (_tracks_3D_stack[a][i][p][z].getUid() == track_id) {

            coords[0] = _tracks_3D_stack[a][i][p][z].getStart()->getX();
            coords[1] = _tracks_3D_stack[a][i][p][z].getStart()->getY();
            coords[2] = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
            coords[3] = _tracks_3D_stack[a][i][p][z].getEnd()->getX();
            coords[4] = _tracks_3D_stack[a][i][p][z].getEnd()->getY();
            coords[5] = _tracks_3D_stack[a][i][p][z].getEnd()->getZ();
            return;
          }
  log_printf(ERROR, "Unable to find a 3D track associated with the given track"
                    "ID during coordinate retrieval");
  return;
}


/**
 * @brief Computes and returns an array of volumes indexed by FSR for
          on-the-fly computation.
 * @details Segment lengths are computed on-the-fly and subsequently used to
            tally FSR volumes. Note: It is the function caller's responsibility
            to deallocate the memory reserved for the FSR volume array.
 * @return a pointer to the array of FSR volumes
 */
FP_PRECISION* TrackGenerator::get3DFSRVolumesOTF() {

  int num_FSRs = _geometry->getNumFSRs();
  FP_PRECISION *FSR_volumes = new FP_PRECISION[num_FSRs];
  memset(FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));

  VolumeKernel kernel;
  kernel.setBuffer(FSR_volumes);

  /* Calculate each FSR's "volume" by accumulating the total length of
   * all Track segments multiplied by the Track "widths" for each FSR.  */
  for (int ext_id=0; ext_id < _num_2D_tracks; ext_id++) {

    /* Extract indices of 3D tracks associated with the extruded track */
    Track* flattened_track = _flattened_tracks[ext_id];
    int a = flattened_track->getAzimIndex();
    int azim_index = _quadrature->getFirstOctantAzim(a);
    int i = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < _num_polar; p++) {

      /* Extract polar angle */
      int polar_index = _quadrature->getFirstOctantPolar(p);

      /* Calculate the weight of the track */
      FP_PRECISION weight = _quadrature->getAzimWeight(azim_index)
          * _quadrature->getPolarWeight(azim_index, polar_index)
          * getAzimSpacing(azim_index)
          * getPolarSpacing(azim_index, polar_index);

      kernel.setWeight(weight);

      /* Loop over z-stacked rays */
      for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

        /* Extract track and starting point */
        Track3D* curr_track = &_tracks_3D_stack[a][i][p][z];
        Point* start = curr_track->getStart();
        double theta = curr_track->getTheta();

        traceSegmentsOTF(flattened_track, start, theta, &kernel);

      }
    }
  }

  for (int i=0; i < num_FSRs; i++)
    if (FSR_volumes[i] == 0)
      log_printf(ERROR, "Zero volume calculated in an FSR region since no "
               "track traversed the FSR. Use a finer track laydown to ensure "
               "every FSR is traversed.");

  return FSR_volumes;
}


/**
 * @brief Computes 3D segment lengths for a given associated 2D Track with a
 *        starting point and an angle on-the-fly and stores the lengths in the
 *        kernel passed by the user.
 * @details Segment lengths are computed on-the-fly using 2D segment lengths
 *          stored in a 2D Track object and 1D meshes from the extruded
 *          FSRs. Note: before calling this funciton with a SegmentationKernel,
 *          the memory for the segments should be allocated and referenced by
 *          the kernel using the setSegments routine in the kernels.
 * @param flattened_track the 2D track associated with the 3D track for which
 *        3D segments are computed
 * @param start the starting coordinates of the 3D track
 * @param theta the polar angle of the 3D track
 * @param kernel An MOCKernel object to apply to the calculated 3D segments
 */
void TrackGenerator::traceSegmentsOTF(Track* flattened_track, Point* start, 
                                      double theta, MOCKernel* kernel) {

  /* Create unit vector */
  double phi = flattened_track->getPhi();
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  int sign = (cos_theta > 0) - (cos_theta < 0);

  /* Extract starting coordinates */
  double x_start_3D = start->getX();
  double x_start_2D = flattened_track->getStart()->getX();
  double z_coord = start->getZ();

  /* Find 2D distance from 2D edge to start of track */
  double start_dist_2D = (x_start_3D - x_start_2D) / cos(phi);

  /* Find starting 2D segment */
  int seg_start = 0;
  segment* segments_2D = flattened_track->getSegments();
  for (int s=0; s < flattened_track->getNumSegments(); s++) {

    /* Determine if start point of track is beyond current 2D segment */
    double seg_len_2D = segments_2D[s]._length;
    if (start_dist_2D > seg_len_2D) {
      start_dist_2D -= seg_len_2D;
      seg_start++;
    }
    else {
      break;
    }
  }

  /* Track current location in root universe */
  Cmfd* cmfd = _geometry->getCmfd();
  LocalCoords curr_coords(start->getX(), start->getY(), z_coord);
  curr_coords.setUniverse(_geometry->getRootUniverse());

  FP_PRECISION tiny_delta_x = sin_theta * cos(phi) * TINY_MOVE;
  FP_PRECISION tiny_delta_y = sin_theta * sin(phi) * TINY_MOVE;
  FP_PRECISION tiny_delta_z = cos_theta * TINY_MOVE;

  /* Extract the appropriate starting mesh */
  int num_fsrs;
  FP_PRECISION* axial_mesh;
  if (_contains_global_z_mesh) {
    num_fsrs = _global_z_mesh.size() - 1;
    axial_mesh = &_global_z_mesh[0];
  }
  else {
    int extruded_fsr_id = segments_2D[seg_start]._region_id;
    ExtrudedFSR* extruded_FSR = _geometry->getExtrudedFSR(extruded_fsr_id);
    num_fsrs = extruded_FSR->_num_fsrs;
    axial_mesh = extruded_FSR->_mesh;
  }

  /* Get the starting z index */
  int z_ind = binarySearch(axial_mesh, num_fsrs+1, z_coord, sign);

  /* Loop over 2D segments */
  bool first_segment = true;
  bool segments_complete = false;
  for (int s=seg_start; s < flattened_track->getNumSegments(); s++) {

    /* Extract extruded FSR */
    int extruded_fsr_id = segments_2D[s]._region_id;
    ExtrudedFSR* extruded_FSR = _geometry->getExtrudedFSR(extruded_fsr_id);

    /* Determine new mesh and z index */
    if (first_segment || _contains_global_z_mesh) {
      first_segment = false;
    }
    else {
      /* Determine the axial region */
      num_fsrs = extruded_FSR->_num_fsrs;
      axial_mesh = extruded_FSR->_mesh;
      z_ind = binarySearch(axial_mesh, num_fsrs+1, z_coord, sign);
    }

    /* Extract 2D segment length */
    double remaining_length_2D = segments_2D[s]._length - start_dist_2D;
    start_dist_2D = 0;

    /* Transport along the 2D segment until it is completed */
    while (remaining_length_2D > 0) {

      /* Calculate 3D distance to z intersection */
      double z_dist_3D;
      if (sign > 0)
        z_dist_3D = (axial_mesh[z_ind+1] - z_coord) / cos_theta;
      else
        z_dist_3D = (axial_mesh[z_ind] - z_coord) / cos_theta;

      /* Calculate 3D distance to end of segment */
      double seg_dist_3D = remaining_length_2D / sin_theta;

      /* Calcualte shortest distance to intersection */
      double dist_2D;
      double dist_3D;
      int z_move;
      if (z_dist_3D <= seg_dist_3D) {
        dist_2D = z_dist_3D * sin_theta;
        dist_3D = z_dist_3D;
        z_move = sign;
      }
      else {
        dist_2D = remaining_length_2D;
        dist_3D = seg_dist_3D;
        z_move = 0;
      }

      /* Calculate CMFD surface */
      int cmfd_surface_bwd = -1;
      int cmfd_surface_fwd = -1;
      if (cmfd != NULL && dist_3D > TINY_MOVE) {

        /* Determine if this is the first 3D segment handled for the flattened
           2D segment. If so, get the 2D cmfd surface. */
        if (segments_2D[s]._length - remaining_length_2D <= TINY_MOVE)
          cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;

        /* Determine if this is the last 3D segment handled for the flattened
           2D segment. If so, get the 2D cmfd surface. */
        double next_dist_3D = (remaining_length_2D - dist_2D) / sin_theta;
        if (z_move == 0 || next_dist_3D <= TINY_MOVE)
          cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;

        /* Adjust coordinats forward to find cell */
        curr_coords.adjustCoords(tiny_delta_x, tiny_delta_y, tiny_delta_z);
        int cmfd_cell = cmfd->findCmfdCell(&curr_coords);

        /* Adjust coordinates back to find the backwards surface */
        curr_coords.adjustCoords(-tiny_delta_x, -tiny_delta_y, -tiny_delta_z);
        cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, &curr_coords, 
            cmfd_surface_bwd);

        /* Move coordinates to end of segment */
        FP_PRECISION delta_x = sin_theta * cos(phi) * dist_3D;
        FP_PRECISION delta_y = sin_theta * sin(phi) * dist_3D;
        FP_PRECISION delta_z = cos_theta * dist_3D;
        curr_coords.adjustCoords(delta_x, delta_y, delta_z);

        /* Find forward surface */
        cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, &curr_coords,
            cmfd_surface_fwd);
      }

      /* Operate on segment */
      if (dist_3D > TINY_MOVE)
        kernel->execute(dist_3D, extruded_FSR->_materials[z_ind],
                        extruded_FSR->_fsr_ids[z_ind], cmfd_surface_fwd,
                        cmfd_surface_bwd);

      /* Shorten remaining 2D segment length and move axial level */
      remaining_length_2D -= dist_2D;
      z_coord += dist_3D * cos_theta;
      z_ind += z_move;

      /* Check if the track has crossed a Z boundary */
      if (z_ind < 0 or z_ind >= num_fsrs) {

        /* Reset z index */
        if (z_ind < 0)
          z_ind = 0;
        else
          z_ind = num_fsrs - 1;

        /* Mark the 2D segment as complete */
        segments_complete = true;
        break;
      }
    }

    /* Check if the track is completed due to an axial boundary */
    if (segments_complete)
      break;
  }
}


// TODO FIXME
/**
 * @brief Computes 3D segment lengths for a given associated 2D Track with a
 *        starting point and an angle on-the-fly and stores the lengths in the
 *        kernel passed by the user.
 * @details Segment lengths are computed on-the-fly using 2D segment lengths
 *          stored in a 2D Track object and 1D meshes from the extruded
 *          FSRs. Note: before calling this funciton with a SegmentationKernel,
 *          the memory for the segments should be allocated and referenced by
 *          the kernel using the setSegments routine in the kernels.
 * @param flattened_track the 2D track associated with the 3D track for which
 *        3D segments are computed
 * @param start the starting coordinates of the 3D track
 * @param theta the polar angle of the 3D track
 * @param kernel An MOCKernel object to apply to the calculated 3D segments
 */
// TODO FIXME
void TrackGenerator::traceStackOTF(Track* flattened_track, int polar_index,
                                    MOCKernel** kernels) {

  /* Extract information about the z-stack */
  int azim_index = flattened_track->getAzimIndex();
  int ai = _quadrature->getFirstOctantAzim(azim_index);
  int pi = _quadrature->getFirstOctantPolar(polar_index);
  double z_spacing = _dz_eff[ai][pi];
  double theta = _quadrature->getTheta(azim_index, polar_index);
  int track_index = flattened_track->getXYIndex();

  /* Create unit vector */
  double phi = flattened_track->getPhi();
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double tan_theta = sin_theta / cos_theta;
  int sign = (cos_theta > 0) - (cos_theta < 0);

  /* Find 2D distance from 2D edge to start of track */
  Track3D* first = &_tracks_3D_stack[azim_index][track_index][polar_index][0];
  double x_start_3D = first->getStart()->getX();
  double x_start_2D = flattened_track->getStart()->getX();
  double start_dist_2D = (x_start_3D - x_start_2D) / cos(phi);

  /* Calculate starting intersection of lowest track with z-axis */
  double z0 = first->getStart()->getZ();
  double start_z = z0 - start_dist_2D / tan_theta;

  /* Track current location in root universe FIXME */
  Cmfd* cmfd = _geometry->getCmfd();

  /* Extract the appropriate starting mesh */
  int num_fsrs;
  FP_PRECISION* axial_mesh;
  if (_contains_global_z_mesh) {
    num_fsrs = _global_z_mesh.size() - 1;
    axial_mesh = &_global_z_mesh[0];
  }

  /* Loop over 2D segments */
  double first_start_z = start_z;
  segment* segments_2D = flattened_track->getSegments();
  for (int s=0; s < flattened_track->getNumSegments(); s++) {
  
    /* Get segment length and extruded FSR */
    FP_PRECISION seg_length_2D = segments_2D[s]._length;
    int extruded_fsr_id = segments_2D[s]._region_id;
    ExtrudedFSR* extruded_FSR = _geometry->getExtrudedFSR(extruded_fsr_id);

    /* Determine new mesh and z index */
    if (!_contains_global_z_mesh) {
      num_fsrs = extruded_FSR->_num_fsrs;
      axial_mesh = extruded_FSR->_mesh;
    }

    /* Calculate the end z coordinate of the first track */
    double first_end_z = first_start_z + seg_length_2D / tan_theta;

    /* Find the upper and lower z coordinates of the first track */
    double first_track_lower_z;
    double first_track_upper_z;
    if (sign > 0) {
      first_track_lower_z = first_start_z;
      first_track_upper_z = first_end_z;
    }
    else {
      first_track_lower_z = first_end_z;
      first_track_upper_z = first_start_z;
    }
    
    /* Loop over all 3D FSRs in the Extruded FSR to find intersections */
    for (int z_ind = 0; z_ind < num_fsrs; z_ind++) {
  
      /* Extract the FSR ID and Material ID of this 3D FSR */
      int fsr_id = extruded_FSR->_fsr_ids[z_ind];
      int material_id = extruded_FSR->_fsr_ids[z_ind];

      /* Get boundaries of the current mesh cell */
      double z_min = axial_mesh[z_ind];
      double z_max = axial_mesh[z_ind+1];

      /* Calculate z-stack track indexes that cross the 3D FSR */
      int start_track = (z_min - first_track_upper_z) / z_spacing + 1;
      int start_full = (z_min - first_track_lower_z) / z_spacing + 1;
      int end_full = (z_max - first_track_upper_z) / z_spacing + 1;
      int end_track = (z_max - first_track_lower_z) / z_spacing + 1;
      
      /* Check track bounds */
      start_track = std::max(start_track, 0);
      end_track = std::min(end_track, num_z_stack);
      
      /* Treat lower tracks that do not cross the entire 2D length */
      int min_lower = std::min(start_full, end_full);
      for (int i = start_track; i < min_lower; i++) {

        /* Calculate distance traveled in 3D FSR */
        double end_z = first_track_upper_z + i * z_spacing;
        double seg_len_3D = (end_z - z_min) / std::abs(cos_theta);

        /* Operate on segment */
        kernels[i]->execute(seg_len_3D, material_id, fsr_id, -1, -1);
      }

      /* Treat tracks that do cross the entire 2D length */
      for (int i = start_full; i < end_full; i++) {

        /* Calculate distance traveled in 3D FSR */
        double seg_len_3D = seg_length_2D / sin_theta;

        /* Operate on segment */
        kernels[i]->execute(seg_len_3D, material_id, fsr_id, -1, -1);
      }

      /* Treat tracks that cross through both the upper and lower boundary
         NOTE: this loop will only execute if there are no tracks that cross
         the entire 2D length in the FSR */
      for (int i = end_full; i < start_full; i++) {
        
        /* Calculate distance traveled in 3D FSR */
        double seg_len_3D = (z_max - z_min) / std::abs(cos_theta);

        /* Operate on segment */
        kernels[i]->execute(seg_len_3D, material_id, fsr_id, -1, -1);
      }

      /* Treat upper tracks that do not cross the entire 2D length */
      int min_upper = std::max(start_full, end_full);
      for (int i = min_upper; i < end_track; i++) {

        /* Calculate distance traveled in 3D FSR */
        double start_z = first_track_lower_z + i * z_spacing;
        double seg_len_3D = (z_max - start_z) / std::abs(cos_theta);

        /* Operate on segment */
        kernels[i]->execute(seg_len_3D, material_id, fsr_id, -1, -1);
      }
    }
    /* Traverse segment on first track */
    first_start_z = first_end_z;
  }
}


/**
 * @brief Counts the number of 3D segments in the Geomtry
 * @details All 3D segments are computed on-the-fly subject to the max optical
 *          path length to determine the number of 3D segments in the Geometry
 */
void TrackGenerator::countSegments() {

  /* Calculate each FSR's "volume" by accumulating the total length of
   * all Track segments multiplied by the Track "widths" for each FSR.  */
  #pragma omp parallel for
  for (int ext_id=0; ext_id < _num_2D_tracks; ext_id++) {

    CounterKernel counter;
    counter.setMaxVal(_max_optical_length);

    /* Extract indices of 3D tracks associated with the extruded track */
    Track* flattened_track = _flattened_tracks[ext_id];
    int a = flattened_track->getAzimIndex();
    int azim_index = _quadrature->getFirstOctantAzim(a);
    int i = flattened_track->getXYIndex();

    /* Loop over polar angles */
    for (int p=0; p < _num_polar; p++) {

      /* Extract polar angle */
      int polar_index = _quadrature->getFirstOctantPolar(p);

      /* Loop over z-stacked rays */
      for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

        /* Extract track and starting point */
        Track3D* curr_track = &_tracks_3D_stack[a][i][p][z];
        double theta = curr_track->getTheta();
        Point* start = curr_track->getStart();

        /* Trace 3D segments */
        traceSegmentsOTF(flattened_track, start, theta, &counter);

        /* Set the number of segments for the track */
        int num_segments = counter.getCount();
        curr_track->setNumSegments(num_segments);
        if (num_segments > _max_num_segments)
          _max_num_segments = num_segments;
        counter.resetCount();
      }
    }
  }
}


/**
 * @brief Sets the track periodic indices of all 2D Tracks
 * @details Periodic cylces are traversed until all 2D Tracks are visited and
 *          their periodic indices are set
 */
void TrackGenerator::initialize2DTrackPeriodicIndices() {

  log_printf(NORMAL, "Initializing track periodic indices...");

  if (!_periodic)
    return;

  Track* track;
  int track_index;

  /* Set the track periodic cycle indices for 2D tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

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
 * @brief Sets the track periodic indices of all 3D Tracks
 * @details Periodic cylces are traversed until all 3D Tracks are visited and
 *          their periodic indices are set
 */
void TrackGenerator::initialize3DTrackPeriodicIndices() {

  if (!_periodic)
    return;

  Track* track;
  int track_index;

  /* Set the track periodic cycle indices for 3D tracks */
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      for (int p=0; p < _num_polar; p++) {
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
          track = &_tracks_3D_stack[a][i][p][z];

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
 * @brief Creates a Track array by increasing uid 
 * @details An array is created which indexes Tracks by increasing uid.
 *          Parallel groups are also initialized -- groups of Tracks that can
 *          be computed in parallel without the potential of overwriting 
 *          angular fluxes of connecting tracks prematurely.
 */
void TrackGenerator::initializeTracksArray() {

  log_printf(NORMAL, "Initializing tracks array...");

  Track* track;
  int uid = 0;
  int num_tracks;
  int azim_group_id, periodic_group_id, polar_group_id;
  int track_azim_group_id, track_periodic_group_id, track_polar_group_id;
  int track_periodic_index;

  /* Set the number of parallel track groups */
  if (_solve_3D) {
    if (_periodic)
      _num_parallel_track_groups = 12;
    else
      _num_parallel_track_groups = 4;
  }
  else {
    if (_periodic)
      _num_parallel_track_groups = 6;
    else
      _num_parallel_track_groups = 2;
  }

  /* Create the array of track ids separating the parallel groups */
  _num_tracks_by_parallel_group = new int[_num_parallel_track_groups + 1];

  /* Set the first index in the num tracks by parallel group array to 0 */
  _num_tracks_by_parallel_group[0] = 0;

  /* Set the number of tracks to allocate pointers for in tracks array */
  if (_solve_3D)
    num_tracks = getNum3DTracks();
  else
    num_tracks = getNum2DTracks();

  /* Allocate memory for tracks array */
  _tracks = new Track*[num_tracks];

  /* Recalibrate the tracks to the origin and set the uid. Note that the
   * loop structure is unconventional in order to preserve a monotonically
   * increasing track uid value in the Solver's tracks array. The tracks array
   * is oriented such the tracks can be broken up into 4 sub arrays that are
   * guaranteed to contain tracks that do not transport into other tracks both
   * reflectively and periodically. This is done to guarantee reproducability
   * in parallel runs. */
  if (!_solve_3D || _OTF) {
    for (int g = 0; g < _num_parallel_track_groups; g++) {

      /* Set the azimuthal and periodic group ids */
      azim_group_id = g % 2;
      periodic_group_id = g / 2;

      /* Loop over all 2D tracks */
      for (int a = 0; a < _num_azim / 2; a++) {
        for (int i=0; i < getNumX(a) + getNumY(a); i++) {

          /* Get current track and azim group ids */
          track = &_tracks_2D[a][i];

          /* Get the track azim group id */
          track_azim_group_id = a / (_num_azim / 4);

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

          /* Check if track has current azim_group_id and periodic_group_id */
          if (azim_group_id == track_azim_group_id &&
              periodic_group_id == track_periodic_group_id) {
            track->setUid(uid);
            _tracks[uid] = track;
            uid++;
          }
        }
      }

      /* Set the track index boundary for this parallel group */
      _num_tracks_by_parallel_group[g + 1] = uid;
    }
  }
  if(_solve_3D) {

    /* Reset UID in case 2D tracks were intiialized */
    uid = 0;

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
              track = &_tracks_3D_stack[a][i][p][z];

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
                track->setUid(uid);
                _tracks[uid] = track;
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
}


/**
 * @brief Return an array with the Track uid separating the azimuthal and
 * periodic halfspaces
 * @return array with the Track uid separating the azimuthal and periodic
 *         halfspaces
 */
int* TrackGenerator::getNumTracksByParallelGroupArray() {
  if (!_contains_2D_tracks && !_contains_3D_tracks)
    log_printf(ERROR, "Unable to return the array with the Track uid "
               "separating the azimuthal and periodic halspaces since "
               "Tracks have not yet been generated.");

  return _num_tracks_by_parallel_group;
}


/**
 * @brief Returns the number of Track groups which can be executed in parallel
 *        without conflicts
 * @return the number of parallel grooups
 */
int TrackGenerator::getNumParallelTrackGroups() {
  return _num_parallel_track_groups;
}


/**
 * @brief returns whether periodic boundaries are present in Track generation
 * @return a boolean value - true if periodic; false otherwise
 */
bool TrackGenerator::getPeriodic() {
  return _periodic;
}


/**
 * @brief allocates memory for 3D Tracks
 * @details Before calling this function, the number of tracks per z-stack
 *          should be known and initialized in the _tracks_per_stack 3D array
 */
void TrackGenerator::create3DTracksArrays() {

  _tracks_3D_stack = new Track3D***[_num_azim/2];
  for (int a=0; a < _num_azim/2; a++) {
    _tracks_3D_stack[a] = new Track3D**[getNumX(a) + getNumY(a)];
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {
      _tracks_3D_stack[a][i] = new Track3D*[_num_polar];
      for (int p=0; p < _num_polar; p++) {
        _tracks_3D_stack[a][i][p] = new Track3D[_tracks_per_stack[a][i][p]];
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
 * @brief Sets a flag to record all segment information in the tracking file
 * @param A boolean value to determine whether or not to record segment
 *        information in the tracking file: true to record, false not to record
 */
void TrackGenerator::setDumpSegments(bool dump_segments) {
  _dump_segments = dump_segments;
}
