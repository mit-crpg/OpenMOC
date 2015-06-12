#include "TrackGenerator.h"


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
  _num_2D_tracks = 0;
  _num_3D_tracks = 0;
  _num_2D_segments = 0;
  _num_3D_segments = 0;
  _contains_2D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_tracks = false;
  _contains_3D_segments = false;
  _quadrature = NULL;
  _z_level = 0.0;
  _solve_3D = true;
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator::~TrackGenerator() {

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_3D_tracks){    
    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++){
        for (int p=0; p < _num_polar; p++){
          delete [] _tracks_per_plane[a][c][p];
          for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
            delete [] _tracks_3D[a][c][p][i];
          }
        }
      }
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

    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++)
        delete [] _tracks_2D[a][c];
      delete [] _tracks_2D[a];
    }

    delete [] _num_x;
    delete [] _num_y;
    delete [] _dx_eff;
    delete [] _dy_eff;
    delete [] _azim_spacings;
    delete [] _cycles_per_azim;
    delete [] _tracks_per_cycle;
    delete [] _cycle_length;
  }
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
 * @details This will return the user-specified track spacing and NOT the
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
 * @brief Return the total number of Tracks across the Geometry.
 * @return the total number of Tracks
 */
int TrackGenerator::getNum2DTracks() {
  return _num_2D_tracks;
}


/**
 * @brief Return the total number of Tracks across the Geometry.
 * @return the total number of Tracks
 */
int TrackGenerator::getNum3DTracks() {
  return _num_3D_tracks;
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator::getNum2DSegments() {
  return _num_2D_segments;
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator::getNum3DSegments() {
  return _num_3D_segments;
}


/**
 * @brief Returns a 3D jagged array of the Tracks.
 * @details The first index into the array is the azimuthal angle and the
 *          second index is polar angle and the third index is the 
 *          Track number.
 * @return the 3D jagged array of Tracks
 */
Track2D*** TrackGenerator::get2DTracks() {

  if (!_contains_2D_tracks)
    log_printf(ERROR, "Unable to return the 3D ragged array of the 2D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_2D;
}


Track3D***** TrackGenerator::get3DTracks() {

  if (!_contains_3D_tracks)
    log_printf(ERROR, "Unable to return the 3D ragged array of the 3D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_3D;
}


double* TrackGenerator::getAzimSpacings(){
  return _azim_spacings;
}


double TrackGenerator::getAzimSpacing(int azim){
  azim = _quadrature->getFirstOctantAzim(azim);
  return _azim_spacings[azim];
}


double** TrackGenerator::getPolarSpacings(){
  return _polar_spacings;
}


double TrackGenerator::getPolarSpacing(int azim, int polar){
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

  if (_solve_3D){
    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++){
        for (int p=0; p < _num_polar; p++){
          for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
            for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
              for (int s=0; s < _tracks_3D[a][c][p][i][t].getNumSegments(); s++){
                curr_segment = _tracks_3D[a][c][p][i][t].getSegment(s);
                length = curr_segment->_length;
                material = curr_segment->_material;
                sigma_t = material->getSigmaT();
                
                for (int e=0; e < material->getNumEnergyGroups(); e++)
                  max_optical_length = std::max(max_optical_length, length*sigma_t[e]);
              }
            }
          }
        }
      }
    }
  }
  else{
    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++){
        for (int t=0; t < _tracks_per_cycle[a]; t++){
          for (int s=0; s < _tracks_2D[a][c][t].getNumSegments(); s++){
            curr_segment = _tracks_2D[a][c][t].getSegment(s);
            length = curr_segment->_length;
            material = curr_segment->_material;
            sigma_t = material->getSigmaT();
            
            for (int e=0; e < material->getNumEnergyGroups(); e++)
              max_optical_length = std::max(max_optical_length, length*sigma_t[e]);
          }
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


int* TrackGenerator::getTracksPerCycle(){
  return _tracks_per_cycle;
}


int**** TrackGenerator::getTracksPerPlane(){
  return _tracks_per_plane;
}


int* TrackGenerator::getCyclesPerAzim(){
  return _cycles_per_azim;
}


double TrackGenerator::getCycleLength(int azim){
  azim = _quadrature->getFirstOctantAzim(azim);
  return _cycle_length[azim];
}


int TrackGenerator::getNumX(int azim){
  azim = _quadrature->getFirstOctantAzim(azim);
  return _num_x[azim];
}


int TrackGenerator::getNumY(int azim){
  azim = _quadrature->getFirstOctantAzim(azim);
  return _num_y[azim];
}


int TrackGenerator::getNumZ(int azim, int polar){
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
  return _num_z[azim][polar];
}


int TrackGenerator::getNumL(int azim, int polar){
  azim = _quadrature->getFirstOctantAzim(azim);
  polar = _quadrature->getFirstOctantPolar(polar);
  return _num_l[azim][polar];
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
   * all Track segments multipled by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int t=0; t < _tracks_per_cycle[a]; t++){
        for (int s=0; s < _tracks_2D[a][c][t].getNumSegments(); s++){
          segment = _tracks_2D[a][c][t].getSegment(s);
          volume = segment->_length * _quadrature->getAzimWeight(a)
            * getAzimSpacing(a);
          FSR_volumes[segment->_region_id] += volume;
        }
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

  if (!contains3DSegments())
    log_printf(ERROR, "Unable to get the FSR volumes since 3D tracks "
               "have not yet been generated");

  int num_FSRs = _geometry->getNumFSRs();
  FP_PRECISION *FSR_volumes = new FP_PRECISION[num_FSRs];
  memset(FSR_volumes, 0., num_FSRs*sizeof(FP_PRECISION));

  segment* segment;
  FP_PRECISION volume;

  /* Calculate each FSR's "volume" by accumulating the total length of * 
   * all Track segments multipled by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
            for (int s=0; s < _tracks_3D[a][c][p][i][t].getNumSegments(); s++){
              segment = _tracks_3D[a][c][p][i][t].getSegment(s);
              volume = segment->_length * _quadrature->getAzimWeight(a)
                * _quadrature->getPolarWeight(a, p) * getAzimSpacing(a)
                * getPolarSpacing(a,p);
              FSR_volumes[segment->_region_id] += volume;
            }
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
   * all Track segments multipled by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int t=0; t < _tracks_per_cycle[a]; t++){
        for (int s=0; s < _tracks_2D[a][c][t].getNumSegments(); s++){
          segment = _tracks_2D[a][c][t].getSegment(s);
          if (segment->_region_id == fsr_id)
            volume += segment->_length * _quadrature->getAzimWeight(a)
              * getAzimSpacing(a);
        }
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
   * all Track segments multipled by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
            for (int s=0; s < _tracks_3D[a][c][p][i][t].getNumSegments(); s++){
              segment = _tracks_3D[a][c][p][i][t].getSegment(s);
              if (segment->_region_id == fsr_id)
                volume += segment->_length * _quadrature->getAzimWeight(a)
                * _quadrature->getPolarWeight(a, p) * getAzimSpacing(a)
                * getPolarSpacing(a,p);
            }
          }
        }
      }
    }
  }

  return volume;
}


double TrackGenerator::getZLevel(){
  return _z_level;
}


Quadrature* TrackGenerator::getQuadrature(){
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
    log_printf(ERROR, "Unable to set the number of polar angles to %d for "
               "the TrackGenerator since it is not a multiple of 2", num_polar);

  _num_polar = num_polar;
  _contains_2D_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
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
  _num_2D_tracks = 0;
  _num_3D_tracks = 0;
  _num_2D_segments = 0;
  _num_3D_segments = 0;
  _contains_2D_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
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
  _num_2D_tracks = 0;
  _num_3D_tracks = 0;
  _num_2D_segments = 0;
  _num_3D_segments = 0;
  _contains_2D_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
}


/**
 * @brief Set a pointer to the Geometry to use for track generation.
 * @param geometry a pointer to the Geometry
 */
void TrackGenerator::setGeometry(Geometry* geometry) {
  _geometry = geometry;
  _num_2D_tracks = 0;
  _num_3D_tracks = 0;
  _num_2D_segments = 0;
  _num_3D_segments = 0;
  _contains_2D_tracks = false;
  _contains_3D_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
}


void TrackGenerator::setSolve2D(){
  _solve_3D = false;
}


void TrackGenerator::setSolve3D(){
  _solve_3D = true;
}


void TrackGenerator::setZLevel(double z_level){
  _z_level = z_level;
}


void TrackGenerator::setQuadrature(Quadrature* quadrature){
  _quadrature = quadrature;
}


/**
 * @brief Returns whether or not the TrackGenerator contains Track that are
 *        for its current number of azimuthal angles, track spacing and
 *        geometry.
 * @return true if the TrackGenerator conatains Tracks; false otherwise
 */
bool TrackGenerator::contains2DTracks() {
  return _contains_2D_tracks;
}


bool TrackGenerator::contains3DTracks() {
  return _contains_3D_tracks;
}


bool TrackGenerator::contains2DSegments() {
  return _contains_2D_segments;
}


bool TrackGenerator::contains3DSegments() {
  return _contains_3D_segments;
}


/**
 * @brief Fills an array with the x,y coordinates for each Track.
 * @details This class method is intended to be called by the OpenMOC
 *          Python "plotter" module as a utility to assist in plotting
 *          tracks. Although this method appears to require two arguments,
 *          in reality it only requires on due to SWIG and would be called
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
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int t = 0; t < _tracks_per_cycle[a]; t++){
        coords[counter]   = _tracks_2D[a][c][t].getStart()->getX();
        coords[counter+1] = _tracks_2D[a][c][t].getStart()->getY();
        coords[counter+2] = _tracks_2D[a][c][t].getEnd()->getX();
        coords[counter+3] = _tracks_2D[a][c][t].getEnd()->getY();
        
        counter += 4;
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
 *          in reality it only requires on due to SWIG and would be called
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
void TrackGenerator::retrieve3DTrackCoords(double* coords, int num_tracks) {

  if (num_tracks != 6*getNum3DTracks())
    log_printf(ERROR, "Unable to retrieve the Track coordinates since the "
               "TrackGenerator contains %d Tracks with %d coordinates but an "
               "array of length %d was input",
               getNum3DTracks(), 6*getNum3DTracks(), num_tracks);
  
  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
            coords[counter]   = _tracks_3D[a][c][p][i][t].getStart()->getX();
            coords[counter+1] = _tracks_3D[a][c][p][i][t].getStart()->getY();
            coords[counter+2] = _tracks_3D[a][c][p][i][t].getStart()->getZ();
            coords[counter+3] = _tracks_3D[a][c][p][i][t].getEnd()->getX();
            coords[counter+4] = _tracks_3D[a][c][p][i][t].getEnd()->getY();
            coords[counter+5] = _tracks_3D[a][c][p][i][t].getEnd()->getZ();            
            counter += 6;
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
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
              
            x0    = _tracks_3D[a][c][p][i][t].getStart()->getX();
            y0    = _tracks_3D[a][c][p][i][t].getStart()->getY();
            z0    = _tracks_3D[a][c][p][i][t].getStart()->getZ();
            phi   = _tracks_3D[a][c][p][i][t].getPhi();
            theta = _tracks_3D[a][c][p][i][t].getTheta();
            
            segments = _tracks_3D[a][c][p][i][t].getSegments();
            
            for (int s=0; s < _tracks_3D[a][c][p][i][t].getNumSegments(); s++) {
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
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int t = 0; t < _tracks_per_cycle[a]; t++){
              
        x0    = _tracks_2D[a][c][t].getStart()->getX();
        y0    = _tracks_2D[a][c][t].getStart()->getY();
        phi   = _tracks_2D[a][c][t].getPhi();
        
        segments = _tracks_2D[a][c][t].getSegments();
        
        for (int s=0; s < _tracks_2D[a][c][t].getNumSegments(); s++) {
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
  }
  
  return;
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
    log_printf(ERROR, "The total height, width, and depth of the Geometry must be "
               "nonzero for Track generation. Create a CellFill which "
               "is filled by the entire geometry and bounded by XPlanes, "
               "YPlanes, and ZPlanes to enable the Geometry to determine the total "
               "width, height, and depth of the model.");

  /* Generate Tracks, perform ray tracing across the geometry, and store
   * the data to a Track file */
  try {

    /* Create default quadrature set if user one has not been set */
    if (_quadrature == NULL){
      if (_solve_3D)
        _quadrature = new EqualWeightPolarQuad();
      else
        _quadrature = new TYPolarQuad();
    }
    
    _quadrature->setNumPolarAngles(_num_polar);
    _quadrature->setNumAzimAngles(_num_azim);
    
    /* Initialize the quadrature set */
    _quadrature->initialize();

    /* Initialize the 2D tracks */
    initialize2DTracks();
    initialize2DTrackReflections();

    /* If 3D problem, initialize the 3D tracks */
    if (_solve_3D){
      initialize3DTracks();
      initialize3DTrackReflections();
    }
    
    /* Recalibrate the 2D tracks back to the geometry origin */
    recalibrate2DTracksToOrigin();
    
    /* If 3D problem, recalibrate the 3D tracks back to the geometry origin */
    if (_solve_3D)
      recalibrate3DTracksToOrigin();

    initializeTrackFileDirectory();

    if (_use_input_file == false){

      /* Segmentize the tracks */
      if (_solve_3D){
        segmentize3D();
        dump3DSegmentsToFile();
      }
      else{
        segmentize2D();
        dump2DSegmentsToFile();
      }
    }
    
    /* Precompute the quadrature weights */
    _quadrature->precomputeWeights();
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to allocate memory needed to generate "
               "Tracks. Backtrace:\n%s", e.what());
  }

  return;
}


double TrackGenerator::leastCommonMultiple(double a, double b){

  bool _found = false;
  int lcm_a = 1;
  int lcm_b;
  double residual;

  /* For efficiency, make a the longer length */
  if (a < b){
    double a_temp = a;
    a = b;
    b = a_temp;
  }
  
  while (!_found){
    
    lcm_b = (int) round((lcm_a * a) / b);
    residual = fabs(lcm_a * a - lcm_b * b);

    if (residual < LCM_TOLERANCE)
      _found = true;
    else
      lcm_a++;
  }

  return lcm_a * a;
}


bool TrackGenerator::isSolve2D(){
  return !_solve_3D;
}


bool TrackGenerator::isSolve3D(){
  return _solve_3D;
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
  _tracks_2D        = new Track2D**[_num_azim/4];
  _num_x            = new int[_num_azim/4];
  _num_y            = new int[_num_azim/4];
  _cycle_length     = new double[_num_azim/4];
  _azim_spacings    = new double[_num_azim/4];

  double x1, x2, y1, y2;
  double phi;
  double width  = _geometry->getWidth();
  double height = _geometry->getHeight();
  
  /* Determine angular quadrature and track spacing */
  for (int i = 0; i < _num_azim/4; i++) {

    /* Get the desired azimuthal angle */
    phi = _quadrature->getPhi(i);

    /* The number of intersections with x,y-axes */
    _num_x[i] = (int) (fabs(width / _azim_spacing * sin(phi))) + 1;
    _num_y[i] = (int) (fabs(height / _azim_spacing * cos(phi))) + 1;

    /* Effective/actual angle (not the angle we desire, but close) */
    _quadrature->setPhi(atan((height * _num_x[i]) / (width * _num_y[i])), i);

    /* Effective Track spacing (not spacing we desire, but close) */
    _dx_eff[i]   = (width / _num_x[i]);
    _dy_eff[i]   = (height / _num_y[i]);
    _azim_spacings[i] = (_dx_eff[i] * sin(_quadrature->getPhi(i)));

    /* The length of all tracks in a 2D cycle */
    _cycle_length[i] = _dx_eff[i] / cos(_quadrature->getPhi(i)) * 
      leastCommonMultiple(2*_num_x[i], 2 * height /
                          (tan(_quadrature->getPhi(i)) * _dx_eff[i]));

    /* Get the number of tracks per cycle */
    _tracks_per_cycle[i] =
      (int) (round(_cycle_length[i] * sin(_quadrature->getPhi(i)) / width) +
             round(_cycle_length[i] * cos(_quadrature->getPhi(i)) / height));

    /* Compute the number of cycles */
    _cycles_per_azim[i] = (_num_x[i] + _num_y[i]) * 2 / _tracks_per_cycle[i];
  }

  /* Generate the 2D track cycles */
  for (int i = 0; i < _num_azim/4; i++) {
    
    /* Allocate memory for the 2D tracks array */
    _tracks_2D[i] = new Track2D*[_cycles_per_azim[i]];
    
    /* Start making the track cycles */
    for (int c = 0; c < _cycles_per_azim[i]; c++){

      /* Set the initial azimuthal angle */
      phi = _quadrature->getPhi(i);
      
      /* Get pointer to first track */
      _tracks_2D[i][c] = new Track2D[_tracks_per_cycle[i]];
      Track2D* track = &_tracks_2D[i][c][0];
      Track2D* track_prev = track;

      /* Create first track in cycle */
      track->getStart()->setCoords(width - _dx_eff[i] * (c + 0.5), 0.0);
      track->setPhi(phi);
      phi = findTrackEndPoint(track, phi, i);
      track->setBCIn(_geometry->getMinYBoundaryType());
      
      /* Generate the 2D track start and end points in this cycle */
      for (int t = 1; t < _tracks_per_cycle[i]; t++){

        /* Set start point and other parameters for the track */
        track = &_tracks_2D[i][c][t];
        track->getStart()->setCoords(track_prev->getEnd()->getX(),
                                     track_prev->getEnd()->getY());
        track->setPhi(phi);
        track->setBCIn(track_prev->getBCOut());
        
        /* Get the endpoint of the track and return the complementary azimuthal angle */
        phi = findTrackEndPoint(track, phi, i);
        
        /* Set previous track to current track */
        track_prev = track;
      }
    }
  }

  _contains_2D_tracks = true;
}


void TrackGenerator::initialize2DTrackReflections(){

  /* Generate the 2D track cycles */
  for (int a = 0; a < _num_azim/4; a++) {
    for (int c = 0; c < _cycles_per_azim[a]; c++){
      for (int t = 0; t < _tracks_per_cycle[a]; t++){      
        
        Track2D* track = &_tracks_2D[a][c][t];

        if (t == 0)
          track->setTrackIn(&_tracks_2D[a][c][_tracks_per_cycle[a]-1]);
        else
          track->setTrackIn(&_tracks_2D[a][c][t-1]);
        
        if (t == _tracks_per_cycle[a] - 1)
          track->setTrackOut(&_tracks_2D[a][c][0]);
        else
          track->setTrackOut(&_tracks_2D[a][c][t+1]);
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
  _tracks_3D        = new Track3D****[_num_azim/4];
  _tracks_per_plane = new int***[_num_azim/4];
  _num_z            = new int*[_num_azim/4];
  _num_l            = new int*[_num_azim/4];
  _polar_spacings   = new double*[_num_azim/4];
    
  for (int i=0; i < _num_azim/4; i++){
    _dz_eff[i]         = new double[_num_polar/2];
    _dl_eff[i]         = new double[_num_polar/2];
    _num_z[i]          = new int[_num_polar/2];
    _num_l[i]          = new int[_num_polar/2];
    _polar_spacings[i] = new double[_num_polar/2];
  }

  double x1, x2, y1, y2, z1, z2;
  double theta;
  double depth  = _geometry->getDepth();

  /* Determine angular quadrature and track spacing */
  for (int i = 0; i < _num_azim/4; i++) {
    
    /* Determine the polar angles and spacing for this azimuthal angle */
    for (int j=0; j < _num_polar/2; j++){

      /* Compute the cosine weighted average angle */
      theta = _quadrature->getTheta(i, j);

      /* The number of intersections with xy (denoted "l") and z-axes */
      _num_l[i][j] = (int) (fabs(_cycle_length[i] * tan(M_PI_2 - theta) * sin(theta)
                                 / _polar_spacing)) + 1;
      _num_z[i][j] = (int) (fabs(depth * _num_l[i][j] * tan(theta)
                                 / _cycle_length[i])) + 1;

      /* Effective track spacing */
      _dl_eff[i][j]          = _cycle_length[i] / _num_l[i][j];
      _dz_eff[i][j]          = depth            / _num_z[i][j];
      _quadrature->setTheta(atan(_dl_eff[i][j] / _dz_eff[i][j]), i, j);
      _polar_spacings[i][j] = _dz_eff[i][j] * sin(_quadrature->getTheta(i, j));
    }
  }

  /* Allocate memory for 3D tracks */
  for (int a = 0; a < _num_azim/4; a++) {
    _tracks_3D[a] = new Track3D***[_cycles_per_azim[a]];
    for (int c = 0; c < _cycles_per_azim[a]; c++){
      _tracks_3D[a][c] = new Track3D**[_num_polar];
      for (int p=0; p < _num_polar; p++)
        _tracks_3D[a][c][p] = new Track3D*[getNumZ(a,p) + getNumL(a,p)];
    }
  }

  Track3D track_3D;
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
  for (int a = 0; a < _num_azim/4; a++) {
    _tracks_per_plane[a] = new int**[_cycles_per_azim[a]];
    for (int c = 0; c < _cycles_per_azim[a]; c++){
      _tracks_per_plane[a][c] = new int*[_num_polar];

      /* Loop over polar angles < PI/2 */
      for (int p=0; p < _num_polar/2; p++){

        /* Allocate the memory for the number of tracks per lz track */
        _tracks_per_plane[a][c][p] = new int[_num_l[a][p] + _num_z[a][p]];
        
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
        for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++){
          
          /* Get the starting point */
          if (i < _num_l[a][p]){
            l_start = _cycle_length[a] - (i + 0.5) * _dl_eff[a][p];
            x1 = convertLtoX(l_start, a, c);
            y1 = convertLtoY(l_start, a, c);
            z1 = 0.0;
            track_3D.setBCIn(_geometry->getMinZBoundaryType());
          }
          else{
            l_start = 0.0;
            x1 = _tracks_2D[a][c][0].getStart()->getX();
            y1 = _tracks_2D[a][c][0].getStart()->getY();
            z1 = _dz_eff[a][p] * (i - _num_l[a][p] + 0.5);
            track_3D.setBCIn(_geometry->getMinYBoundaryType());
          }
          
          /* Get the end point */
          if (i < _num_z[a][p]){
            l_end = _cycle_length[a];
            x2 = _tracks_2D[a][c][0].getStart()->getX();
            y2 = _tracks_2D[a][c][0].getStart()->getY();
            z2 = _dz_eff[a][p] * (i + 0.5);
            track_3D.setBCOut(_geometry->getMinYBoundaryType());
          }
          else{
            l_end = _cycle_length[a] - _dl_eff[a][p] * (i - _num_z[a][p] + 0.5);
            x2 = convertLtoX(l_end, a, c);
            y2 = convertLtoY(l_end, a, c);
            z2 = depth;
            track_3D.setBCOut(_geometry->getMaxZBoundaryType());
          }

          /* Set start and end points and save polar angle */
          track_3D.getStart()->setCoords(x1, y1, z1);
          track_3D.getEnd()->setCoords(x2, y2, z2);
          track_3D.setTheta(_quadrature->getTheta(a,p));

          /* Decompose the track in the LZ plane by splitting it
           * based on the x and y geometry boundaries */
          decomposeLZTrack(&track_3D, l_start, l_end, a, c, p, i);
        }
      }

      /* Create tracks for polar angles [PI/2,PI] */
      for (int p=0; p < _num_polar/2; p++){
        
        pc = _num_polar-p-1;
        _tracks_per_plane[a][c][pc] = new int[_num_l[a][p] + _num_z[a][p]];
        
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
        for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++){
          
          /* Get the starting point */
          if (i < _num_z[a][p]){
            l_start = 0.0;
            x1 = _tracks_2D[a][c][0].getStart()->getX();
            y1 = _tracks_2D[a][c][0].getStart()->getY();
            z1 = _dz_eff[a][p] * (i + 0.5);
            track_3D.setBCIn(_geometry->getMinYBoundaryType());
          }
          else{
            l_start = _dl_eff[a][p] * (i - _num_z[a][p] + 0.5);
            x1 = convertLtoX(l_start, a, c);
            y1 = convertLtoY(l_start, a, c);
            z1 = depth;
            track_3D.setBCIn(_geometry->getMaxZBoundaryType());
          }
          
          /* Get the end point */
          if (i < _num_l[a][p]){
            l_end = _dl_eff[a][p] * (i + 0.5);
            x2 = convertLtoX(l_end, a, c);
            y2 = convertLtoY(l_end, a, c);
            z2 = 0.0;
            track_3D.setBCOut(_geometry->getMinZBoundaryType());
          }
          else{
            l_end = _cycle_length[a];
            x2 = _tracks_2D[a][c][0].getStart()->getX();
            y2 = _tracks_2D[a][c][0].getStart()->getY();
            z2 = _dz_eff[a][p] * (i - _num_l[a][p] + 0.5);
            track_3D.setBCOut(_geometry->getMinYBoundaryType());
          }

          /* Set start and end points and save polar angle */
          track_3D.getStart()->setCoords(x1, y1, z1);
          track_3D.getEnd()->setCoords(x2, y2, z2);
          track_3D.setTheta(M_PI - _quadrature->getTheta(a,p));
          
          /* Decompose the track in the LZ plane by splitting it
           * based on the x and y geometry boundaries */
          decomposeLZTrack(&track_3D, l_start, l_end, a, c, pc, i);
        }
      }
    }
  }

  _contains_3D_tracks = true;
}


void TrackGenerator::initialize3DTrackReflections(){

  int pc;
  
  /* Set reflective tracks indices */
  for (int a = 0; a < _num_azim/4; a++) {
    for (int c = 0; c < _cycles_per_azim[a]; c++){

      /* Loop over polar angles < PI/2 */
      for (int p=0; p < _num_polar/2; p++){

        /* Set the complementary polar angle */
        pc = _num_polar-p-1;

        for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){

            /* Set outgoing track */
            if (t == _tracks_per_plane[a][c][p][i]-1){
              if (i < _num_z[a][p]){
                _tracks_3D[a][c][p][i][t].setTrackOut
                  (&_tracks_3D[a][c][p][_num_l[a][p] + i][0]);
              }
              else{
                _tracks_3D[a][c][p][i][t].setTrackOut
                  (&_tracks_3D[a][c][pc][_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]);
              }
            }
            else{
              _tracks_3D[a][c][p][i][t].setTrackOut
                (&_tracks_3D[a][c][p][i][t+1]);
            }
              
            /* Set the incoming track */
            if (t == 0){
              if (i < _num_l[a][p]){
                _tracks_3D[a][c][p][i][t].setTrackIn
                  (&_tracks_3D[a][c][pc][_num_l[a][p] - i - 1]
                   [_tracks_per_plane[a][c][pc][_num_l[a][p] - i - 1] - 1]);
              }
              else{
                _tracks_3D[a][c][p][i][t].setTrackIn
                  (&_tracks_3D[a][c][p][i - _num_l[a][p]]
                   [_tracks_per_plane[a][c][p]
                    [i - _num_l[a][p]] - 1]);
              }
            }
            else{
              _tracks_3D[a][c][p][i][t].setTrackIn
                (&_tracks_3D[a][c][p][i][t-1]);
            }
          }
        }
      }

      /* Loop over polar angles > PI/2 */
      for (int p = _num_polar/2; p < _num_polar; p++){
        pc = _num_polar-p-1;
        for (int i=0; i < _num_l[a][pc] + _num_z[a][pc]; i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){

            /* Set outgoing track */
            if (t == _tracks_per_plane[a][c][p][i]-1){
              if (i < _num_l[a][pc]){
                _tracks_3D[a][c][p][i][t].setTrackOut
                  (&_tracks_3D[a][c][pc][_num_l[a][pc] - i - 1][0]);
              }
              else{
                _tracks_3D[a][c][p][i][t].setTrackOut
                  (&_tracks_3D[a][c][p][i - _num_l[a][pc]][0]);
              }
            }
            else{
              _tracks_3D[a][c][p][i][t].setTrackOut
                (&_tracks_3D[a][c][p][i][t+1]);
            }
              
            /* Set the incoming track */
            if (t == 0){
              if (i < _num_z[a][pc]){
              _tracks_3D[a][c][p][i][t].setTrackIn
                (&_tracks_3D[a][c][p][_num_l[a][pc] + _num_z[a][pc] - i - 1]
                 [_tracks_per_plane[a][c][p][_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]);
              }
              else{
              _tracks_3D[a][c][p][i][t].setTrackIn
                (&_tracks_3D[a][c][pc][2*_num_z[a][pc] + _num_l[a][pc] - i - 1]
                 [_tracks_per_plane[a][c][pc]
                  [2*_num_z[a][pc] + _num_l[a][pc] - i - 1] - 1]);
              }
            }
            else{
              _tracks_3D[a][c][p][i][t].setTrackIn
                (&_tracks_3D[a][c][p][i][t-1]);
            }
          }
        }
      }
    }
  }
}

 
void TrackGenerator::decomposeLZTrack(Track3D* track, double l_start, double l_end,
                                      int azim, int cycle, int polar, int lz_index){

  if (l_start > _cycle_length[azim] || l_start < 0.0)
    log_printf(ERROR, "Unable to tally tracks per stack for track with "
               "starting length not within "
               "(0, %f), l_start: %f", _cycle_length[azim], l_start);
  else if  (l_start > _cycle_length[azim] || l_start < 0.0)
    log_printf(ERROR, "Unable to tally tracks per stack for track with "
               "ending length not within "
               "(0, %f), l_end: %f", _cycle_length[azim], l_end);

  double length_sum = 0.0;
  int first_stack = 0;
  int last_stack = _tracks_per_cycle[azim] - 1;
  
  /* Find the last cycle index */
  for (int i=0; i < _tracks_per_cycle[azim]; i++){
    if (l_end < length_sum + _tracks_2D[azim][cycle][i].getLength() + ON_SURFACE_THRESH){
      last_stack = i;
      break;
    }    
    else{
      length_sum += _tracks_2D[azim][cycle][i].getLength();
    }
  }
  
  length_sum = 0.0;
  
  /* Find the first cycle index */
  for (int i=0; i < _tracks_per_cycle[azim]; i++){
    if (l_start < length_sum + _tracks_2D[azim][cycle][i].getLength() - ON_SURFACE_THRESH){
      first_stack = i;
      break;
    }
    else{
      length_sum += _tracks_2D[azim][cycle][i].getLength();
    }
  }
  
  if (last_stack < first_stack)
    log_printf(ERROR, "Unable to tally tracks per stack for track with "
               "last stack less than first stack. "
               "first stack: %i, last stack: %i", first_stack, last_stack);

  /* Set the number of tracks in the lz track */
  _tracks_per_plane[azim][cycle][polar][lz_index] =
    last_stack - first_stack + 1;

  /* Allocate memory for 3D tracks */
  _tracks_3D[azim][cycle][polar][lz_index] =
    new Track3D[_tracks_per_plane[azim][cycle][polar][lz_index]];

  Track3D* track2;
  double x1, y1, z1;
  double x2, y2, z2;
  int t=0;
  double theta = track->getTheta();
  
  /* Set the start and end point for each 3D track */
  for (int i=first_stack; i <= last_stack; i++){

    track2 = &_tracks_3D[azim][cycle][polar][lz_index][t];
    
    /* Find the starting coords */
    if (i == first_stack){
      x1 = track->getStart()->getX();
      y1 = track->getStart()->getY();
      z1 = track->getStart()->getZ();
      track2->setBCIn(track->getBCIn());
    }
    else{
      x1 = x2;
      y1 = y2;
      z1 = z2;
      track2->setBCIn(_tracks_2D[azim][cycle][i].getBCIn());
    }
    
    /* Set the starting point */
    track2->getStart()->setCoords(x1, y1, z1);

    /* Find the ending length */
    if (i == last_stack){
      x2 = track->getEnd()->getX();
      y2 = track->getEnd()->getY();
      z2 = track->getEnd()->getZ();
      track2->setBCOut(track->getBCOut());
    }
    else{
      length_sum += _tracks_2D[azim][cycle][i].getLength();
      x2 = _tracks_2D[azim][cycle][i].getEnd()->getX();
      y2 = _tracks_2D[azim][cycle][i].getEnd()->getY();
      z2 = z1 + (length_sum - l_start) * tan(M_PI_2 - theta);
      l_start = length_sum;
      track2->setBCOut(_tracks_2D[azim][cycle][i].getBCOut());
    }      
    
    /* Set the ending point and angles */
    track2->getEnd()->setCoords(x2, y2, z2);
    track2->setTheta(theta);
    track2->setPhi(_tracks_2D[azim][cycle][i].getPhi());
    track2->setAzimAngleIndex(_tracks_2D[azim][cycle][i].getAzimAngleIndex());
    
    t++;
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

  int uid = 0;
  _num_2D_tracks = 0;
  
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int t=0; t < _tracks_per_cycle[a]; t++){
        _tracks_2D[a][c][t].setUid(uid);
        uid++;
            
        double x0 = _tracks_2D[a][c][t].getStart()->getX();
        double y0 = _tracks_2D[a][c][t].getStart()->getY();
        double x1 = _tracks_2D[a][c][t].getEnd()->getX();
        double y1 = _tracks_2D[a][c][t].getEnd()->getY();
        double new_x0 = x0 + _geometry->getMinX();
        double new_y0 = y0 + _geometry->getMinY();
        double new_x1 = x1 + _geometry->getMinX();
        double new_y1 = y1 + _geometry->getMinY();
                
        _tracks_2D[a][c][t].setCoords(new_x0, new_y0, new_x1, new_y1);
        _num_2D_tracks++;
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
void TrackGenerator::recalibrate3DTracksToOrigin() {

  int uid = 0;
  _num_3D_tracks = 0;
  
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
            _tracks_3D[a][c][p][i][t].setUid(uid);
            _tracks_3D[a][c][p][i][t].setPolarAngleIndex(p);
            uid++;
            
            double x0 = _tracks_3D[a][c][p][i][t].getStart()->getX();
            double y0 = _tracks_3D[a][c][p][i][t].getStart()->getY();
            double z0 = _tracks_3D[a][c][p][i][t].getStart()->getZ();
            double x1 = _tracks_3D[a][c][p][i][t].getEnd()->getX();
            double y1 = _tracks_3D[a][c][p][i][t].getEnd()->getY();
            double z1 = _tracks_3D[a][c][p][i][t].getEnd()->getZ();
            double new_x0 = x0 + _geometry->getMinX();
            double new_y0 = y0 + _geometry->getMinY();
            double new_z0 = z0 + _geometry->getMinZ();
            double new_x1 = x1 + _geometry->getMinX();
            double new_y1 = y1 + _geometry->getMinY();
            double new_z1 = z1 + _geometry->getMinZ();
            
            _tracks_3D[a][c][p][i][t]
              .setCoords(new_x0, new_y0, new_z0, new_x1, new_y1, new_z1);
            _num_3D_tracks++;            
          }
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

  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      
      log_printf(NORMAL, "segmenting 2D cycle %i, %i", a, c);

      #pragma omp parallel for
      for (int t=0; t < _tracks_per_cycle[a]; t++){
        _geometry->segmentize2D(&_tracks_2D[a][c][t], _z_level);
      }
    }
  }

  _num_2D_segments = 0;

  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int t=0; t < _tracks_per_cycle[a]; t++)
        _num_2D_segments += _tracks_2D[a][c][t].getNumSegments();
    }
  }

  log_printf(NORMAL, "num 2d segments: %i", _num_2D_segments);
  
  _contains_2D_segments = true;
  
  return;
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize3D() {

  log_printf(NORMAL, "Ray tracing for 3D track segmentation...");

  int nt;
  
  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      log_printf(NORMAL, "segmenting 3D cycle %i, %i", a, c);
      for (int p=0; p < _num_polar; p++){
        nt = getNumZ(a,p) + getNumL(a,p);
        #pragma omp parallel for
        for (int i=0; i < nt; i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
            _geometry->segmentize3D(&_tracks_3D[a][c][p][i][t]);
          }
        }
      }
    }
  }

  _num_3D_segments = 0;
  
  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++)
            _num_3D_segments += _tracks_3D[a][c][p][i][t].getNumSegments();
        }
      }
    }
  }

  _contains_3D_segments = true;
  
  return;
}


double TrackGenerator::findTrackEndPoint(Track2D* track, double phi, int azim_index){

  /* Get start and end points of new track */  
  double x = track->getStart()->getX();
  double y = track->getStart()->getY();
  double width = _geometry->getWidth();
  double height = _geometry->getHeight();

  if (phi < M_PI_2)
    track->setAzimAngleIndex(azim_index);
  else if (phi < M_PI)
    track->setAzimAngleIndex(_num_azim/2 - azim_index - 1);
  else if (phi < 3.0 * M_PI_2)
    track->setAzimAngleIndex(_num_azim/2 + azim_index);
  else
    track->setAzimAngleIndex(_num_azim - azim_index - 1);

  /* X_MIN side */
  if (x == 0.0){
    if (phi < M_PI_2){
      /* X_MIN to X_MAX */
      if (y + width * tan(phi) < height){
        track->getEnd()->setCoords(width, y + width * tan(phi));
        phi = M_PI - phi;
        track->setBCOut(_geometry->getMaxXBoundaryType());
      }
      
      /* X_MIN to Y_MAX */
      else{
        track->getEnd()->setCoords((height-y) * tan(M_PI_2 - phi), height);
        phi = 2.0 * M_PI - phi;
        track->setBCOut(_geometry->getMaxYBoundaryType());
      }

    }
    else{
      /* X_MIN to X_MAX */
      if (y + width * tan(phi) > 0.0){
        track->getEnd()->setCoords(width, y + width * tan(phi));
        phi = 3.0 * M_PI - phi;
        track->setBCOut(_geometry->getMaxXBoundaryType());
      }
        
      /* X_MIN to Y_MIN */
      else{
        track->getEnd()->setCoords(-y * tan(M_PI_2 - phi), 0.0);
        phi = 2.0 * M_PI - phi;
        track->setBCOut(_geometry->getMinYBoundaryType());
      }
    }
  }
  
  /* Y_MIN side */
  else if (y == 0.0){
    if (phi < M_PI_2){
      /* Y_MIN to X_MAX */
      if ((width-x) * tan(phi) < height){
        track->getEnd()->setCoords(width, (width-x) * tan(phi));
        phi = M_PI - phi;
        track->setBCOut(_geometry->getMaxXBoundaryType());
      }
      
      /* Y_MIN to Y_MAX */
      else{
        track->getEnd()->setCoords(x + height * tan(M_PI_2 - phi), height);
        phi = 2.0 * M_PI - phi;
        track->setBCOut(_geometry->getMaxYBoundaryType());
      }
    }
    else{
      /* Y_MIN to X_MIN */
      if (-x * tan(phi) < height){
        track->getEnd()->setCoords(0.0, -x * tan(phi));
        phi = M_PI - phi;
        track->setBCOut(_geometry->getMinXBoundaryType());
      }
      
      /* Y_MIN to Y_MAX */
      else{
        track->getEnd()->setCoords(x + height * tan(M_PI_2 - phi), height);
        phi = 2.0 * M_PI - phi;
        track->setBCOut(_geometry->getMaxYBoundaryType());
      }
    }
  }
  
  /* X_MAX side */
  else if (x == width){
    if (phi > M_PI){
      /* X_MAX to X_MIN */
      if (y - width * tan(phi) > 0.0){
        track->getEnd()->setCoords(0.0, y - width * tan(phi));
        phi = 3.0 * M_PI - phi;
        track->setBCOut(_geometry->getMinXBoundaryType());
      }
        
      /* X_MAX to Y_MIN */
      else{
        track->getEnd()->setCoords(width - y * tan(M_PI_2 - phi), 0.0);
       phi = 2.0 * M_PI - phi;
       track->setBCOut(_geometry->getMinYBoundaryType());
      }
    }
    else{
      /* X_MAX to X_MIN */
      if (y - width * tan(phi) < height){
        track->getEnd()->setCoords(0.0, y - width * tan(phi));
        phi = M_PI - phi;
        track->setBCOut(_geometry->getMinXBoundaryType());
      }
        
      /* X_MAX to Y_MAX */
      else{
        track->getEnd()->setCoords(width + (height-y) * tan(M_PI_2 - phi), height);
        phi = 2.0 * M_PI - phi;
        track->setBCOut(_geometry->getMaxYBoundaryType());
      }
    }
  }
  
  /* Y_MAX side */
  else {
    if (phi < 3.0 * M_PI_2){
      /* Y_MAX to X_MIN */
      if (height - x * tan(phi) > 0.0){
        track->getEnd()->setCoords(0.0, height - x * tan(phi));
        phi = 3.0 * M_PI - phi;
        track->setBCOut(_geometry->getMinXBoundaryType());
      }
      
      /* Y_MAX to Y_MIN */
      else{
        track->getEnd()->setCoords(x - height * tan(M_PI_2 - phi), 0.0);
        phi = 2.0 * M_PI - phi;
        track->setBCOut(_geometry->getMinYBoundaryType());
      }
    }
    else{
      /* Y_MAX to Y_MIN */
      if (x - height / tan(phi) < width){
        track->getEnd()->setCoords(x - height * tan(M_PI_2 - phi), 0.0);
        phi = 2.0 * M_PI - phi;
        track->setBCOut(_geometry->getMinYBoundaryType());
      }
      
      /* Y_MAX to X_MAX */
      else{
        track->getEnd()->setCoords(width, height + (width-x) * tan(phi));
        phi = 3.0 * M_PI - phi;
        track->setBCOut(_geometry->getMaxXBoundaryType());
      }
    }
  }

  return phi;
}


/* Convert a length traveled along a track cycle to value in x */
double TrackGenerator::convertLtoX(double l, int azim, int cycle){

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to X since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);
  
  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++){
    if (l <= length_sum + _tracks_2D[azim][cycle][i].getLength()){
      track_index = i;
      break;
    }
    else{
      length_sum += _tracks_2D[azim][cycle][i].getLength();
    }
  }

  if (l - length_sum < 0.0)
    log_printf(NORMAL, "found negative length residual in converting l to x");
  
  double x1 = _tracks_2D[azim][cycle][track_index].getStart()->getX();
  double x2 = _tracks_2D[azim][cycle][track_index].getEnd()->getX();
  double l_rel = (l - length_sum) / _tracks_2D[azim][cycle][track_index].getLength();
  double x = l_rel * x2 + (1.0 - l_rel) * x1;

  return x;
}


/* Convert a length traveled along a track cycle to value in y */
double TrackGenerator::convertLtoY(double l, int azim, int cycle){

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to Y since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);
  
  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++){
    if (l <= length_sum + _tracks_2D[azim][cycle][i].getLength()){
      track_index = i;
      break;
    }
    else{
      length_sum += _tracks_2D[azim][cycle][i].getLength();
    }
  }

  if (l - length_sum < 0.0)
    log_printf(NORMAL, "found negative length residual in converting l to y");
  
  double y1 = _tracks_2D[azim][cycle][track_index].getStart()->getY();
  double y2 = _tracks_2D[azim][cycle][track_index].getEnd()->getY();
  double l_rel = (l - length_sum) / _tracks_2D[azim][cycle][track_index].getLength();
  double y = l_rel * y2 + (1.0 - l_rel) * y1;

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

  if (_solve_3D){

    /* Get the quadrature type */
    std::string quad_type;
    
    if (_quadrature->getQuadratureType() == EQUAL_WEIGHT)
      quad_type = "EQ_WGT";
    else if (_quadrature->getQuadratureType() == EQUAL_ANGLE)
      quad_type = "EQ_ANG";
    else
      log_printf(ERROR, "Unable to solve 3D problem with quadrature type"
                 " other than EQUAL_WEIGHT or EQUAL_ANGLE");

    if (_geometry->getCmfd() != NULL){
      test_filename << directory.str() << "/3D_"
                    << _num_azim << "_azim_"
                    << _num_polar << "_polar_"
                    << _azim_spacing << "x" << _polar_spacing << "_cm_spacing_cmfd_"
                    << _geometry->getCmfd()->getNumX()
                    << "x" << _geometry->getCmfd()->getNumY()
                    << "x" << _geometry->getCmfd()->getNumZ()
                    << "_quad_" << quad_type
                    << ".data";
    }
    else{
      test_filename << directory.str() << "/3D_"
                    << _num_azim << "_azim_"
                    << _num_polar << "_polar_"
                    << _azim_spacing << "x" << _polar_spacing << "_cm_spacing_quad_"
                    << quad_type << ".data";
    }
  }
  else{
    if (_geometry->getCmfd() != NULL){
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
    if (_solve_3D){
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

  log_printf(NORMAL, "Dumping 2D segments to file...");

  if (!_contains_2D_segments)
    log_printf(ERROR, "Unable to dump 2D Segments to a file since no Segments have "
               "been generated for %d azimuthal angles and %f track spacing",
               _num_azim, _azim_spacing);

  FILE* out;
  out = fopen(_tracks_filename.c_str(), "w");

  /* Get a string representation of the Geometry's attributes. This is used to
   * check whether or not ray tracing has been performed for this Geometry */
  std::string geometry_to_string = _geometry->toString();
  int string_length = geometry_to_string.length() + 1;

  /* Write geometry metadata to the Track file */
  fwrite(&string_length, sizeof(int), 1, out);
  fwrite(geometry_to_string.c_str(), sizeof(char)*string_length, 1, out);
  fwrite(&_z_level, sizeof(double), 1, out);
  
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
  for (int a=0; a < _num_azim/4; a++) {
    for (int c=0; c < _cycles_per_azim[a]; c++) {
      for (int t=0; t < _tracks_per_cycle[a]; t++) {

        /* Get data for this Track */
        curr_track = &_tracks_2D[a][c][t];
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
          if (cmfd != NULL){
            cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
            cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;
            fwrite(&cmfd_surface_fwd, sizeof(int), 1, out);
            fwrite(&cmfd_surface_bwd, sizeof(int), 1, out);
          }
        }
      }
    }
  }

  /* Get FSR vector maps */
  std::unordered_map<std::size_t, fsr_data> FSR_keys_map = _geometry->getFSRKeysMap();
  std::unordered_map<std::size_t, fsr_data>::iterator iter;
  std::vector<std::size_t> FSRs_to_keys = _geometry->getFSRsToKeys();
  std::vector<int> FSRs_to_material_IDs = _geometry->getFSRsToMaterialIDs();
  std::size_t fsr_key;
  int fsr_id;
  int fsr_counter = 0;
  double x, y, z;
  
  /* Write number of FSRs */
  int num_FSRs = _geometry->getNumFSRs();
  fwrite(&num_FSRs, sizeof(int), 1, out);

  /* Write FSR vector maps to file */
  for (iter = FSR_keys_map.begin(); iter != FSR_keys_map.end(); ++iter){
    
    /* Write data to file from FSR_keys_map */
    fsr_key = iter->first;
    fsr_id = iter->second._fsr_id;
    x = iter->second._point->getX();
    y = iter->second._point->getY();
    z = iter->second._point->getZ();
    fwrite(&fsr_key, sizeof(std::size_t), 1, out);
    fwrite(&fsr_id, sizeof(int), 1, out);
    fwrite(&x, sizeof(double), 1, out);
    fwrite(&y, sizeof(double), 1, out);
    fwrite(&z, sizeof(double), 1, out);

    /* Write data to file from FSRs_to_material_IDs */
    fwrite(&(FSRs_to_material_IDs.at(fsr_counter)), sizeof(int), 1, out);
    
    /* Write data to file from FSRs_to_keys */
    fwrite(&(FSRs_to_keys.at(fsr_counter)), sizeof(std::size_t), 1, out);
    
    /* Increment FSR ID counter */
    fsr_counter++;
  }

  /* Write cmfd_fsrs vector of vectors to file */
  if (cmfd != NULL){
    std::vector< std::vector<int> > cell_fsrs = cmfd->getCellFSRs();
    std::vector<int>::iterator iter;
    int num_cells = cmfd->getNumCells();
    fwrite(&num_cells, sizeof(int), 1, out);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++){
      num_FSRs = cell_fsrs.at(cell).size();
      fwrite(&num_FSRs, sizeof(int), 1, out);

      /* Loop over FSRs within cell */
      for (iter = cell_fsrs.at(cell).begin(); iter != cell_fsrs.at(cell).end();
           ++iter)
        fwrite(&(*iter), sizeof(int), 1, out);
    }
  }

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

  log_printf(NORMAL, "Dumping 3D segments to file...");

  if (!_contains_3D_segments)
    log_printf(ERROR, "Unable to dump 3D Segments to a file since no Segments have "
               "been generated for %d azimuthal angles and %f track spacing",
               _num_azim, _azim_spacing);

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
  for (int a=0; a < _num_azim/4; a++) {
    for (int c=0; c < _cycles_per_azim[a]; c++) {
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
        
            /* Get data for this Track */
            curr_track = &_tracks_3D[a][c][p][i][t];
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
              if (cmfd != NULL){
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
  }

  /* Get FSR vector maps */
  std::unordered_map<std::size_t, fsr_data> FSR_keys_map = _geometry->getFSRKeysMap();
  std::unordered_map<std::size_t, fsr_data>::iterator iter;
  std::vector<std::size_t> FSRs_to_keys = _geometry->getFSRsToKeys();
  std::vector<int> FSRs_to_material_IDs = _geometry->getFSRsToMaterialIDs();
  std::size_t fsr_key;
  int fsr_id;
  int fsr_counter = 0;
  double x, y, z;
  
  /* Write number of FSRs */
  int num_FSRs = _geometry->getNumFSRs();
  fwrite(&num_FSRs, sizeof(int), 1, out);

  /* Write FSR vector maps to file */
  for (iter = FSR_keys_map.begin(); iter != FSR_keys_map.end(); ++iter){
    
    /* Write data to file from FSR_keys_map */
    fsr_key = iter->first;
    fsr_id = iter->second._fsr_id;
    x = iter->second._point->getX();
    y = iter->second._point->getY();
    z = iter->second._point->getZ();
    fwrite(&fsr_key, sizeof(std::size_t), 1, out);
    fwrite(&fsr_id, sizeof(int), 1, out);
    fwrite(&x, sizeof(double), 1, out);
    fwrite(&y, sizeof(double), 1, out);
    fwrite(&z, sizeof(double), 1, out);

    /* Write data to file from FSRs_to_material_IDs */
    fwrite(&(FSRs_to_material_IDs.at(fsr_counter)), sizeof(int), 1, out);
    
    /* Write data to file from FSRs_to_keys */
    fwrite(&(FSRs_to_keys.at(fsr_counter)), sizeof(std::size_t), 1, out);
    
    /* Increment FSR ID counter */
    fsr_counter++;
  }

  /* Write cmfd_fsrs vector of vectors to file */
  if (cmfd != NULL){
    std::vector< std::vector<int> > cell_fsrs = cmfd->getCellFSRs();
    std::vector<int>::iterator iter;
    int num_cells = cmfd->getNumCells();
    fwrite(&num_cells, sizeof(int), 1, out);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++){
      num_FSRs = cell_fsrs.at(cell).size();
      fwrite(&num_FSRs, sizeof(int), 1, out);

      /* Loop over FSRs within cell */
      for (iter = cell_fsrs.at(cell).begin(); iter != cell_fsrs.at(cell).end();
           ++iter)
        fwrite(&(*iter), sizeof(int), 1, out);
    }
  }

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
  double z_level;
  
  /* Import Geometry metadata from the Track file */
  ret = fread(&string_length, sizeof(int), 1, in);
  char* geometry_to_string = new char[string_length];
  ret = fread(geometry_to_string, sizeof(char)*string_length, 1, in);
  ret = fread(&z_level, sizeof(double), 1, in);
  
  /* Check if our Geometry is exactly the same as the Geometry in the
   * Track file for this number of azimuthal angles and track spacing */
  if (_geometry->toString().compare(std::string(geometry_to_string)) != 0 ||
      _z_level != z_level)
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
  //segment curr_segment;

  std::map<int, Material*> materials = _geometry->getAllMaterials();

  int uid = 0;

  /* Loop over Tracks */
  for (int a=0; a < _num_azim/4; a++) {
    for (int c=0; c < _cycles_per_azim[a]; c++) {
      for (int t=0; t < _tracks_per_cycle[a]; t++) {

        /* Import data for this Track from Track file */
        ret = fread(&num_segments, sizeof(int), 1, in);
        
        /* Initialize a Track with this data */
        curr_track = &_tracks_2D[a][c][t];

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
          if (cmfd != NULL){
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

  /* Create FSR vector maps */
  std::unordered_map<std::size_t, fsr_data> FSR_keys_map;
  std::vector<int> FSRs_to_material_IDs;
  std::vector<std::size_t> FSRs_to_keys;
  int num_FSRs;
  std::size_t fsr_key;
  int fsr_key_id;
  double x, y, z;
  
  /* Get number of FSRs */
  ret = fread(&num_FSRs, sizeof(int), 1, in);
  _geometry->setNumFSRs(num_FSRs);
  
  /* Read FSR vector maps from file */
  for (int fsr_id=0; fsr_id < num_FSRs; fsr_id++){

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
    FSR_keys_map[fsr_key] = *fsr;

    /* Read data from file for FSR_to_materials_IDs */
    ret = fread(&material_id, sizeof(int), 1, in);
    FSRs_to_material_IDs.push_back(material_id);

    /* Read data from file for FSR_to_keys */
    ret = fread(&fsr_key, sizeof(std::size_t), 1, in);
    FSRs_to_keys.push_back(fsr_key);
  }

  /* Set FSR vector maps */
  _geometry->setFSRKeysMap(FSR_keys_map);
  _geometry->setFSRsToMaterialIDs(FSRs_to_material_IDs);
  _geometry->setFSRsToKeys(FSRs_to_keys);

  /* Read cmfd cell_fsrs vector of vectors from file */
  if (cmfd != NULL){
    std::vector< std::vector<int> > cell_fsrs;
    int num_cells, fsr_id;
    ret = fread(&num_cells, sizeof(int), 1, in);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++){
      std::vector<int> *fsrs = new std::vector<int>;
      cell_fsrs.push_back(*fsrs);
      ret = fread(&num_FSRs, sizeof(int), 1, in);

      /* Loop over FRSs within cell */
      for (int fsr = 0; fsr < num_FSRs; fsr++){
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
  //segment curr_segment;

  std::map<int, Material*> materials = _geometry->getAllMaterials();

  int uid = 0;

  /* Loop over Tracks */
  for (int a=0; a < _num_azim/4; a++) {
    for (int c=0; c < _cycles_per_azim[a]; c++) {
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
            
            /* Import data for this Track from Track file */
            ret = fread(&num_segments, sizeof(int), 1, in);

            /* Get data for this Track */
            curr_track = &_tracks_3D[a][c][p][i][t];
                        
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
              if (cmfd != NULL){
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
  }

  /* Create FSR vector maps */
  std::unordered_map<std::size_t, fsr_data> FSR_keys_map;
  std::vector<int> FSRs_to_material_IDs;
  std::vector<std::size_t> FSRs_to_keys;
  int num_FSRs;
  std::size_t fsr_key;
  int fsr_key_id;
  double x, y, z;
  
  /* Get number of FSRs */
  ret = fread(&num_FSRs, sizeof(int), 1, in);
  _geometry->setNumFSRs(num_FSRs);
  
  /* Read FSR vector maps from file */
  for (int fsr_id=0; fsr_id < num_FSRs; fsr_id++){

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
    FSR_keys_map[fsr_key] = *fsr;

    /* Read data from file for FSR_to_materials_IDs */
    ret = fread(&material_id, sizeof(int), 1, in);
    FSRs_to_material_IDs.push_back(material_id);

    /* Read data from file for FSR_to_keys */
    ret = fread(&fsr_key, sizeof(std::size_t), 1, in);
    FSRs_to_keys.push_back(fsr_key);
  }

  /* Set FSR vector maps */
  _geometry->setFSRKeysMap(FSR_keys_map);
  _geometry->setFSRsToMaterialIDs(FSRs_to_material_IDs);
  _geometry->setFSRsToKeys(FSRs_to_keys);

  /* Read cmfd cell_fsrs vector of vectors from file */
  if (cmfd != NULL){
    std::vector< std::vector<int> > cell_fsrs;
    int num_cells, fsr_id;
    ret = fread(&num_cells, sizeof(int), 1, in);

    /* Loop over CMFD cells */
    for (int cell=0; cell < num_cells; cell++){
      std::vector<int> *fsrs = new std::vector<int>;
      cell_fsrs.push_back(*fsrs);
      ret = fread(&num_FSRs, sizeof(int), 1, in);

      /* Loop over FRSs within cell */
      for (int fsr = 0; fsr < num_FSRs; fsr++){
        ret = fread(&fsr_id, sizeof(int), 1, in);
        cell_fsrs.at(cell).push_back(fsr_id);
      }
    }

    /* Set CMFD cell_fsrs vector of vectors */
    cmfd->setCellFSRs(cell_fsrs);
  }

  _num_3D_segments = 0;
  
  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/4; a++){
    for (int c=0; c < _cycles_per_azim[a]; c++){
      for (int p=0; p < _num_polar; p++){
        for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
          for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++)
            _num_3D_segments += _tracks_3D[a][c][p][i][t].getNumSegments();
        }
      }
    }
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

  int num_cuts, min_num_cuts;
  segment* curr_segment;

  FP_PRECISION length, tau;
  int fsr_id;
  Material* material;
  FP_PRECISION* sigma_t;
  int num_groups;

  if (_solve_3D){
    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++){
        for (int p=0; p < _num_polar; p++){
          for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
            for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
              for (int s=0; s < _tracks_3D[a][c][p][i][t].getNumSegments(); s++){
                
                /* Extract data from this segment to compute it optical length */
                curr_segment = _tracks_3D[a][c][p][i][t].getSegment(s);
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
                    new_segment->_cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;
                  
                  if (k == min_num_cuts-1)
                    new_segment->_cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
                  
                  /* Insert the new segment to the Track */
                  _tracks_3D[a][c][p][i][t].insertSegment(s+k+1, new_segment);
                }
                
                /* Remove the original segment from the Track */
                _tracks_3D[a][c][p][i][t].removeSegment(s);
              }
            }
          }
        }
      }
    }
  }
  else{
    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++){
        for (int t=0; t < _tracks_per_cycle[a]; t++){
          for (int s=0; s < _tracks_2D[a][c][t].getNumSegments(); s++){
            
            /* Extract data from this segment to compute it optical length */
            curr_segment = _tracks_2D[a][c][t].getSegment(s);
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
                new_segment->_cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;
              
              if (k == min_num_cuts-1)
                new_segment->_cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
              
              /* Insert the new segment to the Track */
              _tracks_2D[a][c][t].insertSegment(s+k+1, new_segment);
            }
            
            /* Remove the original segment from the Track */
            _tracks_2D[a][c][t].removeSegment(s);
          }
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
void TrackGenerator::generateFSRCentroids(){

  int num_FSRs = _geometry->getNumFSRs();
  FP_PRECISION* FSR_volumes;
  
  /* Create temporary array of centroids and initialize to origin */
  Point** centroids = new Point*[num_FSRs];
  for (int r=0; r < num_FSRs; r++){
    centroids[r] = new Point();
    centroids[r]->setCoords(0.0, 0.0, 0.0);
  }

  if (_solve_3D){

    FSR_volumes = get3DFSRVolumes();
    
    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++){
        for (int p=0; p < _num_polar; p++){
          for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++){
            for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){

              int num_segments = _tracks_3D[a][c][p][i][t].getNumSegments();
              segment* segments = _tracks_3D[a][c][p][i][t].getSegments();
              double x = _tracks_3D[a][c][p][i][t].getStart()->getX();
              double y = _tracks_3D[a][c][p][i][t].getStart()->getY();
              double z = _tracks_3D[a][c][p][i][t].getStart()->getZ();
              double phi = _tracks_3D[a][c][p][i][t].getPhi();
              double theta = _quadrature->getTheta(a, p);
              double wgt = _quadrature->getAzimWeight(a) *
                _quadrature->getPolarWeight(a, p) * getAzimSpacing(a)
                * getPolarSpacing(a,p);
              
              for (int s=0; s < num_segments; s++) {
                segment* curr_segment = &segments[s];
                int fsr = curr_segment->_region_id;
                double volume = FSR_volumes[fsr];
                centroids[fsr]->
                  setX(centroids[fsr]->getX() + wgt *
                       (x + cos(phi) * sin(theta) * curr_segment->_length / 2.0) *
                       curr_segment->_length / FSR_volumes[fsr]);
                
                centroids[fsr]->
                  setY(centroids[fsr]->getY() + wgt *
                       (y + sin(phi) * sin(theta) * curr_segment->_length / 2.0) *
                       curr_segment->_length / FSR_volumes[fsr]);
                
                centroids[fsr]->
                  setZ(centroids[fsr]->getZ() + wgt *
                       (z + cos(theta) * curr_segment->_length / 2.0) *
                       curr_segment->_length / FSR_volumes[fsr]);
                
                x += cos(phi) * sin(theta) * curr_segment->_length;
                y += sin(phi) * sin(theta) * curr_segment->_length;
                z += cos(theta) * curr_segment->_length;
              }
            }
          }
        }
      }
    }
  }
  else{
    
    FSR_volumes = get2DFSRVolumes();

    for (int a=0; a < _num_azim/4; a++){
      for (int c=0; c < _cycles_per_azim[a]; c++){
        for (int t=0; t < _tracks_per_cycle[a]; t++){

          int num_segments = _tracks_2D[a][c][t].getNumSegments();
          segment* segments = _tracks_2D[a][c][t].getSegments();
          double x = _tracks_2D[a][c][t].getStart()->getX();
          double y = _tracks_2D[a][c][t].getStart()->getY();
          double phi = _tracks_2D[a][c][t].getPhi();
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
  }

  /* Set the centroid for the FSR */
  for (int r=0; r < num_FSRs; r++){
    _geometry->setFSRCentroid(r, centroids[r]);
    delete centroids[r];
  }

  /* Delete temporary array of centroids and FSR volumes */
  delete [] centroids;
  delete [] FSR_volumes;
}
