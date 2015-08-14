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
  _OTF = false;
  _track_generation_method = GLOBAL_TRACKING;
}


/**
 * @brief Destructor frees memory for all Tracks.
 */
TrackGenerator::~TrackGenerator() {

  /* Deletes Tracks arrays if Tracks have been generated */
  if (_contains_3D_tracks){    
    for (int a=0; a < _num_azim/2; a++){
      for (int i=0; i < getNumX(a) + getNumY(a); i++){
        for (int p=0; p < _num_polar; p++){
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

    for (int a = 0; a < _num_azim/4; a++) {
      for (int c = 0; c < _cycles_per_azim[a]; c++){
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
    
    for (int a=0; a < _num_azim/4; a++){
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

    for (int a=0; a < _num_azim/2; a++)
      delete [] _tracks_2D[a];
    delete [] _tracks_2D;

    delete [] _num_x;
    delete [] _num_y;
    delete [] _dx_eff;
    delete [] _dy_eff;
    delete [] _azim_spacings;
    delete [] _cycles_per_azim;
    delete [] _tracks_per_cycle;
    delete [] _cycle_length;
  }

  if (_contains_extruded_tracks)
    delete[] _extruded_tracks;
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

  if (!contains2DSegments())
    log_printf(ERROR, "Cannot get the number of 2D segments since they "
               "have not been generated.");
  
  _num_2D_segments = 0;

  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      _num_2D_segments += _tracks_2D[a][i].getNumSegments();
    }
  }

  return _num_2D_segments;
}


/**
 * @brief Return the total number of Track segments across the Geometry.
 * @return the total number of Track segments
 */
int TrackGenerator::getNum3DSegments() {

  if (!contains3DSegments())
    log_printf(ERROR, "Cannot get the number of 3D segments since they "
               "have not been generated.");

  _num_3D_segments = 0;
  
  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          _num_3D_segments += _tracks_3D_stack[a][i][p][z].getNumSegments();
      }
    }
  }
  
  return _num_3D_segments;
}


/**
 * @brief Returns a 3D jagged array of the Tracks.
 * @details The first index into the array is the azimuthal angle and the
 *          second index is polar angle and the third index is the 
 *          Track number.
 * @return the 3D jagged array of Tracks
 */
Track2D** TrackGenerator::get2DTracks() {

  if (!_contains_2D_tracks)
    log_printf(ERROR, "Unable to return the 3D ragged array of the 2D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_2D;
}


Track3D**** TrackGenerator::get3DTracks() {

  if (!_contains_3D_tracks)
    log_printf(ERROR, "Unable to return the 3D ragged array of the 3D Tracks "
               "since Tracks have not yet been generated.");

  return _tracks_3D_stack;
}

ExtrudedTrack* TrackGenerator::getExtrudedTracks() {
  if (!_contains_extruded_tracks)
    log_printf(ERROR, "Unable to return the array of extruded tracks "
               "since extruded tracks have not yet been generated.");
  return _extruded_tracks;
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

  if (_solve_3D) {
    if (!_OTF) {
      for (int a=0; a < _num_azim/2; a++) {
        for (int i=0; i < getNumX(a) + getNumY(a); i++) {
          for (int p=0; p < _num_polar; p++) {
            for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
              for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments();
                  s++) {
                curr_segment = _tracks_3D_stack[a][i][p][z].getSegment(s);
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
    }
  }
  else{
    for (int a=0; a < _num_azim/2; a++){
      for (int i=0; i < getNumX(a) + getNumY(a); i++){
        for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++){
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


int* TrackGenerator::getTracksPerCycle(){
  return _tracks_per_cycle;
}


int*** TrackGenerator::getTracksPerStack(){
  return _tracks_per_stack;
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


double TrackGenerator::getDxEff(int azim){
  azim = _quadrature->getFirstOctantAzim(azim);
  return _dx_eff[azim];
}


double TrackGenerator::getDyEff(int azim){
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
   * all Track segments multipled by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++){
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
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
          for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments(); s++){
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
   * all Track segments multipled by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++){
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
   * all Track segments multipled by the Track "widths" for each FSR.  */
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
          for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments(); s++){
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
  _contains_extruded_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _contains_extruded_segments = false;
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
  _contains_extruded_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _contains_extruded_segments = false;
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
  _contains_extruded_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _contains_extruded_segments = false;
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
  _contains_extruded_tracks = false;
  _contains_2D_segments = false;
  _contains_3D_segments = false;
  _contains_extruded_segments = false;
}


void TrackGenerator::setSolve2D(){
  _solve_3D = false;
}


void TrackGenerator::setSolve3D(){
  _solve_3D = true;
}

void TrackGenerator::setOTF() {
  _OTF = true;
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
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      coords[counter]   = _tracks_2D[a][i].getStart()->getX();
      coords[counter+1] = _tracks_2D[a][i].getStart()->getY();
      coords[counter+2] = _tracks_2D[a][i].getEnd()->getX();
      coords[counter+3] = _tracks_2D[a][i].getEnd()->getY();
      
      counter += 4;
    }
  }

  return;
}


void TrackGenerator::retrieve2DPeriodicCycleCoords(double* coords, int num_tracks) {

  if (num_tracks != 5*getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the 2D Track periodic cycle coordinates "
               "since the TrackGenerator contains %d Tracks with %d coordinates"
               " but an array of length %d was input",
               getNum2DTracks(), 5*getNum2DTracks(), num_tracks);
  
  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
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


void TrackGenerator::retrieve2DReflectiveCycleCoords(double* coords, int num_tracks) {

  if (num_tracks != 5*getNum2DTracks())
    log_printf(ERROR, "Unable to retrieve the 2D Track reflective cycle coordinates "
               "since the TrackGenerator contains %d Tracks with %d coordinates"
               " but an array of length %d was input",
               getNum2DTracks(), 5*getNum2DTracks(), num_tracks);
  
  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
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


void TrackGenerator::retrieve3DPeriodicCycleCoords(double* coords, int num_tracks) {

  if (num_tracks != 7*getNum3DTracks())
    log_printf(ERROR, "Unable to retrieve the 3D Track periodic cycle coordinates "
               "since the TrackGenerator contains %d Tracks with %d coordinates"
               " but an array of length %d was input",
               getNum3DTracks(), 7*getNum3DTracks(), num_tracks);
  
  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
          coords[counter]   = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D_stack[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D_stack[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D_stack[a][i][p][z].getEnd()->getZ();
          coords[counter+6] = _tracks_3D_stack[a][i][p][z].getPeriodicCycleId();
          counter += 7;
        }
      }
    }
  }

  return;
}


void TrackGenerator::retrieve3DReflectiveCycleCoords(double* coords, int num_tracks) {

  if (num_tracks != 7*getNum3DTracks())
    log_printf(ERROR, "Unable to retrieve the 3D Track reflective cycle coordinates "
               "since the TrackGenerator contains %d Tracks with %d coordinates"
               " but an array of length %d was input",
               getNum3DTracks(), 7*getNum3DTracks(), num_tracks);
  
  /* Fill the array of coordinates with the Track start and end points */
  int counter = 0;
  
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
          coords[counter]   = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          coords[counter+1] = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          coords[counter+2] = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          coords[counter+3] = _tracks_3D_stack[a][i][p][z].getEnd()->getX();
          coords[counter+4] = _tracks_3D_stack[a][i][p][z].getEnd()->getY();
          coords[counter+5] = _tracks_3D_stack[a][i][p][z].getEnd()->getZ();
          coords[counter+6] = _tracks_3D_stack[a][i][p][z].getReflectiveCycleId();
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
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
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
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){

          x0    = _tracks_3D_stack[a][i][p][z].getStart()->getX();
          y0    = _tracks_3D_stack[a][i][p][z].getStart()->getY();
          z0    = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
          phi   = _tracks_3D_stack[a][i][p][z].getPhi();
          theta = _tracks_3D_stack[a][i][p][z].getTheta();
          
          segments = _tracks_3D_stack[a][i][p][z].getSegments();
          
          for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments(); s++) {
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
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      
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
    log_printf(ERROR, "The total height, width, and depth of the Geometry must "
               "be nonzero for Track generation. Create a CellFill which "
               "is filled by the entire geometry and bounded by XPlanes, "
               "YPlanes, and ZPlanes to enable the Geometry to determine the "
               "total width, height, and depth of the model.");

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

    /* Check periodic BCs for symmetry */
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
    if ((_geometry->getMinXBoundaryType() == PERIODIC ||
         _geometry->getMinYBoundaryType() == PERIODIC ||
         _geometry->getMinZBoundaryType() == PERIODIC ||
         _geometry->getMaxXBoundaryType() == PERIODIC ||
         _geometry->getMaxYBoundaryType() == PERIODIC ||
         _geometry->getMaxZBoundaryType() == PERIODIC) &&
        _track_generation_method != MODULAR_RAY_TRACING &&
        _track_generation_method != SIMPLIFIED_MODULAR_RAY_TRACING &&
        _solve_3D)
      log_printf(ERROR, "Cannot create tracks for a geometry containing a"
                 " periodic BC with a track generation method that is not"
                 " modular");

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
  _tracks_2D        = new Track2D*[_num_azim/2];
  _num_x            = new int[_num_azim/4];
  _num_y            = new int[_num_azim/4];
  _cycle_length     = new double[_num_azim/4];
  _azim_spacings    = new double[_num_azim/4];

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
  for (int a=0; a < _num_azim/2; a++){

    /* Allocate memory for the 2D tracks array */
    _tracks_2D[a] = new Track2D[getNumX(a) + getNumY(a)];

    /* Set the azimuthal angle */
    phi = _quadrature->getPhi(a);

    for (int i=0; i < getNumX(a) + getNumY(a); i++){

      Track2D* track = (&_tracks_2D[a][i]);
      track->setPhi(phi);
      track->setAzimIndex(a);
      track->setXYIndex(i);

      /* Set start point */
      if (a < _num_azim/4){
        if (i < getNumX(a))
          track->getStart()->setCoords(width - getDxEff(a) * (i + 0.5), 0.0);
        else
          track->getStart()->setCoords(0.0, getDyEff(a) * ((i-getNumX(a)) + 0.5));
      }
      else{
        if (i < getNumX(a))
          track->getStart()->setCoords(getDxEff(a) * (i + 0.5), 0.0);
        else
          track->getStart()->setCoords(width, getDyEff(a) *
                                       ((i-getNumX(a)) + 0.5));
      }

      /* Set end point */
      if (a < _num_azim/4){
        if (i < getNumY(a))
          track->getEnd()->setCoords(width, getDyEff(a) * (i + 0.5));
        else
          track->getEnd()->setCoords(width - getDxEff(a) * ((i-getNumY(a)) + 0.5),
                                     height);
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

  _contains_2D_tracks = true;
}


void TrackGenerator::initialize2DTrackReflections(){

  Track* track;
  int ac;
  
  /* Generate the 2D track cycles */
  for (int a=0; a < _num_azim/2; a++){
    ac = _num_azim/2 - a - 1;
    for (int i=0; i < getNumX(a) + getNumY(a); i++){

      /* Get current track */
      track = &_tracks_2D[a][i];

      /* Set the direction of periodic tracks */
      track->setPrdcFwdFwd(true);
      track->setPrdcBwdFwd(false);

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
      if (a < _num_azim/4){
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

  int id = 0;
  bool fwd;
  
  /* Set the periodic track cycle ids */
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      
      track = &_tracks_2D[a][i];
      
      if (track->getPeriodicCycleId() == -1) {
        while (track->getPeriodicCycleId() == -1){
          
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
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){

      track = &_tracks_2D[a][i];
      fwd = true;

      if (track->getReflectiveCycleId() == -1) {
        while (track->getReflectiveCycleId() == -1){

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
  
  for (int i=0; i < _num_azim/4; i++){
    _dz_eff[i]         = new double[_num_polar/2];
    _dl_eff[i]         = new double[_num_polar/2];
    _num_z[i]          = new int[_num_polar/2];
    _num_l[i]          = new int[_num_polar/2];
    _polar_spacings[i] = new double[_num_polar/2];
  }

  /* Allocate memory for tracks per stack */
  _tracks_per_stack = new int**[_num_azim/2];
  for (int a=0; a < _num_azim/2; a++){
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
    for (int j=0; j < _num_polar/2; j++){

      /* Compute the cosine weighted average angle */
      theta = _quadrature->getTheta(i, j);

      /* The number of intersections with xy (denoted "l") plane */
      if (_track_generation_method == GLOBAL_TRACKING){
        _num_l[i][j] = (int) (fabs(_cycle_length[i] * tan(M_PI_2 - theta)
                                   * sin(theta) / _polar_spacing)) + 1;
      }
      else if (_track_generation_method == MODULAR_RAY_TRACING){
        _num_l[i][j] = 2 * (int)
          (fabs(_cycle_length[i] * 0.5 * tan(M_PI_2 - theta)
                * sin(theta) / _polar_spacing) + 1);
      }
      else if (_track_generation_method == SIMPLIFIED_MODULAR_RAY_TRACING){
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
    for (int c = 0; c < _cycles_per_azim[a]; c++){
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
    if (create_tracks) {
      _tracks_3D_stack = new Track3D***[_num_azim/2];
      for (int a=0; a < _num_azim/2; a++){
        _tracks_3D_stack[a] = new Track3D**[getNumX(a) + getNumY(a)];
        for (int i=0; i < getNumX(a) + getNumY(a); i++){
          _tracks_3D_stack[a][i] = new Track3D*[_num_polar];
          for (int p=0; p < _num_polar; p++) {
            _tracks_3D_stack[a][i][p] = new Track3D[_tracks_per_stack[a][i][p]];
            _tracks_per_stack[a][i][p] = 0;
          }
        }
      }

      /* Allocate memory for 3D tracks cycles */
      for (int a = 0; a < _num_azim/4; a++) {
        _tracks_3D_cycle[a] = new Track3D****[_cycles_per_azim[a]];
        for (int c = 0; c < _cycles_per_azim[a]; c++){
          _tracks_3D_cycle[a][c] = new Track3D***[_num_polar];
          for (int p=0; p < _num_polar; p++){
            _tracks_3D_cycle[a][c][p] =
              new Track3D**[getNumZ(a,p) + getNumL(a,p)];
            for (int i=0; i < getNumZ(a,p) + getNumL(a,p); i++)
              _tracks_3D_cycle[a][c][p][i] =
                new Track3D*[_tracks_per_train[a][c][p][i]];
          }
        }
      }
    }

    /* Loop over 3D track cycles */
    for (int a = 0; a < _num_azim/4; a++) {
      for (int c = 0; c < _cycles_per_azim[a]; c++){
        
        /* Loop over polar angles < PI/2 */
        for (int p=0; p < _num_polar/2; p++){
          
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
            }
            else{
              l_start = 0.0;
              track_2D = getTrack2DByCycle(a, c, 0);
              x1 = track_2D->getStart()->getX();
              y1 = track_2D->getStart()->getY();
              z1 = _dz_eff[a][p] * (i - _num_l[a][p] + 0.5);
            }
            
            /* Get the end point */
            if (i < _num_z[a][p]){
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
        for (int p=0; p < _num_polar/2; p++){
          
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
          for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++){

            /* Get the starting point */
            if (i < _num_z[a][p]){
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
            if (i < _num_l[a][p]){
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

}


void TrackGenerator::initialize3DTrackReflections(){

  int pc;
  int ai, xi, pi, zi, xp, zp, pp, ci;
  bool periodic;

  /* Check to see if periodic connecting tracks need to be set */ 
  if (_geometry->getMinXBoundaryType() == PERIODIC ||
      _geometry->getMinYBoundaryType() == PERIODIC ||
      _geometry->getMinZBoundaryType() == PERIODIC ||
      _geometry->getMaxXBoundaryType() == PERIODIC ||
      _geometry->getMaxYBoundaryType() == PERIODIC ||
      _geometry->getMaxZBoundaryType() == PERIODIC)
    periodic = true;
  else
    periodic = false;

  /* Set reflective tracks and periodic top and bottom indices */
  for (int a = 0; a < _num_azim/4; a++) {
    for (int c = 0; c < _cycles_per_azim[a]; c++){

      /* Loop over polar angles < PI/2 */
      for (int p=0; p < _num_polar/2; p++){

        /* Set the complementary polar angle */
        pc = _num_polar-p-1;

        for (int i=0; i < _num_l[a][p] + _num_z[a][p]; i++){
          for (int t=0; t < _tracks_per_train[a][c][p][i]; t++){

            if (t == _tracks_per_train[a][c][p][i]-1){

              /* SURFACE_Y_MIN */
              if (i < _num_z[a][p]){
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMinYBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][p][_num_l[a][p] + i][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][p][_num_l[a][p] + i][0]);

                  /* PERIODIC */
                  if (periodic) {
                    ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][_num_l[a][p] + i][0]->getZIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
                else {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (_tracks_3D_cycle[a][c][p][_num_l[a][p] + i][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][p][_num_l[a][p] + i][0]);
                                    
                  /* PERIODIC */
                  if (periodic) {
                    ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][_num_l[a][p] + i][0]->getZIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
              }

              /* SURFACE_Z_MAX */
              else{
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMaxZBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]);

                  /* PERIODIC */
                  if (periodic) {
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                    _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcFwd
                      (_tracks_3D_cycle[a][c][p][i - _num_z[a][p]][0]);
                  }
                }
                else{

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMaxZBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc]
                     [_num_l[a][p] + 2*_num_z[a][p] - i - 1][0]);

                  /* PERIODIC */
                  if (periodic) {
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcBwd
                      (_tracks_3D_cycle[a][c][p][i - _num_z[a][p]][0]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURAFCE_Y_MAX */
            else{
              if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCFwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t+1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCBwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t+1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
            }
            
            if (t == 0){
              
              /* SURFACE_Z_MIN */
              if (i < _num_l[a][p]){
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMinZBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (!_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]
                     ->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]);

                  /* PERIODIC */
                  if (periodic) {
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcBwd
                      (_tracks_3D_cycle[a][c][p][i + _num_z[a][p]]
                       [_tracks_per_train[a][c][p][i + _num_z[a][p]]-1]);
                  }
                }
                else {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMinZBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (!_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]
                     ->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][p] - i - 1]
                     [_tracks_per_train[a][c][pc][_num_l[a][p] - i - 1] - 1]);

                  /* PERIODIC */
                  if (periodic) {
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                    _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcFwd
                      (_tracks_3D_cycle[a][c][p][i + _num_z[a][p]]
                       [_tracks_per_train[a][c][p][i + _num_z[a][p]]-1]);
                  }
                }
              }

              /* SURFACE_Y_MIN */
              else{
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMinYBoundaryType());
                  
                  /* REFFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (!_tracks_3D_cycle[a][c][p][i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][p][i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]);

                  /* PERIODIC */
                  if (periodic) {
                    ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][i - _num_l[a][p]]
                      [_tracks_per_train[a][c][p]
                       [i - _num_l[a][p]] - 1]->getZIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
                else {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMinYBoundaryType());
                  
                  /* REFFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][p][i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][p][i - _num_l[a][p]]
                     [_tracks_per_train[a][c][p]
                      [i - _num_l[a][p]] - 1]);

                  /* PERIODIC */
                  if (periodic) {
                    ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    xp = getTrack2DByCycle(a, c, xi)->getTrackPrdcFwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][i - _num_l[a][p]]
                      [_tracks_per_train[a][c][p]
                       [i - _num_l[a][p]] - 1]->getZIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, SURAFCE_Y_MAX */
            else{
              if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCBwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                  (!_tracks_3D_cycle[a][c][p][i][t-1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                  (_tracks_3D_cycle[a][c][p][i][t-1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t-1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCFwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                  (!_tracks_3D_cycle[a][c][p][i][t-1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                  (_tracks_3D_cycle[a][c][p][i][t-1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t-1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
            }
          }
        }
      }

      /* Loop over polar angles > PI/2 */
      for (int p = _num_polar/2; p < _num_polar; p++){
        pc = _num_polar-p-1;
        for (int i=0; i < _num_l[a][pc] + _num_z[a][pc]; i++){
          for (int t=0; t < _tracks_per_train[a][c][p][i]; t++){
            
            if (t == _tracks_per_train[a][c][p][i]-1){
              
              /* SURFACE_Z_MIN */
              if (i < _num_l[a][pc]){
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMinZBoundaryType());

                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]);
                  
                  /* PERIODIC */
                 if (periodic) {
                   _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                   _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcFwd
                     (_tracks_3D_cycle[a][c][p][i + _num_z[a][pc]][0]);
                 }
                }
                else {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMinZBoundaryType());

                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc][_num_l[a][pc] - i - 1][0]);
                  
                  /* PERIODIC */
                  if (periodic) {
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcBwd
                      (_tracks_3D_cycle[a][c][p][i + _num_z[a][pc]][0]);
                  }
                }
              }

              /* SURFACE_Y_MIN */
              else{
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMinYBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][p][i - _num_l[a][pc]][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][p][i - _num_l[a][pc]][0]);

                  /* PERIODIC */
                  if (periodic) {
                    ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][i - _num_l[a][pc]][0]->getZIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
                else {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMinYBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (_tracks_3D_cycle[a][c][p][i - _num_l[a][pc]][0]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][p][i - _num_l[a][pc]][0]);

                  /* PERIODIC */
                  if (periodic) {
                    ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][i - _num_l[a][pc]][0]->getZIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
              } 
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, or SURFACE_Y_MAX */
            else{
              if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {
                
                _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCFwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t+1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCBwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                  (_tracks_3D_cycle[a][c][p][i][t+1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t+1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
            }
            
            if (t == 0){

              /* SURFACE_Y_MIN */
              if (i < _num_z[a][pc]){
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMinYBoundaryType());

                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (!_tracks_3D_cycle[a][c][p][_num_l[a][pc] + i]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] + i] - 1]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][p][_num_l[a][pc] + i]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] + i] - 1]);

                  /* PERIODIC */
                  if (periodic) {
                    ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][_num_l[a][pc] + i]
                      [_tracks_per_train[a][c][p]
                       [_num_l[a][pc] + i] - 1]->getZIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
                else {
                  
                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMinYBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (_tracks_3D_cycle[a][c][p][_num_l[a][pc] + _num_z[a][pc] - i - 1]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][p][_num_l[a][pc] + _num_z[a][pc] - i - 1]
                     [_tracks_per_train[a][c][p]
                      [_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]);
                  
                  /* PERIODIC */
                  if (periodic) {ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                    xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                    pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                    zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                    ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                    xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()->getXYIndex();
                    zp = _tracks_3D_cycle[a][c][p][_num_l[a][pc] + _num_z[a][pc] - i - 1]
                      [_tracks_per_train[a][c][p]
                       [_num_l[a][pc] +_num_z[a][pc] - i - 1] - 1]->getZIndex();
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                    _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcFwd
                      (&_tracks_3D_stack[ai][xp][pi][zp]);
                  }
                }
              }

              /* SURFACE_Z_MAX */
              else{
                if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                    (_geometry->getMaxZBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                    (!_tracks_3D_cycle[a][c][pc][2*_num_z[a][pc] + _num_l[a][pc] - i - 1]
                     [_tracks_per_train[a][c][pc]
                      [2*_num_z[a][pc] + _num_l[a][pc] - i - 1] - 1]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                    (_tracks_3D_cycle[a][c][pc][2*_num_z[a][pc] + _num_l[a][pc] - i - 1]
                     [_tracks_per_train[a][c][pc]
                      [2*_num_z[a][pc] + _num_l[a][pc] - i - 1] - 1]);

                  /* PERIODIC */
                  if (periodic) {
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                    _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcBwd
                      (_tracks_3D_cycle[a][c][p][i - _num_z[a][pc]]
                       [_tracks_per_train[a][c][p][i - _num_z[a][pc]]-1]);
                  }
                }
                else {

                  _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                    (_geometry->getMaxZBoundaryType());
                  
                  /* REFLECTIVE */
                  _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                    (!_tracks_3D_cycle[a][c][pc][2*_num_z[a][pc] + _num_l[a][pc] - i - 1]
                     [_tracks_per_train[a][c][pc]
                      [2*_num_z[a][pc] + _num_l[a][pc] - i - 1] - 1]->getCycleFwd());
                  _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                    (_tracks_3D_cycle[a][c][pc][2*_num_z[a][pc] + _num_l[a][pc] - i - 1]
                     [_tracks_per_train[a][c][pc]
                      [2*_num_z[a][pc] + _num_l[a][pc] - i - 1] - 1]);

                  /* PERIODIC */
                  if (periodic) {
                    _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
                    _tracks_3D_cycle[a][c][p][i][t]->setTrackPrdcFwd
                      (_tracks_3D_cycle[a][c][p][i - _num_z[a][pc]]
                       [_tracks_per_train[a][c][p][i - _num_z[a][pc]]-1]);
                  }
                }
              }
            }

            /* SURFACE_X_MIN, SURFACE_X_MAX, SURFACE_Y_MIN, or SURFACE_Y_MAX */
            else{
              if (_tracks_3D_cycle[a][c][p][i][t]->getCycleFwd()) {

                _tracks_3D_cycle[a][c][p][i][t]->setBCBwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCBwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflBwdFwd
                  (!_tracks_3D_cycle[a][c][p][i][t-1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflBwd
                  (_tracks_3D_cycle[a][c][p][i][t-1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcBwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t-1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcBwdFwd(false);
                  _tracks_3D_stack[ai][xi][pi][zi].setTrackPrdcBwd
                    (&_tracks_3D_stack[ai][xp][pi][zp]);
                }
              }
              else {

                _tracks_3D_cycle[a][c][p][i][t]->setBCFwd
                  (getTrack2DByCycle(a, c, _tracks_3D_cycle[a][c][p][i][t]
                                     ->getCycleTrackIndex())->getBCFwd());
                
                /* REFLECTIVE */
                _tracks_3D_cycle[a][c][p][i][t]->setReflFwdFwd
                  (!_tracks_3D_cycle[a][c][p][i][t-1]->getCycleFwd());
                _tracks_3D_cycle[a][c][p][i][t]->setTrackReflFwd
                  (_tracks_3D_cycle[a][c][p][i][t-1]);

                /* PERIODIC */
                if (periodic) {
                  ai = _tracks_3D_cycle[a][c][p][i][t]->getAzimIndex();
                  xi = _tracks_3D_cycle[a][c][p][i][t]->getXYIndex();
                  pi = _tracks_3D_cycle[a][c][p][i][t]->getPolarIndex();
                  zi = _tracks_3D_cycle[a][c][p][i][t]->getZIndex();
                  ci = _tracks_3D_cycle[a][c][p][i][t]->getCycleTrackIndex();
                  xp = getTrack2DByCycle(a, c, ci)->getTrackPrdcFwd()->getXYIndex();
                  zp = _tracks_3D_cycle[a][c][p][i][t-1]->getZIndex();
                  _tracks_3D_cycle[a][c][p][i][t]->setPrdcFwdFwd(true);
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

  int id = 0;
  Track* track;
  bool fwd;
  
  /* Set the periodic track cycle ids */
  if (periodic) {
    for (int a=0; a < _num_azim/2; a++){
      for (int i=0; i < getNumX(a) + getNumY(a); i++){
        for (int p=0; p < _num_polar; p++){
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
            
            track = &_tracks_3D_stack[a][i][p][z];
            
            if (track->getPeriodicCycleId() == -1) {
              while (track->getPeriodicCycleId() == -1){
                
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
  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++){
          
          track = &_tracks_3D_stack[a][i][p][z];
          fwd = true;
          
          if (track->getReflectiveCycleId() == -1) {
            while (track->getReflectiveCycleId() == -1){
              
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


 
void TrackGenerator::decomposeLZTrack(Track3D* track, double l_start,
                                      double l_end, int azim, int cycle,
                                      int polar, int lz_index,
                                      bool create_tracks){

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
  for (int i=0; i < _tracks_per_cycle[azim]; i++){
    
    track_2d = getTrack2DByCycle(azim, cycle, i);
    
    if (l_end < length_sum + track_2d->getLength() + 1.05 * nudge){
      last_stack = i;
      break;
    }    
    else{
      length_sum += track_2d->getLength();
    }
  }
  
  length_sum = 0.0;
  
  /* Find the first cycle index */
  for (int i=0; i < _tracks_per_cycle[azim]; i++){

    track_2d = getTrack2DByCycle(azim, cycle, i);
    
    if (l_start < length_sum + track_2d->getLength() - 1.05 * nudge){
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
    for (int i=first_stack; i <= last_stack; i++){
      
      /* Get the 2D track associated with this 3D track */
      track_2d = getTrack2DByCycle(azim, cycle, i);
      fwd = getCycleDirection(azim, cycle, i);
      ai = track_2d->getAzimIndex();
      ti = track_2d->getXYIndex();

      /* Set the start and end points */
      if (i == first_stack){
        x1 = track->getStart()->getX();
        y1 = track->getStart()->getY();
        z1 = track->getStart()->getZ();
      }
      else{
        x1 = x2;
        y1 = y2;
        z1 = z2;
      }

      if (i == last_stack){
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
  for (int i=first_stack; i <= last_stack; i++){
    track_2d = getTrack2DByCycle(azim, cycle, i);
    fwd = getCycleDirection(azim, cycle, i);
    ti = track_2d->getXYIndex();
    ai = track_2d->getAzimIndex();

    /* Computing the endpoint z value */
    if (i == first_stack)
      z1 = track->getStart()->getZ();
    else
      z1 = z2;
    
    if (i == last_stack){
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
 * @brief TODO
 */
void TrackGenerator::initializeExtrudedTracks() {

  ExtrudedTrack* _extruded_tracks = new ExtrudedTrack[_num_2D_tracks];

  size_t index = 0;
  for (int a=0; a < _num_azim/2; a++) {
    for (int i=0; i < getNumX(a) + getNumY(a); i++) {

      _extruded_tracks[index]._azim_index = a;
      _extruded_tracks[index]._track_index = i;
      _extruded_tracks[index]._num_segments = 0;
      _extruded_tracks[index]._track_2D = &_tracks_2D[a][i];
      
      index++;
    }
  }
  _contains_extruded_tracks = true;
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
  int num_x, num_y, azim_period;

  /* Recalibrate the tracks to the origin and set the uid. Note that the 
   * loop structure is unconventional in order to preserve an increasing 
   * track uid value in the Solver's tracks array. The tracks array is 
   * oriented with this loop structure in order to maintain reproducability
   * for parallel runs.
   */

  /* Loop over azim reflective halfspaces */
  for (int azim_refl_half=0; azim_refl_half < 2; azim_refl_half++) {

    /* Loop over azim periodic halfspaces */
    for (int azim_prdc_half=0; azim_prdc_half < 2; azim_prdc_half++) {

      /* Loop over all tracks in current azim reflective halfspace */
      for (int a=azim_refl_half*_num_azim/4;
           a < (azim_refl_half+1)*_num_azim/4; a++) {
        num_x = getNumX(a);
        num_y = getNumY(a);
        azim_period = std::min(num_x, num_y);

        /* Loop over all xy tracks */
        for (int i=0; i < num_x + num_y; i++) {

          /* Check if track is in current azim periodic halfspace */
          if ((i / azim_period) % 2 == azim_prdc_half) {
            _tracks_2D[a][i].setUid(uid);
            uid++;
            
            double x0 = _tracks_2D[a][i].getStart()->getX();
            double y0 = _tracks_2D[a][i].getStart()->getY();
            double x1 = _tracks_2D[a][i].getEnd()->getX();
            double y1 = _tracks_2D[a][i].getEnd()->getY();
            double new_x0 = x0 + _geometry->getMinX();
            double new_y0 = y0 + _geometry->getMinY();
            double new_x1 = x1 + _geometry->getMinX();
            double new_y1 = y1 + _geometry->getMinY();
            
            _tracks_2D[a][i].setCoords(new_x0, new_y0, new_x1, new_y1);
            _num_2D_tracks++;
          }
        }
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
  int num_x, num_y, num_l, num_z;
  int polar_period, azim_period;
  
  /* Recalibrate the tracks to the origin and set the uid. Note that the 
   * loop structure is unconventional in order to preserve an increasing 
   * track uid value in the Solver's tracks array. The tracks array is 
   * oriented with this loop structure in order to maintain reproducability
   * for parallel runs.
   */

  /* Loop over azim reflective halfspaces */
  for (int azim_refl_half=0; azim_refl_half < 2; azim_refl_half++) {

    /* Loop over azim periodic halfspaces */
    for (int azim_prdc_half=0; azim_prdc_half < 2; azim_prdc_half++) {

      /* Loop over polar reflective halfspaces */
      for (int polar_refl_half=0; polar_refl_half < 2; polar_refl_half++) {

        /* Loop over polar periodic halfspaces */
        for (int polar_prdc_half=0; polar_prdc_half < 2; polar_prdc_half++) {

          /* Loop over azim angles in current azim reflective halfspace */
          for (int a=azim_refl_half*_num_azim/4;
               a < (azim_refl_half+1)*_num_azim/4; a++) {
            num_x = getNumX(a);
            num_y = getNumY(a);
            azim_period = std::min(num_x, num_y);

            /* Loop over all xy tracks */
            for (int i=0; i < num_x + num_y; i++) {

              /* Check if track is in current azim periodic halfspace */
              if ((i / azim_period) % 2 == azim_prdc_half) {

                /* Loop over polar angles in current polar reflective
                 * halfspace */
                for (int p=polar_refl_half*_num_polar/2;
                     p < (polar_refl_half+1)*_num_polar/2; p++) {
                  
                  num_l = getNumL(a, p);
                  num_z = getNumZ(a, p);
                  polar_period = std::min(num_l, num_z);
                  
                  /* Loop over all tracks in z stack */
                  for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

                    /* Check if track is in current polar periodic halfspace */
                    if ((_tracks_3D_stack[a][i][p][z].getLZIndex()
                         / polar_period) % 2 == polar_prdc_half) {

                      _tracks_3D_stack[a][i][p][z].setUid(uid);
                      uid++;
                      
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
                      _num_3D_tracks++;            
                    }
                  }
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
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize2D() {

  log_printf(NORMAL, "Ray tracing for 2D track segmentation...");

  int tracks_segmented = 0;
  
  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++){
    log_printf(NORMAL, "segmenting 2D tracks - Percent complete: %5.2f %%",
               double(tracks_segmented) / _num_2D_tracks * 100.0);
    #pragma omp parallel for
    for (int i=0; i < getNumX(a) + getNumY(a); i++)
      _geometry->segmentize2D(&_tracks_2D[a][i], _z_level);

    tracks_segmented += getNumX(a) + getNumY(a);
  }

  _num_2D_segments = 0;

  for (int a=0; a < _num_azim/2; a++){
    for (int i=0; i < getNumX(a) + getNumY(a); i++)
      _num_2D_segments += _tracks_2D[a][i].getNumSegments();
  }

  _geometry->initializeFSRVectors();
  _contains_2D_segments = true;
  
  return;
}

/**
 * @brief TODO YADA YADA YADA segmentize
 */
void TrackGenerator::segmentizeExtruded() {

  log_printf(NORMAL, "Ray tracing for axially extruded track segmentation...");

  /* Loop over all extruded Tracks */
  #pragma omp parallel for
  for (int index=0; index < _num_2D_tracks; index++) {

    log_printf(NORMAL, "segmenting axially extruded tracks - Percent complete:"
       " %5.2f %%", double(index) / _num_2D_tracks * 100.0);
    
    _geometry->segmentizeExtruded(&_extruded_tracks[index]);
  }

  _geometry->initializeFSRVectors();
  _contains_extruded_segments = true;
  
  return;
}


/**
 * @brief Generate segments for each Track across the Geometry.
 */
void TrackGenerator::segmentize3D() {

  log_printf(NORMAL, "Ray tracing for 3D track segmentation...");

  int tracks_segmented = 0;
  
  /* Loop over all Tracks */
  for (int a=0; a < _num_azim/2; a++){
    
    log_printf(NORMAL, "segmenting 3D tracks - Percent complete: %5.2f %%",
               double(tracks_segmented) / _num_3D_tracks * 100.0);
    
    #pragma omp parallel for
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          _geometry->segmentize3D(&_tracks_3D_stack[a][i][p][z]);
      }
    }
    
    for (int i=0; i < getNumX(a) + getNumY(a); i++){
      for (int p=0; p < _num_polar; p++){
        for (int z=0; z < _tracks_per_stack[a][i][p]; z++)
          tracks_segmented++;
      }
    }
  }

  _geometry->initializeFSRVectors();
  _contains_3D_segments = true;
  
  return;
}


/* Convert a length traveled along a track cycle to value in x */
double TrackGenerator::convertLtoX(double l, int azim, int cycle){

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to X since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);
  
  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++){
    if (l <= length_sum + getTrack2DByCycle(azim, cycle, i)->getLength()){
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


/* Convert a length traveled along a track cycle to value in y */
double TrackGenerator::convertLtoY(double l, int azim, int cycle){

  if (l > _cycle_length[azim] || l < 0.0)
    log_printf(ERROR, "Unable to convert L to Y since l is not within "
               "(0, %f), l: %f", _cycle_length[azim], l);
  
  double length_sum = 0.0;
  int track_index = _tracks_per_cycle[azim] - 1;

  for (int i=0; i < _tracks_per_cycle[azim]; i++){
    if (l <= length_sum + getTrack2DByCycle(azim, cycle, i)->getLength()){
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

  if (_solve_3D){

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

    if (_geometry->getCmfd() != NULL){
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
      //FIXME
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
        if (cmfd != NULL){
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

  int num_cuts, min_num_cuts;
  segment* curr_segment;

  FP_PRECISION length, tau;
  int fsr_id;
  Material* material;
  FP_PRECISION* sigma_t;
  int num_groups;

  // FIXME
  if (_solve_3D){
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {
            for (int s=0; s < _tracks_3D_stack[a][i][p][z].getNumSegments(); s++) {
                
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
                _tracks_3D_stack[a][i][p][z].insertSegment(s+k+1, new_segment);
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
        for (int s=0; s < _tracks_2D[a][i].getNumSegments(); s++){
          
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
    
    for (int a=0; a < _num_azim/2; a++) {
      for (int i=0; i < getNumX(a) + getNumY(a); i++) {
        for (int p=0; p < _num_polar; p++) {
          for (int z=0; z < _tracks_per_stack[a][i][p]; z++) {

            int num_segments = _tracks_3D_stack[a][i][p][z].getNumSegments();
            segment* segments = _tracks_3D_stack[a][i][p][z].getSegments();
            double xx = _tracks_3D_stack[a][i][p][z].getStart()->getX();
            double yy = _tracks_3D_stack[a][i][p][z].getStart()->getY();
            double zz = _tracks_3D_stack[a][i][p][z].getStart()->getZ();
            double phi = _tracks_3D_stack[a][i][p][z].getPhi();
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
  for (int r=0; r < num_FSRs; r++){
    _geometry->setFSRCentroid(r, centroids[r]);
  }
  
  /* Delete temporary array of centroids and FSR volumes */
  delete [] FSR_volumes;
}


void TrackGenerator::setTrackGenerationMethod(int method) {

  if (method != GLOBAL_TRACKING && method != MODULAR_RAY_TRACING &&
      method != SIMPLIFIED_MODULAR_RAY_TRACING)
    log_printf(ERROR, "Unable to set Track Generation Method to %i. Valid"
               " methods include GLOBAL_TRACKING, MODULAR_RAY_TRACING, "
               "and SIMPLIFIED_MODULAR_RAY_TRACING", method);

  _track_generation_method = method;
}


int TrackGenerator::getTrackGenerationMethod() {
  return _track_generation_method;
}


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
