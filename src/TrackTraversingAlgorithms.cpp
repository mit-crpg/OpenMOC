#include "TrackTraversingAlgorithms.h"
#include "CPUSolver.h"
#include "Quadrature.h"

/**
 * @brief Constructor for MaxOpticalLength calls the TraverseSegments
 *        constructor and sets the max optical path length to zero
 * @param track_generator The TrackGenerator to pull tracking information from
 */
MaxOpticalLength::MaxOpticalLength(TrackGenerator* track_generator)
                                 : TraverseSegments(track_generator) {
  _max_tau = 0;
}


/**
 * @brief Determines the maximum optical path length for the TrackGenerator
 *        provided during construction
 * @details The maximum optical path length is initialized to infinity for
 *          segmentation within the TrackGenerator and then SegmentationKernels
 *          are allocated to store temporary segmnents. Tracks are traversed
 *          and onTrack(...) is applied, calculating a maximum optical length
 *          and setting it on the TrackGenerator.
*/
void MaxOpticalLength::execute() {
  FP_PRECISION infinity = std::numeric_limits<FP_PRECISION>::max();
  _track_generator->setMaxOpticalLength(infinity);
#pragma omp parallel
  {
    MOCKernel* kernel = getKernel<SegmentationKernel>();
    loopOverTracks(kernel);
  }
  _track_generator->setMaxOpticalLength(_max_tau);
}


/**
 * @brief Calculates the optical path length for the provided segments and
 *        updates the maximum optical path length if necessary
 * @param track The track associated with the segments
 * @param segments The segments for which the optical path length is calculated
 */
void MaxOpticalLength::onTrack(Track* track, segment* segments) {
  for (int s=0; s < track->getNumSegments(); s++) {
    FP_PRECISION length = segments[s]._length;
    Material* material = segments[s]._material;
    FP_PRECISION* sigma_t = material->getSigmaT();

    for (int e=0; e < material->getNumEnergyGroups(); e++) {
      FP_PRECISION tau = length*sigma_t[e];
      if (tau > _max_tau) {
#pragma omp critical
        _max_tau = std::max(_max_tau, tau);
      }
    }
  }
}


/**
 * @brief Constructor for SegmentCounter calls the TraverseSegments
 *        constructor and sets the max number of segments per Track to zero
 * @param track_generator The TrackGenerator to pull tracking information from
 */
SegmentCounter::SegmentCounter(TrackGenerator* track_generator)
                               : TraverseSegments(track_generator) {
  _max_num_segments = 0;
  _total_num_segments = 0;
  _count_total_segments = false;
  _total_segments_counted = false;
}


/**
 * @brief Determines the maximum number of segments per Track
 * @details CounterKernels are initialized to count segments along each Track.
 *          Then Tracks are traversed, saving the maximum number of segments
 *          and setting the corresponding parameter on the TrackGenerator.
*/
void SegmentCounter::execute() {
  _total_segments_counted = false;
#pragma omp parallel
  {
    MOCKernel* kernel = getKernel<CounterKernel>();
    loopOverTracks(kernel);
  }
  _track_generator->setMaxNumSegments(_max_num_segments);
  if (_count_total_segments)
    _total_segments_counted = true;
}


//FIXME
void SegmentCounter::countTotalNumSegments() {
  _count_total_segments = true;
}


//FIXME
long SegmentCounter::getTotalNumSegments() {
  if (!_total_segments_counted)
    log_printf(ERROR, "The total number of segments have not been counted. "
               "The SegmentCounter was not instructed to count segments "
               "before execution");
  return _total_num_segments;
}


/**
 * @brief Updates the maximum number of segments per Track if a Track with a
 *        larger number of segments is observed
 * @param track The Track whose segments are counted
 * @param segments The segments associated with the Track
 */
void SegmentCounter::onTrack(Track* track, segment* segments) {
  if (track->getNumSegments() > _max_num_segments) {
#pragma omp critical
    _max_num_segments = std::max(_max_num_segments, track->getNumSegments());
  }
  if (_count_total_segments) {
#pragma omp critical
    _total_num_segments += track->getNumSegments();
  }
}


/**
 * @brief Constructor for SegmentSplitter calls the TraverseSegments
 *        constructor
 * @param track_generator The TrackGenerator to pull tracking information from
 */
SegmentSplitter::SegmentSplitter(TrackGenerator* track_generator)
                               : TraverseSegments(track_generator) {
}


/**
 * @brief Splits segments stored explicity along each Track
 * @details No MOCKernels are initialized for this function.
 */
void SegmentSplitter::execute() {
#pragma omp parallel
  {
    loopOverTracks(NULL);
  }
}


/**
 * @brief Segments for the provided Track are split so that no segment has a
 *        larger optical path length than the maximum optical path length
 * @param track The Track whose segments are potentially split
 * @param segments The segments associated with the Track
 */
void SegmentSplitter::onTrack(Track* track, segment* segments) {

  /* Get the max optical length from the TrackGenerator */
  FP_PRECISION max_optical_length =
    _track_generator->retrieveMaxOpticalLength();

  /* Extract data from this segment to compute its optical
   * length */
  for (int s = 0; s < track->getNumSegments(); s++) {
    segment* curr_segment = track->getSegment(s);
    Material* material = curr_segment->_material;
    FP_PRECISION length = curr_segment->_length;
    int fsr_id = curr_segment->_region_id;

    /* Compute number of segments to split this segment into */
    int min_num_cuts = 1;
    int num_groups = material->getNumEnergyGroups();
    FP_PRECISION* sigma_t = material->getSigmaT();

    for (int g=0; g < num_groups; g++) {
      FP_PRECISION tau = length * sigma_t[g];
      int num_cuts = ceil(tau / max_optical_length);
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
      track->insertSegment(s+k+1, new_segment);
    }

    /* Remove the original segment from the Track */
    track->removeSegment(s);
  }
}


/**
 * @brief Constructor for SegmentSplitter calls the TraverseSegments
 *        constructor
 * @param track_generator The TrackGenerator to pull tracking information from
 */
VolumeCalculator::VolumeCalculator(TrackGenerator* track_generator)
                                  : TraverseSegments(track_generator) {
}


/**
 * @brief FSR volumes are calculated and saved in the TrackGenerator's FSR
 *        volumes buffer
 * @details VolumeKernels are created and used to loop over all segments and
 *          tally each segments contribution to FSR volumes.
 */
void VolumeCalculator::execute() {
#pragma omp parallel
  {
    MOCKernel* kernel = getKernel<VolumeKernel>();
    loopOverTracks(kernel);
  }
}


/**
 * @brief No functionality is applied for each Track during the execution of
 *        the VolumeCalculator
 * @param track The current Track
 * @param segments The segments associated with the Track
 */
void VolumeCalculator::onTrack(Track* track, segment* segments) {
}


/**
 * @brief Constructor for CentroidGenerator calls the TraverseSegments
 *        constructor and imports refernces to both the FSR volumes and FSR
 *        locks arrays
 * @param track_generator The TrackGenerator to pull tracking information from
 */
CentroidGenerator::CentroidGenerator(TrackGenerator* track_generator)
                                  : TraverseSegments(track_generator) {

  _FSR_volumes = track_generator->getFSRVolumesBuffer();
  _FSR_locks = track_generator->getFSRLocks();
  _quadrature = track_generator->getQuadrature();

  int num_threads = omp_get_max_threads();
  _starting_points = new Point*[num_threads];

  int num_rows = 1;
  if (_track_generator_3D != NULL)
    num_rows = _track_generator_3D->getMaxNumTracksPerStack();
  for (int i=0; i < num_threads; i++)
    _starting_points[i] = new Point[num_rows];
}


//FIXME
CentroidGenerator::~CentroidGenerator() {
  int num_threads = omp_get_max_threads();
  for (int i=0; i < num_threads; i++)
    delete [] _starting_points[i];
  delete [] _starting_points;
}


/**
 * @brief Calculates the centroid of every FSR
 * @details SegmentationKernels are created to temporarily save segments for
 *          on-the-fly methods. Then on each segment, onTrack(...) calculates
 *          the contribution to each FSR centroid and saves the centroids in
 *          the centroids array provided by setCentroids(...).
 */
void CentroidGenerator::execute() {
#pragma omp parallel
  {
    MOCKernel* kernel = getKernel<SegmentationKernel>();
    loopOverTracks(kernel);
  }
}


/**
 * @brief Specifies an array to save calculated FSR centroids
 * @brief centroids The array of FSR centroids pointers
 */
void CentroidGenerator::setCentroids(Point** centroids) {
  _centroids = centroids;
}


/**
 * @brief Centroid contributions are calculated for every segment in the Track
 * @param track The Track associated with the segments
 * @param segments The segments whose contributions are added to the centroids
 */
void CentroidGenerator::onTrack(Track* track, segment* segments) {

  /* Extract common information from the Track */
  Point* start = track->getStart();
  int azim_index = track->getAzimIndex();
  double phi = track->getPhi();

  /* Compute the Track azimuthal weight */
  double wgt = _quadrature->getAzimSpacing(azim_index)
      * _quadrature->getAzimWeight(azim_index);

  //FIXME
  int tid = omp_get_thread_num();
  _starting_points[tid][0].copyCoords(track->getStart());
  Track3D* current_stack = NULL;

  /* Get polar angles depending on the dimensionality */
  double sin_theta = 1;
  double cos_theta = 0;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL) {
    double theta = track_3D->getTheta();
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    int polar_index = track_3D->getPolarIndex();
    wgt *= _quadrature->getPolarSpacing(azim_index, polar_index)
        *_quadrature->getPolarWeight(azim_index, polar_index);

    if (_segment_formation == OTF_STACKS) {
      int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
      current_stack = _track_generator_3D->getTemporary3DTracks(tid);
      int xy_index = track->getXYIndex();
      int stack_size = tracks_per_stack[azim_index][xy_index][polar_index];
      for (int i=1; i < stack_size; i++)
        _starting_points[tid][i].copyCoords(current_stack[i].getStart());
    }
  }

  /* Pre-compute azimuthal angles for efficiency */
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);

  /* Loop over segments to accumlate contribution to centroids */
  for (int s=0; s < track->getNumSegments(); s++) {

    segment* curr_segment = &segments[s];
    int fsr = curr_segment->_region_id;
    int track_idx = curr_segment->_track_idx;

    /* Extract information */
    double volume = _FSR_volumes[fsr];
    double x = _starting_points[tid][track_idx].getX();
    double y = _starting_points[tid][track_idx].getY();
    double z = _starting_points[tid][track_idx].getZ();

    /* Set the lock for this FSR */
    omp_set_lock(&_FSR_locks[fsr]);

    _centroids[fsr]->
        setX(_centroids[fsr]->getX() + wgt *
        (x + cos_phi * sin_theta * curr_segment->_length / 2.0)
        * curr_segment->_length / _FSR_volumes[fsr]);

    _centroids[fsr]->
        setY(_centroids[fsr]->getY() + wgt *
        (y + sin_phi * sin_theta * curr_segment->_length / 2.0)
        * curr_segment->_length / _FSR_volumes[fsr]);

    _centroids[fsr]->
        setZ(_centroids[fsr]->getZ() + wgt *
        (z + cos_theta * curr_segment->_length / 2.0)
        * curr_segment->_length / _FSR_volumes[fsr]);

    /* Unset the lock for this FSR */
    omp_unset_lock(&_FSR_locks[fsr]);

    x += cos_phi * sin_theta * curr_segment->_length;
    y += sin_phi * sin_theta * curr_segment->_length;
    z += cos_theta * curr_segment->_length;

    _starting_points[tid][track_idx].setX(x);
    _starting_points[tid][track_idx].setY(y);
    _starting_points[tid][track_idx].setZ(z);
  }
}


/**
 * @brief Constructor for TransportSweep calls the TraverseSegments
 *        constructor and initializes the associated CPUSolver to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
TransportSweep::TransportSweep(TrackGenerator* track_generator)
                              : TraverseSegments(track_generator) {
  _cpu_solver = NULL;
  _tracks_per_stack = NULL;

  /* Allocate temporary storage of FSR fluxes */
  int num_threads = omp_get_max_threads();
  int num_groups = track_generator->getGeometry()->getNumEnergyGroups();
  _thread_fsr_fluxes = new FP_PRECISION*[num_threads];
  for (int i=0; i < num_threads; i++)
    _thread_fsr_fluxes[i] = new FP_PRECISION[2*num_groups];

  /* Get the number of tracks per stack if 3D calculation */
  TrackGenerator3D* TG_3D = dynamic_cast<TrackGenerator3D*>(track_generator);
  if (TG_3D != NULL)
    _tracks_per_stack = TG_3D->getTracksPerStack();
}


//FIXME
TransportSweep::~TransportSweep() {
  int num_threads = omp_get_max_threads();
  for (int i=0; i < num_threads; i++)
    delete [] _thread_fsr_fluxes[i];
  delete [] _thread_fsr_fluxes;
}



/**
 * @brief MOC equations are applied to every segment in the TrackGenerator
 * @details SegmntationKernels are allocated to temporarily save segments. Then
 *          onTrack(...) applies the MOC equations to each segment and
 *          transfers boundary fluxes for the corresponding Track.
 */
void TransportSweep::execute() {
#pragma omp parallel
  {
    MOCKernel* kernel = getKernel<SegmentationKernel>();
    loopOverTracks(kernel);
  }
}


/**
 * @brief Sets the CPUSolver so that TransportSweep can apply MOC equations
 * @details This allows TransportSweep to transfer boundary fluxes from the
 *          CPUSolver and tally scalar fluxes
 * @param cpu_solver The CPUSolver which applies the MOC equations
 */
void TransportSweep::setCPUSolver(CPUSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


/**
 * @brief Applies the MOC equations the Track and segments
 * @details The MOC equations are applied to each segment, attenuating the
 *          Track's angular flux and tallying FSR contributions. Finally,
 *          Track boundary fluxes are transferred.
 * @param track The Track for which the angular flux is attenuated and
 *        transferred
 * @param segments The segments over which the MOC equations are applied
 */
void TransportSweep::onTrack(Track* track, segment* segments) {

  /* Get the temporary FSR flux */
  int tid = omp_get_thread_num();
  FP_PRECISION* thread_fsr_flux = _thread_fsr_fluxes[tid];

  /* Extract Track information */
  int track_id = track->getUid();
  int azim_index = track->getAzimIndex();
  int num_segments = track->getNumSegments();
  FP_PRECISION* track_flux;

  /* Extract the polar index if a 3D track */
  int polar_index = 0;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL)
    polar_index = track_3D->getPolarIndex();

  /* Extract the maximum track index and get Track data */
  TrackStackIndexes tsi;
  Track** tracks_array;
  int max_track_index = 0;
  if (_segment_formation == OTF_STACKS) {

    /* Extract indexes */
    int xy_index = track->getXYIndex();
    max_track_index = _tracks_per_stack[azim_index][xy_index][polar_index] - 1;
    tsi._azim = azim_index;
    tsi._xy = xy_index;
    tsi._polar = polar_index;

    /* Get Track data for the entire z-stack */
    Track3D* tracks_3D = _track_generator_3D->getTemporary3DTracks(tid);
    for (int z=0; z < max_track_index+1; z++) {
      tsi._z = z;
      _track_generator_3D->getTrackOTF(&tracks_3D[z], &tsi);
    }
    tracks_array = _track_generator_3D->getTemporaryTracksArray(tid);
  }
  else {
    tracks_array = &track;
  }

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {

    /* Get the forward track flux */
    segment* curr_segment = &segments[s];
    int curr_track_id = track_id + curr_segment->_track_idx;
    track_flux = _cpu_solver->getBoundaryFlux(curr_track_id, true);

    /* Apply MOC equations */
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, polar_index,
                                 track_flux, thread_fsr_flux);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, polar_index,
                                     track_flux, true);
  }

  /* Transfer boundary angular flux to outgoing Track */
  for (int i=0; i <= max_track_index; i++) {
    track_flux = _cpu_solver->getBoundaryFlux(track_id+i, true);
    _cpu_solver->transferBoundaryFlux(tracks_array[i], azim_index, polar_index,
                                      true, track_flux);
  }

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {

    /* Get the backward track flux */
    segment* curr_segment = &segments[s];
    int curr_track_id = track_id + curr_segment->_track_idx;
    track_flux = _cpu_solver->getBoundaryFlux(curr_track_id, false);

    /* Apply MOC equations */
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, polar_index,
                                 track_flux, thread_fsr_flux);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, polar_index,
                                     track_flux, false);
  }

  /* Transfer boundary angular flux to outgoing Track */
  for (int i=0; i <= max_track_index; i++) {
    track_flux = _cpu_solver->getBoundaryFlux(track_id+i, false);
    _cpu_solver->transferBoundaryFlux(tracks_array[i], azim_index, polar_index,
                                      false, track_flux);
  }
}


/**
 * @brief Constructor for DumpSegments calls the TraverseSegments
 *        constructor and initializes the output FILE to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
DumpSegments::DumpSegments(TrackGenerator* track_generator)
                           : TraverseSegments(track_generator) {
  _out = NULL;
}


/**
 * @brief Wrties all tracking information to file
 * @details SegmentationKernels are created to temporarily store segments for
 *          on-the-fly method. For each Track, onTrack(...) writes the tracking
 *          information to file.
 */
void DumpSegments::execute() {
  MOCKernel* kernel = getKernel<SegmentationKernel>();
  loopOverTracks(kernel);
}


/**
 * @brief Sets the file which to write tracking information
 * @param out the file which to write tracking infmormation
 */
void DumpSegments::setOutputFile(FILE* out) {
  _out = out;
}


/**
 * @brief Writes tracking information to file for a Track and associated
 *        segments
 * @param track The Track whose information is written to file
 * @param segments The segments associated with the Track whose information is
 *        written to file
 */
void DumpSegments::onTrack(Track* track, segment* segments) {

  /* Write data for this Track to the Track file */
  int num_segments = track->getNumSegments();
  fwrite(&num_segments, sizeof(int), 1, _out);

  /* Get CMFD mesh object */
  Cmfd* cmfd = _track_generator->getGeometry()->getCmfd();

  /* Loop over all segments for this Track */
  for (int s=0; s < num_segments; s++) {

    /* Get data for this segment */
    segment* curr_segment = &segments[s];
    FP_PRECISION length = curr_segment->_length;
    int material_id = curr_segment->_material->getId();
    int region_id = curr_segment->_region_id;
    int track_idx = curr_segment->_track_idx;

    /* Write data for this segment to the Track file */
    fwrite(&length, sizeof(double), 1, _out);
    fwrite(&material_id, sizeof(int), 1, _out);
    fwrite(&region_id, sizeof(int), 1, _out);
    fwrite(&track_idx, sizeof(int), 1, _out);

    /* Write CMFD-related data for the Track if needed */
    if (cmfd != NULL) {
      int cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
      int cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;
      fwrite(&cmfd_surface_fwd, sizeof(int), 1, _out);
      fwrite(&cmfd_surface_bwd, sizeof(int), 1, _out);
    }
  }
}


/**
 * @brief Constructor for ReadSegments calls the TraverseSegments
 *        constructor and initializes the input FILE to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
ReadSegments::ReadSegments(TrackGenerator* track_generator)
                           : TraverseSegments(track_generator) {
  _in = NULL;
}


/**
 * @brief Reads a tracking file and saves the information explicitly for every
 *        Track in the TrackGenerator.
 * @details The tracking file is set by setInputFile(...)
 */
void ReadSegments::execute() {
  loopOverTracks(NULL);
}


/**
 * @brief Sets the input file to read in tracking information
 * @param in The input tracking file
 */
void ReadSegments::setInputFile(FILE* input) {
  _in = input;
}


/**
 * @brief Saves tracking information to the corresponding Track explicity
 * @param track The track for which all tracking information is explicitly
 *        saved (including segments)
 * @param segments The segments associated with the Track
 */
void ReadSegments::onTrack(Track* track, segment* segments) {

  /* Get CMFD mesh object */
  Geometry* geometry = _track_generator->getGeometry();
  Cmfd* cmfd = geometry->getCmfd();
  int ret;

  /* Get materials map */
  std::map<int, Material*> materials = geometry->getAllMaterials();

  /* Import data for this Track from Track file */
  int num_segments;
  ret = fread(&num_segments, sizeof(int), 1, _in);

  /* Loop over all segments in this Track */
  for (int s=0; s < num_segments; s++) {

    /* Import data for this segment from Track file */
    double length;
    ret = fread(&length, sizeof(double), 1, _in);
    int material_id;
    ret = fread(&material_id, sizeof(int), 1, _in);
    int region_id;
    ret = fread(&region_id, sizeof(int), 1, _in);
    int  track_idx;
    ret = fread(&track_idx, sizeof(int), 1, _in);

    /* Initialize segment with the data */
    segment* curr_segment = new segment;
    curr_segment->_length = length;
    curr_segment->_material = materials[material_id];
    curr_segment->_region_id = region_id;
    curr_segment->_track_idx = track_idx;

    /* Import CMFD-related data if needed */
    if (cmfd != NULL) {
      int cmfd_surface_fwd;
      ret = fread(&cmfd_surface_fwd, sizeof(int), 1, _in);
      curr_segment->_cmfd_surface_fwd = cmfd_surface_fwd;
      int cmfd_surface_bwd;
      ret = fread(&cmfd_surface_bwd, sizeof(int), 1, _in);
      curr_segment->_cmfd_surface_bwd = cmfd_surface_bwd;
    }

    /* Add this segment to the Track */
    track->addSegment(curr_segment);
  }
}


//FIXME
TransportSweepOTF::TransportSweepOTF(TrackGenerator* track_generator)
                                   : TraverseSegments(track_generator) {
  _cpu_solver = NULL;
}


//FIXME
void TransportSweepOTF::execute() {
#pragma omp parallel
  {
    TransportKernel kernel(_track_generator, 0);
    kernel.setCPUSolver(_cpu_solver);
    loopOverTracksByStackTwoWay(&kernel);
  }
}


//FIXME
void TransportSweepOTF::setCPUSolver(CPUSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}

//FIXME
void TransportSweepOTF::onTrack(Track* track, segment* segments) {
}
