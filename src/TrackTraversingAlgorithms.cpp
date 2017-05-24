#include "TrackTraversingAlgorithms.h"
#include "CPUSolver.h"
#include "CPULSSolver.h"
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
  FP_PRECISION infinity = std::numeric_limits<NEW_PRECISION>::max();
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
    NEW_PRECISION length = segments[s]._length;
    Material* material = segments[s]._material;
    NEW_PRECISION* sigma_t = material->getSigmaT();

    for (int e=0; e < material->getNumEnergyGroups(); e++) {
      NEW_PRECISION tau = length*sigma_t[e];
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

  /* Get the direction of travel */
  double phi = track->getPhi();
  double theta = M_PI_2;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != 0)
    theta = track_3D->getTheta();
  double xdir = cos(phi) * sin(theta);
  double ydir = sin(phi) * sin(theta);
  double zdir = cos(theta);

  /* Extract data from this segment to compute its optical
   * length */
  for (int s = 0; s < track->getNumSegments(); s++) {
    segment* curr_segment = track->getSegment(s);
    Material* material = curr_segment->_material;
    double length = curr_segment->_length;
    long fsr_id = curr_segment->_region_id;

    /* Compute number of segments to split this segment into */
    int min_num_cuts = 1;
    int num_groups = material->getNumEnergyGroups();
    NEW_PRECISION* sigma_t = material->getSigmaT();

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

    /* Extract the current starting points */
    double x_curr = curr_segment->_starting_position[0];
    double y_curr = curr_segment->_starting_position[1];
    double z_curr = curr_segment->_starting_position[2];

    /* Split the segment into sub-segments */
    for (int k=0; k < min_num_cuts; k++) {

      /* Create a new Track segment */
      segment* new_segment = new segment;
      new_segment->_material = material;
      new_segment->_length = length / min_num_cuts;
      new_segment->_region_id = fsr_id;

      /* Assign CMFD surface boundaries */
      if (k == 0)
        new_segment->_cmfd_surface_bwd = cmfd_surface_bwd;

      if (k == min_num_cuts-1)
        new_segment->_cmfd_surface_fwd = cmfd_surface_fwd;

      /* Set the starting position */
      new_segment->_starting_position[0] = x_curr;
      new_segment->_starting_position[1] = y_curr;
      new_segment->_starting_position[2] = z_curr;
      x_curr += new_segment->_length * xdir;
      y_curr += new_segment->_length * ydir;
      z_curr += new_segment->_length * zdir;

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
      int xy_index = track->getXYIndex();
      int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
      Track3D* current_stack = _track_generator_3D->getTemporary3DTracks(tid);
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
    long fsr = curr_segment->_region_id;
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


//FIXME
LinearExpansionGenerator::LinearExpansionGenerator(CPULSSolver* solver)
    : TraverseSegments(solver->getTrackGenerator()) {

  _lin_exp_coeffs = solver->getLinearExpansionCoeffsBuffer();
  _src_constants = solver->getSourceConstantsBuffer();
  TrackGenerator* track_generator = solver->getTrackGenerator();
  _FSR_volumes = track_generator->getFSRVolumesBuffer();
  _FSR_locks = track_generator->getFSRLocks();
  _quadrature = track_generator->getQuadrature();
  _num_groups = track_generator->getGeometry()->getNumEnergyGroups();

  _num_flat = 0;
  int num_rows = 1;
  _num_coeffs = 3;
  if (_track_generator_3D != NULL) {
    num_rows = _track_generator_3D->getMaxNumTracksPerStack();
    _num_coeffs = 6;
  }

  int num_threads = omp_get_max_threads();
  _starting_points = new Point*[num_threads];
  _thread_source_constants = new FP_PRECISION*[num_threads];

  for (int i=0; i < num_threads; i++) {
    _thread_source_constants[i] = new FP_PRECISION[_num_coeffs * _num_groups];
    _starting_points[i] = new Point[num_rows];
  }

  _exp_evaluator = new ExpEvaluator();

  std::string msg = "Initializing track linear source componenets";
  _progress = new Progress(_track_generator->getNumTracks(), msg, 0.1,
                    track_generator->getGeometry(), true);
}

//FIXME destructor
LinearExpansionGenerator::~LinearExpansionGenerator() {
  int num_threads = omp_get_max_threads();
  for (int i=0; i < num_threads; i++) {
    delete [] _thread_source_constants[i];
    delete [] _starting_points[i];
  }
  delete [] _thread_source_constants;
  delete [] _starting_points;
  delete _exp_evaluator;
  delete _progress;
}


//FIXME
void LinearExpansionGenerator::execute() {
#pragma omp parallel
  {
    MOCKernel* kernel = getKernel<SegmentationKernel>();
    loopOverTracks(kernel);
  }

  Geometry* geometry = _track_generator->getGeometry();
  long num_FSRs = geometry->getNumFSRs();
  //FIXME
  FP_PRECISION* inv_lin_exp_coeffs = new FP_PRECISION[num_FSRs*_num_coeffs];
  memset(inv_lin_exp_coeffs, 0., num_FSRs*_num_coeffs*sizeof(FP_PRECISION));

  FP_PRECISION* lem = _lin_exp_coeffs;
  FP_PRECISION* ilem = inv_lin_exp_coeffs;
  int nc = _num_coeffs;

  /* Invert the expansion coefficient matrix */
  if (_track_generator_3D != NULL) {
#pragma omp parallel for
    for (int r=0; r < num_FSRs; r++) {
      FP_PRECISION det;
      det = lem[r*nc + 0] * lem[r*nc + 1] * lem[r*nc + 5] +
            lem[r*nc + 2] * lem[r*nc + 4] * lem[r*nc + 3] +
            lem[r*nc + 3] * lem[r*nc + 2] * lem[r*nc + 4] -
            lem[r*nc + 0] * lem[r*nc + 4] * lem[r*nc + 4] -
            lem[r*nc + 3] * lem[r*nc + 1] * lem[r*nc + 3] -
            lem[r*nc + 2] * lem[r*nc + 2] * lem[r*nc + 5];

      double volume = _FSR_volumes[r];
      if (std::abs(det) < MIN_DET || volume < 1e-6) {
        log_printf(INFO, "Unable to form linear source components in "
                   "source region %d. Switching to flat source in that "
                   "source region.", r);
        _num_flat++;
        ilem[r*nc + 0] = 0.0;
        ilem[r*nc + 1] = 0.0;
        ilem[r*nc + 2] = 0.0;
        ilem[r*nc + 3] = 0.0;
        ilem[r*nc + 4] = 0.0;
        ilem[r*nc + 5] = 0.0;
      }
      else {

        ilem[r*nc + 0] = (lem[r*nc + 1] * lem[r*nc + 5] -
                          lem[r*nc + 4] * lem[r*nc + 4]) / det;
        ilem[r*nc + 1] = (lem[r*nc + 0] * lem[r*nc + 5] -
                          lem[r*nc + 3] * lem[r*nc + 3]) / det;
        ilem[r*nc + 2] = (lem[r*nc + 3] * lem[r*nc + 4] -
                          lem[r*nc + 2] * lem[r*nc + 5]) / det;
        ilem[r*nc + 3] = (lem[r*nc + 2] * lem[r*nc + 4] -
                          lem[r*nc + 3] * lem[r*nc + 1]) / det;
        ilem[r*nc + 4] = (lem[r*nc + 3] * lem[r*nc + 2] -
                          lem[r*nc + 0] * lem[r*nc + 4]) / det;
        ilem[r*nc + 5] = (lem[r*nc + 0] * lem[r*nc + 1] -
                          lem[r*nc + 2] * lem[r*nc + 2]) / det;

      }
    }
  }
  else {
#pragma omp parallel for
    for (long r=0; r < num_FSRs; r++) {

      FP_PRECISION det;
      det = lem[r*nc  ] * lem[r*nc + 1] - lem[r*nc + 2] * lem[r*nc + 2];

      if (std::abs(det) < MIN_DET) {
        ilem[r*nc + 0] = 0.0;
        ilem[r*nc + 1] = 0.0;
        ilem[r*nc + 2] = 0.0;
        log_printf(INFO, "Unable to form linear source components in "
                   "source region %d. Switching to flat source in that "
                   "source region.", r);
        _num_flat++;
      }
      else {
        ilem[r*nc + 0] =  lem[r*nc + 1] / det;
        ilem[r*nc + 1] =  lem[r*nc + 0] / det;
        ilem[r*nc + 2] = -lem[r*nc + 2] / det;
      }
    }
  }

  //FIXME
  if (false) {
#pragma omp parallel for
    for (long r = 0; r < num_FSRs; r++) {
      if (_track_generator->getFSRVolume(r) < 1e-3)
        for (int i=0; i < nc; i++)
          ilem[r*nc+i] = 0.0;
    }
  }

  //FIXME
  if (true) {
      //FIXME -2 -> +4
    double max_linear_radius = (8.5*15-2)*1.26; // (8.5*15+4) * 1.26 || * -2 * || 10*17*1.26
    Universe* root_universe = geometry->getRootUniverse();
    double center_x = (root_universe->getMinX() + root_universe->getMaxX()) / 2;
    double center_y = (root_universe->getMinY() + root_universe->getMaxY()) / 2;
    double center_z = (root_universe->getMinZ() + root_universe->getMaxZ()) / 2;
#pragma omp parallel for
    for (long r = 0; r < num_FSRs; r++) {
      Point* centroid = geometry->getFSRCentroid(r);
      double dist = centroid->distance(center_x, center_y, centroid->getZ());
      if (dist > max_linear_radius) {
        for (int i=0; i < nc; i++)
          ilem[r*nc+i] = 0.0;
      }
    }
  }

  memcpy(_lin_exp_coeffs, inv_lin_exp_coeffs,
         num_FSRs*_num_coeffs*sizeof(FP_PRECISION));
  delete [] inv_lin_exp_coeffs;

  /* Notify user of any regions needing to use a flat source approximation */
  int total_num_flat = _num_flat;
  long total_num_FSRs = num_FSRs;
#ifdef MPIx
  if (geometry->isDomainDecomposed()) {
    MPI_Allreduce(&_num_flat, &total_num_flat, 1, MPI_INT, MPI_SUM,
                  geometry->getMPICart());
    MPI_Allreduce(&num_FSRs, &total_num_FSRs, 1, MPI_LONG, MPI_SUM,
                  geometry->getMPICart());
    }
#endif
  if (total_num_flat > 0)
    if (geometry->isRootDomain())
      log_printf(WARNING, "Unable to form linear source components in %d / %d "
                          "source regions. Switching to flat source in those "
                          "source regions.", total_num_flat, total_num_FSRs);
}


//FIXME
void LinearExpansionGenerator::onTrack(Track* track, segment* segments) {

  /* Extract track information */
  Point* start = track->getStart();
  double phi = track->getPhi();
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);
  int azim_index = track->getAzimIndex();
  int xy_index = track->getXYIndex();

  /* Use local array accumulator to prevent false sharing */
  int tid = omp_get_thread_num();
  FP_PRECISION* thread_src_constants = _thread_source_constants[tid];

  /* Calculate the azimuthal weight */
  double wgt = _quadrature->getAzimSpacing(azim_index)
      * _quadrature->getAzimWeight(azim_index);

  /* Get polar angles and weight depending on the dimensionality */
  double sin_theta = 1;
  double cos_theta = 0;
  int polar_index = 0;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL) {
    double theta = track_3D->getTheta();
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    int polar_index = track_3D->getPolarIndex();
    wgt *= _quadrature->getPolarSpacing(azim_index, polar_index)
        *_quadrature->getPolarWeight(azim_index, polar_index);
  }

  /* Loop over segments to accumlate contribution to centroids */
  Geometry* geometry = _track_generator->getGeometry();
  for (int s=0; s < track->getNumSegments(); s++) {

    /* Extract segment information */
    segment* curr_segment = &segments[s];
    long fsr = curr_segment->_region_id;
    int track_idx = curr_segment->_track_idx;
    NEW_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
    FP_PRECISION length = curr_segment->_length;
    FP_PRECISION length_2 = length * length;

    /* Extract FSR information */
    double volume = _FSR_volumes[fsr];

    /* Extract the starting points of the segment */
    double x = curr_segment->_starting_position[0];
    double y = curr_segment->_starting_position[1];
    double z = curr_segment->_starting_position[2];

    /* Get the centroid of the segment in the local coordinate system */
    double xc = x + length * 0.5 * cos_phi * sin_theta;
    double yc = y + length * 0.5 * sin_phi * sin_theta;
    double zc = z + length * 0.5 * cos_theta;

    /* Set the FSR src constants buffer to zero */
    memset(thread_src_constants, 0.0, _num_groups * _num_coeffs *
           sizeof(FP_PRECISION));

    FP_PRECISION vol_impact = wgt * length / volume;
    for (int g=0; g < _num_groups; g++) {

      thread_src_constants[g*_num_coeffs] += vol_impact * xc * xc;
      thread_src_constants[g*_num_coeffs + 1] += vol_impact * yc * yc;
      thread_src_constants[g*_num_coeffs + 2] += vol_impact * xc * yc;

      if (track_3D != NULL) {
        thread_src_constants[g*_num_coeffs + 3] += vol_impact * xc * zc;
        thread_src_constants[g*_num_coeffs + 4] += vol_impact * yc * zc;
        thread_src_constants[g*_num_coeffs + 5] += vol_impact * zc * zc;
      }

      /* Calculate the optical path length and source contribution */
      FP_PRECISION tau = length * sigma_t[g];
      FP_PRECISION src_constant = vol_impact * length / 2.0;

      if (track_3D == NULL) {
        for (int p=0; p < _quadrature->getNumPolarAngles()/2; p++) {

          FP_PRECISION sin_theta = _quadrature->getSinTheta(azim_index, p);
          FP_PRECISION G2_src =
              length * _exp_evaluator->computeExponentialG2(tau / sin_theta)
              * src_constant * 2 * _quadrature->getPolarWeight(azim_index, p)
              * sin_theta;

          thread_src_constants[g*_num_coeffs] += cos_phi * cos_phi * G2_src;
          thread_src_constants[g*_num_coeffs + 1] += sin_phi * sin_phi
              * G2_src;
          thread_src_constants[g*_num_coeffs + 2] += sin_phi * cos_phi
              * G2_src;
        }
      }
      else {

        FP_PRECISION G2_src = _exp_evaluator->computeExponentialG2(tau) *
            length * src_constant;

        thread_src_constants[g*_num_coeffs] += cos_phi * cos_phi * G2_src
            * sin_theta * sin_theta;
        thread_src_constants[g*_num_coeffs + 1] += sin_phi * sin_phi * G2_src
            * sin_theta * sin_theta;
        thread_src_constants[g*_num_coeffs + 2] += sin_phi * cos_phi * G2_src
            * sin_theta * sin_theta;
        thread_src_constants[g*_num_coeffs + 3] += cos_phi * cos_theta * G2_src
            * sin_theta;
        thread_src_constants[g*_num_coeffs + 4] += sin_phi * cos_theta * G2_src
            * sin_theta;
        thread_src_constants[g*_num_coeffs + 5] += cos_theta * cos_theta * G2_src;
      }
    }

    /* Set the lock for this FSR */
    omp_set_lock(&_FSR_locks[fsr]);

    _lin_exp_coeffs[fsr*_num_coeffs] += wgt * length / volume *
        (xc * xc + pow(cos_phi * sin_theta * length, 2) / 12.0);
    _lin_exp_coeffs[fsr*_num_coeffs + 1] += wgt * length / volume *
        (yc * yc + pow(sin_phi * sin_theta * length, 2) / 12.0);
    _lin_exp_coeffs[fsr*_num_coeffs + 2] += wgt * length / volume *
        (xc * yc + sin_phi * cos_phi * pow(sin_theta * length, 2) / 12.0);

    if (track_3D != NULL) {
      _lin_exp_coeffs[fsr*_num_coeffs + 3] += wgt * length / volume *
          (xc * zc + cos_phi * cos_theta * sin_theta * pow(length, 2) / 12.0);
      _lin_exp_coeffs[fsr*_num_coeffs + 4] += wgt * length / volume *
          (yc * zc + sin_phi * cos_theta * sin_theta * pow(length, 2) / 12.0);
      _lin_exp_coeffs[fsr*_num_coeffs + 5] += wgt * length / volume *
          (zc * zc + pow(cos_theta * length, 2) / 12.0);
    }

	/* Set the source constants for all groups and coefficients */
    for (int g=0; g < _num_groups; g++) {
      for (int i=0; i < _num_coeffs; i++)
        _src_constants[fsr*_num_groups*_num_coeffs + g*_num_coeffs + i] +=
          thread_src_constants[g*_num_coeffs + i];
    }

    /* Unset the lock for this FSR */
    omp_unset_lock(&_FSR_locks[fsr]);
  }

  /* Determine progress */
  int max_track_index = 0;
  if (_segment_formation == OTF_STACKS) {
    int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
    max_track_index = tracks_per_stack[azim_index][xy_index][polar_index] - 1;
  }
  for (int i=0; i <= max_track_index; i++)
    _progress->incrementCounter();
}


/**
 * @brief Constructor for TransportSweep calls the TraverseSegments
 *        constructor and initializes the associated CPUSolver to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
TransportSweep::TransportSweep(CPUSolver* cpu_solver)
    : TraverseSegments(cpu_solver->getTrackGenerator()) {

  _cpu_solver = cpu_solver;
  _ls_solver = dynamic_cast<CPULSSolver*>(cpu_solver);
  TrackGenerator* track_generator = cpu_solver->getTrackGenerator();
  _geometry = _track_generator->getGeometry();

  /* Determine size of temporary storage for FSR fluxes */
  int num_threads = omp_get_max_threads();
  int num_groups = _geometry->getNumEnergyGroups();
  int size;
  if (_ls_solver != NULL)
    size = 4 * num_groups;
  else
    size = num_groups;

  /* Pad the buffer to prevent false sharing */
  size += 8;

  /* Allocate temporary storage of FSR fluxes */
  _thread_fsr_fluxes = new FP_PRECISION*[num_threads];
  for (int i=0; i < num_threads; i++)
    _thread_fsr_fluxes[i] = new FP_PRECISION[size];
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
  long track_id = track->getUid();
  int azim_index = track->getAzimIndex();
  int xy_index = track->getXYIndex();
  int num_segments = track->getNumSegments();
    //FIXME MEM : float / FP_PRECISION
  float* track_flux;

  /* Extract the polar index if a 3D track */
  int polar_index = 0;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL)
    polar_index = track_3D->getPolarIndex();

  /* Compute unit vector if necessary */
  NEW_PRECISION direction[3];
  if (_ls_solver != NULL) {
    double phi = track->getPhi();
    double cos_theta = 0.0;
    double sin_theta = 1.0;
    if (track_3D != NULL) {
      double theta = track_3D->getTheta();
      cos_theta = cos(theta);
      sin_theta = sin(theta);
    }
    direction[0] = cos(phi) * sin_theta;
    direction[1] = sin(phi) * sin_theta;
    direction[2] = cos_theta;
  }

  /* Extract the maximum track index */
  Track** tracks_array = &track;
  int max_track_index = 0;
  if (_segment_formation == OTF_STACKS) {
    int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
    max_track_index = tracks_per_stack[azim_index][xy_index][polar_index] - 1;
    tracks_array = _track_generator_3D->getTemporaryTracksArray(tid);
  }

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {

    /* Get the forward track flux */
    segment* curr_segment = &segments[s];
    long curr_track_id = track_id + curr_segment->_track_idx;
    track_flux = _cpu_solver->getBoundaryFlux(curr_track_id, true);
   
    /* Apply MOC equations */
    if (_ls_solver == NULL)
      _cpu_solver->tallyScalarFlux(curr_segment, azim_index, polar_index,
                                   track_flux, thread_fsr_flux);
    else
      _ls_solver->tallyLSScalarFlux(curr_segment, azim_index, polar_index,
                                    track_flux, thread_fsr_flux, direction);

    /* Tally the current for CMFD */
    _cpu_solver->tallyCurrent(curr_segment, azim_index, polar_index,
                              track_flux, true);
  }

  /* Transfer boundary angular flux to outgoing Track */
  for (int i=0; i <= max_track_index; i++) {
    track_flux = _cpu_solver->getBoundaryFlux(track_id+i, true);
    _cpu_solver->transferBoundaryFlux(tracks_array[i], azim_index, polar_index,
                                      true, track_flux);
  }

  /* Reverse the direction */
  for (int i=0; i<3; i++)
    direction[i] *= -1;

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {

    /* Get the backward track flux */
    segment* curr_segment = &segments[s];
    long curr_track_id = track_id + curr_segment->_track_idx;
    track_flux = _cpu_solver->getBoundaryFlux(curr_track_id, false);
    
    /* Apply MOC equations */
    if (_ls_solver == NULL)
      _cpu_solver->tallyScalarFlux(curr_segment, azim_index, polar_index,
                                   track_flux, thread_fsr_flux);
    else
      _ls_solver->tallyLSScalarFlux(curr_segment, azim_index, polar_index,
                                    track_flux, thread_fsr_flux, direction);


    /* Tally the current for CMFD */
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
    int material_id;
    if (curr_segment->_material != NULL)
      material_id = curr_segment->_material->getId();
    else
      material_id = -1;
    long region_id = curr_segment->_region_id;
    int track_idx = curr_segment->_track_idx;
    double start_x = curr_segment->_starting_position[0];
    double start_y = curr_segment->_starting_position[1];
    double start_z = curr_segment->_starting_position[2];

    /* Write data for this segment to the Track file */
    fwrite(&length, sizeof(double), 1, _out);
    fwrite(&material_id, sizeof(int), 1, _out);
    fwrite(&region_id, sizeof(long), 1, _out);
    fwrite(&track_idx, sizeof(int), 1, _out);
    fwrite(&start_x, sizeof(double), 1, _out);
    fwrite(&start_y, sizeof(double), 1, _out);
    fwrite(&start_z, sizeof(double), 1, _out);

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
  ret = geometry->twiddleRead(&num_segments, sizeof(int), 1, _in);

  /* Loop over all segments in this Track */
  for (int s=0; s < num_segments; s++) {

    /* Import data for this segment from Track file */
    double length;
    ret = geometry->twiddleRead(&length, sizeof(double), 1, _in);
    int material_id;
    ret = geometry->twiddleRead(&material_id, sizeof(int), 1, _in);
    long region_id;
    ret = geometry->twiddleRead(&region_id, sizeof(long), 1, _in);
    int  track_idx;
    ret = geometry->twiddleRead(&track_idx, sizeof(int), 1, _in);
    double start_x;
    ret = geometry->twiddleRead(&start_x, sizeof(double), 1, _in);
    double start_y;
    ret = geometry->twiddleRead(&start_y, sizeof(double), 1, _in);
    double start_z;
    ret = geometry->twiddleRead(&start_z, sizeof(double), 1, _in);

    /* Initialize segment with the data */
    segment* curr_segment = new segment;
    curr_segment->_length = length;
    curr_segment->_material = materials[material_id];
    curr_segment->_region_id = region_id;
    curr_segment->_track_idx = track_idx;
    curr_segment->_starting_position[0] = start_x;
    curr_segment->_starting_position[1] = start_y;
    curr_segment->_starting_position[2] = start_z;

    /* Import CMFD-related data if needed */
    if (cmfd != NULL) {
      int cmfd_surface_fwd;
      ret = geometry->twiddleRead(&cmfd_surface_fwd, sizeof(int), 1, _in);
      curr_segment->_cmfd_surface_fwd = cmfd_surface_fwd;
      int cmfd_surface_bwd;
      ret = geometry->twiddleRead(&cmfd_surface_bwd, sizeof(int), 1, _in);
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


//FIXME
RecenterSegments::RecenterSegments(TrackGenerator* track_generator)
                                   : TraverseSegments(track_generator) {
  _geometry = _track_generator->getGeometry();
}


//FIXME
void RecenterSegments::execute() {
#pragma omp parallel
  {
    MOCKernel* kernel = getKernel<SegmentationKernel>();
    loopOverTracks(kernel);
  }
}


//FIXME
void RecenterSegments::onTrack(Track* track, segment* segments) {
  if (_geometry->containsFSRCentroids()) {
    for (int s=0; s < track->getNumSegments(); s++) {
      long fsr_id = segments[s]._region_id;
      Point* centroid = _geometry->getFSRCentroid(fsr_id);
      segments[s]._starting_position[0] -= centroid->getX();
      segments[s]._starting_position[1] -= centroid->getY();
      segments[s]._starting_position[2] -= centroid->getZ();
    }
  }
}


/**
 * @brief Constructor for DumpSegments calls the TraverseSegments
 *        constructor and initializes the output FILE to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
//FIXME
 */
PrintSegments::PrintSegments(TrackGenerator* track_generator)
                           : TraverseSegments(track_generator) {
  _out = NULL;
}


/**
 * @brief Wrties all tracking information to file
 * @details SegmentationKernels are created to temporarily store segments for
 *          on-the-fly method. For each Track, onTrack(...) writes the tracking
 *          information to file.
 //FIXME
 */
void PrintSegments::execute() {
  MOCKernel* kernel = getKernel<SegmentationKernel>();
  loopOverTracks(kernel);
}


/**
 * @brief Sets the file which to write tracking information
 * @param out the file which to write tracking infmormation
 //FIXME
 */
void PrintSegments::setOutputFile(FILE* out) {
  _out = out;
}


/**
 * @brief Writes tracking information to file for a Track and associated
 *        segments
 * @param track The Track whose information is written to file
 * @param segments The segments associated with the Track whose information is
 *        written to file
 //FIXME
 */
void PrintSegments::onTrack(Track* track, segment* segments) {

  /* Write data for this Track to the Track file */
  int num_segments = track->getNumSegments();

  /* Get CMFD mesh object */
  Cmfd* cmfd = _track_generator->getGeometry()->getCmfd();

  /* Loop over all segments for this Track */
  for (int s=0; s < num_segments; s++) {

    /* Get data for this segment */
    segment* curr_segment = &segments[s];
    FP_PRECISION length = curr_segment->_length;
    int material_id = curr_segment->_material->getId();
    long region_id = curr_segment->_region_id;
    int track_idx = curr_segment->_track_idx;
    double start_x = curr_segment->_starting_position[0];
    double start_y = curr_segment->_starting_position[1];
    double start_z = curr_segment->_starting_position[2];

    /* Write data for this segment to the Track file */
    fprintf(_out, "%6.4f %d %ld %6.4f %6.4f %6.4f", length, material_id,
            region_id, start_x, start_y, start_z);

    /* Write CMFD-related data for the Track if needed */
    if (cmfd != NULL) {
      int cmfd_surface_fwd = curr_segment->_cmfd_surface_fwd;
      int cmfd_surface_bwd = curr_segment->_cmfd_surface_bwd;
      fprintf(_out, " %d %d", cmfd_surface_fwd, cmfd_surface_bwd);
    }
    fprintf(_out, "\n");
  }
}


