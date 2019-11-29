#include "TrackTraversingAlgorithms.h"
#include "CPUSolver.h"
#include "CPULSSolver.h"
#include "Quadrature.h"

/**
 * @brief Constructor for MaxOpticalLength calls the TraverseSegments
 *        constructor and sets the max optical path length to zero.
 * @param track_generator The TrackGenerator to pull tracking information from
 */
MaxOpticalLength::MaxOpticalLength(TrackGenerator* track_generator)
                                 : TraverseSegments(track_generator) {
  _min_tau = 1e6;
  _max_tau = 0;
}


/**
 * @brief Determines the maximum optical path length for the TrackGenerator
 *        provided during construction.
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
    // OTF ray tracing requires segmentation of tracks
    if (_segment_formation != EXPLICIT_2D &&
        _segment_formation != EXPLICIT_3D) {
      MOCKernel* kernel = getKernel<SegmentationKernel>();
      loopOverTracks(kernel);
    }
    else
      loopOverTracks(NULL);
  }
  _track_generator->setMaxOpticalLength(_max_tau);

  /* Notify user of the range of optical lengths present in the geometry */
  Geometry* geometry = _track_generator->getGeometry();
  float global_min_tau = _min_tau;
  float global_max_tau = _max_tau;
#ifdef MPIx
  if (geometry->isDomainDecomposed()) {
    MPI_Allreduce(&_min_tau, &global_min_tau, 1, MPI_FLOAT, MPI_MIN,
                  geometry->getMPICart());
    MPI_Allreduce(&_max_tau, &global_max_tau, 1, MPI_FLOAT, MPI_MAX,
                  geometry->getMPICart());
  }
#endif
  log_printf(INFO_ONCE, "Min/max optical lengths in geometry %.2e / %.2e",
             global_min_tau, global_max_tau);
}


/**
 * @brief Calculates the optical path length for the provided segments and
 *        updates the maximum optical path length if necessary.
 * @param track The track associated with the segments
 * @param segments The segments for which the optical path length is calculated
 */
void MaxOpticalLength::onTrack(Track* track, segment* segments) {

  FP_PRECISION sin_theta = 1.0;
  Track3D* track_3d = dynamic_cast<Track3D*>(track);
  if (track_3d != NULL)
    sin_theta = sin(track_3d->getTheta());

  for (int s=0; s < track->getNumSegments(); s++) {
    FP_PRECISION length = segments[s]._length * sin_theta;
    Material* material = segments[s]._material;
    FP_PRECISION _max_sigma_t = material->getMaxSigmaT();
    FP_PRECISION tau = length * _max_sigma_t;

    //FIXME Potential but inexistent race condition
    if (tau > _max_tau) {
#pragma omp critical
      _max_tau = std::max(_max_tau, tau);
    }
    if (tau < _min_tau) {
#pragma omp critical
      _min_tau = std::min(_min_tau, tau);
    }
  }
}


/**
 * @brief Constructor for SegmentCounter calls the TraverseSegments
 *        constructor and sets the max number of segments per Track to zero.
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
 * @brief Determines the maximum number of segments per Track.
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


/**
 * @brief Turn on counting segments.
*/
void SegmentCounter::countTotalNumSegments() {
  _count_total_segments = true;
}


/**
 * @brief Get the total number of segments.
 * @return The total number of segments
*/
long SegmentCounter::getTotalNumSegments() {
  if (!_total_segments_counted)
    log_printf(ERROR, "The total number of segments have not been counted. "
               "The SegmentCounter was not instructed to count segments "
               "before execution");
  return _total_num_segments;
}


/**
 * @brief Updates the maximum number of segments per Track if a Track with a
 *        larger number of segments is observed.
 * @param track The Track whose segments are counted
 * @param segments The segments associated with the Track
 */
void SegmentCounter::onTrack(Track* track, segment* segments) {
#pragma omp critical //TODO Optimize by having a _max_num_segments per thread
  {
    if (track->getNumSegments() > _max_num_segments)
      _max_num_segments = std::max(_max_num_segments, track->getNumSegments());
  }

  if (_count_total_segments) {
#pragma omp atomic update
    _total_num_segments += track->getNumSegments();
  }
}


/**
 * @brief Constructor for SegmentSplitter calls the TraverseSegments
 *        constructor.
 * @param track_generator The TrackGenerator to pull tracking information from
 */
SegmentSplitter::SegmentSplitter(TrackGenerator* track_generator)
                               : TraverseSegments(track_generator) {
}


/**
 * @brief Splits segments stored explicity along each Track.
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
 *        larger optical path length than the maximum optical path length.
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

  /* Extract data from this segment to compute its optical length */
  for (int s = 0; s < track->getNumSegments(); s++) {
    segment* curr_segment = track->getSegment(s);
    Material* material = curr_segment->_material;
    double length = curr_segment->_length;
    long fsr_id = curr_segment->_region_id;

    /* Compute number of segments to split this segment into */
    int num_groups = material->getNumEnergyGroups();
    FP_PRECISION max_sigma_t = material->getMaxSigmaT();
    FP_PRECISION max_tau = length * max_sigma_t;
    int min_num_cuts = ceil(max_tau / max_optical_length);

    /* If the segment does not need subdivisions, go to next segment */
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
      segment new_segment;
      new_segment._material = material;
      new_segment._length = length / min_num_cuts;
      new_segment._region_id = fsr_id;

      /* Assign CMFD surface boundaries */
      if (k == 0)
        new_segment._cmfd_surface_bwd = cmfd_surface_bwd;

      if (k == min_num_cuts-1)
        new_segment._cmfd_surface_fwd = cmfd_surface_fwd;

      /* Set the starting position */
      new_segment._starting_position[0] = x_curr;
      new_segment._starting_position[1] = y_curr;
      new_segment._starting_position[2] = z_curr;
      x_curr += new_segment._length * xdir;
      y_curr += new_segment._length * ydir;
      z_curr += new_segment._length * zdir;

      /* Insert the new segment to the Track */
      track->insertSegment(s+k+1, &new_segment);
    }

    /* Remove the original segment from the Track */
    track->removeSegment(s);
  }
}


/**
 * @brief Constructor for SegmentSplitter calls the TraverseSegments
 *        constructor.
 * @param track_generator The TrackGenerator to pull tracking information from
 */
VolumeCalculator::VolumeCalculator(TrackGenerator* track_generator)
                                  : TraverseSegments(track_generator) {
}


/**
 * @brief FSR volumes are calculated and saved in the TrackGenerator's FSR
 *        volumes buffer.
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
 *        the VolumeCalculator.
 * @param track The current Track
 * @param segments The segments associated with the Track
 */
void VolumeCalculator::onTrack(Track* track, segment* segments) {
}


/**
 * @brief Constructor for CentroidGenerator calls the TraverseSegments
 *        constructor and imports references to both the FSR volumes and FSR
 *        locks arrays.
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


/**
 * @brief Destructor for the CentroidGenerator.
*/
CentroidGenerator::~CentroidGenerator() {
  int num_threads = omp_get_max_threads();
  for (int i=0; i < num_threads; i++)
    delete [] _starting_points[i];
  delete [] _starting_points;
}


/**
 * @brief Calculates the centroid of every FSR.
 * @details SegmentationKernels are created to temporarily save segments for
 *          on-the-fly methods. Then on each segment, onTrack(...) calculates
 *          the contribution to each FSR centroid and saves the centroids in
 *          the centroids array provided by setCentroids(...).
 */
void CentroidGenerator::execute() {
#pragma omp parallel
  {
    // OTF ray tracing requires segmentation of tracks
    if (_segment_formation != EXPLICIT_2D &&
        _segment_formation != EXPLICIT_3D) {
      MOCKernel* kernel = getKernel<SegmentationKernel>();
      loopOverTracks(kernel);
    }
    else
      loopOverTracks(NULL);
  }
}


/**
 * @brief Specifies an array to save calculated FSR centroids
 * @param centroids The array of pointer to FSR centroids
 */
void CentroidGenerator::setCentroids(Point** centroids) {
  _centroids = centroids;
}


/**
 * @brief Centroid contributions are calculated for every segment in the Track.
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

  /* Use local array accumulator to prevent false sharing */
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

  /* Loop over segments to accumulate contribution to centroids */
  for (int s=0; s < track->getNumSegments(); s++) {

    segment* curr_segment = &segments[s];
    long fsr = curr_segment->_region_id;
    int track_idx = curr_segment->_track_idx;

    /* Extract information */
    FP_PRECISION volume = _FSR_volumes[fsr];
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
#ifdef INTEL
#pragma omp flush
#endif

    x += cos_phi * sin_theta * curr_segment->_length;
    y += sin_phi * sin_theta * curr_segment->_length;
    z += cos_theta * curr_segment->_length;

    _starting_points[tid][track_idx].setX(x);
    _starting_points[tid][track_idx].setY(y);
    _starting_points[tid][track_idx].setZ(z);
  }
}


/**
 * @brief Constructor for LinearExpansionGenerator calls the TraverseSegments
 *        constructor, allocates memory for the linear expansion terms and
 *        initializes its own exponential evaluator.
 * @param solver the linear source solver used
 */
LinearExpansionGenerator::LinearExpansionGenerator(CPULSSolver* solver)
    : TraverseSegments(solver->getTrackGenerator()) {

  /* Import data from the Solver and TrackGenerator */
  TrackGenerator* track_generator = solver->getTrackGenerator();
  _solver = solver;
  _FSR_volumes = track_generator->getFSRVolumesBuffer();
  _FSR_locks = track_generator->getFSRLocks();
  _quadrature = track_generator->getQuadrature();
#ifndef NGROUPS
  _NUM_GROUPS = solver->getNumEnergyGroups();
#endif

  /* Determine the number of linear coefficients */
  _num_flat = 0;
  int num_rows = 1;
  if (_track_generator_3D != NULL)
    num_rows = _track_generator_3D->getMaxNumTracksPerStack();
#ifndef THREED
  _NUM_COEFFS = 3;
  if (_track_generator_3D != NULL)
    _NUM_COEFFS = 6;
#endif

  /* Reset linear source coefficients to zero */
  long size = track_generator->getGeometry()->getNumFSRs() * _NUM_COEFFS;
  _lin_exp_coeffs = new double[size]();
  size *= _NUM_GROUPS;
  _src_constants = new double[size]();

  /* Create local thread tallies */
  int num_threads = omp_get_max_threads();
  _starting_points = new Point*[num_threads];

  for (int i=0; i < num_threads; i++)
    _starting_points[i] = new Point[num_rows];

  _exp_evaluator = new ExpEvaluator();

  std::string msg = "Initializing linear source constant components";
  _progress = new Progress(_track_generator->getNumTracks(), msg, 0.1,
                    track_generator->getGeometry(), true);
}

/**
 * @brief Destructor for the LinearExpansionGenerator.
 */
LinearExpansionGenerator::~LinearExpansionGenerator() {
  int num_threads = omp_get_max_threads();
  for (int i=0; i < num_threads; i++)
    delete [] _starting_points[i];

  delete [] _starting_points;
  delete [] _src_constants;
  delete [] _lin_exp_coeffs;
  delete _exp_evaluator;
  delete _progress;
}


/**
 * @brief When executed, the LinearExpansionGenerator Kernel loops over all 
 *        Tracks to compute constant terms used to compute the linear source.
 */
void LinearExpansionGenerator::execute() {
#pragma omp parallel
  {
    // OTF ray tracing requires segmentation of tracks
    if (_segment_formation != EXPLICIT_2D &&
        _segment_formation != EXPLICIT_3D) {
      MOCKernel* kernel = getKernel<SegmentationKernel>();
      loopOverTracks(kernel);
    }
    else
      loopOverTracks(NULL);
  }

  Geometry* geometry = _track_generator->getGeometry();
  long num_FSRs = geometry->getNumFSRs();

  double* lem = _lin_exp_coeffs;
  double* ilem = _solver->getLinearExpansionCoeffsBuffer();
  int nc = _NUM_COEFFS;
  double max_ilem = 0;

  /* Invert the expansion coefficient matrix */
  if (_track_generator_3D != NULL) {
#pragma omp parallel for
    for (long r=0; r < num_FSRs; r++) {
      double det;
      det = lem[r*nc + 0] * lem[r*nc + 1] * lem[r*nc + 5] +
            lem[r*nc + 2] * lem[r*nc + 4] * lem[r*nc + 3] +
            lem[r*nc + 3] * lem[r*nc + 2] * lem[r*nc + 4] -
            lem[r*nc + 0] * lem[r*nc + 4] * lem[r*nc + 4] -
            lem[r*nc + 3] * lem[r*nc + 1] * lem[r*nc + 3] -
            lem[r*nc + 2] * lem[r*nc + 2] * lem[r*nc + 5];

      double volume = _FSR_volumes[r];
      if (std::abs(det) < MIN_DET || volume < 1e-6) {
        if (volume > 0)
          log_printf(DEBUG, "Unable to form linear source components in "
                     "source region %d : determinant %.2e volume %.2e", r, det,
                     volume);

#pragma omp atomic update
        _num_flat++;
        ilem[r*nc + 0] = 0.0;
        ilem[r*nc + 1] = 0.0;
        ilem[r*nc + 2] = 0.0;
        ilem[r*nc + 3] = 0.0;
        ilem[r*nc + 4] = 0.0;
        ilem[r*nc + 5] = 0.0;
      }
      else {

        double curr_ilem[6];

        curr_ilem[0] = (lem[r*nc + 1] * lem[r*nc + 5] -
                        lem[r*nc + 4] * lem[r*nc + 4]) / det;
        curr_ilem[1] = (lem[r*nc + 0] * lem[r*nc + 5] -
                        lem[r*nc + 3] * lem[r*nc + 3]) / det;
        curr_ilem[2] = (lem[r*nc + 3] * lem[r*nc + 4] -
                        lem[r*nc + 2] * lem[r*nc + 5]) / det;
        curr_ilem[3] = (lem[r*nc + 2] * lem[r*nc + 4] -
                        lem[r*nc + 3] * lem[r*nc + 1]) / det;
        curr_ilem[4] = (lem[r*nc + 3] * lem[r*nc + 2] -
                        lem[r*nc + 0] * lem[r*nc + 4]) / det;
        curr_ilem[5] = (lem[r*nc + 0] * lem[r*nc + 1] -
                        lem[r*nc + 2] * lem[r*nc + 2]) / det;

        /* Copy inverses */
        long ind = r*nc;
        for (int i=0; i < 6; i++) {
          ilem[ind+i] = curr_ilem[i];
          if (curr_ilem[i] > max_ilem)
#pragma omp critical
            max_ilem = std::max(max_ilem, curr_ilem[i]);
        }
      }
    }
  }
  else {
#pragma omp parallel for
    for (long r=0; r < num_FSRs; r++) {

      double det;
      det = lem[r*nc  ] * lem[r*nc + 1] - lem[r*nc + 2] * lem[r*nc + 2];

      if (std::abs(det) < MIN_DET) {
        ilem[r*nc + 0] = 0.0;
        ilem[r*nc + 1] = 0.0;
        ilem[r*nc + 2] = 0.0;
        log_printf(DEBUG, "Unable to form linear source components in "
                   "source region %d : determinant = %.2e", r, det);
#pragma omp atomic update
        _num_flat++;
      }
      else {
        ilem[r*nc + 0] =  lem[r*nc + 1] / det;
        ilem[r*nc + 1] =  lem[r*nc + 0] / det;
        ilem[r*nc + 2] = -lem[r*nc + 2] / det;
      }
    }
  }

  /* Copy the source constants to buffer */
  FP_PRECISION* src_constants_buffer = _solver->getSourceConstantsBuffer();
  long size = num_FSRs * _NUM_COEFFS * _NUM_GROUPS;
#pragma omp parallel for
  for (long i=0; i < size; i++)
    src_constants_buffer[i] = _src_constants[i];

  /* Notify user of any very large linear expansion matrix coefficient */
  if (max_ilem > 1e10)
    log_printf(INFO, "Max inverse linear expansion matrix coeff %e", max_ilem);

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


/**
 * @brief Contributions to the linear source expansion terms and constant terms
 *        are calculated for every segment in the Track.
 * @param track Track the Kernel is acting on
 * @param segments segments on that Track
 */
void LinearExpansionGenerator::onTrack(Track* track, segment* segments) {

  /* Extract track information */
  Point* start = track->getStart();
  double phi = track->getPhi();
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);
  int azim_index = track->getAzimIndex();
  int xy_index = track->getXYIndex();

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
    polar_index = track_3D->getPolarIndex();
    wgt *= _quadrature->getPolarSpacing(azim_index, polar_index)
        *_quadrature->getPolarWeight(azim_index, polar_index);
  }

  /* Loop over segments to accumulate contribution to centroids */
  Geometry* geometry = _track_generator->getGeometry();
  for (int s=0; s < track->getNumSegments(); s++) {

    /* Extract segment information */
    segment* curr_segment = &segments[s];
    long fsr = curr_segment->_region_id;
    int track_idx = curr_segment->_track_idx;
    FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
    double length = curr_segment->_length;
    double length_2 = length * length;

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

    /* Allocate a buffer for the FSR source constants on the stack */
    double thread_src_constants[_NUM_GROUPS * _NUM_COEFFS]  __attribute__
       ((aligned (VEC_ALIGNMENT)));

    /* Pre-compute non-energy dependent source constant terms */
    double vol_impact = wgt * length / volume;
    double src_constant = vol_impact * length / 2.0;

#pragma omp simd aligned(sigma_t)
    for (int g=0; g < _NUM_GROUPS; g++) {

      thread_src_constants[g] = vol_impact * xc * xc;
      thread_src_constants[_NUM_GROUPS + g] = vol_impact * yc * yc;
      thread_src_constants[2*_NUM_GROUPS + g] = vol_impact * xc * yc;

#ifndef THREED
      if (track_3D != NULL) {
#endif
        thread_src_constants[3*_NUM_GROUPS + g] = vol_impact * xc * zc;
        thread_src_constants[4*_NUM_GROUPS + g] = vol_impact * yc * zc;
        thread_src_constants[5*_NUM_GROUPS + g] = vol_impact * zc * zc;
#ifndef THREED
      }
#endif

      double tau = length * sigma_t[g];

#ifndef THREED
      if (track_3D == NULL) {
        for (int p=0; p < _quadrature->getNumPolarAngles()/2; p++) {

          double sin_theta = _quadrature->getSinTheta(azim_index, p);
          double G2_src =
              length * _exp_evaluator->computeExponentialG2(tau / sin_theta)
              * src_constant * 2 * _quadrature->getPolarWeight(azim_index, p)
              * sin_theta;

          thread_src_constants[g] += cos_phi * cos_phi * G2_src;
          thread_src_constants[_NUM_GROUPS + g] += sin_phi * sin_phi
              * G2_src;
          thread_src_constants[2*_NUM_GROUPS + g] += sin_phi * cos_phi
              * G2_src;
        }
      }
      else {
#endif
        double G2_src = _exp_evaluator->computeExponentialG2(tau) *
            length * src_constant;

        thread_src_constants[g] += cos_phi * cos_phi * G2_src * sin_theta
             * sin_theta;
        thread_src_constants[_NUM_GROUPS + g] += sin_phi * sin_phi * G2_src
             * sin_theta * sin_theta;
        thread_src_constants[2*_NUM_GROUPS + g] += sin_phi * cos_phi * G2_src
             * sin_theta * sin_theta;
        thread_src_constants[3*_NUM_GROUPS + g] += cos_phi * cos_theta * G2_src
             * sin_theta;
        thread_src_constants[4*_NUM_GROUPS + g] += sin_phi * cos_theta * G2_src
             * sin_theta;
        thread_src_constants[5*_NUM_GROUPS + g] += cos_theta * cos_theta * G2_src;
#ifndef THREED
      }
#endif
    }

    /* Set the lock for this FSR */
    omp_set_lock(&_FSR_locks[fsr]);

    _lin_exp_coeffs[fsr*_NUM_COEFFS] += wgt * length / volume *
        (xc * xc + pow(cos_phi * sin_theta * length, 2) / 12.0);
    _lin_exp_coeffs[fsr*_NUM_COEFFS + 1] += wgt * length / volume *
        (yc * yc + pow(sin_phi * sin_theta * length, 2) / 12.0);
    _lin_exp_coeffs[fsr*_NUM_COEFFS + 2] += wgt * length / volume *
        (xc * yc + sin_phi * cos_phi * pow(sin_theta * length, 2) / 12.0);

    if (track_3D != NULL) {
      _lin_exp_coeffs[fsr*_NUM_COEFFS + 3] += wgt * length / volume *
          (xc * zc + cos_phi * cos_theta * sin_theta * pow(length, 2) / 12.0);
      _lin_exp_coeffs[fsr*_NUM_COEFFS + 4] += wgt * length / volume *
          (yc * zc + sin_phi * cos_theta * sin_theta * pow(length, 2) / 12.0);
      _lin_exp_coeffs[fsr*_NUM_COEFFS + 5] += wgt * length / volume *
          (zc * zc + pow(cos_theta * length, 2) / 12.0);
    }

    /* Set the source constants for all groups and coefficients */
#pragma omp simd
    for (int g=0; g < _NUM_GROUPS; g++) {
      for (int i=0; i < _NUM_COEFFS; i++)
        _src_constants[fsr*_NUM_GROUPS*_NUM_COEFFS + i*_NUM_GROUPS + g] +=
          thread_src_constants[i*_NUM_GROUPS + g];
    }

    /* Unset the lock for this FSR */
    omp_unset_lock(&_FSR_locks[fsr]);
#ifdef INTEL
#pragma omp flush
#endif
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
 *        constructor and initializes the associated CPUSolver to NULL.
 * @param cpu_solver The CPUSolver to use for propagating angular fluxes
 */
TransportSweep::TransportSweep(CPUSolver* cpu_solver)
    : TraverseSegments(cpu_solver->getTrackGenerator()) {

  _cpu_solver = cpu_solver;
  _ls_solver = dynamic_cast<CPULSSolver*>(cpu_solver);
  TrackGenerator* track_generator = cpu_solver->getTrackGenerator();
  _geometry = _track_generator->getGeometry();
#ifndef NGROUPS
  _NUM_GROUPS = _geometry->getNumEnergyGroups();
#endif
}


/**
 * @brief Destructor for the TransportSweep.
 */
TransportSweep::~TransportSweep() {
}



/**
 * @brief MOC equations are applied to every segment in the TrackGenerator
 * @details SegmentationKernels are allocated to temporarily save segments. Then
 *          onTrack(...) applies the MOC equations to each segment and
 *          transfers boundary fluxes for the corresponding Track.
 */
void TransportSweep::execute() {
#pragma omp parallel
  {
    // OTF ray tracing requires segmentation of tracks
    if (_segment_formation != EXPLICIT_2D &&
        _segment_formation != EXPLICIT_3D) {
      MOCKernel* kernel = getKernel<SegmentationKernel>();
      loopOverTracks(kernel);
    }
    else
      loopOverTracks(NULL);
  }
}


/**
 * @brief Applies the MOC equations the Track and segments.
 * @details The MOC equations are applied to each segment, attenuating the
 *          Track's angular flux and tallying FSR contributions. Finally,
 *          Track boundary fluxes are transferred.
 * @param track The Track for which the angular flux is attenuated and
 *        transferred
 * @param segments The segments over which the MOC equations are applied
 */
void TransportSweep::onTrack(Track* track, segment* segments) {

  /* Get the thread number */
  int tid = omp_get_thread_num();

  /* Extract Track information */
  long track_id = track->getUid();
  int azim_index = track->getAzimIndex();
  int xy_index = track->getXYIndex();
  int num_segments = track->getNumSegments();
  float* track_flux;

  /* Extract the polar index and quadrature weight if a 3D track */
  int polar_index = 0;
  FP_PRECISION weight = 1;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL) {
    polar_index = track_3D->getPolarIndex();
    weight = _track_generator->getQuadrature()->getWeightInline(azim_index,
                                                                polar_index);
  }

  /* Compute unit vector if necessary */
  FP_PRECISION direction[3];
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

  /* Allocate a temporary flux buffer on the stack (free) and initialize it */
  /* Select right size for buffer */
#ifndef NGROUPS
  int _NUM_GROUPS = _cpu_solver->getNumEnergyGroups();
#endif
#ifndef LINEARSOURCE
  int num_moments = 1;
  if (_ls_solver != NULL)
    num_moments = 4;
#else
  const int num_moments = 4;
#endif
  int vec_alignment = VEC_ALIGNMENT / sizeof(FP_PRECISION);
  int num_groups_aligned = (_NUM_GROUPS / vec_alignment +
                            (_NUM_GROUPS % vec_alignment != 0)) * vec_alignment;

  /* Allocate an aligned buffer on the stack */
  FP_PRECISION fsr_flux[num_moments * num_groups_aligned] __attribute__
       ((aligned (VEC_ALIGNMENT)));
  memset(fsr_flux, 0, num_moments * num_groups_aligned * sizeof(FP_PRECISION));
  FP_PRECISION* fsr_flux_x = &fsr_flux[num_groups_aligned];
  FP_PRECISION* fsr_flux_y = &fsr_flux[2*num_groups_aligned];
  FP_PRECISION* fsr_flux_z = &fsr_flux[3*num_groups_aligned];

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {

    /* Get the forward track flux */
    segment* curr_segment = &segments[s];
    long curr_track_id = track_id + curr_segment->_track_idx;
    track_flux = _cpu_solver->getBoundaryFlux(curr_track_id, true);
    long fsr_id = curr_segment->_region_id;

    /* Apply MOC equations */
#ifndef LINEARSOURCE
    if (_ls_solver == NULL)
      _cpu_solver->tallyScalarFlux(curr_segment, azim_index, fsr_flux,
                                   track_flux);
    else
#endif
      _ls_solver->tallyLSScalarFlux(curr_segment, azim_index, polar_index,
                                    fsr_flux, fsr_flux_x, fsr_flux_y, 
                                    fsr_flux_z, track_flux, direction);

    /* Accumulate contribution of segments to scalar flux before changing fsr */
    if (s < num_segments - 1 && fsr_id != (&segments[s+1])->_region_id) {
#ifndef LINEARSOURCE
      if (_ls_solver == NULL)
        _cpu_solver->accumulateScalarFluxContribution(fsr_id, weight, fsr_flux);
      else
#endif
        _ls_solver->accumulateLinearFluxContribution(fsr_id, weight, fsr_flux);
    }

    /* Tally the current for CMFD */
    _cpu_solver->tallyCurrent(curr_segment, azim_index, polar_index,
                              track_flux, true);
  }

#ifndef ONLYVACUUMBC
  /* Transfer boundary angular flux to outgoing Track */
  for (int i=0; i <= max_track_index; i++) {
    track_flux = _cpu_solver->getBoundaryFlux(track_id+i, true);
    _cpu_solver->transferBoundaryFlux(tracks_array[i], azim_index, polar_index,
                                      true, track_flux);
  }
#endif

  /* Reverse the direction */
  for (int i=0; i<3; i++)
    direction[i] *= -1;

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {

    /* Get the backward track flux */
    segment* curr_segment = &segments[s];
    long curr_track_id = track_id + curr_segment->_track_idx;
    track_flux = _cpu_solver->getBoundaryFlux(curr_track_id, false);
    long fsr_id = curr_segment->_region_id;

    /* Apply MOC equations */
#ifndef LINEARSOURCE
    if (_ls_solver == NULL)
      _cpu_solver->tallyScalarFlux(curr_segment, azim_index, fsr_flux,
                                   track_flux);
    else
#endif
      _ls_solver->tallyLSScalarFlux(curr_segment, azim_index, polar_index,
                                    fsr_flux, fsr_flux_x, fsr_flux_y, 
                                    fsr_flux_z, track_flux, direction);

    /* Accumulate contribution of segments to scalar flux before changing fsr */
    if (s == 0 || fsr_id != (&segments[s-1])->_region_id) {
#ifndef LINEARSOURCE
      if (_ls_solver == NULL)
        _cpu_solver->accumulateScalarFluxContribution(fsr_id, weight, fsr_flux);
      else
#endif
        _ls_solver->accumulateLinearFluxContribution(fsr_id, weight, fsr_flux);
    }

    /* Tally the current for CMFD */
    _cpu_solver->tallyCurrent(curr_segment, azim_index, polar_index,
                              track_flux, false);
  }

#ifndef ONLYVACUUMBC
  /* Transfer boundary angular flux to outgoing Track */
  for (int i=0; i <= max_track_index; i++) {
    track_flux = _cpu_solver->getBoundaryFlux(track_id+i, false);
    _cpu_solver->transferBoundaryFlux(tracks_array[i], azim_index, polar_index,
                                      false, track_flux);
  }
#endif
}


/**
 * @brief Constructor for DumpSegments calls the TraverseSegments
 *        constructor and initializes the output FILE to NULL.
 * @param track_generator The TrackGenerator to pull tracking information from
 */
DumpSegments::DumpSegments(TrackGenerator* track_generator)
                           : TraverseSegments(track_generator) {
  _out = NULL;
}


/**
 * @brief Writes all tracking information to file.
 * @details SegmentationKernels are created to temporarily store segments for
 *          on-the-fly method. For each Track, onTrack(...) writes the tracking
 *          information to file.
 */
void DumpSegments::execute() {

  // OTF ray tracing requires segmentation of tracks
  if (_segment_formation != EXPLICIT_2D &&
      _segment_formation != EXPLICIT_3D) {
    MOCKernel* kernel = getKernel<SegmentationKernel>();
    loopOverTracks(kernel);
  }
  else
    loopOverTracks(NULL);
}


/**
 * @brief Sets the file in which to write tracking information.
 * @param out the file in which to write tracking information
 */
void DumpSegments::setOutputFile(FILE* out) {
  _out = out;
}


/**
 * @brief Writes tracking information to file for a Track and associated
 *        segments.
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
 *        constructor and initializes the input FILE to NULL.
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
 * @brief Sets the input file to read in tracking information.
 * @param input the tracking file
 */
void ReadSegments::setInputFile(FILE* input) {
  _in = input;
}


/**
 * @brief Saves tracking information to the corresponding Track explicitly.
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
    segment curr_segment;
    curr_segment._length = length;
    curr_segment._material = materials[material_id];
    curr_segment._region_id = region_id;
    curr_segment._track_idx = track_idx;
    curr_segment._starting_position[0] = start_x;
    curr_segment._starting_position[1] = start_y;
    curr_segment._starting_position[2] = start_z;

    /* Import CMFD-related data if needed */
    if (cmfd != NULL) {
      int cmfd_surface_fwd;
      ret = geometry->twiddleRead(&cmfd_surface_fwd, sizeof(int), 1, _in);
      curr_segment._cmfd_surface_fwd = cmfd_surface_fwd;
      int cmfd_surface_bwd;
      ret = geometry->twiddleRead(&cmfd_surface_bwd, sizeof(int), 1, _in);
      curr_segment._cmfd_surface_bwd = cmfd_surface_bwd;
    }

    /* Add this segment to the Track */
    track->addSegment(&curr_segment);
  }
}


/**
 * @brief Constructor for TransportSweepOTF calls the TraverseSegments 
 *        constructor.
 * @param track_generator Track generator to generate the tracks
 */
TransportSweepOTF::TransportSweepOTF(TrackGenerator* track_generator)
                                   : TraverseSegments(track_generator) {
  _cpu_solver = NULL;
}


/**
 * @brief When executed, the Kernel loops over all tracks, both generating them
 *        and solving the MOC equations.
 */
void TransportSweepOTF::execute() {
#pragma omp parallel
  {
    TransportKernel kernel(_track_generator);
    kernel.setCPUSolver(_cpu_solver);
    loopOverTracksByStackTwoWay(&kernel);
  }
}


/**
 * @brief Set the solver for the OTF TransportSweep
 * @param cpu_solver Solver to use to solve MOC equations
 */
void TransportSweepOTF::setCPUSolver(CPUSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


/**
 * @brief Placeholder, an onTrack routine is not required when performing
 *        track generation and transport simultaneously.
 * @param track the Track of interest
 * @param segments array of segments on that track
 */
void TransportSweepOTF::onTrack(Track* track, segment* segments) {
}


/**
 * @brief Constructor for the RecenterSegments calls the TraverseSegments
 *        constructor and sets the track generator.
 * @param track_generator Track generator to obtain and re-center segments.
 */
RecenterSegments::RecenterSegments(TrackGenerator* track_generator)
                                   : TraverseSegments(track_generator) {
  _geometry = _track_generator->getGeometry();
}


/**
 * @brief When executed, the Kernel loops over all Tracks to recenter their
 *        segments.
 */
void RecenterSegments::execute() {
#pragma omp parallel
  {
    // OTF ray tracing requires segmentation of tracks
    if (_segment_formation != EXPLICIT_2D &&
        _segment_formation != EXPLICIT_3D) {
      MOCKernel* kernel = getKernel<SegmentationKernel>();
      loopOverTracks(kernel);
    }
    else
      loopOverTracks(NULL);
  }
}


/**
 * @brief Loops over all segments provided, obtain their region (FSR) centroid
 *        and re-center the segment.
 * @param track Track which contains the segments
 * @param segments array of segments to re-center
 */
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
 * @brief Constructor for PrintSegments calls the TraverseSegments
 *        constructor and initializes the output FILE to NULL.
 * @param track_generator The TrackGenerator to pull tracking information from
 */
PrintSegments::PrintSegments(TrackGenerator* track_generator)
                           : TraverseSegments(track_generator) {
  _out = NULL;
}


/**
 * @brief Writes all tracking information to file.
 * @details SegmentationKernels are created to temporarily store segments for
 *          on-the-fly methods. For each Track, onTrack(...) writes the tracking
 *          information to file.
 //FIXME debug ?
 */
void PrintSegments::execute() {

  // OTF ray tracing requires segmentation of tracks
  if (_segment_formation != EXPLICIT_2D &&
      _segment_formation != EXPLICIT_3D) {
    MOCKernel* kernel = getKernel<SegmentationKernel>();
    loopOverTracks(kernel);
  }
  else
    loopOverTracks(NULL);
}


/**
 * @brief Sets the file which to write tracking information.
 * @param out the file which to write tracking infmormation
 */
void PrintSegments::setOutputFile(FILE* out) {
  _out = out;
}


/**
 * @brief Writes tracking information to file for a Track and associated
 *        segments.
 * @param track The Track whose information is written to file
 * @param segments The segments associated with the Track whose information is
 *        written to file
 //FIXME debug ?
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
