#include "TrackTraversingAlgorithms.h"
#include "CPUSolver.h"
#include "Quadrature.h"


/**
 * @brief Constructor for MaxOpticalLength calls the TraverseTracks
 *        constructor and sets the max optical path length to zero
 * @param track_generator The TrackGenerator to pull tracking information from
 */
MaxOpticalLength::MaxOpticalLength(TrackGenerator* track_generator)
                                 : TraverseTracks(track_generator) {
  _max_tau = 0;
}


/**
 * @brief Determines the maximum optical path length for the TrackGenerator
 *        provided during construction
 * @details The maximum optical path length is initialized to infinity for
 *          segmentation within the TrackGenerator. Then Tracks are traversed
 *          and onTrack(...) is applied, calculating a maximum optical length
 *          and setting it on the TrackGenerator.
*/
void MaxOpticalLength::execute() {
  FP_PRECISION infinity = std::numeric_limits<FP_PRECISION>::max();
  _track_generator->setMaxOpticalLength(infinity);
#pragma omp parallel
  {
    loopOverTracks(NULL);
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
 * @brief Constructor for SegmentSplitter calls the TraverseTracks
 *        constructor
 * @param track_generator The TrackGenerator to pull tracking information from
 */
SegmentSplitter::SegmentSplitter(TrackGenerator* track_generator)
                               : TraverseTracks(track_generator) {
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
 * @brief Constructor for VolumeCalculator calls the TraverseTracks
 *        constructor
 * @param track_generator The TrackGenerator to pull tracking information
 */
VolumeCalculator::VolumeCalculator(TrackGenerator* track_generator)
                                  : TraverseTracks(track_generator) {
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
    VolumeKernel kernel(_track_generator);
    loopOverTracks(&kernel);
  }
}


/**
 * @brief Constructor for CentroidGenerator calls the TraverseTracks
 *        constructor and imports refernces to both the FSR volumes and FSR
 *        locks arrays
 * @param track_generator The TrackGenerator to pull tracking information from
 */
CentroidGenerator::CentroidGenerator(TrackGenerator* track_generator)
                                  : TraverseTracks(track_generator) {

  _FSR_volumes = track_generator->getFSRVolumes();
  _FSR_locks = track_generator->getFSRLocks();
  _quadrature = track_generator->getQuadrature();
}


/**
 * @brief Calculates the centroid of every FSR
 * @details On each segment, onTrack(...) calculates the contribution to each
 *          FSR centroid and saves the centroids in the centroids array
 *          provided by setCentroids(...).
 */
void CentroidGenerator::execute() {
#pragma omp parallel
  {
    loopOverTracks(NULL);
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

  /* Compute the Track cross-sectional area */
  double wgt = _quadrature->getAzimSpacing(azim_index)
      * _quadrature->getAzimWeight(azim_index);

  /* Pre-compute azimuthal angles for efficiency */
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);

  /* Extract starting points */
  double x = track->getStart()->getX();
  double y = track->getStart()->getY();

  /* Loop over segments to accumlate contribution to centroids */
  for (int s=0; s < track->getNumSegments(); s++) {

    segment* curr_segment = &segments[s];
    int fsr = curr_segment->_region_id;

    /* Set the lock for this FSR */
    omp_set_lock(&_FSR_locks[fsr]);

    _centroids[fsr]->
        setX(_centroids[fsr]->getX() + wgt *
        (x + cos_phi * curr_segment->_length / 2.0)
        * curr_segment->_length / _FSR_volumes[fsr]);

    _centroids[fsr]->
        setY(_centroids[fsr]->getY() + wgt *
        (y + sin_phi * curr_segment->_length / 2.0)
        * curr_segment->_length / _FSR_volumes[fsr]);

    /* Unset the lock for this FSR */
    omp_unset_lock(&_FSR_locks[fsr]);

    x += cos_phi * curr_segment->_length;
    y += sin_phi * curr_segment->_length;
  }
}



/**
 * @brief Constructor for TransportSweep calls the TraverseTracks
 *        constructor and initializes the associated CPUSolver to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
TransportSweep::TransportSweep(TrackGenerator* track_generator)
                              : TraverseTracks(track_generator) {
  _cpu_solver = NULL;
}


/**
 * @brief MOC equations are applied to every segment in the TrackGenerator
 * @details onTrack(...) applies the MOC equations to each segment and
 *          transfers boundary fluxes for the corresponding Track.
 */
void TransportSweep::execute() {
#pragma omp parallel
  {
    loopOverTracks(NULL);
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

  /* Allocate temporary FSR flux locally */
  int num_groups = _track_generator->getGeometry()->getNumEnergyGroups();
  FP_PRECISION thread_fsr_flux[num_groups];

  /* Extract Track information */
  int track_id = track->getUid();
  int azim_index = track->getAzimIndex();
  int num_segments = track->getNumSegments();
  FP_PRECISION* track_flux;

  /* Correct azimuthal index to first octant */
  Quadrature* quad = _track_generator->getQuadrature();

  /* Get the forward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, true);

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {
    segment* curr_segment = &segments[s];
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, track_flux,
                                 thread_fsr_flux);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, true);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, true, track_flux);

  /* Get the backward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, false);

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {
    segment* curr_segment = &segments[s];
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, track_flux,
                                 thread_fsr_flux);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, false);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, false, track_flux);
}


/**
 * @brief Constructor for DumpSegments calls the TraverseTracks
 *        constructor and initializes the output FILE to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
DumpSegments::DumpSegments(TrackGenerator* track_generator)
                           : TraverseTracks(track_generator) {
  _out = NULL;
}


/**
 * @brief Wrties all tracking information to file
 * @details For each Track, onTrack(...) writes the tracking information to
 *          file.
 */
void DumpSegments::execute() {
  loopOverTracks(NULL);
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

  /* Get data for this Track */
  double x0 = track->getStart()->getX();
  double y0 = track->getStart()->getY();
  double z0 = track->getStart()->getZ();
  double x1 = track->getEnd()->getX();
  double y1 = track->getEnd()->getY();
  double z1 = track->getEnd()->getZ();
  double phi = track->getPhi();
  int azim_angle_index = track->getAzimIndex();
  int num_segments = track->getNumSegments();

  /* Write data for this Track to the Track file */
  fwrite(&x0, sizeof(double), 1, _out);
  fwrite(&y0, sizeof(double), 1, _out);
  fwrite(&z0, sizeof(double), 1, _out);
  fwrite(&x1, sizeof(double), 1, _out);
  fwrite(&y1, sizeof(double), 1, _out);
  fwrite(&z1, sizeof(double), 1, _out);
  fwrite(&phi, sizeof(double), 1, _out);
  fwrite(&azim_angle_index, sizeof(int), 1, _out);
  fwrite(&num_segments, sizeof(int), 1, _out);

  /* Get CMFD mesh object */
  Cmfd* cmfd = _track_generator->getGeometry()->getCmfd();
  Quadrature* quad = _track_generator->getQuadrature();

  /* Loop over all segments for this Track */
  for (int s=0; s < num_segments; s++) {

    /* Get data for this segment */
    segment* curr_segment = &segments[s];
    FP_PRECISION length = curr_segment->_length;
    int material_id = curr_segment->_material->getId();
    int region_id = curr_segment->_region_id;

    /* Write data for this segment to the Track file */
    fwrite(&length, sizeof(double), 1, _out);
    fwrite(&material_id, sizeof(int), 1, _out);
    fwrite(&region_id, sizeof(int), 1, _out);

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
 * @brief Constructor for ReadSegments calls the TraverseTracks
 *        constructor and initializes the input FILE to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
ReadSegments::ReadSegments(TrackGenerator* track_generator)
                           : TraverseTracks(track_generator) {
  _in = NULL;
  _quadrature = track_generator->getQuadrature();
  _num_azim_2 = track_generator->getNumAzim() / 2;
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
void ReadSegments::setInputFile(FILE* in) {
  _in = in;
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
  double x0, y0, z0, x1, y1, z1;
  ret = fread(&x0, sizeof(double), 1, _in);
  ret = fread(&y0, sizeof(double), 1, _in);
  ret = fread(&z0, sizeof(double), 1, _in);
  ret = fread(&x1, sizeof(double), 1, _in);
  ret = fread(&y1, sizeof(double), 1, _in);
  ret = fread(&z1, sizeof(double), 1, _in);
  double phi;
  ret = fread(&phi, sizeof(double), 1, _in);
  int azim_angle_index;
  ret = fread(&azim_angle_index, sizeof(int), 1, _in);
  int num_segments;
  ret = fread(&num_segments, sizeof(int), 1, _in);

  /* Initialize a Track with this data */
  track->setValues(x0, y0, z0, x1, y1, z1, phi);
  track->setAzimAngleIndex(azim_angle_index);
  if (azim_angle_index < _num_azim_2 / 2)
    _quadrature->setPhi(phi, azim_angle_index);

  /* Loop over all segments in this Track */
  for (int s=0; s < num_segments; s++) {

    /* Import data for this segment from Track file */
    double length;
    ret = fread(&length, sizeof(double), 1, _in);
    int material_id;
    ret = fread(&material_id, sizeof(int), 1, _in);
    int region_id;
    ret = fread(&region_id, sizeof(int), 1, _in);

    /* Initialize segment with the data */
    segment* curr_segment = new segment;
    curr_segment->_length = length;
    curr_segment->_material = materials[material_id];
    curr_segment->_region_id = region_id;

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


