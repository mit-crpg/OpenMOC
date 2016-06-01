#include "TraverseTracks.h"

/**
 * @brief Constructor for the TravseSegments class assigns the TrackGenerator
 *        and pulls relevant information from it.
 */
TraverseTracks::TraverseTracks(TrackGenerator* track_generator) {

  /* Save the track generator */
  _track_generator = track_generator;

  /* Determine the type of segment formation used */
  _segment_formation = track_generator->getSegmentFormation();
}


/**
 * @brief Destructor for TraverseTracks
 */
TraverseTracks::~TraverseTracks() {
}


/**
 * @brief Loops over Tracks, applying the provided kernel to all segments and
 *        the functionality described in onTrack(...) to all Tracks.
 * @details The segment formation method imported from the TrackGenerator
 *          during construction is used to redirect to the appropriate looping
 *          scheme.
 * @param kernel The MOCKernel to apply to all segments
 */
void TraverseTracks::loopOverTracks(MOCKernel* kernel) {
  switch (_segment_formation) {
    case EXPLICIT_2D:
      loopOverTracks2D(kernel);
      break;
    default:
      log_printf(ERROR, "Segment formation type not currently supported");
  }
}


/**
 * @brief Loops over all explicit 2D Tracks
 * @details The onTrack(...) function is applied to all 2D Tracks and the
 *          specified kernel is applied to all segments. If NULL is provided
 *          for the kernel, only the onTrack(...) functionality is applied.
 * @param kernel The MOCKernel to apply to all segments
 */
void TraverseTracks::loopOverTracks2D(MOCKernel* kernel) {

  /* Loop over all parallel tracks for each azimuthal angle */
  Track** tracks_2D = _track_generator->getTracks();
  int num_azim = _track_generator->getNumAzim();
  for (int a=0; a < num_azim/2; a++) {
    int num_xy = _track_generator->getNumX(a) + _track_generator->getNumY(a);
#pragma omp for
    for (int i=0; i < num_xy; i++) {

      Track* track_2D = &tracks_2D[a][i];

      /* Apply the kernel to segments if necessary */
      if (kernel != NULL) {
        kernel->newTrack(track_2D);
        traceSegmentsExplicit(track_2D, kernel);
      }

      /* Operate on the Track */
      segment* segments = track_2D->getSegments();
      onTrack(track_2D, segments);
    }
  }
}


/**
 * @brief Loops over segments in a Track when segments are explicitly generated
 * @details All segments in the provided Track are looped over and the provided
 *          MOCKernel is applied to them.
 * @param track The Track whose segments will be traversed
 * @param kernel The kernel to apply to all segments
 */
void TraverseTracks::traceSegmentsExplicit(Track* track, MOCKernel* kernel) {
  for (int s=0; s < track->getNumSegments(); s++) {
    segment* seg = track->getSegment(s);
    kernel->execute(seg->_length, seg->_material, seg->_region_id,
                    seg->_cmfd_surface_fwd, seg->_cmfd_surface_bwd);
  }
}


/**
 * @brief Dummy function for default onTrack implementation
 */
void TraverseTracks::onTrack(Track* track, segment* segments) {
}
