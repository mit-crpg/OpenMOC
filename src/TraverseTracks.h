/**
 * @file TraverseTracks.h
 * @brief A TraverseTracks object
 * @date February 15, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef TRAVERSE_TRACKS_H_
#define TRAVERSE_TRACKS_H_

#ifdef SWIG
#include "Python.h"
#endif
#include "Track.h"
#include "Geometry.h"
#include "TrackGenerator.h"
#include "MOCKernel.h"


/**
 * @class TraverseTracks TraverseTracks.h "src/TraverseTracks.h"
 * @brief An TraverseTracks object defines how to loop over Tracks given
 *        various different segmentation schemes and how to apply algorithms
 *        to the Tracks and associated segments.
 * @details A TraverseTracks object sketches how to loop over Tracks for
 *          various different segmentation schemes such as 2D explicit, 3D
 *          explicit, and on-the-fly ray tracing. A TraverseTracks object
 *          is an abstract class meant to be extended by classes defined in
 *          TrackTraversingAlgorithms.h. This parent class's main purpose is to
 *          abstract the looping procedure and apply a function onTrack(...) to
 *          each Track and apply supplied MOCKernels to each segment. If NULL
 *          is provided for the MOCKernels, only the functionality defined in
 *          onTrack(...) is applied to each Track.
 */
class TraverseTracks {

private:

  /* Functions defining how to loop over Tracks */
  void loopOverTracks2D(MOCKernel* kernel);

  /* Functions defining how to traverse segments */
  void traceSegmentsExplicit(Track* track, MOCKernel* kernel);

protected:

  /** Pointer to the associated TrackGenerator */
  TrackGenerator* _track_generator;

  /** The type of segmentation used for segment formation */
  segmentationType _segment_formation;

  TraverseTracks(TrackGenerator* track_generator);
  virtual ~TraverseTracks();

  /* Functions defining how to loop over and operate on Tracks */
  void loopOverTracks(MOCKernel* kernel);
  virtual void onTrack(Track* track, segment* segments);

public:

  virtual void execute() = 0;

};

#endif
