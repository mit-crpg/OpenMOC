/**
 * @file TraverseSegments.h
 * @brief A TraverseSegments object
 * @date February 15, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef TRAVERSE_SEGMENTS_H_
#define TRAVERSE_SEGMENTS_H_

#ifdef SWIG
#include "Python.h"
#endif
#include "Track.h"
#include "Track3D.h"
#include "Geometry.h"
#include "TrackGenerator3D.h"


/**
 * @class TraverseSegments TraverseSegments.h "src/TraverseSegments.h"
 * @brief An TraverseSegments object defines how to loop over Tracks given
 *        various different segmentation schemes and how to apply algorithms
 *        to the Tracks and associated segments.
 * @details A TraverseSegments object sketches how to loop over Tracks for
 *          various different segmentation schemes such as 2D explicit, 3D
 *          explicit, and on-the-fly ray tracing. A TraverseSegments object
 *          is an abstract class meant to be extended by classes defined in
 *          TrackTraversingAlgorithms.h. This parent class's main purpose is to
 *          abstract the looping procedure and apply a function onTrack(...) to
 *          each Track and apply the supplied MOCKernel to each segment. If
 *          NULL is provided for the MOCKernel, only the functionality defined
 *          in onTrack(...) is applied to each Track.
 */
class TraverseSegments {

private:

  /* Functions defining how to loop over Tracks */
  void loopOverTracks2D(MOCKernel* kernel);
  void loopOverTracksExplicit(MOCKernel* kernel);
  void loopOverTracksByTrackOTF(MOCKernel* kernel);
  void loopOverTracksByStackOTF(MOCKernel* kernel);

  /* Functions defining how to traverse segments */
  void traceSegmentsExplicit(Track* track, MOCKernel* kernel);
  void traceSegmentsOTF(Track* flattened_track, Point* start,
                        double theta, MOCKernel* kernel);
  void traceStackOTF(Track* flattened_track, int polar_index,
                     MOCKernel* kernel);

  void traceStackTwoWay(Track* flattened_track, int polar_index,
                        TransportKernel* kernel);


  int findMeshIndex(double* values, int size, double val, int sign);

protected:

  /** Pointer to the associated TrackGenerator */
  TrackGenerator* _track_generator;
  TrackGenerator3D* _track_generator_3D;

  /** A pointer to the associated global z-mesh (if applicable) */
  double* _global_z_mesh;

  /** The size of the global z-mesh */
  int _mesh_size;

  /** The type of segmentation used for segment formation */
  segmentationType _segment_formation;

  TraverseSegments(TrackGenerator* track_generator);
  virtual ~TraverseSegments();

  /* Functions defining how to loop over and operate on Tracks */
  void loopOverTracks(MOCKernel* kernel);
  virtual void onTrack(Track* track, segment* segments) = 0;

  //FIXME Rework function calls to make this private
  void loopOverTracksByStackTwoWay(TransportKernel* kernel);

  /* Returns a kernel of the requested type */
  template <class KernelType>
  MOCKernel* getKernel() {

    /* Check for segmentation kernel in explicit methods */
    if ((typeid(KernelType) != typeid(SegmentationKernel)) ||
        ((_segment_formation != EXPLICIT_2D) &&
        (_segment_formation != EXPLICIT_3D))) {

      /* Allocate kernel */
      MOCKernel* kernel = new KernelType(_track_generator);
      return kernel;
    }
    else
      return NULL;
  }

public:
  virtual void execute() = 0;

};

#endif
