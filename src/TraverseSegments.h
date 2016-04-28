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
 *          each Track and apply supplied MOCKernels to each segment. If NULL
 *          is provided for the MOCKernels, only the functionality defined in
 *          onTrack(...) is applied to each Track.
 */
class TraverseSegments {

private:

  /* Functions defining how to loop over Tracks */
  void loopOverTracks2D(MOCKernel** kernels);
  void loopOverTracksExplicit(MOCKernel** kernels);
  void loopOverTracksByTrackOTF(MOCKernel** kernels);
  void loopOverTracksByStackOTF(MOCKernel** kernels);

  /* Functions defining how to traverse segments */
  void traceSegmentsExplicit(Track* track, MOCKernel* kernel);
  void traceSegmentsOTF(Track* flattened_track, Point* start,
                        double theta, MOCKernel* kernel);
  void traceStackOTF(Track* flattened_track, int polar_index,
                     MOCKernel** kernels);


  int findMeshIndex(FP_PRECISION* values, int size, FP_PRECISION val,
                    int sign);

protected:

  /** Pointer to the associated TrackGenerator */
  TrackGenerator* _track_generator;
  TrackGenerator3D* _track_generator_3D;

  /** A pointer to the associated global z-mesh (if applicable) */
  FP_PRECISION* _global_z_mesh;

  /** The size of the global z-mesh */
  int _mesh_size;

  /** The type of segmentation used for segment formation */
  segmentationType _segment_formation;

  TraverseSegments(TrackGenerator* track_generator);
  virtual ~TraverseSegments();

  /* Functions defining how to loop over and operate on Tracks */
  void loopOverTracks(MOCKernel** kernels);
  virtual void onTrack(Track* track, segment* segments) = 0;

  /* Returns a matrix of kernels of the requested type */
  template <class KernelType>
  MOCKernel** getKernels() {

    /* Check for segmentation kernels in explicit methods */
    if ((typeid(KernelType) != typeid(SegmentationKernel)) ||
        ((_segment_formation != EXPLICIT_2D) &&
        (_segment_formation != EXPLICIT_3D))) {

      /* Allocate kernels */
      int num_segment_matrix_rows = 1;
      MOCKernel** kernels = new MOCKernel*[num_segment_matrix_rows];
      for (int z=0; z < num_segment_matrix_rows; z++)
        kernels[z] = new KernelType(_track_generator, z);
      return kernels;
    }
    else
      return NULL;
  }

public:
  virtual void execute() = 0;

};

#endif
