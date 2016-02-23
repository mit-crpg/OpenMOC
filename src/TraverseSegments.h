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
#include "Track2D.h"
#include "Track3D.h"
#include "Geometry.h"
#include "TrackGenerator.h"


class TraverseSegments {

private:

  // descriptions
  void loopOverTracks2D(MOCKernel** kernels);
  void loopOverTracksExplicit(MOCKernel** kernels);
  void loopOverTracksByTrackOTF(MOCKernel** kernels);
  void loopOverTracksByStackOTF(MOCKernel** kernels);

  void traceSegmentsExplicit(Track* track, MOCKernel* kernel);
  void traceSegmentsOTF(Track* flattened_track, Point* start,
                        double theta, MOCKernel* kernel);
  void traceStackOTF(Track* flattened_track, int polar_index,
                     MOCKernel** kernels);


  int binarySearch(FP_PRECISION* values, int size, FP_PRECISION val, int sign);

protected:

  /** Pointer to the associated TrackGenerator */
  TrackGenerator* _track_generator;

  //descriptions
  segmentationType _segment_formation;
  FP_PRECISION* _global_z_mesh;
  int _mesh_size;

  // descriptions
  TraverseSegments(TrackGenerator* track_generator);
  virtual ~TraverseSegments();

  void loopOverTracks(MOCKernel** kernels);
  virtual void onTrack(Track* track, segment* segments) = 0;

  // description
  template <class KernelType>
  MOCKernel** getKernels() {
    int num_rows = _track_generator->getNumRows();
    MOCKernel** kernels = new MOCKernel*[num_rows];
    for (int z=0; z < num_rows; z++)
      kernels[z] = new KernelType(_track_generator, z);
    return kernels;
  }

public:
  virtual void execute() = 0;

};

#endif
