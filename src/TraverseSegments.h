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
  CounterKernel* initializeKernel<CounterKernel>();
  VolumeKernel* initializeKernel<VolumeKernel>();
  SegmentationKernel* initializeKernel<SegmentationKernel>();

  void loopOverTracks2D();
  void loopOverTracksExplicit();
  void loopOverTracksByTrackOTF();
  void loopOverTracksByStackOTF();

protected:

  /** Pointer to the associated TrackGenerator */
  TrackGenerator* _track_generator;

  // descriptions
  MOCKernel** _kernels;
  segment** _temporary_segments;
  FP_PRECISION* _FSR_volumes;
  omp_lock_t* _FSR_locks;

  //descriptions
  segmentationType _segment_formation;
  FP_PRECISION _max_optical_length;
  int _curr_z_index;

  // descriptions
  TraverseSegments(TrackGenerator* track_generator);
  virtual ~TraverseSegments();

  void allocateTemporarySegmentStorage();
  void deallocateTemporarySegmentStorage();

  //TODO TEMPLATE
  template <class KernelType>
  void allocateKernels<KernelType>();
  template <class KernelType>
  void deallocateKernels<KernelType>();
  void loopOverTracks();

public:
  virtual void execute() = 0;
  virtual void onTrack(Track* track, segment* segments) = 0;

};

#endif
