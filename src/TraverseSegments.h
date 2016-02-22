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
  void loopOverTracks2D();
  void loopOverTracksExplicit();
  void loopOverTracksByTrackOTF();
  void loopOverTracksByStackOTF();

protected:

  /** Pointer to the associated TrackGenerator */
  TrackGenerator* _track_generator;

  // descriptions
  MOCKernel** _kernels;

  //descriptions
  segmentationType _segment_formation;

  // descriptions
  TraverseSegments(TrackGenerator* track_generator);
  virtual ~TraverseSegments();

  void loopOverTracks();
  virtual void onTrack(Track* track, segment* segments) = 0;

  // description
  template <class KernelType>
  void setKernel() {
    int num_rows = _track_generator->getNumRows();
    _kernels = new MOCKernel*[num_rows];
    for (int z=0; z < num_rows; z++)
      _kernels[z] = new KernelType(_track_generator, z);
  }

public:
  virtual void execute() = 0;
};

#endif
