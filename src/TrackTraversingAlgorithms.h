#ifndef TRACK_TRAVERSING_ALGORITHMS_H_
#define TRACK_TRAVERSING_ALGORITHMS_H_

#include "TraverseSegments.h"

//TODO: description
class MaxOpticalLength: public TraverseSegments {
private:
  FP_PRECISION _max_tau;

public:

  MaxOpticalLength(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
};


//TODO: description
class SegmentCounter: public TraverseSegments {
private:
  int _max_num_segments;

public:

  SegmentCounter(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
};


// TODO: description
class VolumeCalculator: public TraverseSegments {

public:

  VolumeCalculator(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
};


/*
class TransportSweep: public TraverseSegments {

public:

  TransportSweep(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
};
*/


#endif
