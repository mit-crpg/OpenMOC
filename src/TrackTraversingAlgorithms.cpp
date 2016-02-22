#include "TrackTraversingAlgorithms.h"


//TODO: description
MaxOpticalLength::MaxOpticalLength(TrackGenerator* track_generator)
                                 : TraverseSegments(track_generator) {
  _max_tau = 0;
}


//TODO: description
void MaxOpticalLength::execute() {
  FP_PRECISION infinity = std::numeric_limits<FP_PRECISION>::max();
  _track_generator->setMaxOpticalLength(infinity);
#pragma omp parallel
  {
    setKernel<SegmentationKernel>();
    loopOverTracks();
  }
  _track_generator->setMaxOpticalLength(_max_tau);
}


//TODO: description
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


/*
  TODO: class description
*/

//TODO: description
SegmentCounter::SegmentCounter(TrackGenerator* track_generator)
                               : TraverseSegments(track_generator) {
  _max_num_segments = 0;
}


//TODO: description
void SegmentCounter::execute() {
#pragma omp parallel
  {
    setKernel<CounterKernel>();
    loopOverTracks();
  }
  _track_generator->setMaxNumSegments(_max_num_segments);
}


//TODO: description
void SegmentCounter::onTrack(Track* track, segment* segments) {
  if (track->getNumSegments() > _max_num_segments) {
#pragma omp critical
    _max_num_segments = std::max(_max_num_segments, track->getNumSegments());
  }
}


/*
   TODO: class description
*/

//TODO: description
VolumeCalculator::VolumeCalculator(TrackGenerator* track_generator)
                                  : TraverseSegments(track_generator) {
}


//TODO: description
void VolumeCalculator::execute() {
#pragma omp parallel
  {
    setKernel<VolumeKernel>();
    loopOverTracks();
  }
}


//TODO: description
void VolumeCalculator::onTrack(Track* track, segment* segments) {
}


