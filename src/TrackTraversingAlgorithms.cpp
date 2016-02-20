#include "TrackTraversingAlgorithms.h"


//TODO: description
MaxOpticalLength::MaxOpticalLength(TrackGenerator* track_generator)
                                 : TraverseSegments(track_generator) {
  _max_tau = 0;
};
  

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

