#include "TrackTraversingAlgorithms.h"
#include "CPUSolver.h"


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
    MOCKernel** kernels = getKernels<SegmentationKernel>();
    loopOverTracks(kernels);
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
    MOCKernel** kernels = getKernels<CounterKernel>();
    loopOverTracks(kernels);
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
    MOCKernel** kernels = getKernels<VolumeKernel>();
    loopOverTracks(kernels);
  }
}


//TODO: description
void VolumeCalculator::onTrack(Track* track, segment* segments) {
}

/*
   TODO: class description
*/

//TODO: description
TransportSweep::TransportSweep(TrackGenerator* track_generator)
                              : TraverseSegments(track_generator) {
  _cpu_solver = NULL;
}


//TODO: description
void TransportSweep::execute() {
#pragma omp parallel
  {
    MOCKernel** kernels = getKernels<SegmentationKernel>();
    loopOverTracks(kernels);
  }
}


//TODO: description
void TransportSweep::setCPUSolver(CPUSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


//TODO: description
void TransportSweep::onTrack(Track* track, segment* segments) {

  /* Allocate temporary FSR flux locally */
  int num_groups = _track_generator->getGeometry()->getNumEnergyGroups();
  FP_PRECISION thread_fsr_flux[num_groups];

  /* Extract Track information */
  int track_id = track->getUid();
  int azim_index = track->getAzimIndex();
  int num_segments = track->getNumSegments();
  FP_PRECISION* track_flux;

  /* Extract the polar index if a 3D track */
  int polar_index = 0;
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL)
    polar_index = track_3D->getPolarIndex();

  /* Correct azimuthal index to first octant */
  Quadrature* quad = _track_generator->getQuadrature();
  azim_index = quad->getFirstOctantAzim(azim_index);

  /* Get the forward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, true);

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {
    segment* curr_segment = &segments[s];
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, polar_index,
                                 track_flux, thread_fsr_flux);
    _cpu_solver->tallySurfaceCurrent(curr_segment, azim_index, polar_index,
                                     track_flux, true);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, polar_index, true,
                                    track_flux);

  /* Get the backward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, false);

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {
    segment* curr_segment = &segments[s];
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, polar_index,
                                 track_flux, thread_fsr_flux);
    _cpu_solver->tallySurfaceCurrent(curr_segment, azim_index, polar_index,
                                     track_flux, false);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, polar_index, false,
                                    track_flux);
}


