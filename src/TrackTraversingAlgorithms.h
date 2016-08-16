/**
 * @file TrackTraversingAlgorithms.h
 * @brief Contains classes which extend the TraverseTracks class to apply
 *        algorithms to Tracks and possibly their segments
 * @details The classes defined within this file extend the TraverseTracks
 *          class so that they are capable of using the abstract looping
 *          defined in TraverseTracks::loopOverTracks(...). Each class
 *          contains a constructor which pulls data from a provided
 *          TrackGenerator, an onTrack(...) function which specifies what to do
 *          on each Track, and an execute() function which applies the
 *          algorithm. The execute() function should contain a call to
 *          TraverseTracks::loopOverTracks(...). To specify a behavior to
 *          be applied once for each segment, a kernel should be passed to
 *          TraverseTracks::loopOverTracks().
 * @date February 23, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef TRACK_TRAVERSING_ALGORITHMS_H_
#define TRACK_TRAVERSING_ALGORITHMS_H_

#include "TraverseTracks.h"


/** Forward declaration of CPUSolver class */
class CPUSolver;


/**
 * @class VolumeCalculator TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to calculate FSR volumes
 * @details A VolumeCalculator imports a buffer to store FSR volumes from the
 *          provided TrackGenerator and the allocates VolumeKernels to
 *          calculate and update the volumes in each FSR, implicitly writing
 *          the calculated volumes back to the TrackGenerator.
 */
class VolumeCalculator: public TraverseTracks {

public:

  VolumeCalculator(TrackGenerator* track_generator);
  void execute();
};


/**
 * @class TransportSweep TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to apply the MOC transport equations to all segments
 * @details TransportSweep imports data from the provided TrackGenerator and
 *          using a provided CPUSolver, it applies the MOC equations to each
 *          segment, tallying the contributions to each FSR. At the end of each
 *          Track, boundary fluxes are exchanged based on boundary conditions.
 */
class TransportSweep: public TraverseTracks {

private:

  CPUSolver* _cpu_solver;
  FP_PRECISION** _thread_fsr_fluxes;

public:

  TransportSweep(TrackGenerator* track_generator);
  virtual ~TransportSweep();
  void setCPUSolver(CPUSolver* cpu_solver);
  void execute();
  void onTrack(Track* track, segment* segments);
};


#endif
