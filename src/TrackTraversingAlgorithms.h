/**
 * @file TrackTraversingAlgorithms.h
 * @brief Contains classes which extend the TraverseSegments class to apply
 *        algorithms to Tracks and possibly their segments.
 * @details The classes defined within this file extend the TraverseSegments
 *          class so that they are capable of using the abstract looping
 *          defined in TraverseSegments::loopOverTracks(...). Each class
 *          contains a constructor which pulls data from a provided
 *          TrackGenerator, an onTrack(...) function which specifies what to do
 *          on each Track, and an execute() function which applies the
 *          algorithm. The execute() function should contain a call to
 *          TraverseSegments::loopOverTracks(...). To specify a behavior to
 *          be applied once for each segment, the kernel should be initialized
 *          using the TraverseSegments::getKernel<KernelType>() function and
 *          passed to TraverseSegments::loopOverTracks().
 * @date February 23, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef TRACK_TRAVERSING_ALGORITHMS_H_
#define TRACK_TRAVERSING_ALGORITHMS_H_

#include "TraverseSegments.h"
#include "ExpEvaluator.h"

/** Optimization macro to help vectorize computing the LS source contants */
#ifdef THREED
#define _NUM_COEFFS 6
#endif

/** Forward declaration of CPUSolver class */
class CPUSolver;
class CPULSSolver;


/**
 * @class MaxOpticalLength TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to calculate the maximum optical path length across
 *        all segments in the Geometry.
 * @details A MaxOpticalLength allocates SegmentationKernels to temporarily
 *          store segment data. The segments are then traversed afterwards and
 *          the maximum optical path length is calculated.
 */
class MaxOpticalLength: public TraverseSegments {
private:
  FP_PRECISION _min_tau;
  FP_PRECISION _max_tau;

public:

  MaxOpticalLength(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
};


/**
 * @class SegmentCounter TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to count the number of segments on each Track and the
 *        maximum number of segments per Track.
 * @details A SegmentCounter allocates CounterKernels to count segments along
 *          each Track given a maximum optical path length imported from the
 *          provided TrackGenerator and the maximum number of segments per
 *          Track in the TrackGenerator is updated.
 */
class SegmentCounter: public TraverseSegments {
private:
  int _max_num_segments;
  long _total_num_segments;
  bool _count_total_segments;
  bool _total_segments_counted;

public:

  SegmentCounter(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
  void countTotalNumSegments();
  long getTotalNumSegments();
};


/**
 * @class SegmentSplitter TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to split explicit segments along Tracks.
 * @details A SegmentSplitter imports a maximum optical path length from the
 *          provided TrackGenerator and then ensures all segments have an
 *          optical path length less than the maximum optical path length by
 *          splitting segments stored in the Tracks.
 */
class SegmentSplitter: public TraverseSegments {

public:

  SegmentSplitter(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
};


/**
 * @class VolumeCalculator TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to calculate FSR volumes.
 * @details A VolumeCalculator imports a buffer to store FSR volumes from the
 *          provided TrackGenerator and the allocates VolumeKernels to
 *          calculate and update the volumes in each FSR, implicitly writing
 *          the calculated volumes back to the TrackGenerator.
 */
class VolumeCalculator: public TraverseSegments {

public:

  VolumeCalculator(TrackGenerator* track_generator);
  void execute();
  void onTrack(Track* track, segment* segments);
};


/**
 * @class CentroidGenerator TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to calculate the centroids of each FSR.
 * @details A CentroidGenerator imports FSR Volumes and associated locks form
 *          the provided TrackGenerator, then centroids are calculated and
 *          stored in the provided buffer by first allocating
 *          SegmentationKernels to temporarily store segments and then looping
 *          over all segments and adding their contribution to each FSR
 *          centroid.
 */
class CentroidGenerator: public TraverseSegments {

private:

  Point** _centroids;
  FP_PRECISION* _FSR_volumes;
  omp_lock_t* _FSR_locks;
  Quadrature* _quadrature;
  Point** _starting_points;

public:

  CentroidGenerator(TrackGenerator* track_generator);
  virtual ~CentroidGenerator();
  void setCentroids(Point** centroids);
  void execute();
  void onTrack(Track* track, segment* segments);
};


/**
 * @class LinearExpansionGenerator TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to calculate the linear expansion coeffs of each FSR.
 * @details A LinearExpansionGenerator loops through all tracks and computes :
 *          a linear source term that is constant through all iterations, matrix
 *          coefficients in each FSR that help compute the linear source.
 */  //TODO Improve
class LinearExpansionGenerator: public TraverseSegments {

private:

  double* _lin_exp_coeffs;
  FP_PRECISION* _FSR_volumes;
  omp_lock_t* _FSR_locks;
  double* _src_constants;
  Quadrature* _quadrature;
  CPULSSolver* _solver;
  int _num_groups;
#ifndef THREED
  int _NUM_COEFFS;
#endif
  int _num_flat;
  Point** _starting_points;
  double** _thread_source_constants;
  ExpEvaluator* _exp_evaluator;
  Progress* _progress;

public:

  LinearExpansionGenerator(CPULSSolver* solver);
  virtual ~LinearExpansionGenerator();
  void execute();
  void onTrack(Track* track, segment* segments);
};


/**
 * @class TransportSweep TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to apply the MOC transport equations to all segments.
 * @details TransportSweep imports data from the provided TrackGenerator and
 *          using a provided CPUSolver, it applies the MOC equations to each
 *          segment, tallying the contributions to each FSR. At the end of each
 *          Track, boundary fluxes are exchanged based on boundary conditions.
 */
class TransportSweep: public TraverseSegments {

private:

  CPUSolver* _cpu_solver;
  CPULSSolver* _ls_solver;
  Geometry* _geometry;
  int _num_groups;
  FP_PRECISION** _thread_fsr_fluxes;
  FP_PRECISION** _thread_scratch_pads;

public:

  TransportSweep(CPUSolver* cpu_solver);
  virtual ~TransportSweep();
  void execute();
  void onTrack(Track* track, segment* segments);

};


/**
 * @class DumpSegments TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to write tracking data to a file.
 * @details DumpSegments imports Track data from the provided TrackGenerator
 *          and writes the tracking data to the provided file.
 */
class DumpSegments: public TraverseSegments {

private:

  FILE* _out;

public:

  DumpSegments(TrackGenerator* track_generator);
  void setOutputFile(FILE* out);
  void execute();
  void onTrack(Track* track, segment* segments);
};


/**
 * @class ReadSegments TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to read tracking data from a file.
 * @details ReadSegments imports Track data from the provided file and writes
 *          the tracking data to the Tracks in the provided TrackGenerator.
 */
class ReadSegments: public TraverseSegments {

private:

  FILE* _in;

public:

  ReadSegments(TrackGenerator* track_generator);
  void setInputFile(FILE* input);
  void execute();
  void onTrack(Track* track, segment* segments);
};


/**
 * @class TransportSweepOTF TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to ray trace and apply the MOC transport equations 
 *        simultaneously On-The-Fly.
 * @details TransportSweepOTF uses a track generator and a CPUSolver to solve
 *          the ray-trace and solve the MOC equations at the same time. 
 *          NOT FULLY IMPLEMENTED //TODO
 */
class TransportSweepOTF: public TraverseSegments {

private:

  CPUSolver* _cpu_solver;

public:

  TransportSweepOTF(TrackGenerator* track_generator);
  void setCPUSolver(CPUSolver* cpu_solver);
  void onTrack(Track* track, segment* segments);
  void execute();
};


/**
 * @class RecenterSegments TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to recenter Segments once their related FSR centroid is 
 *        known.
 * @details RecenterSegments imports Track data from the provided TrackGenerator
 *          and loops through segments to re-center them once their FSR centroid
 *          has been computed. This is important to use the linear source with 
 *          explicit ray tracing.
 */
class RecenterSegments: public TraverseSegments {

private:

  Geometry* _geometry;

public:

  RecenterSegments(TrackGenerator* track_generator);
  void onTrack(Track* track, segment* segments);
  void execute();
};


/**
 * @class PrintSegments TrackTraversingAlgorithms.h
 *        "src/TrackTraversingAlgorithms.h"
 * @brief A class used to write tracking data to a file.
 * @details PrintSegments imports Track data from the provided TrackGenerator
 *          and writes the tracking data to the provided file.
 */
class PrintSegments: public TraverseSegments {

private:

  FILE* _out;

public:

  PrintSegments(TrackGenerator* track_generator);
  void setOutputFile(FILE* out);
  void execute();
  void onTrack(Track* track, segment* segments);
};

#endif
