/**
 * @file MOCKernel.h
 * @brief An MOCKernel object
 * @date May 5, 2015
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef MOCKERNEL_H_
#define MOCKERNEL_H_

#ifdef SWIG
#include "Python.h"
#endif
#include "Track.h"
#include "Track3D.h"
#include "Geometry.h"
#include "Quadrature.h"


/* Forward declaration of TrackGenerator */
class TrackGenerator;

/* Forward declaration of CPUSolver */
class CPUSolver;

/**
 * @class MOCKernel MOCKernel.h "src/MOCKernel.h"
 * @brief An MOCKernel object specifies a functionality to apply to MOC
 *        segments.
 * @details An MOCKernel is an object that owns some specified data and
 *          contains an "execute" function which applies some functionality
 *          to the data. This is useful in MOC where it is very common to
 *          apply some function to segment data in either a nested loop
 *          structure or from on-the-fly calculations. Kernels specify the
 *          actions applied to the segments, reducing the need for repeated
 *          code. This class is the parent class of CounterKernel,
 *          VolumeKernel, and SegmentationKernel. A generic MOCKernel should
 *          not be explicity instantiated. Instead, an inhereting class should
 *          be instantiated which describes the "execute" function.
 */
class MOCKernel {

protected:

  /** Count referring to the segment number */
  int _count;

  /** Maximum optical path length when forming segments */
  FP_PRECISION _max_tau;

  /** Number of energy groups in the current problem */
  int _num_groups;

public:

  MOCKernel(TrackGenerator* track_generator);
  virtual ~MOCKernel();

  /* Function to get the current segment count */
  int getCount();

  /* Sets the max optical path length to a different value */
  void setMaxOpticalLength(FP_PRECISION max_tau);

  /* Prepare MOCKernel for handling a new track */
  virtual void newTrack(Track* track);

  /* Executing function describes kernel behavior */
  virtual void execute(FP_PRECISION length, Material* mat, long fsr_id,
                       int track_idx, int cmfd_surface_fwd,
                       int cmfd_surface_bwd, FP_PRECISION x_start, FP_PRECISION y_start,
                       FP_PRECISION z_start, FP_PRECISION phi, FP_PRECISION theta)=0;

};


/**
 * @class CounterKernel MOCKernel.h "src/MOCKernel.h"
 * @brief Counts the number of segments of a track
 * @details A CounterKernel inherets from MOCKernel and is a kernel which
 *          counts the number of segments in a track by incrementing the
 *          _count variable by the number of legitimate segment lengths
 *          (less than the max optical path length) in the input length.
 */
class CounterKernel: public MOCKernel {

public:
  CounterKernel(TrackGenerator* track_generator);
  void execute(FP_PRECISION length, Material* mat, long fsr_id,
               int track_idx, int cmfd_surface_fwd, int cmfd_surface_bwd,
               FP_PRECISION x_start, FP_PRECISION y_start, FP_PRECISION z_start,
               FP_PRECISION phi, FP_PRECISION theta);
};


/**
 * @class VolumeKernel MOCKernel.h "src/MOCKernel.h"
 * @brief Calculates the volume in FSRs by adding weighted segment lengths
 * @details A VolumeKernel inherets from MOCKernel and is a kernel which
 *          is initialized with a pointer to floating point data and adds
 *          the product of the length and the weight to the floating point data
 *          at an input index. The weight corresponds to the weight of the
 *          track associated with the segments.
 */
class VolumeKernel: public MOCKernel {

private:

  /** Array of FSR locks */
  omp_lock_t* _FSR_locks;

  /** Pointer to array of FSR volumes */
  FP_PRECISION* _FSR_volumes;

  /** The Track's volume weight */
  FP_PRECISION _weight;

  /** The associated quadrature from which weights are derived */
  Quadrature* _quadrature;

public:

  VolumeKernel(TrackGenerator* track_generator);
  void newTrack(Track* track);
  void execute(FP_PRECISION length, Material* mat, long fsr_id,
               int track_idx, int cmfd_surface_fwd, int cmfd_surface_bwd,
               FP_PRECISION x_start, FP_PRECISION y_start, FP_PRECISION z_start,
               FP_PRECISION phi, FP_PRECISION theta);
};


/**
 * @class SegmentationKernel MOCKernel.h "src/MOCKernel.h"
 * @brief Forms segment data associated with a 3D track
 * @details A SegmentationKernel inherets from MOCKernel and is a kernel which
 *          is initialized with a pointer to segment data. Input data of the
 *          "execute" function is saved to the segment data, forming explicit
 *          segments.
 */
class SegmentationKernel: public MOCKernel {

private:

  /** Pointer to segment data */
  segment* _segments;

public:
  SegmentationKernel(TrackGenerator* track_generator);
  void execute(FP_PRECISION length, Material* mat, long fsr_id,
               int track_idx, int cmfd_surface_fwd, int cmfd_surface_bwd,
               FP_PRECISION x_start, FP_PRECISION y_start, FP_PRECISION z_start,
               FP_PRECISION phi, FP_PRECISION theta);
};

/**
 * @class TransportKernel MOCKernel.h "src/MOCKernel.h"
 * @brief Applies transport equations to segment data
 * @details A TransportKernel inherets from MOCKernel and is a kernel which
 *          is initialized with a pointer to a CPU Solver. Input data of the
 *          "execute" function is used to apply the MOC equations in CPUSolver.
 */
class TransportKernel: public MOCKernel {
private:

  /** Pointer to CPUSolver enabling use of transport functions */
  CPUSolver* _cpu_solver;

  /** Pointer to TrackGenerator enabling use of ray tracer on the fly */
  TrackGenerator* _track_generator;

  /** Pointer to angular flux data in the current direction */
  FP_PRECISION* _thread_fsr_flux;

  /** Azimuthal index of the current track */
  int _azim_index;

  /** XY index of the current track within tracks of that azimuthal angle */
  int _xy_index;

  /** Polar index of the current track */
  int _polar_index;

  /** Unique ID of the current track */
  int _track_id;

  /** Direction of the current track (true = Forward / false = Backward) */
  bool _direction;

  int _min_track_idx;
  int _max_track_idx;

public:
  TransportKernel(TrackGenerator* track_generator);
  virtual ~TransportKernel();
  void newTrack(Track* track);
  void setCPUSolver(CPUSolver* cpu_solver);
  void setTrackFlux(FP_PRECISION* fwd_flux, FP_PRECISION* bwd_flux,
                    int track_id);
  void setTrackIndexes(int azim_index, int polar_index);
  void setDirection(bool direction);
  bool getDirection();
  void execute(FP_PRECISION length, Material* mat, long fsr_id,
               int track_idx, int cmfd_surface_fwd, int cmfd_surface_bwd,
               FP_PRECISION x_start, FP_PRECISION y_start, FP_PRECISION z_start,
               FP_PRECISION phi, FP_PRECISION theta);
  void post();
};


#endif /* MOCKERNEL_H_ */
