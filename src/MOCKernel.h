/**
 * @file MOCKernel.h
 * @brief An MOCKernel object
 * @date May 5, 2015
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef MOCKERNEL_H_
#define MOCKERNEL_H_

#include "Python.h"
#include "Track2D.h"
#include "Track3D.h"
#include "Geometry.h"


/* Forward declaration of TrackGenerator */
class TrackGenerator;


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

public:

  MOCKernel(TrackGenerator* track_generator, int row_num);
  virtual ~MOCKernel();

  /* Function to get the current segment count */
  int getCount();

  /* Sets the max optical path length to a different value */
  void setMaxOpticalLength(FP_PRECISION max_tau);

  /* Prepare MOCKernel for handling a new track */
  void newTrack(Track* track);

  /* Executing function describes kernel behavior */
  virtual void execute(FP_PRECISION length, Material* mat, int id,
      int cmfd_surface_fwd, int cmfd_surface_bwd)=0;

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
  CounterKernel(TrackGenerator* track_generator, int row_num);
  void execute(FP_PRECISION length, Material* mat, int id,
      int cmfd_surface_fwd, int cmfd_surface_bwd);
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

  /** A weight to be applied to segmentation data, if applicable */
  FP_PRECISION _weight;

public:

  VolumeKernel(TrackGenerator* track_generator, int row_num);
  void newTrack(Track* track);
  void execute(FP_PRECISION length, Material* mat, int id,
      int cmfd_surface_fwd, int cmfd_surface_bwd);
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
  SegmentationKernel(TrackGenerator* track_generator, int row_num);
  void execute(FP_PRECISION length, Material* mat, int id,
      int cmfd_surface_fwd, int cmfd_surface_bwd);
};

#endif /* MOCKERNEL_H_ */
