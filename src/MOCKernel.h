/**
 * @file MOCKernel.h
 * @brief An MOCKernel object
 * @date May 5, 2015
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef MOCKERNEL_H_
#define MOCKERNEL_H_

#include "Track.h"
#include "Geometry.h"
#include "Quadrature.h"
#ifdef SWIG
#include "Python.h"
#endif


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
 *          code. This class is the parent class of CounterKernel and
 *          VolumeKernel. A generic MOCKernel should not be explicity
 *          instantiated. Instead, an inheriting class should
 *          be instantiated which implements the "execute" function.
 */
class MOCKernel {

protected:

  /** Count referring to the number of segments the MOCKernel has handled since
   *  the creation or the last reset */
  int _count;

  /** Maximum optical path length when forming segments */
  FP_PRECISION _max_tau;

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
  virtual void execute(FP_PRECISION length, Material* mat, int id,
                       int cmfd_surface_fwd, int cmfd_surface_bwd)=0;

};


/**
 * @class CounterKernel MOCKernel.h "src/MOCKernel.h"
 * @brief Counts the number of segments of a track
 * @details A CounterKernel inherets from MOCKernel and is a kernel which
 *          tallies the number of legitimate segment lengths (less than the max
 *          optical path length) encountered. This is useful for determining
 *          the number of segments in a Track if the number is not yet known.
 */
class CounterKernel: public MOCKernel {

public:

  CounterKernel(TrackGenerator* track_generator);
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

  /** The cross-sectional area of the Track used to weight segment length
   *  contributions to the volume */
  FP_PRECISION _weight;

  /** The associated quadrature from which weights are derived */
  Quadrature* _quadrature;

public:

  VolumeKernel(TrackGenerator* track_generator);
  void newTrack(Track* track);
  void execute(FP_PRECISION length, Material* mat, int id,
               int cmfd_surface_fwd, int cmfd_surface_bwd);
};


#endif /* MOCKERNEL_H_ */
