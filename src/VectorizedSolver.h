/**
 * @file VectorizedSolver.h
 * @brief The VectorizedSolver class.
 * @date May 28, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef VECTORIZEDSOLVER_H_
#define VECTORIZEDSOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <mkl.h>
#include "CPUSolver.h"
#endif

/** Indexing scheme for the optical length (\f$ l\Sigma_t \f$) for a 
 *  given segment for each polar angle and energy group */
#define taus(p,e) (taus[(p)*_num_groups + (e)])

/** Indexing scheme for the exponentials in the neutron transport equation
 *  (\f$ 1 - exp(-\frac{l\Sigma_t}{sin(\theta_p)}) \f$) for a given
 *  segment for each polar angle and energy group */
#define exponentials(p,e) (exponentials[(p)*_num_groups + (e)])

/**
 * @class VectorizedSolver VectorizedSolver.h "openmoc/src/host/VectorizedSolver.h"
 * @brief This is a subclass of the CPUSolver class which uses memory-aligned
 *        data structures and Intel's auto-vectorization.
 */
class VectorizedSolver : public CPUSolver {

protected:

    /** Number of energy groups divided by vector widths (VEC_LENGTH) */
    int _num_vector_lengths;
    
    /** An array for the optical length for each thread in each energy group */
    FP_PRECISION* _thread_taus;

    /** An array for the exponential terms in the transport equation for *
     *  each thread in each energy group and polar angle */
    FP_PRECISION* _thread_exponentials;

    void precomputePrefactors();
    void initializeFluxArrays();
    void initializeSourceArrays();

    void normalizeFluxes();
    FP_PRECISION computeFSRSources();
    void scalarFluxTally(segment* curr_segment, 
			 FP_PRECISION* track_flux,
			 FP_PRECISION* fsr_flux);
    void transferBoundaryFlux(int track_id, bool direction,
			      FP_PRECISION* track_flux);
    void addSourceToScalarFlux();
    void computeKeff();
    

    /**
     * @brief Computes an array of the exponentials in the transport equation,
     *        \f$ exp(-\frac{\Sigma_t * l}{sin(\theta)}) \f$, for each 
     *        energy group and polar angle for a given segment.
     * @param curr_segment pointer to the segment of interest
     * @param exponentials the array to store the exponential values
     */
    virtual void computeExponentials(segment* curr_segment, 
				     FP_PRECISION* exponentials);

public:
    VectorizedSolver(Geometry* geometry=NULL, 
		     TrackGenerator* track_generator=NULL);
    virtual ~VectorizedSolver();
 
    int getNumVectorWidths();

    void setGeometry(Geometry* geometry);
};


#endif /* VECTORIZEDSOLVER_H_ */
