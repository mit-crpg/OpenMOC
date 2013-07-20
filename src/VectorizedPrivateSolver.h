/**
 * @file VectorizedPrivateSolver.h
 * @brief The VectorizedPrivateSolver class.
 * @date July 18, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef VECTORIZEDPRIVATESOLVER_H_
#define VECTORIZEDPRIVATESOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include "VectorizedSolver.h"
#endif


/** Indexing scheme for the thread private scalar flux for each thread in
 *  each flat source region and in each energy group */
#define _thread_flux(tid,r,e) (_thread_flux[(tid)*_num_FSRs*_num_groups+(r)*_num_groups+(e)])


/**
 * @class VectorizedPrivateSolver VectorizedPrivateSolver.h "openmoc/src/host/VectorizedPrivateSolver.h"
 * @brief This is a subclass of the VectorizedSolver class. This class 
 *        uses a thread private array for flat source region scalar 
 *        fluxes during each transport sweep to avoid the use of OpenMP 
 *        atomics. It also uses memory-aligned data structures and 
 *        Intel's auto-vectorization.
 */
class VectorizedPrivateSolver : public VectorizedSolver {

private:

    /** An array for the flat source region scalar fluxes for each thread */
    FP_PRECISION* _thread_flux;

    void initializeFluxArrays();

    void flattenFSRFluxes(FP_PRECISION value);

    void scalarFluxTally(segment* curr_segment, 
			 FP_PRECISION* track_flux,
			 FP_PRECISION* fsr_flux);

    void transportSweep();
    void reduceThreadScalarFluxes();

public:
    VectorizedPrivateSolver(Geometry* geometry=NULL, 
			    TrackGenerator* track_generator=NULL);
    virtual ~VectorizedPrivateSolver();
};


#endif /* VECTORIZEDPRIVATESOLVER_H_ */
