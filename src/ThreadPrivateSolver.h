/**
 * @file ThreadPrivateSolver.h
 * @brief The ThreadPrivateSolver class.
 * @date May 28, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef THREADPRIVATESOLVER_H_
#define THREADPRIVATESOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "CPUSolver.h"
#endif

/** Indexing scheme for the thread private scalar flux for each thread in
 *  each flat source region and in each energy group */
#define _thread_flux(tid,r,e) (_thread_flux[(tid)*_num_FSRs*_num_groups+(r)*_num_groups+(e)])

/**
 * @class ThreadPrivateSolver ThreadPrivateSolver.h "openmoc/src/ThreadPrivateSolver.h"
 * @brief This is a subclass of the CPUSolver which uses thread private 
 *        arrays for the flat source region scalar fluxes to minimize OMP atomics.
 * @details Since this class stores a separate copy of the flat source region scalar
 *          fluxes for each OMP thread, the memory requirements are greater than for
 *          the CPUSolver.
 */
class ThreadPrivateSolver : public CPUSolver {

protected:

    /** An array for the flat source region scalar fluxes for each thread */
    FP_PRECISION* _thread_flux;

    void initializeFluxArrays();

    void flattenFSRFluxes(FP_PRECISION value);
    void scalarFluxTally(segment* curr_segment, 
			 FP_PRECISION* track_flux,
			 FP_PRECISION* fsr_flux);
    void reduceThreadScalarFluxes();
    void transportSweep();

public:
    ThreadPrivateSolver(Geometry* geometry=NULL, 
			TrackGenerator* track_generator=NULL);
    virtual ~ThreadPrivateSolver();
};


#endif /* THREADPRIVATESOLVER_H_ */
