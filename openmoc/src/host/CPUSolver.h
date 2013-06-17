/**
 * @file CPUSolver.h
 * @brief The CPUSolver class.
 * @date February 7, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef CPUSOLVER_H_
#define CPUSOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include "Solver.h"
#endif

#define thread_flux(t,r,e) (thread_flux[(t)*_num_FSRs*_num_groups + (r)*_num_groups + (e)])


/**
 * @class CPUSolver CPUSolver.h "openmoc/src/host/CPUSolver.h"
 * @brief
 */
class CPUSolver : public Solver {

private:
    /** The number of shared memory OpenMP threads */
    int _num_threads;

    void initializeFluxArrays();
    void initializeSourceArrays();
    void initializePowerArrays();
    void initializePolarQuadrature();
    void precomputePrefactors();
    void initializeFSRs();

    void zeroTrackFluxes();
    void flattenFSRFluxes(FP_PRECISION value);
    void flattenFSRSources(FP_PRECISION value);
    void normalizeFluxes();
    FP_PRECISION computeFSRSources();
    void computeKeff();
    bool isScalarFluxConverged();
    void transportSweep(int max_iterations);

public:
    CPUSolver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL);
    virtual ~CPUSolver();

    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();

    int getNumThreads();
    void setNumThreads(int num_threads);
    void computePinPowers();
};


#endif /* CPUSOLVER_H_ */
