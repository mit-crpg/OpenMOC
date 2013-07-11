/**
 * @file CPUSolver.h
 * @brief The CPUSolver class.
 * @date February 7, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef CPUSOLVER_H_
#define CPUSOLVER_H_


#define VEC_LENGTH 8
#define VEC_ALIGNMENT 16

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <mm_malloc.h>
#include "Solver.h"
#endif

#define track_flux(p,e) (track_flux[(e)*_num_polar + (p)])
#define track_out_flux(p,e) (track_out_flux[(e)*_num_polar + (p)])
#define track_leakage(p,e) (track_leakage[(e)*_num_polar + (p)])


/**
 * @class CPUSolver CPUSolver.h "openmoc/src/host/CPUSolver.h"
 * @brief
 */
class CPUSolver : public Solver {

private:

    /** OpenMP locks for atomic scalar flux updates */
    omp_lock_t* _FSR_locks;

    /** The number of shared memory OpenMP threads */
    int _num_threads;

    /** Number of energy groups divided by vector widths (VEC_LENGTH) */
    int _num_groups_vec;

    void initializeFluxArrays();
    void initializeSourceArrays();
    void initializePowerArrays();
    void initializePolarQuadrature();
    void precomputePrefactors();
    void initializeFSRs();

    void zeroTrackFluxes();
    void zeroTrackLeakages();
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
 
    int getNumThreads();
    int getNumGroupVectorWidths();
    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();

    void scalarFluxTally(int fsr_id, double* sigma_t, 
			 FP_PRECISION length, 
			 FP_PRECISION* track_flux,
			 FP_PRECISION* fsr_flux);
    void transferBoundaryFlux(int j, int k, int track_id, 
			      bool direction,
			      FP_PRECISION* track_flux);
    void normalizeFluxToVolume();

    void setNumThreads(int num_threads);
    void setGeometry(Geometry* geometry);
    void computePinPowers();
};


#endif /* CPUSOLVER_H_ */
