/**
 * @file CPUSolver.h
 * @brief The CPUSolver class.
 * @date May 28, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef CPUSOLVER_H_
#define CPUSOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include "Solver.h"
#endif

#define INTERP_EXPONENT


#define _thread_fsr_flux(tid) (_thread_fsr_flux[tid*_num_groups])

#define track_flux(p,e) (track_flux[(p)*_num_groups + (e)])

#define track_out_flux(p,e) (track_out_flux[(p)*_num_groups + (e)])

#define track_leakage(p,e) (track_leakage[(p)*_num_groups + (e)])

#define exponentials(p,e) (exponentials[(p)*_num_groups + (e)])


/**
 * @class CPUSolver CPUSolver.h "openmoc/src/host/CPUSolver.h"
 * @brief This a subclass of the Solver class for multi-core CPUs using 
 *        OpenMP multi-threading.
 */
class CPUSolver : public Solver {

protected:

    /** The number of shared memory OpenMP threads */
    int _num_threads;

    /** OpenMP locks for atomic scalar flux updates */
    omp_lock_t* _FSR_locks;

    /** A buffer for temporary scalar flux updates for each thread */
    FP_PRECISION* _thread_fsr_flux;
    
    /** A boolean indicating whether or not to use linear interpolation 
     *  to comptue the exponential in the transport equation */
    bool _interpolate_exponent;

    FP_PRECISION* _exponentials;

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
    virtual void scalarFluxTally(segment* curr_segment, 
				 FP_PRECISION* track_flux,
				 FP_PRECISION* fsr_flux);
    virtual void transferBoundaryFlux(int track_id, bool direction,
				      FP_PRECISION* track_flux);
    void addSourceToScalarFlux();
    void computeKeff();
    void transportSweep();

    void computeExponentials(segment* curr_segment, 
			     FP_PRECISION* exponentials);

public:
    CPUSolver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL);
    virtual ~CPUSolver();
 
    int getNumThreads();
    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();

    void setNumThreads(int num_threads);
    
    void useExponentialInterpolation();
    void useExponentialIntrinsic();

    void computePinPowers();
};


#endif /* CPUSOLVER_H_ */
