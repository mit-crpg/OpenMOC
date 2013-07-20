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

/** Indexing macro for the thread private flat source region scalar fluxes */
#define _thread_fsr_flux(tid) (_thread_fsr_flux[tid*_num_groups])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for either the forward or reverse direction for a given track */ 
#define track_flux(p,e) (track_flux[(p)*_num_groups + (e)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for the outgoing reflective track from a given track */
#define track_out_flux(p,e) (track_out_flux[(p)*_num_groups + (e)])

/** Indexing macro for the leakage for each polar angle and energy group
 *  for either the forward or reverse direction for a given track */
#define track_leakage(p,e) (track_leakage[(p)*_num_groups + (e)])


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

    /**
     * @brief Computes the contribution to the flat source region scalar flux
     *        from a single track segment.
     * @param curr_segment a pointer to the segment of interest
     * @param track_flux a pointer to the track's angular flux
     * @param fsr_flux a pointer to the temporary flat source region flux buffer
     */
    virtual void scalarFluxTally(segment* curr_segment, 
				 FP_PRECISION* track_flux,
				 FP_PRECISION* fsr_flux);

    /**
     * @brief Updates the boundary flux for a track given boundary conditions.
     * @param track_id the ID number for the track of interest
     * @param direction the track direction (forward - true, reverse - false)
     * @param track_flux a pointer to the track's outgoing angular flux
     */
    virtual void transferBoundaryFlux(int track_id, bool direction,
				      FP_PRECISION* track_flux);

    void addSourceToScalarFlux();
    void computeKeff();
    void transportSweep();

    /**
     * @brief Computes the exponential term in the transport equation for a
     *        track segment.
     * @param sigma_t the total group cross-section at this energy
     * @param length the length of the line segment projected in the xy-plane
     * @param p the polar angle index
     * @return the evaluated exponential
     */
    virtual FP_PRECISION computeExponential(FP_PRECISION sigma_t, 
					    FP_PRECISION length, int p); 
public:
    CPUSolver(Geometry* geometry=NULL, TrackGenerator* track_generator=NULL);
    virtual ~CPUSolver();
 
    int getNumThreads();
    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();

    void setNumThreads(int num_threads);
    
    void computePinPowers();
};


#endif /* CPUSOLVER_H_ */
