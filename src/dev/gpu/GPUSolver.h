/**
 * @file GPUSolver.h
 * @brief The GPUSolver class and GPU solver routines.
 * @date August 5, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GPUSOLVER_H_
#define GPUSOLVER_H_

#ifdef __cplusplus
#include "../../Solver.h"
#endif

#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <sm_13_double_functions.h>
#include <sm_20_atomic_functions.h>
#include "clone.h"

#define scalar_flux(tid,e) (scalar_flux[(tid)*(*num_groups) + (e)])

#define source(tid,e) (source[(tid)*(*num_groups) + (e)])

#define old_source(tid,e) (old_source[(tid)*(*num_groups) + (e)])

#define reduced_source(tid,e) (reduced_source[(tid)*(*num_groups) + (e)])

#define polar_weights(i,p) (polar_weights[(i)*(*num_polar) + (p)])

#define boundary_flux(tid,pe2) (boundary_flux[2*(tid)*(*polar_times_groups)+(pe2)])

#define prefactorindex(tau) (int(tau * _inverse_prefactor_spacing) * _two_times_num_polar)

#define prefactor(index,p,tau) (1. - (_prefactor_array[index+2 * p] * tau + _prefactor_array[index + 2 * p +1]))


/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.5663706143

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.0795774715

/** The maximum number of polar angles to reserve constant memory on GPU */
#define MAX_POLAR_ANGLES 3

/** The maximum number of azimuthal angles to reserve constant memory on GPU */
#define MAX_AZIM_ANGLES 256


/**
 * @class GPUSolver GPUSolver.h "openmoc/src/dev/GPUSolver.h"
 * @brief
 */
class GPUSolver : public Solver {

private:

    /**************************************************************************/
    /*                             Data on the host                           */
    /**************************************************************************/

    /** The number of threadblocks */
    int _B;
    
    /** The number of threads per threadblock */
    int _T;


    /**************************************************************************/
    /*                           Data on the device                           */
    /**************************************************************************/

    /** The flat source region material pointers index by FSR UID */
    int* _FSR_materials;

    /** A pointer to an array of the materials on the device */
    dev_material* _materials;

    /** A pointer to the array of tracks on the device */
    dev_track* _dev_tracks;

    /** An array of the cumulative number of tracks for each azimuthal angle */
    int* _track_index_offsets;

    /** A pointer to the Thrust vector of fission sources in each FSR */
    FP_PRECISION* _fission_source;

    /** A pointer to the Thrust vector of absorption rates in each FSR */
    FP_PRECISION* _tot_absorption;

    /** A pointer to the Thrust vector of fission rates in each FSR */
    FP_PRECISION* _tot_fission;

    /** A pointer to the Thrust vector of leakages for each track */
    FP_PRECISION* _leakage;

    /** Thrust vector of fission sources in each FSR */
    thrust::device_vector<FP_PRECISION> _fission_source_vec;

    /** Thrust vector of fission rates in each FSR */
    thrust::device_vector<FP_PRECISION> _tot_fission_vec;

    /** Thrust vector of absorption rates in each FSR */
    thrust::device_vector<FP_PRECISION> _tot_absorption_vec;

    /** Thrust vector of source residuals in each FSR */
    thrust::device_vector<FP_PRECISION> _source_residuals_vec;

    /** Thrust vector of leakages for each track */
    thrust::device_vector<FP_PRECISION> _leakage_vec;

    void initializePolarQuadrature();
    void initializeFSRs();
    void initializeMaterials();
    void initializeTracks();
    void initializeFluxArrays();
    void initializeSourceArrays();
    void initializePowerArrays();
    void initializeThrustVectors();
    void precomputePrefactors();

    void zeroTrackFluxes();
    void flattenFSRFluxes(FP_PRECISION value);
    void flattenFSRSources(FP_PRECISION value);
    void normalizeFluxes();
    FP_PRECISION computeFSRSources();
    void addSourceToScalarFlux();
    void computeKeff();
    void transportSweep();

public:

    GPUSolver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL);
    virtual ~GPUSolver();
    
    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();

    void setNumThreadBlocks(int num_blocks);
    void setNumThreadsPerBlock(int num_threads);
    void setGeometry(Geometry* geometry);
    void setTrackGenerator(TrackGenerator* track_generator);

    int computeScalarTrackIndex(int i, int j);
    void computePinPowers();
};


#endif /* GPUSOLVER_H_ */
