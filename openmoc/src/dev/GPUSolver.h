/**
 * @file GPUSolver.h
 * @brief The GPUSolver class and GPU solver routines.
 * @date August 5, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GPUSOLVER_H_
#define GPUSOLVER_H_

#ifdef __cplusplus
#include "../host/Solver.h"
#endif

#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <sm_13_double_functions.h>
#include <sm_20_atomic_functions.h>
#include "DeviceFlatSourceRegion.h"
#include "DeviceTrack.h"
#include "DeviceMaterial.h"


#define scalar_flux(tid,e) (scalar_flux[(tid)*(*num_groups) + (e)])
#define old_scalar_flux(tid,e) (old_scalar_flux[(tid)*(*num_groups) + (e)])
#define source(tid,e) (source[(tid)*(*num_groups) + (e)])
#define old_source(tid,e) (old_source[(tid)*(*num_groups) + (e)])
#define ratios(tid,e) (ratios[(tid)*(*num_groups) + (e)])
#define polar_weights(i,p) (polar_weights[(i)*(*num_polar) + (p)])
#define boundary_flux(tid,pe2) (boundary_flux[2*(tid)*(*polar_times_groups)+(pe2)])

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

    /** A pointer to an array of the flat source regions on the device */
    dev_flatsourceregion* _FSRs;

    /** A pointer to an array of the materials on the device */
    dev_material* _materials;

    /** A pointer to the array of tracks on the device */
    dev_track* _dev_tracks;

    /** An array of the cumulative number of tracks for each azimuthal angle */
    int* _track_index_offsets;

    FP_PRECISION* _fission_source;
    FP_PRECISION* _tot_absorption;
    FP_PRECISION* _tot_fission;
    FP_PRECISION* _source_residual;
    FP_PRECISION* _leakage;
    thrust::device_vector<FP_PRECISION> _fission_source_vec;
    thrust::device_vector<FP_PRECISION> _tot_absorption_vec;
    thrust::device_vector<FP_PRECISION> _tot_fission_vec;
    thrust::device_vector<FP_PRECISION> _source_residual_vec;
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
    void computeKeff();
    bool isScalarFluxConverged();
    void transportSweep(int max_iterations);

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
    void setFluxConvergenceThreshold(FP_PRECISION flux_thresh);

    int computeScalarTrackIndex(int i, int j);
    void computePinPowers();
};


#endif /* GPUSOLVER_H_ */
