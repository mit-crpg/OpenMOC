/**
 * @file DeviceSolver.h
 * @brief The DeviceSolver class and GPU solver routines.
 * @date August 5, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef DEVICESOLVER_H_
#define DEVICESOLVER_H_

#ifdef __cplusplus
//#include "../host/Solver.h"
#include "../host/Geometry.h"
#include "../host/TrackGenerator.h"
#include "../host/Quadrature.h"
#endif

#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include "DeviceFlatSourceRegion.h"
#include "DeviceTrack.h"
#include "DeviceMaterial.h"

#define scalar_flux(tid,e) (scalar_flux[(tid)*(*_num_groups_devc) + (e)])
#define old_scalar_flux(tid,e) (old_scalar_flux[(tid)*(*_num_groups_devc) + (e)])
#define source(tid,e) (source[(tid)*(*_num_groups_devc) + (e)])
#define old_source(tid,e) (old_source[(tid)*(*_num_groups_devc) + (e)])
#define ratios(tid,e) (ratios[(tid)*(*_num_groups_devc) + (e)])
#define boundary_flux(tid,pe2) (boundary_flux[2*(tid)*(*_polar_times_groups_devc)+(pe2)])

/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.5663706143

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.0795774715

/** The maximum number of polar angles to reserve constant memory on GPU */
#define MAX_POLAR_ANGLES 3

/** The maximum number of azimuthal angles to reserve constant memory on GPU */
#define MAX_AZIM_ANGLES 256


/**
 * @class DeviceSolver DeviceSolver.h "openmoc/src/dev/DeviceSolver.h"
 * @brief
 */
class DeviceSolver {

private:

    /**************************************************************************/
    /*                             Data on the host                           */
    /**************************************************************************/

    /** The number of threadblocks */
    int _num_blocks;
    
    /** The number of threads per threadblock */
    int _num_threads;

    /** The number of azimuthal angles */
    int _num_azim;

    /** The number of energy groups */
    int _num_groups;

    /** The number of flat source regions */
    int _num_FSRs;

    /** The number of materials */
    int _num_materials;

    /** The number of polar angles */
    int _num_polar;

    /** Twice the number of polar angles */
    int _two_times_num_polar;

    /** The number of polar angles times energy groups */
    int _polar_times_groups;

    /** The type of polar quadrature (TABUCHI or LEONARD) */
    quadratureType _quadrature_type;

    /** A pointer to a trackgenerator which contains tracks */
    TrackGenerator* _track_generator;

    /** A pointer to a geometry with initialized flat source region maps */
    Geometry* _geometry;

    /** A pointer to a polar quadrature */
    Quadrature* _quad;

    /** A pointer to the 2D ragged array of tracks on the host */
    Track** _host_tracks;

    /** A pointer to an array with the number of tracks per azimuthal angle */
    int* _num_tracks;

    /** The total number of tracks */
    int _tot_num_tracks;

    /** The total leakage across vacuum boundaries */
    FP_PRECISION _leakage;

    /** The number of transport sweeps to convergence */
    int _num_iterations;

    /** The current iteration's keff */
    FP_PRECISION _k_eff;

    /** Whether or not the Solver has converged the source */
    bool _converged_source;

    /** The tolerance for converging the source */
    FP_PRECISION _source_convergence_thresh;

    /** The tolerance for converging the flux given a fixed source */
    FP_PRECISION _flux_convergence_thresh;

    //TODO: What are these guys for???
    FP_PRECISION* _FSRs_to_powers;
    FP_PRECISION* _FSRs_to_pin_powers;


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

    /** The angular fluxes for each track for all energy groups, polar angles,
     *  and azimuthal angles. This array stores the boundary fluxes for a
     *  a track along both "forward" and "reverse" directions. */
    FP_PRECISION* _boundary_flux;

    /** The scalar flux for each energy group in each flat source region */
    FP_PRECISION* _scalar_flux;

    /** The scalar flux for each energy group in each flat source region from
     *  the previous iteration */
    FP_PRECISION* _old_scalar_flux;

    /** The source in each energy group in each flat source region */
    FP_PRECISION* _source;

    /** The source in each energy group in each flat source region from the 
     *  previous iteration */
    FP_PRECISION* _old_source;

    /** Pre-computed Ratio of source / sigma_t for each energy group in each
     *  flat source region */
    FP_PRECISION* _ratios;

    FP_PRECISION* _fission_source;
    FP_PRECISION* _tot_abs;
    FP_PRECISION* _tot_fission;
    FP_PRECISION* _source_residual;
    thrust::device_vector<FP_PRECISION> _fission_source_vec;
    thrust::device_vector<FP_PRECISION> _tot_abs_vec;
    thrust::device_vector<FP_PRECISION> _tot_fission_vec;
    thrust::device_vector<FP_PRECISION> _source_residual_vec;

    /** The hashtable of exponential prefactors from the transport equation */
    FP_PRECISION* _prefactor_array;

    void initializeHostMemory();
    void initializePolarQuadrature();
    void initializeGlobalMemory();
    void initializeFSRs();
    void initializeMaterials();
    void initializeTracks();
    void initializeFluxArrays();
    void initializeSourceArrays();
    void initializePowerArrays();
    void initializeThrustVectors();
    void initializeConstantMemory();
    void precomputePrefactors();
    void checkTrackSpacing();

    //FP_PRECISION computeFSRSources();
    //void computeKeff();
    //bool isScalarFluxConverged();
    //int computePrefactorIndex(FP_PRECISION sigma_t_l);
    //void transportSweep(int max_iterations);

public:

    DeviceSolver(Geometry* geometry=NULL, TrackGenerator* track_generator=NULL);
    virtual ~DeviceSolver();
    
    Geometry* getGeometry();
    TrackGenerator* getTrackGenerator();
    int getNumPolarAngles();
    quadratureType getPolarQuadratureType();
    int getNumIterations();
    FP_PRECISION getSourceConvergenceThreshold();
    FP_PRECISION getFluxConvergenceThreshold();
    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();

    void setGeometry(Geometry* geometry);
    void setTrackGenerator(TrackGenerator* track_generator);
    void setPolarQuadratureType(quadratureType quadrature_type);
    void setNumPolarAngles(int num_polar);
    void setSourceConvergenceThreshold(FP_PRECISION source_thresh);
    void setFluxConvergenceThreshold(FP_PRECISION flux_thresh);
    void setNumThreadBlocks(int num_blocks);
    void setNumThreadsPerBlock(int num_threads);

    void allocateDeviceData();
    //    void validateDeviceData();
    int computeScalarTrackIndex(int i, int j);

    FP_PRECISION convergeSource(int max_iterations, int B=64, int T=64);
    //    void computePinPowers();
};


#endif /* DEVICESOLVER_H_ */
