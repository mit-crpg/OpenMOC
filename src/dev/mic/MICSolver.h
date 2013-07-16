#ifndef MICSOLVER_H_
#define MICSOLVER_H_

#ifdef SWIG
#define MIC_ATTRIBUTE
#else
#define MIC_ATTRIBUTE __attribute__((target(mic)))
#endif


#ifdef __cplusplus
#include <stdio.h>
#include <omp.h>
#include <offload.h>
#include "../../Solver.h"
#include "clone.h"
#include "MICQuery.h"
#endif


#define boundary_flux(i,pe2) (boundary_flux[2*(i)*_polar_times_groups+(pe2)])
#define scalar_flux(r,e) (scalar_flux[(r)*_num_groups + (e)])
#define old_scalar_flux(r,e) (old_scalar_flux[(r)*_num_groups + (e)])
#define source(r,e) (source[(r)*_num_groups + (e)])
#define old_source(r,e) (old_source[(r)*_num_groups + (e)])

#define thread_flux(t,r,e) (thread_flux[(t)*_num_FSRs*_num_groups + (r)*_num_groups + (e)])



#define _scalar_flux(r,e) (_scalar_flux[(r)*_num_groups + (e)])
#define _old_scalar_flux(r,e) (_old_scalar_flux[(r)*_num_groups + (e)])
#define _source(r,e) (_source[(r)*_num_groups + (e)])
#define _old_source(r,e) (_old_source[(r)*_num_groups + (e)])
#define _ratios(r,e) (_ratios[(r)*_num_groups + (e)])
#define _polar_weights(i,p) (_polar_weights[(i)*_num_polar + (p)])
#define _boundary_flux(i,pe2) (_boundary_flux[2*(i)*_polar_times_groups+(pe2)])
#define _fission_source(r,e) (_fission_source[(r)*_num_groups + (e)])
#define source_residuals(r,e) (source_residuals[(r)*_num_groups + (e)])

#define prefactor(index,p,sigma_t_l) (1. - (_prefactor_array[index+2 * p] * sigma_t_l + _prefactor_array[index + 2 * p +1]))

/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.5663706143

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.0795774715



class MICSolver : public Solver {

private:

    /** The flat source region material uids */
    int* _FSR_materials;

    size_t FSR_materials;

    size_t FSR_volumes;

    /** The number of shared memory OpenMP threads */
    int _num_threads;

    /** A pointer to an array of the materials on the device */
    dev_material* _materials;

    size_t materials;

    /** A pointer to the array of tracks on the device */
    dev_track* _dev_tracks;

    size_t dev_tracks;

    /** An array of the cumulative number of tracks for each azimuthal angle */
    int* _track_index_offsets;

    size_t track_index_offsets;


    size_t boundary_flux;
    size_t scalar_flux;
    size_t old_scalar_flux;

    size_t fission_source;
    size_t source;
    size_t old_source;
    size_t ratios;

    size_t prefactor_array;
    size_t polar_weights;

    void initializeFluxArrays();
    void initializeSourceArrays();
    void initializePowerArrays();
    void initializePolarQuadrature();
    void precomputePrefactors();
    void initializeFSRs();
    void initializeMaterials();
    void initializeTracks();

    int computeScalarTrackIndex(int i, int j);

    MIC_ATTRIBUTE void zeroTrackFluxes();
    void flattenFSRFluxes(FP_PRECISION value);
    void flattenFSRSources(FP_PRECISION value);
    void normalizeFluxes();
    FP_PRECISION computeFSRSources();
    void computeKeff();
    bool isScalarFluxConverged();
    void transportSweep(int max_iterations);

    MIC_ATTRIBUTE int computePrefactorIndex(FP_PRECISION sigma_t_l);

public:
    MICSolver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL);
    virtual ~MICSolver();

    int getNumThreads();
    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();

    void setNumThreads(int num_threads);
    void setGeometry(Geometry* geometry);
    void setTrackGenerator(TrackGenerator* track_generator);

    FP_PRECISION convergeSource(int max_iterations);

    void computePinPowers();
};


#endif /* MICSOLVER_H_ */
