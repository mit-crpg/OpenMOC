/**
 * @file MICSolver.h
 * @brief The MIColver class.
 * @date June 18, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef MICSOLVER_H_
#define MICSOLVER_H_

#ifdef SWIG
#define MIC_ATTRIBUTE
#else
#define MIC_ATTRIBUTE __attribute__((target(mic)))
#endif


#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <offload.h>
#include "../../Solver.h"
#include "clone.h"
#endif

#define thread_flux(t,r,e) (thread_flux[(t)*_num_FSRs*_num_groups + (r)*_num_groups + (e)])


/**
 * @class MICSolver MICSolver.h "openmoc/src/host/MICSolver.h"
 * @brief
 */
class MICSolver : public Solver {

private:

    /** The angular fluxes for each track for all energy groups, polar angles,
     *  and azimuthal angles. This array stores the boundary fluxes for a
     *  a track along both "forward" and "reverse" directions. */
  //    FP_PRECISION* _boundary_flux;

    /** The number of shared memory OpenMP threads */
    int _num_threads;

    /** A pointer to an array of the materials on the device */
    dev_material* _materials;

    /** A pointer to the array of tracks on the device */
    dev_track* _dev_tracks;

    /** An array of the cumulative number of tracks for each azimuthal angle */
    int* _track_index_offsets;

    void initializeFluxArrays();
    void initializeSourceArrays();
    void initializePowerArrays();
    void initializePolarQuadrature();
    void precomputePrefactors();
    void initializeFSRs();

    void initializeMaterials();
    void initializeTracks();

    void zeroTrackFluxes();
    void flattenFSRFluxes(FP_PRECISION value);
    void flattenFSRSources(FP_PRECISION value);
    void normalizeFluxes();
    FP_PRECISION computeFSRSources();
    void computeKeff();
    bool isScalarFluxConverged();
    void transportSweep(int max_iterations);

public:
    MICSolver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL);
    virtual ~MICSolver();

    FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
    FP_PRECISION* getFSRScalarFluxes();
    FP_PRECISION* getFSRPowers();
    FP_PRECISION* getFSRPinPowers();
    int getNumThreads();

    void setNumThreads(int num_threads);
    void setGeometry(Geometry* geometry);
    void setTrackGenerator(TrackGenerator* track_generator);

    int computeScalarTrackIndex(int i, int j);
    void computePinPowers();
};


#endif /* MICSOLVER_H_ */
