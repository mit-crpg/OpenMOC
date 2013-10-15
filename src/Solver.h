/**
 * @file Solver.h
 * @brief The Solver class.
 * @date February 7, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include "Timer.h"
#include "Quadrature.h"
#include "TrackGenerator.h"
#include "pairwise_sum.h"
#include "Cmfd.h"
#endif

/** Indexing macro for the scalar flux in each flat source region and in
 *  each energy group */
#define _scalar_flux(r,e) (_scalar_flux[(r)*_num_groups + (e)])

/** Indexing macro for the surface currents for each mesh surface and in
 *  each energy group */
#define _surface_currents(r,e) (_surface_currents[(r % _geometry->getMesh()->getNumCurrents())*_num_groups + (e)])

/** Indexing macro for the total source in each flat source region and in
 *  each energy group */
#define _source(r,e) (_source[(r)*_num_groups + (e)])

/** Indexing macro for the total source from the previous source iteration
 *  in each flat source region and in each energy group */
#define _old_source(r,e) (_old_source[(r)*_num_groups + (e)])

/** Indexing macro for the total source divided by the total cross-section
 *  (\f$ \frac{Q}{\Sigma_t} \f$) in each flat source region and in each
 *  energy group */
#define _reduced_source(r,e) (_reduced_source[(r)*_num_groups + (e)])

/** Indexing macro for the polar quadrature weights multiplied by the 
 *  azimuthal angle quadrature weights */
#define _polar_weights(i,p) (_polar_weights[(i)*_num_polar + (p)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for the outgoing reflective track for both the forward and
 *  reverse direction for a given track */
#define _boundary_flux(i,j,p,e) (_boundary_flux[(i)*2*_polar_times_groups + (j)*_polar_times_groups + (p)*_num_groups + (e)])

/** Indexing macro for the leakage for each polar angle and energy group
 *  for both the forward and reverse direction for each track */
#define _boundary_leakage(i,pe2) (_boundary_leakage[2*(i)*_polar_times_groups+(pe2)])

/** Indexing scheme for the total fission source (\f$ \nu\Sigma_f\Phi \f$) 
 *  for each flat source region in each energy group */
#define _fission_sources(r,e) (_fission_sources[(r)*_num_groups + (e)])

/** Indexing scheme for the total in-scatter source (\f$ \Sigma_s\Phi \f$) 
 *  for each flat source region in each energy group */
#define _scatter_sources(r,e) (_scatter_sources[(r)*_num_groups + (e)])

/** Indexing scheme for the residual between sources from this iteration
 *  and the previous iteration in each flat source region and energy group */
#define _source_residuals(r,e) (_source_residuals[(r)*_num_groups + (e)])

/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.5663706143

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.0795774715


/**
 * @class Solver Solver.h "openmoc/src/host/Solver.h"
 * @brief This is an abstract base class from which different types of Solvers subclass for
 *        different architectures or using different algorithms.
 */
class Solver {

protected:

    /** The number of azimuthal angles */
    int _num_azim;

    /** The number of energy groups */
    int _num_groups;

    /** The number of flat source regions */
    int _num_FSRs;

    /** The number of mesh cells */
    int _num_mesh_cells;

    /** The flat source region "volumes" (ie, areas) index by FSR UID */
    FP_PRECISION* _FSR_volumes;

    /** The flat source region material pointers index by FSR UID */
    Material** _FSR_materials;

    /** A pointer to a trackgenerator which contains tracks */
    TrackGenerator* _track_generator;

    /** A pointer to a geometry with initialized flat source region maps */
    Geometry* _geometry;

    /** The number of materials */
    int _num_materials;

    /** A pointer to a polar quadrature */
    Quadrature* _quad;

    /** The number of polar angles */
    int _num_polar;

    /** Twice the number of polar angles */
    int _two_times_num_polar;

    /** The number of polar angles times energy groups */
    int _polar_times_groups;

    /** The type of polar quadrature (TABUCHI or LEONARD) */
    quadratureType _quadrature_type;

    /** A pointer to the 2D ragged array of tracks */
    Track** _tracks;

    /** A pointer to an array with the number of tracks per azimuthal angle */
    int* _num_tracks;

    /** The total number of tracks */
    int _tot_num_tracks;

    /** The weights for each azimuthal angle */
    double* _azim_weights;

    /** The weights for each polar angle in the polar angle quadrature */
    FP_PRECISION* _polar_weights;

    /** The angular fluxes for each track for all energy groups, polar angles,
     *  and azimuthal angles. This array stores the boundary fluxes for a
     *  a track along both "forward" and "reverse" directions. */
    FP_PRECISION* _boundary_flux;

    /** The angular leakages for each track for all energy groups, polar angles,
     *  and azimuthal angles. This array stores the weighted outgoing fluxes 
     *  for a track along both "forward" and "reverse" directions. */
    FP_PRECISION* _boundary_leakage;

    /* Flat source regions */
    /** The scalar flux for each energy group in each flat source region */
    FP_PRECISION* _scalar_flux;

    /* mesh surfaces */
    /** The surface currents for each energy group for each mesh surface */
    double* _surface_currents;

    /** The fission source in each energy group in each flat source region */
    FP_PRECISION* _fission_sources;

    /** The in-scatter source in each energy group in each flat source region */
    FP_PRECISION* _scatter_sources;

    /** The source in each energy group in each flat source region */
    FP_PRECISION* _source;

    /** The source in each energy group in each flat source region from the 
     *  previous iteration */
    FP_PRECISION* _old_source;

    /** Pre-computed ratio of source / sigma_t for each energy group in each
     *  flat source region */
    FP_PRECISION* _reduced_source;

    /** An array of the residuals between the old source and the new source
     *  on each iteration in each flat source region and energy group */
    FP_PRECISION* _source_residuals;

    /** The current iteration's approximation to k-effective */
    FP_PRECISION _k_eff; 

    /** The total leakage across vacuum boundaries */
    FP_PRECISION _leakage;

    /** The number of transport sweeps to convergence */
    int _num_iterations;

    /** Whether or not the Solver has converged the source */
    bool _converged_source;

    /** The tolerance for converging the source */
    FP_PRECISION _source_convergence_thresh;

    /** A boolean indicating whether or not to use linear interpolation 
     *  to comptue the exponential in the transport equation */
    bool _interpolate_exponential;

    /* Exponential pre-factor hash table */
    /** The hashtable of exponential prefactors from the transport equation */
    FP_PRECISION* _prefactor_array;

    /** The size of the exponential prefactor array */
    int _prefactor_array_size;

    /** The maximum index of the exponential prefactor array */
    int _prefactor_max_index;

    /** The spacing for the exponential prefactor array */
    FP_PRECISION _prefactor_spacing;

    /** The inverse spacing for the exponential prefactor array */
    FP_PRECISION _inverse_prefactor_spacing;

    /** A timer to record timing data for a simulation */
    Timer* _timer;

    Cmfd* _cmfd;

    /**
     * @brief Creates a polar quadrature object for the solver.
     */
    virtual void initializePolarQuadrature() =0;

    /**
     * @brief Initializes the volumes and material arrays for each flat source 
     *        region. 
     */
    virtual void initializeFSRs() =0;

    /**
     * @brief Allocates memory for track boundary angular fluxes and leakages
     *        flat source region scalar fluxes.
     */
    virtual void initializeFluxArrays() =0;

    /**
     * @brief Allocates memory for flat source region source arrays.
     */
    virtual void initializeSourceArrays() =0;

    /**
     * @brief Builds an interpolation table for the exponential prefactors 
     *        referenced for each segment in the transport equation.
     */
    virtual void precomputePrefactors() =0;

    virtual void checkTrackSpacing();

    /**
     * @brief Zero each track's boundary fluxes for each energy group and polar
     *        angle in the "forward" and "reverse" directions.
     */
    virtual void zeroTrackFluxes() =0;

    /**
     * @brief Set the scalar flux for each energy group inside each flat source 
     *        region to a constant value.
     * @param value the value to assign to each flat source region flux
     */
    virtual void flattenFSRFluxes(FP_PRECISION value) =0;

    /**
     * @brief Set the source for each energy group in each flat source region
     *        to a constant value.
     * @param value the value to assign to each flat source region source
     */
    virtual void flattenFSRSources(FP_PRECISION value) =0;

    /**
     * @brief Normalizes all flat source region scalar fluxes and track boundary
     *        angular fluxes to the total fission source (times \f$ \nu \f$).
     */
    virtual void normalizeFluxes() =0;

    /**
     * @brief Computes the total source (fission and scattering) in each flat 
     *        source region.
     *
     * @return the residual between this source and the previous source
     */
    virtual FP_PRECISION computeFSRSources() =0;

    /**
     * @brief Compute \f$ k_{eff} \f$ from total fission and absorption rates.
     */
    virtual void computeKeff() =0;

    /**
     * @brief Add the source term contribution in the transport equation to 
     *        the flat source region scalar flux
     */
    virtual void addSourceToScalarFlux() =0;

    /**
     * @brief This method performs one transport sweep of all azimuthal angles, 
     *        tracks, segments, polar angles and energy groups.
     */
    virtual void transportSweep() =0;

    void clearTimerSplits();

public:
    Solver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL, Cmfd* cmfd=NULL);
    virtual ~Solver();

    Geometry* getGeometry();
    TrackGenerator* getTrackGenerator();
    int getNumPolarAngles();
    quadratureType getPolarQuadratureType();
    int getNumIterations();
    FP_PRECISION getSourceConvergenceThreshold();

    /**
     * @brief Returns the scalar flux for a flat source region
     * @param fsr_id the ID for the FSR of interest
     * @param energy_group the energy group of interest
     * @return the flat source region scalar flux
     */
    virtual FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group) =0;

    /**
     * @brief Returns an array of the scalar flux in each flat source
     *        region in each energy group.
     * @return an array of flat source region scalar fluxes
     */
    virtual FP_PRECISION* getFSRScalarFluxes() =0;

    virtual void setGeometry(Geometry* geometry);
    virtual void setTrackGenerator(TrackGenerator* track_generator);
    virtual void setPolarQuadratureType(quadratureType quadrature_type);
    virtual void setNumPolarAngles(int num_polar);
    virtual void setSourceConvergenceThreshold(FP_PRECISION source_thresh);

    void useExponentialInterpolation();
    void useExponentialIntrinsic();

    virtual FP_PRECISION convergeSource(int max_iterations);
    
    virtual void computeFSRFissionRates(double* fission_rates, int num_FSRs) =0;
    void printTimerReport();
    void initializeCmfd();
};


#endif /* SOLVER_H_ */
