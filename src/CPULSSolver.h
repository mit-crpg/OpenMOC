/**
 * @file CPULSSolver.h
 * @brief The CPULSSolver class.
 * @date February 19, 2016
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */


#ifndef CPULSSOLVER_H_
#define CPULSSOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include "CPUSolver.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#endif


/** Indexing macro for the scalar flux in each FSR and energy group */
#define _scalar_flux_x(r,e) (_scalar_flux_x[(r)*_num_groups + (e)])

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _scalar_flux_y(r,e) (_scalar_flux_y[(r)*_num_groups + (e)])

/** Indexing macro for the total source divided by the total cross-section
 *  (\f$ \frac{Q}{\Sigma_t} \f$) in each FSR and energy group */
#define _reduced_sources_x(r,e) (_reduced_sources_x[(r)*_num_groups + (e)])

/** Indexing macro for the total source divided by the total cross-section
 *  (\f$ \frac{Q}{\Sigma_t} \f$) in each FSR and energy group */
#define _reduced_sources_y(r,e) (_reduced_sources_y[(r)*_num_groups + (e)])


/**
 * @class CPULSSolver CPULSSolver.h "src/CPULSSolver.h"
 * @brief This a subclass of the CPUSolver class for using the linear source
 *        approximation.
 */
class CPULSSolver : public CPUSolver {

protected:

  /** The scalar flux moments for each energy group in each FSR */
  FP_PRECISION* _scalar_flux_x;
  FP_PRECISION* _scalar_flux_y;

  /** Ratios of source moments to total cross-section for each FSR and energy
   *  group */
  FP_PRECISION* _reduced_sources_x;
  FP_PRECISION* _reduced_sources_y;

  /** The FSR linear expansion matrix values for each FSR */
  FP_PRECISION* _FSR_lin_exp_matrix;

  /** The FSR source constants for each FSR and energy group */
  FP_PRECISION* _FSR_source_constants;

  /** Array of sin(phi) */
  double* _sin_phi;

  /** Array of cos(phi) */
  double* _cos_phi;

  /**
   * @brief Computes the contribution to the FSR flux from a Track segment.
   * @param curr_segment a pointer to the Track segment of interest
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param track_flux a pointer to the Track's angular flux
   * @param fsr_flux a pointer to the temporary FSR scalar flux buffer
   * @param x the x-coord of the segment starting point
   * @param y the y-coord of the segment starting point
   * @param fwd bool indicating whether the segment is pointing forward or
   *            backwards
   */
  void tallyLSScalarFlux(segment* curr_segment, int azim_index,
                         FP_PRECISION* track_flux,
                         FP_PRECISION* fsr_flux, double x, double y,
                         bool fwd);

public:
  CPULSSolver(TrackGenerator* track_generator=NULL);
  virtual ~CPULSSolver();

  void initializeFluxArrays();
  void initializeSourceArrays();

  void flattenFSRFluxes(FP_PRECISION value);
  void normalizeFluxes();
  void computeFSRSources();
  void transportSweep();
  void addSourceToScalarFlux();

  FP_PRECISION getFluxByCoords(LocalCoords* coords, int group);
};


#endif /* CPULSSOLVER_H_ */
