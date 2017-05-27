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
#define _scalar_flux_xyz(r,e,x) (_scalar_flux_xyz[(r)*_num_groups*3 + (e)*3 + (x)])

/** Indexing macro for the total source divided by the total cross-section
 *  (\f$ \frac{Q}{\Sigma_t} \f$) in each FSR and energy group */
#define _reduced_sources_xyz(r,e,x) (_reduced_sources_xyz[(r)*_num_groups*3 + (e)*3 + (x)])


/**
 * @class CPULSSolver CPULSSolver.h "src/CPULSSolver.h"
 * @brief This a subclass of the CPUSolver class for using the linear source
 *        approximation.
 */
class CPULSSolver : public CPUSolver {

protected:

  /** The FSR linear expansion matrix values for each FSR */
  NEW_PRECISION* _FSR_lin_exp_matrix;

  /** The FSR source constants for each FSR and energy group */
  NEW_PRECISION* _FSR_source_constants;

  /** An array of the scalar flux x, y, and z terms */
  NEW_PRECISION* _scalar_flux_xyz;

  /** An array of the reduced source x, y, and z terms */
  NEW_PRECISION* _reduced_sources_xyz;

public:
  CPULSSolver(TrackGenerator* track_generator=NULL);
  virtual ~CPULSSolver();

  void initializeFluxArrays();
  void initializeSourceArrays();
  void initializeCmfd();
  void initializeExpEvaluators();
  void initializeFSRs();

  void flattenFSRFluxes(NEW_PRECISION value);
  double normalizeFluxes();
  void computeFSRSources(int iteration);
  void addSourceToScalarFlux();

    //FIXME MEM : float / FP_PRECISION
  void tallyLSScalarFlux(segment* curr_segment, int azim_index, int polar_index,
                         float* track_flux, NEW_PRECISION* fsr_flux,
                         NEW_PRECISION direction[3]);

  FP_PRECISION getFluxByCoords(LocalCoords* coords, int group);
  NEW_PRECISION* getLinearExpansionCoeffsBuffer();
  NEW_PRECISION* getSourceConstantsBuffer();
};


#endif /* CPULSSOLVER_H_ */
