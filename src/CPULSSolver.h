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
#define _scalar_flux_xyz(r,e,x) (_scalar_flux_xyz[(r)*_num_groups*3 + (x)*_num_groups + (e)])

/** Indexing macro for the total source divided by the total cross-section
 *  (\f$ \frac{Q}{\Sigma_t} \f$) in each FSR and energy group */
#define _reduced_sources_xyz(r,e,x) (_reduced_sources_xyz[(r)*_num_groups*3 + (x)*_num_groups + (e)])

/** Indexing macro for the stabilizing scalar flux moments in each FSR and
 *  energy group */
#define _stabilizing_flux_xyz(r,e,x) (_stabilizing_flux_xyz[(r)*_num_groups*3 + (e)*3 + (x)])

/** Indexing scheme for fixed source moments for each FSR and energy group */
#define _fixed_sources_xyz(r,e,i) (_fixed_sources_xyz[(r)*_num_groups + (e)][i])

/**
 * @class CPULSSolver CPULSSolver.h "src/CPULSSolver.h"
 * @brief This a subclass of the CPUSolver class for using the linear source
 *        approximation.
 */
class CPULSSolver : public CPUSolver {

protected:

  /** The FSR linear expansion matrix values for each FSR */
  double* _FSR_lin_exp_matrix;

  /** The FSR source constants for each FSR and energy group */
  FP_PRECISION* _FSR_source_constants;

  /** An array of the scalar flux x, y, and z terms */
  FP_PRECISION* _scalar_flux_xyz;

  /** An array of the reduced source x, y, and z terms */
  FP_PRECISION* _reduced_sources_xyz;

  /** The stabilizing flux for each energy group in each FSR */
  FP_PRECISION* _stabilizing_flux_xyz;

  /** Optional user-specified fixed linear source in each FSR & energy group */
  std::vector<std::vector<double> > _fixed_sources_xyz;

  /** A map of fixed sources moments keyed by the pair (FSR ID, energy group) */
  std::map< std::pair<int, int>, std::vector<double> > _fix_src_xyz_FSR_map;

  /** A map of fixed sources moments keyed by the pair (Cell*, energy group) */
  std::map< std::pair<Cell*, int>, std::vector<double> > _fix_src_xyz_cell_map;

  /** Whether to stabilize the flux moments */
  bool _stabilize_moments;

  /** Whether fixed linear source moments have been provided */
  bool _fixed_source_moments_on;

public:
  CPULSSolver(TrackGenerator* track_generator=NULL);
  virtual ~CPULSSolver();

  /* Initialization routines */
  void initializeFluxArrays();
  void initializeSourceArrays();
  void initializeCmfd();
  void initializeExpEvaluators();
  void initializeFSRs();

  /* Routines to handle fixed source moments */
  void initializeFixedSources();
  void setFixedSourceMomentsByCell(Cell* cell, int group, double source_x,
                                   double source_y, double source_z);
  void setFixedSourceMomentByFSR(long fsr_id, int group, double source_x,
                                 double source_y, double source_z);
  void resetFixedSources();

  /* Worker routines */
  void flattenFSRFluxes(FP_PRECISION value);
  double normalizeFluxes();
  void computeFSRSources(int iteration);
  void tallyLSScalarFlux(segment* curr_segment, int azim_index,
                                    int polar_index,
                                    FP_PRECISION* fsr_flux,
                                    FP_PRECISION* fsr_flux_x,
                                    FP_PRECISION* fsr_flux_y,
                                    FP_PRECISION* fsr_flux_z,
                                    float* track_flux,
                                    FP_PRECISION direction[3]);
  void accumulateLinearFluxContribution(long fsr_id, FP_PRECISION weight,
                                        FP_PRECISION* fsr_flux);
  void addSourceToScalarFlux();

  /* Transport stabilization routines */
  void computeStabilizingFlux();
  void stabilizeFlux();
  void checkLimitXS(int iteration);

  /* Routines to handle constant part of linear source */
  void initializeLinearSourceConstants();
  double* getLinearExpansionCoeffsBuffer();
  FP_PRECISION* getSourceConstantsBuffer();

  /* Getter routine */
  FP_PRECISION getFluxByCoords(LocalCoords* coords, int group);
};


#endif /* CPULSSOLVER_H_ */
