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

/** Indexing macro for the stabilizing scalar flux moments in each FSR and
 *  energy group */
#define _stabilizing_flux_xyz(r,e,x) (_stabilizing_flux_xyz[(r)*_num_groups*3 + (e)*3 + (x)])

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

  /** Whether to stabilize the flux moments */
  bool _stabilize_moments;

public:
  CPULSSolver(TrackGenerator* track_generator=NULL);
  virtual ~CPULSSolver();

  void initializeFluxArrays();
  void initializeSourceArrays();
  void initializeCmfd();
  void initializeExpEvaluators();
  void initializeFSRs();

  void flattenFSRFluxes(FP_PRECISION value);
  double normalizeFluxes();
  void computeFSRSources(int iteration);
  void addSourceToScalarFlux();
  
  /* Transport stabilization routines */
  void computeStabilizingFlux();
  void stabilizeFlux();
  void checkLimitXS(int iteration);

  FP_PRECISION getFluxByCoords(LocalCoords* coords, int group);
  void initializeLinearSourceConstants();
  double* getLinearExpansionCoeffsBuffer();
  FP_PRECISION* getSourceConstantsBuffer();

/**
 * @brief Computes the contribution to the LSR scalar flux from a Track segment.
 * @details This method integrates the angular flux for a Track segment across
 *          energy groups and polar angles, and tallies it into the LSR
 *          scalar flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index azimuthal angle index for this 3D Track
 * @param polar_index polar angle index for this 3D Track
 * @param track_flux a pointer to the Track's angular flux
 * @param direction the segment's direction
 */
inline void tallyLSScalarFlux(segment* curr_segment, int azim_index,
                                    int polar_index,
                                    float* track_flux,
                                    FP_PRECISION direction[3]) {

  long fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION* position = curr_segment->_starting_position;
  ExpEvaluator* exp_evaluator = _exp_evaluators[azim_index][polar_index];

  long start_scalar_idx = fsr_id * _num_groups;
  long start_linear_idx = start_scalar_idx * 3;

  /* Allocate a temporary flux buffer on the stack (free) and initialize it */
  FP_PRECISION fsr_flux[_num_groups] = {0.0};
  FP_PRECISION fsr_flux_x[_num_groups] = {0.0};
  FP_PRECISION fsr_flux_y[_num_groups] = {0.0};
  FP_PRECISION fsr_flux_z[_num_groups] = {0.0};

  if (_solve_3D) {

    /* Compute the segment midpoint */
    FP_PRECISION center_x2[3];
    for (int i=0; i<3; i++)
      center_x2[i] = 2 * (position[i] + 0.5 * length * direction[i]);

    FP_PRECISION wgt = _quad->getWeightInline(azim_index, polar_index);
    FP_PRECISION length_2D = exp_evaluator->convertDistance3Dto2D(length);

    // Compute the exponential terms
    FP_PRECISION exp_F1[_num_groups];
    FP_PRECISION exp_F2[_num_groups];
    FP_PRECISION exp_H[_num_groups];
#pragma omp simd
    for (int e=0; e < _num_groups; e++) {
      FP_PRECISION tau = sigma_t[e] * length_2D;
      exp_evaluator->retrieveExponentialComponents(tau, 0, &exp_F1[e],
                                                   &exp_F2[e],
                                                   &exp_H[e]);
    }

    // Compute the sources
    FP_PRECISION src_flat[_num_groups];
    FP_PRECISION src_linear[_num_groups];
#pragma omp simd
    for (int e=0; e < _num_groups; e++) {
      src_flat[e] = _reduced_sources[start_scalar_idx+e];
      for (int i=0; i<3; i++)
        src_flat[e] += _reduced_sources_xyz[start_linear_idx+3*e+i] *
                       center_x2[i];
      src_linear[e] = _reduced_sources_xyz[start_linear_idx+3*e] *
                         direction[0];
      src_linear[e] += _reduced_sources_xyz[start_linear_idx+3*e+1] *
                         direction[1];
      src_linear[e] += _reduced_sources_xyz[start_linear_idx+3*e+2] *
                         direction[2];
    }

    // Compute the flux attenuation and tally contribution
#pragma omp simd
    for (int e=0; e < _num_groups; e++) {

      FP_PRECISION tau = sigma_t[e] * length_2D;

      // Compute the change in flux across the segment
      exp_H[e] *= length * track_flux[e] * tau * wgt;
      FP_PRECISION delta_psi = (tau * track_flux[e] - length_2D * src_flat[e]) *
          exp_F1[e] - src_linear[e] * length_2D * length_2D *
          exp_F2[e];
      track_flux[e] -= delta_psi;

      // Increment the fsr scalar flux and scalar flux moments
      fsr_flux[e] += wgt * delta_psi;
      fsr_flux_x[e] += exp_H[e] * direction[0] + wgt * 
                                 delta_psi * position[0];
      fsr_flux_y[e] += exp_H[e] * direction[1] + wgt * 
                                 delta_psi * position[1];
      fsr_flux_z[e] += exp_H[e] * direction[2] + wgt * 
                                 delta_psi * position[2];
    }
  }
  else {

    int num_polar_2 = _num_polar / 2;

    /* Compute the segment midpoint */
    FP_PRECISION center[2];
    for (int i=0; i<2; i++)
      center[i] = position[i] + 0.5 * length * direction[i];

    /* Compute exponentials */
    FP_PRECISION exp_F1[num_polar_2*_num_groups];
    FP_PRECISION exp_F2[num_polar_2*_num_groups];
    FP_PRECISION exp_H[num_polar_2*_num_groups];
    for (int p=0; p < num_polar_2; p++) {
#pragma omp simd
      for (int e=0; e < _num_groups; e++) {
        FP_PRECISION tau = sigma_t[e] * length;
        exp_evaluator->retrieveExponentialComponents(tau, p, 
                                                     &exp_F1[p*_num_groups+e],
                                                     &exp_F2[p*_num_groups+e],
                                                     &exp_H[p*_num_groups+e]);
      }
    }

    /* Compute flat part of source */
    FP_PRECISION src_flat[_num_groups] = {0.0};
#pragma omp simd
    for (int e=0; e < _num_groups; e++) {
      for (int i=0; i<2; i++)
        src_flat[e] += _reduced_sources_xyz(fsr_id, e, i) * center[i];
      src_flat[e] *= 2;
      src_flat[e] += _reduced_sources(fsr_id, e);
    }

    /* Compute linear part of source */
    FP_PRECISION src_linear[num_polar_2 * _num_groups] = {0.0};
    for (int p=0; p < num_polar_2; p++) {
      FP_PRECISION sin_theta = _quad->getSinTheta(azim_index, p);
#pragma omp simd
      for (int e=0; e < _num_groups; e++) {
        for (int i=0; i<2; i++)
          src_linear[p*_num_groups+e] += direction[i] * sin_theta *
              _reduced_sources_xyz(fsr_id, e, i);
      }
    }

    /* Compute attenuation and tally flux */
    for (int p=0; p < num_polar_2; p++) {
      FP_PRECISION wgt = _quad->getWeightInline(azim_index, p);
#pragma omp simd
      for (int e=0; e < _num_groups; e++) {
        FP_PRECISION tau = sigma_t[e] * length;
        exp_H[p*_num_groups+e] *=  tau * length * track_flux[p*_num_groups+e];

        // Compute the change in flux across the segment
        FP_PRECISION delta_psi = (tau * track_flux[p*_num_groups+e] - length
              * src_flat[e]) * exp_F1[p*_num_groups+e] - length * length 
              * src_linear[p*_num_groups+e] * exp_F2[p*_num_groups+e];
        track_flux[p*_num_groups+e] -= delta_psi;

        // Increment the fsr scalar flux and scalar flux moments
        fsr_flux[e] += wgt * delta_psi;
        fsr_flux_x[e] += wgt * (exp_H[p*_num_groups+e] * direction[0] +
              delta_psi * position[0]);
        fsr_flux_y[e] += wgt * (exp_H[p*_num_groups+e] * direction[1] +
              delta_psi * position[1]);
      }
    }
  }

  // Atomically increment the FSR scalar flux from the temporary array
  omp_set_lock(&_FSR_locks[fsr_id]);
#pragma omp simd
  for (int e=0; e < _num_groups; e++) {

    // Add to global scalar flux vector
    _scalar_flux[start_scalar_idx + e] += fsr_flux[e];
    _scalar_flux_xyz[start_linear_idx + 3*e] += fsr_flux_x[e];
    _scalar_flux_xyz[start_linear_idx + 3*e + 1] += fsr_flux_y[e];
    _scalar_flux_xyz[start_linear_idx + 3*e + 2] += fsr_flux_z[e];
  }

  omp_unset_lock(&_FSR_locks[fsr_id]);

  for (int i=0; i < 3; i++)
    position[i] += direction[i] * length;
}

};


#endif /* CPULSSOLVER_H_ */
