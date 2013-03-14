/*
 * configurations.h
 *
 *  Created on: Jan 16, 2012
 *      Author: will
 */

#ifndef CONFIGURATIONS_H_
#define CONFIGURATIONS_H_


/******************************************************************************
 ****************************** USER DEFINED **********************************
 *****************************************************************************/

/* Define number or polar angles (2 or 3) and the number of energy groups */
#define NUM_POLAR_ANGLES 3
#define NUM_ENERGY_GROUPS 7
#define GRP_TIMES_ANG NUM_POLAR_ANGLES*NUM_ENERGY_GROUPS


/* Floating point precision (DOUBLE or SINGLE) */
#define SINGLE

#ifdef DOUBLE
#define FP_PRECISION double
#else
#define FP_PRECISION float
#endif



/* Convergence threshold for computing k_eff */
#define SOURCE_CONVERG_THRESH 1E-3


/* Convergence threshold for scalar flux in each region during fixed source
 * iteration */
#define FLUX_CONVERGENCE_THRESH 1E-5


/** Maximum number of fixed source iterations allowed */
#define MAX_ITERATIONS 10E3


/* Precompute and store exponential pre-factors in transport equation */
#define STORE_PREFACTORS false


/******************************************************************************
 *********************** PHYSICAL CONSTANTS ***********************************
 *****************************************************************************/

#define PI 3.1415926536
#define FOUR_PI 12.5663706143
#define ONE_OVER_FOUR_PI 0.0795774715


/******************************************************************************
 *************************** ERROR THRESHOLDS *********************************
 *****************************************************************************/

/* Error threshold for determining how close the sum of sigma_a and sigma_s
 * must match that of sigma_t for each energy group */
#define SIGMA_T_THRESH 1E-3


/* Error threshold for determining how close a point needs to be to a surface
 * to be considered on it */
#define ON_SURFACE_THRESH 1E-12


/* Error threshold for determining how close to the boundary of a lattice cell
 * a point needs to be to be considered on it */
#define ON_LATTICE_CELL_THRESH 1E-12


/* Distance a point is moved to cross over a surface into a new cell during
 * track segmentation */
#define TINY_MOVE 1E-10


#endif /* CONFIGURATIONS_H_ */
