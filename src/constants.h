/**
 * @file constants.h
 * @brief Math constants and comparision tolerances.
 * @date April 9, 2015.
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_


/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.5663706143

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.0795774715

/** A negligible cross-section value to over-ride user-defined
 *  cross-sections very near zero (e.g., within (-1E-10, 1E-10)) */
#define ZERO_SIGMA_T 1E-10

/** Threshold to determine how close the sum of \f$ \Sigma_a \f$ and 
 *  \f$ \Sigma_s \f$ must match \f$ \Sigma_t \f$ for each energy group */
#define SIGMA_T_THRESH 1E-3

/** Distance a Point is moved to cross over a Surface into a new Cell */
#define TINY_MOVE 1E-10

/** Threshold to determine if a Point is on the boundary of a Lattice cell */
#define ON_LATTICE_CELL_THRESH 1E-12

/** Error threshold to determine if a point is to be considered on a Surface */
#define ON_SURFACE_THRESH 1E-12

/** Tolerance for difference of the sum of polar weights with respect to 1.0 */
#define POLAR_WEIGHT_SUM_TOL 1E-5

/** The default maximum optical path length */
#define MAX_OPTICAL_LENGTH FP_PRECISION(10.)

/** The minimum acceptable precision for exponential evaluations from
 *  the ExpEvaluator's linear interpolation table. This default precision
 *  was selected based on analysis by Yamamoto's 2004 paper on the topic. */
#define EXP_PRECISION FP_PRECISION(1E-5)

/** The maximum number of polar angles to reserve constant memory on GPU */
#define MAX_POLAR_ANGLES_GPU 10

/** The maximum number of azimuthal angles to reserve constant memory on GPU */
#define MAX_AZIM_ANGLES_GPU 256


#endif /* CONSTANTS_H_ */
