/**
 * @file constants.h
 * @brief Math constants and comparision tolerances.
 * @date April 9, 2015.
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_


/* The single line width permissible for timer / memory logger reports */
#define REPORT_WIDTH 53

/** The minimum auto ID used for Surfaces, Cells, Materials and Universes */
#define DEFAULT_INIT_ID 10000

/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.566370614359172

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.07957747154594767

/** A negligible cross-section value to over-ride user-defined
 *  cross-sections very near zero (e.g., within (-1E-10, 1E-10)) */
#define ZERO_SIGMA_T 1E-10

/** Threshold to determine how close the sum of \f$ \Sigma_a \f$ and
 *  \f$ \Sigma_s \f$ must match \f$ \Sigma_t \f$ for each energy group */
#define SIGMA_T_THRESH 1E-3

/** Distance a Point is moved to cross over a Surface into a new Cell */
#define TINY_MOVE 1E-8

/** Threshold to determine if a Point is on the boundary of a Lattice cell */
#define ON_LATTICE_CELL_THRESH 1E-12

/** Error threshold to determine if a point is to be considered on a Surface */
#define ON_SURFACE_THRESH 1E-12

/** Tolerance for difference of the sum of polar weights with respect to 1.0 */
#define POLAR_WEIGHT_SUM_TOL 1E-5

/** The default maximum optical path length */
#define MAX_OPTICAL_LENGTH FP_PRECISION(10.)

/** A small amount to increment the tau, the max optical path length, to ensure
 *  that tracks with the max optical path length are not split. */
#define TAU_NUDGE 1E-12

/** The minimum acceptable precision for exponential evaluations from
 *  the ExpEvaluator's linear interpolation table. This default precision
 *  was selected based on analysis by Yamamoto's 2004 paper on the topic. */
#define EXP_PRECISION FP_PRECISION(1E-5)

/** The maximum number of iterations allowed for a power method eigenvalue
 *  solve in linalg.cpp */
#define MIN_LINALG_POWER_ITERATIONS 10
#define MAX_LINALG_POWER_ITERATIONS 25000

/** The maximum number of iterations allowed for a linear solve in linalg.cpp */
#define MIN_LINEAR_SOLVE_ITERATIONS 10
#define MAX_LINEAR_SOLVE_ITERATIONS 1000

/** The faces and edges that collectively make up the surfaces of a
 *  horizontal slice of a rectangular prism. The faces are denoted
 *  as "f" and edges denoted as "e" on the illustration below:
 *
 *          e            f            e
 *           +----------------------+
 *           |                      |
 *           |                      |
 *           |                      |
 *         f |                      | f
 *           |                      |
 *           |                      |
 *           |                      |
 *           +----------------------+
 *          e            f            e
 *
 */
#define NUM_FACES 4
#define NUM_EDGES 4
#define NUM_SURFACES 8
#define SURFACE_X_MIN 0
#define SURFACE_Y_MIN 1
#define SURFACE_X_MAX 2
#define SURFACE_Y_MAX 3
#define SURFACE_X_MIN_Y_MIN 4
#define SURFACE_X_MAX_Y_MIN 5
#define SURFACE_X_MAX_Y_MAX 6
#define SURFACE_X_MIN_Y_MAX 7

/** The number of values used in representing a Track when the Tracks are
 *  retrieved from the TrackGenerator. */
#define NUM_VALUES_PER_RETRIEVED_TRACK 6

/** The number of values used in representing a Segment when the Segments are
 *  retrieved from the TrackGenerator. */
#define NUM_VALUES_PER_RETRIEVED_SEGMENT 7

#ifdef NVCC

/** The maximum number of polar angles to reserve constant memory on GPU */
#define MAX_POLAR_ANGLES_GPU 10

/** The maximum number of azimuthal angles to reserve constant memory on GPU */
#define MAX_AZIM_ANGLES_GPU 256

#endif

#endif /* CONSTANTS_H_ */
