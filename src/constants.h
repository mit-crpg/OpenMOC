/**
 * @file constants.h
 * @brief Math constants and comparision tolerances.
 * @date April 9, 2015.
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

/** Threshold to determine if a float equals to 0.0 */
#define FLT_EPSILON 1.0E-12

/** Threshold to determine if a flux is equal to 0.0 */
#define FLUX_EPSILON FP_PRECISION(1.0E-25)

/** Threshold to determine if a float is equal to infinity */
#define FLT_INFINITY 1.0E300

/* The single line width permissible for timer / memory logger reports */
#define REPORT_WIDTH 53

/** The minimum auto ID used for Surfaces, Cells, Materials and Universes */
#define DEFAULT_INIT_ID 1000000

/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.566370614359172

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.07957747154594767

/** A negligible cross-section value to over-ride user-defined
 *  cross-sections very near zero (e.g., within (-1E-10, 1E-10)) */
#define ZERO_SIGMA_T 1E-6

/** Threshold to determine how close the sum of \f$ \Sigma_a \f$ and
 *  \f$ \Sigma_s \f$ must match \f$ \Sigma_t \f$ for each energy group */
#define SIGMA_T_THRESH 1E-10

/** Distance a Point is moved to cross over a Surface into a new Cell */
#define TINY_MOVE 1E-8

/** Threshold to determine if a Point is on the boundary of a Lattice cell */
#define ON_LATTICE_CELL_THRESH 1E-12

/** Error threshold to determine if a point is to be considered on a Surface */
#define ON_SURFACE_THRESH 1E-12

/** Tolerance for difference of the sum of polar weights with respect to 1.0 */
#define POLAR_WEIGHT_SUM_TOL 1E-5

/** The default maximum optical path length */
#define MAX_OPTICAL_LENGTH FP_PRECISION(100.)

/** A small amount to increment the tau, the max optical path length, to ensure
 *  that tracks with the max optical path length are not split. */
#define TAU_NUDGE 1E-12

/** The minimum acceptable precision for exponential evaluations from
 *  the ExpEvaluator's linear interpolation table. This default precision
 *  was selected based on analysis by Yamamoto's 2004 paper on the topic. */
#define EXP_PRECISION FP_PRECISION(1E-5)

/** The minimum number of interpolation points to be used in an exponential
 *  lookup table */
#define MIN_EXP_INTERP_POINTS 100

/** The minimum calculated determinant to allow for the calculation of a matrix
  * inverse. */
#define MIN_DET 1E-10

/** The maximum number of iterations allowed for a power method eigenvalue
 *  solve in linalg.cpp */
#define MIN_LINALG_POWER_ITERATIONS 25
#define MAX_LINALG_POWER_ITERATIONS 25000

#define MIN_LINALG_TOLERANCE LINALG_TOL

/** The maximum number of iterations allowed for a linear solve in linalg.cpp */
#define MIN_LINEAR_SOLVE_ITERATIONS 25
#define MAX_LINEAR_SOLVE_ITERATIONS 10000

#ifdef MPIx
//TODO Make tracks per buffer dependent on number of processes, and groups
#define TRACKS_PER_BUFFER 2000
#define CMFD_BUFFER_SIZE 10000
#endif

#define LOCAL_COORDS_LEN 16
#define MAX_VERSION_NUM 50

/** The faces, edges, and vertices that collectively make up the surfaces of a
 *  rectangular prism. The edges denoted as "e" and vertices as "v" on the
 *  illustration below:
 *
 *                                   e
 *                       v +--------------------+ v
 *                        /|                   /|
 *                       / |                  / |
 *                      /  |                 /  |
 *                    e/  e|                /e  |e
 *                    /    |               /    |
 *                   /     |         e    /     |
 *                  /    v +-------------/------+ v
 *               v +--------------------+ v    /
 *                 |     /              |     /
 *                 |    /               |    /
 *                e|  e/                |e e/
 *                 |  /                 |  /
 *                 | /                  | /
 *                 |/                   |/
 *               v +--------------------+ v
 *                           e
 *
 */
#define NUM_FACES 6
#define NUM_EDGES 12
#define NUM_VERTICES 8
#define NUM_SURFACES 26
#define SURFACE_X_MIN 0
#define SURFACE_Y_MIN 1
#define SURFACE_Z_MIN 2
#define SURFACE_X_MAX 3
#define SURFACE_Y_MAX 4
#define SURFACE_Z_MAX 5
#define SURFACE_X_MIN_Y_MIN 6
#define SURFACE_X_MAX_Y_MIN 7
#define SURFACE_X_MIN_Y_MAX 8
#define SURFACE_X_MAX_Y_MAX 9
#define SURFACE_X_MIN_Z_MIN 10
#define SURFACE_X_MAX_Z_MIN 11
#define SURFACE_X_MIN_Z_MAX 12
#define SURFACE_X_MAX_Z_MAX 13
#define SURFACE_Y_MIN_Z_MIN 14
#define SURFACE_Y_MAX_Z_MIN 15
#define SURFACE_Y_MIN_Z_MAX 16
#define SURFACE_Y_MAX_Z_MAX 17
#define SURFACE_X_MIN_Y_MIN_Z_MIN 18
#define SURFACE_X_MIN_Y_MIN_Z_MAX 19
#define SURFACE_X_MIN_Y_MAX_Z_MIN 20
#define SURFACE_X_MIN_Y_MAX_Z_MAX 21
#define SURFACE_X_MAX_Y_MIN_Z_MIN 22
#define SURFACE_X_MAX_Y_MIN_Z_MAX 23
#define SURFACE_X_MAX_Y_MAX_Z_MIN 24
#define SURFACE_X_MAX_Y_MAX_Z_MAX 25

/** The number of values used in representing a Track when the Tracks are
   *  retrieved from the TrackGenerator. */
#define NUM_VALUES_PER_RETRIEVED_TRACK 6

/** The number of values used in representing a Segment when the Segments are
   *  retrieved from the TrackGenerator. */
#define NUM_VALUES_PER_RETRIEVED_SEGMENT 7

/** Least common multiple tolerance */
#define LCM_TOLERANCE 1.e-8

#ifdef NVCC

/** The maximum number of polar angles to reserve constant memory on GPU */
#define MAX_POLAR_ANGLES_GPU 10

/** The maximum number of azimuthal angles to reserve constant memory on GPU */
#define MAX_AZIM_ANGLES_GPU 256

#endif

#endif /* CONSTANTS_H_ */
