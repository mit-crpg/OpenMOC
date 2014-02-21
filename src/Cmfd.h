/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CMFD_H_
#define CMFD_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <utility>
#include <math.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <fstream>
#include "Quadrature.h"
#include "log.h"
#include "Mesh.h"
#include "Material.h"
#include "Surface.h"
#include "Geometry.h"
#include "Timer.h"
#endif


/**
 * @enum cellType
 * @brief Eigenvalue solution methods.
*/
enum eigenMethod {

  /** The power iteration method */
  POWER,

  /** The Wielandt iteration method */
  WIELANDT
};



/**
 * @class Cell Cell.h "src/host/Cell.h"
 * @brief A class for Coarse Mesh Finite Difference (CMFD) acceleration.
 */
class Cmfd {

private:

  /** Pointer to polar Quadrature object */
  Quadrature* _quad;

  /** Pointer to Geometry object */
  Geometry* _geometry;

  /** Pointer to Mesh object */
  Mesh* _mesh;

  /** Th keff eigenvalue */
  double _k_eff;

  /** The A matrix */
  double** _A;

  /** The M matrix */
  double** _M;

  /** The AM matrix */
  double** _AM;

  /** The old source vector */
  double* _old_source;

  /** The new source vector */
  double* _new_source;

  /** A temporary scalar flux vector */
  double* _phi_temp;

  /** Gauss-Seidel SOR factor */
  double _omega;

  /** L2 norm */
  double _l2_norm;

  /** Keff convergence criteria */
  double _conv_criteria;

  /** Number of cells in x direction */
  int _cx;

  /** Number of cells in y direction */
  int _cy;

  /** Pointer to Timer object */
  Timer* _timer;

  /** Flux type (PRIMAL or ADJOINT) */
  fluxType _flux_type;

  /** Number of energy groups */
  int _num_groups;

  /** Number of coarse energy groups */
  int _num_cmfd_groups;

  /** Coarse energy indices for fine energy groups */
  int* _group_indices;

  /** Number of fine energy groups per coarse energy group */
  int _group_width;

  /** Number of FSRs */
  int _num_FSRs;

  /** Solution method (DIFFUSION or MOC) */
  solveType _solve_method;

  /** The volumes (areas) for each FSR */
  FP_PRECISION* _FSR_volumes;

  /** Pointers to Materials for each FSR */
  Material** _FSR_materials;

  /** The FSR scalar flux in each energy group */
  FP_PRECISION* _FSR_fluxes;

  /* Eigenvalue method */
  eigenMethod _eigen_method;

public:

  Cmfd(Geometry* geometry, double criteria=1e-8);
  virtual ~Cmfd();

  /* Worker functions */
  void constructMatrices();
  void computeDs();
  void computeXS();
  void updateMOCFlux();
  double computeDiffCorrect(double d, double h);
  double computeKeff();
  void initializeFSRs();
  void rescaleFlux();
  void linearSolve(double** mat, double* vec_x, double* vec_b,
                   double conv, int max_iter=10000);

  /* Matrix and Vector functions */
  void dumpVec(double* vec, int length);
  void matZero(double** mat, int width);
  void vecCopy(double* vec_from, double* vec_to);
  double vecSum(double* vec);
  void matMultM(double** mat, double* vec_x, double* vec_y);
  void matMultA(double** mat, double* vec_x, double* vec_y);
  void vecNormal(double** mat, double* vec);
  void vecSet(double* vec, double val);
  void vecScale(double* vec, double scale_val);
  void matSubtract(double** AM, double** A, double omega, double** M);
  double vecMax(double* vec);
  double rayleighQuotient(double* x, double* snew, double* sold);
  void createGroupStructure();

  /* Get parameters */
  Mesh* getMesh();
  double getKeff();
  int getNumCmfdGroups();
  int getCmfdGroupWidth();

  /* Set parameters */
  void setOmega(double omega);
  void setFluxType(const char* flux_type);
  void setEigenMethod(const char* eigen_method);
  void setNumCmfdGroups(int num_cmfd_groups);

  /* Set FSR parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
};

#endif /* CMFD_H_ */
