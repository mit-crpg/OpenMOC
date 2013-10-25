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
#include <unordered_map>
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

class Cmfd {

private:

  /* pointer to quadrature object */
  Quadrature* _quad;

  /* pointer to geometry object */
  Geometry* _geometry;

  /* pointer to mesh object */
  Mesh* _mesh;

  /* keff */
  double _k_eff;

  /* matrix and vector objects */
  double* _A;
  double* _M;
  double* _sold;
  double* _snew;
  double* _phi_temp;

  /* Gauss Seidel SOR factor */
  double _omega;

  /* l2 norm */
  double _l2_norm;

  /* keff convergence criteria */
  double _conv_criteria;

  /* num cells in x and y direction */
  int _cells_x;
  int _cells_y;

  /* pointer to timer object */
  Timer* _timer;

  /* flux type (PRIMAL or ADJOINT) */
  fluxType _flux_type;

  /* number of groups */
  int _num_groups;

  /* number of fsrs */
  int _num_fsrs;

  /* solve method (DIFFUSION or MOC) */
  solveType _solve_method;
  
  /* arrays for fsr parameters */
  FP_PRECISION* _FSR_volumes;
  Material** _FSR_materials;
  FP_PRECISION* _FSR_fluxes;
  
public:
	
  Cmfd(Geometry* geometry, double criteria=1e-8);
  virtual ~Cmfd();

  /* worker functions */
  void constructMatrices();
  void computeDs();
  void computeXS();
  void updateMOCFlux();
  double computeDiffCorrect(double d, double h);
  double computeKeff();
  void initializeFSRs();
  void rescaleFlux();
  void linearSolve(double* mat, double* vec_x, double* vec_b, double conv);

  /* matrix and vector functions */
  void dumpVec(double* vec, int length);
  void matZero(double* mat, int width);
  void vecCopy(double* vec_from, double* vec_to);
  double vecSum(double* vec);
  void matMult(double* mat, double* vec_x, double* vec_y);
  void vecNormal(double* mat, double* vec);
  void vecSet(double* vec, double val);
  void vecScale(double* vec, double scale_val);
      
  /* get parameters */
  Mesh* getMesh();
  double getKeff();

  /* set parameters */
  void setOmega(double omega);
  void setFluxType(const char* flux_type);
  
  /* set fsr parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
};

#endif /* CMFD_H_ */
