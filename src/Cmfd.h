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

/* petsc input files */
#ifdef CMFD
#ifdef _mm_malloc
#undef _mm_malloc
#undef _mm_free
#endif
#include "petsc.h"
#include <petscmat.h>
#endif
#endif

/**
 * Flux types
 */
enum fluxType {
	PRIMAL,
	ADJOINT
};


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

#ifdef CMFD
  /* loss matrix */
  Mat _A;

  /* source matrix */
  Mat _M;

  /* row insertion arrays */
  PetscScalar* _A_array;
  PetscScalar* _M_array;
  PetscInt* _indices_y_M;
  PetscInt* _indices_y_A;

  /* flux, source, and residual vectors */
  Vec _phi_new;
  Vec _phi_old;
  Vec _source_old;
  Vec _sold;
  Vec _snew;
  Vec _res;

#endif

  /* l2 norm */
  double _l2_norm;

  /* keff convergence criteria */
  double _conv_criteria;

  /* num cells in x and y direction */
  int _cells_x;
  int _cells_y;

  /* pointer to timer object */
  Timer* _timer;

  /* flag to determine whether we need to assemble M */
  bool _assemble_M;
  
  /* flux type (PRIMAL or ADJOINT) */
  fluxType _flux_method;

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
  int constructMatrices();
  void computeDs();
  void computeXS();
  void updateMOCFlux();
  double computeDiffCorrect(double d, double h);
  double computeKeff();
  int createAMPhi();
  void initializeFSRs();
  int rescaleFlux();

  /* get arrays, values, and objects */
#ifdef CMFD
  Mat getA();
  Mat getM();
#endif
  Mesh* getMesh();
  double getKeff();

  /* set parameters */
  int setMeshCellFlux();
  void assembleM(bool assembleM);
  void toggleFluxType(fluxType flux_method);

  /* set fsr parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
};

#endif /* CMFD_H_ */
