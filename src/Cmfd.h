/*
 * Cmfd.h
 *
 *  Created on: September 13, 2012
 *      Author: Sam Shaner
 *				MIT, Course 22
 *              shaner@mit.edu
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
#ifdef _mm_malloc
#undef _mm_malloc
#undef _mm_free
#endif
#include "petsc.h"
#include <petscmat.h>
#endif

/**
 * Solve types
 */
enum solveType {
	DIFFUSION,
	MOC
};

/**
 * Flux types
 */
enum fluxType {
	PRIMAL,
	ADJOINT
};


class Cmfd {
private:
	Quadrature* _quad;
	Geometry* _geometry;
	Mesh* _mesh;
	double _k_eff;
	Mat _A;
	Mat _M;
	PetscScalar* _A_array;
	PetscScalar* _M_array;
	PetscInt* _indices_y_M;
	PetscInt* _indices_y_A;
	Vec _phi_new;
	Vec _phi_old;
	Vec _source_old;
	Vec _sold;
	Vec _snew;
	Vec _res;
	double _keff;
	double _l2_norm;
	double _l2_norm_old;
	double _keff_0;
	double _keff_conv;
	int _cells_x;
	int _cells_y;
	Timer* _timer;
	bool _assemble_M;
	double _relax_factor;
	solveType _solve_method;
	fluxType _flux_method;
	int _num_groups;
	int _num_fsrs;
	bool _optically_thick;

	FP_PRECISION* _FSR_volumes;
	Material** _FSR_materials;
	FP_PRECISION* _FSR_fluxes;


public:
	Cmfd(Geometry* geometry, solveType solve_method, double relax_factor=0.6);
	virtual ~Cmfd();
 	void computeDs();
 	void computeXS();
 	double computeDiffCorrect(double d, double h);
 	void updateMOCFlux();
 	int constructMatrices();
 	double computeKeff();
 	Mat getA();
 	Mat getM();
 	Vec getPhiNew();
 	int createAMPhi();
 	double getKeff();
	void initializeFSRs();
	int fisSourceNorm(Vec snew);
	double getL2Norm();
	int rescaleFlux();
	int setMeshCellFlux();
	void assembleM(bool assembleM);
	solveType getSolveType();
	void setRelaxFactor(double relax_factor);
	void setFSRMaterials(Material** FSR_materials);
	void setFSRVolumes(FP_PRECISION* FSR_volumes);
	void setFSRFluxes(FP_PRECISION* scalar_flux);
	Mesh* getMesh();
	void toggleFluxType(fluxType flux_method);
	void opticallyThick(bool thick);
};

#endif /* CMFD_H_ */
