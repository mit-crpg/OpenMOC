/**
 * @file linalg.h
 * @details This file contains a library of functions for performing linear
 *          algebra operations.
 * @date July 3, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "log.h"
#include "Matrix.h"
#include "Vector.h"
#include "constants.h"
#include <math.h>
#include <vector>
#include <omp.h>
#endif


//FIXME
struct ConvergenceData {
  double pf;
  double cmfd_res_1;
  double cmfd_res_end;
  double linear_res_1;
  double linear_res_end;
  int cmfd_iters;
  int linear_iters_1;
  int linear_iters_end;
};


//FIXME
struct DomainCommunicator {
  int _num_domains_x;
  int _num_domains_y;
  int _num_domains_z;
  int _domain_idx_x;
  int _domain_idx_y;
  int _domain_idx_z;
  int _local_num_x;
  int _local_num_y;
  int _local_num_z;
  int _offset;
  int** num_connections;
  int*** indexes;
  int*** domains;
  CMFD_PRECISION*** coupling_coeffs;
  CMFD_PRECISION*** fluxes;
  CMFD_PRECISION** buffer;
  int num_groups;
  bool stop;
#ifdef MPIx
  MPI_Comm _MPI_cart;
#endif
};

//FIXME
#ifdef MPIx
void getCouplingTerms(DomainCommunicator* comm, int color, int*& coupling_sizes,
                      int**& coupling_indexes, CMFD_PRECISION**& coupling_coeffs,
                      CMFD_PRECISION**& coupling_fluxes,
                      CMFD_PRECISION* curr_fluxes, int& offset);
#endif

double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double k_eff,
                       double tol, double SOR_factor=1.5,
                       ConvergenceData* convergence_data = NULL,
                       DomainCommunicator* comm = NULL);
bool linearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                 double SOR_factor=1.5,
                 ConvergenceData* convergence_data = NULL,
                 DomainCommunicator* comm = NULL);
void matrixMultiplication(Matrix* A, Vector* X, Vector* B);
double computeRMSE(Vector* x, Vector* y, bool integrated, int it,
                         DomainCommunicator* comm = NULL);
void oldLinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                    double SOR_factor=1.5,
                    ConvergenceData* convergence_data = NULL);

/**
 * @brief Transpose a 2D matrix.
 * @param matrix array to transpose
 * @param dim1 first dimension length
 * @param dim2 second dimension length
 */
template<typename T>
inline void matrix_transpose(T* matrix, int dim1, int dim2) {

  std::vector<T> temp(dim1 * dim2);

  for (int i=0; i < dim1; i++) {
    for (int j=0; j < dim2; j++)
      temp[i * dim1 + j] = matrix[j * dim1 + i];
  }

  std::copy(temp.begin(), temp.end(), matrix);
}

#endif /* LINALG_H_ */
