
/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#ifdef __cplusplus
#include<math.h>
#include "log.h"
#include<omp.h>
#include "Matrix.h"
#include "Vector.h"
#endif

double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double tol, double SOR_factor=1.5);
void linearSolve(Matrix* A, Vector* X, Vector* B, double tol, double SOR_factor=1.5);
void matrixMultiplication(Matrix* A, Vector* X, Vector* B);

#endif /* LINALG_H_ */
