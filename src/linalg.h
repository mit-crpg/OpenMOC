
/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#ifdef __cplusplus
#include <math.h>
#include "log.h"
#include <omp.h>
#include "Matrix.h"
#include "Vector.h"
#endif

FP_PRECISION eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, FP_PRECISION tol,
                             FP_PRECISION SOR_factor=1.5);
void linearSolve(Matrix* A, Vector* X, Vector* B, FP_PRECISION tol,
                 FP_PRECISION SOR_factor=1.5);
void matrixMultiplication(Matrix* A, Vector* X, Vector* B);

#endif /* LINALG_H_ */
