/* File: linalg.cpp */

#include "linalg.h"

FP_PRECISION eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, FP_PRECISION tol,
                             FP_PRECISION SOR_factor) {

  log_printf(INFO, "computing eigenvalue...");

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != M->getNumX() || A->getNumX() != X->getNumX())
    log_printf(NORMAL, "Cannot compute eigenvalue with different x dimensions"
               " for the A matrix, M matrix, and X vector: (%i, %i, %i)",
               A->getNumX(), M->getNumX(), X->getNumX());
  else if (A->getNumX() != M->getNumY() || A->getNumY() != X->getNumY())
    log_printf(NORMAL, "Cannot compute eigenvalue with different y dimensions"
               " for the A matrix, M matrix, and X vector: (%i, %i, %i)",
               A->getNumY(), M->getNumY(), X->getNumY());
  else if (A->getNumGroups() != M->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log_printf(NORMAL, "Cannot compute eigenvalue with different num groups"
               " for the A matrix, M matrix, and X vector: (%i, %i, %i)",
               A->getNumGroups(), M->getNumGroups(), X->getNumGroups());
  
  /* Initialize variables */
  int num_rows = X->getNumRows();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_groups = X->getNumGroups();
  Vector old_source(num_x, num_y, num_groups);
  Vector new_source(num_x, num_y, num_groups);
  FP_PRECISION residual, _k_eff;

  /* Compute and normalize the initial source */
  matrixMultiplication(M, X, &old_source);
  old_source.scaleByValue(num_rows / old_source.getSum());
  X->scaleByValue(num_rows / old_source.getSum());

  /* Power iteration diffusion solver */
  for (int iter = 0; iter < 25000; iter++) {

    /* Solve phi = A^-1 * old_source */
    linearSolve(A, X, &old_source, tol*1e1, SOR_factor);

    /* Compute the new source */
    matrixMultiplication(M, X, &new_source);

    /* Compute and set keff */
    _k_eff = new_source.getSum() / num_rows;

    /* Scale the old source by keff */
    old_source.scaleByValue(_k_eff);

    /* Compute the L2 norm of source error */
    #pragma omp parallel for
    for (int i = 0; i < num_x*num_y; i++) {
      for (int g = 0; g < num_groups; g++) {
        if (new_source.getValue(i, g) != 0.0)
          old_source.setValue
            (i, g, pow((new_source.getValue(i, g) - old_source.getValue(i, g))
                       / new_source.getValue(i, g), 2));
      }
    }

    /* Compute the source RMS error */
    residual = sqrt(old_source.getSum() / num_rows);

    /* Normalize the new source to have an average value of 1.0 */
    new_source.scaleByValue(num_rows / new_source.getSum());
    new_source.copyTo(&old_source);

    log_printf(INFO, "CMFD iter: %i, keff: %f, error: %g",
               iter, _k_eff, residual);

    /* Check for convergence */
    if (residual < tol && iter > 10)
      break;
  }

  return _k_eff;
}


/**
 * @brief Solve the linear system AX=B using Gauss Seidel with SOR.
 * @param pointer to A matrix
 * @param pointer to X vector
 * @param pointer to B vector
 * @param flux convergence criteria
 * @param the maximum number of iterations
 */
void linearSolve(Matrix* A, Vector* X, Vector* B, FP_PRECISION tol,
                 FP_PRECISION SOR_factor) {

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != B->getNumX() || A->getNumX() != X->getNumX())
    log_printf(NORMAL, "Cannot perform linear solve with different x dimensions"
               " for the A matrix, B vector, and X vector: (%i, %i, %i)",
               A->getNumX(), B->getNumX(), X->getNumX());
  else if (A->getNumX() != B->getNumY() || A->getNumY() != X->getNumY())
    log_printf(NORMAL, "Cannot perform linear solve with different y dimensions"
               " for the A matrix, B vector, and X vector: (%i, %i, %i)",
               A->getNumY(), B->getNumY(), X->getNumY());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log_printf(NORMAL, "Cannot perform linear solve with different num groups"
               " for the A matrix, B vector, and X vector: (%i, %i, %i)",
               A->getNumGroups(), B->getNumGroups(), X->getNumGroups());
  
  FP_PRECISION residual = 1E10;
  int iter = 0;
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_groups = X->getNumGroups();
  int num_rows = X->getNumRows();
  Vector X_old(num_rows);
  FP_PRECISION* x_old = X_old.getArray();
  int* IA = A->getIA();
  int* JA = A->getJA();
  FP_PRECISION* DIAG = A->getDiag();
  FP_PRECISION* a = A->getA();
  FP_PRECISION* x = X->getArray();
  FP_PRECISION* b = B->getArray();
  int row, col;

  while (iter < 1000) {

    /* Pass new flux to old flux */
    X->copyTo(&X_old);

    /* Iteration over red/black cells */
    for (int color = 0; color < 2; color++) {
      #pragma omp parallel for private(row, col)
      for (int cy = 0; cy < num_y; cy++) {
        for (int cx = 0; cx < num_x; cx++) {
          
          /* check for correct color */
          if (((cx % 2)+(cy % 2)) % 2 == color) {

            for (int g = 0; g < num_groups; g++) {
              
              row = (cy*num_x + cx)*num_groups + g;
              
              /* Over-relax the xx array */
              x[row] = (1.0 - SOR_factor) * x[row];
              
              for (int i = IA[row]; i < IA[row+1]; i++) {
                
                col = JA[i];
                
                if (row == col)
                  x[row] += SOR_factor * b[row] / DIAG[row];
                else
                  x[row] -= SOR_factor * a[i] * x[col] / DIAG[row];
              }
            }

            for (int g = num_groups-1; g >= 0; g--) {
              
              row = (cy*num_x + cx)*num_groups + g;
              
              /* Over-relax the xx array */
              x[row] = (1.0 - SOR_factor) * x[row];
              
              for (int i = IA[row]; i < IA[row+1]; i++) {
                
                col = JA[i];
                
                if (row == col)
                  x[row] += SOR_factor * b[row] / DIAG[row];
                else
                  x[row] -= SOR_factor * a[i] * x[col] / DIAG[row];
              }
            }
          }
        }
      }
    }
  
    /* Compute the average residual */
    for (int i = 0; i < num_rows; i++) {
      if (x[i] != 0.0)
        x_old[i] = pow((x[i] - x_old[i]) / x[i], 2);
      else
        x_old[i] = 0.0;
    }
    
    residual = sqrt(pairwise_sum(x_old, num_rows) / num_rows);
    
    /* Increment the interations counter */
    iter++;
    
    log_printf(INFO, "SOR iter: %i, res: %g", iter, residual);
    
    if (residual < tol && iter > 10)
      break;
  }

  log_printf(DEBUG, "linear solver iterations: %i", iter);
}


/**
 * @brief Multiply matrix by vector (i.e., B = M*X).
 * @param matrix source matrix
 * @param vector_x x vector
 * @param vector_y y vector
 * @param num_blocks number of cell blocks in M matrix.
 * @param block_width number of elements in cell blocks in M matrix.
 */
void matrixMultiplication(Matrix* A, Vector* X, Vector* B) {
  
  int* IA = A->getIA();
  int* JA = A->getJA();
  FP_PRECISION* a = A->getA();
  FP_PRECISION* x = X->getArray();
  FP_PRECISION* b = B->getArray();
  int num_rows = X->getNumRows();
  std::fill_n(b, num_rows, 0.0);
  
  #pragma omp parallel for
  for (int row = 0; row < num_rows; row++) {
    for (int i = IA[row]; i < IA[row+1]; i++)
      b[row] += a[i] * x[JA[i]];
  }
}
