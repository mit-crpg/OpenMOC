/* File: linalg.cpp */

#include "linalg.h"

double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double tol, double SOR_factor){

  log_printf(NORMAL, "computing eigenvalue...");

  /* Initialize variables */
  int num_rows = X->getNumRows();
  Vector old_source(num_rows);
  Vector new_source(num_rows);
  double sum_new, val, residual, scale_val, _k_eff;
  int row;

  /* Compute and normalize the initial source */
  matrixMultiplication(M, X, &old_source);
  sum_new = pairwise_sum(old_source.getArray(), num_rows);
  scale_val = double(num_rows) / sum_new;
  old_source.scaleByValue(scale_val);
  X->scaleByValue(scale_val);

  /* Power iteration diffusion solver */
  for (int iter = 0; iter < 25000; iter++){

    /* Solve phi = A^-1 * old_source */
    linearSolve(A, X, &old_source, tol*1e2, SOR_factor);

    /* Compute the new source */
    matrixMultiplication(M, X, &new_source);
    sum_new = pairwise_sum(new_source.getArray(), num_rows);

    /* Compute and set keff */
    _k_eff = sum_new / double(num_rows);

    /* Scale the old source by keff */
    old_source.scaleByValue(_k_eff);

    /* Compute the L2 norm of source error */
    #pragma omp parallel for
    for (int i = 0; i < num_rows; i++){
      if (new_source.getValue(i) != 0.0)
        old_source.setValue(i, pow((new_source.getValue(i) - old_source.getValue(i))
                                    / new_source.getValue(i), 2));
    }

    /* Compute the source RMS error */
    residual = sqrt(pairwise_sum(old_source.getArray(), num_rows) / num_rows);

    /* Normalize the new source to have an average value of 1.0 */
    new_source.scaleByValue(num_rows / sum_new);
    new_source.copyTo(&old_source);

    log_printf(NORMAL, "CMFD iter: %i, keff: %f, error: %f",
               iter, _k_eff, residual);

    /* Check for convergence */
    if (residual < tol && iter > 10)
      break;
  }

  return _k_eff;
}


/**
 * @brief Solve the linear system Ax=b using Gauss Seidel with SOR.
 * @param pointer to A matrix
 * @param pointer to x vector
 * @param pointer to b vector
 * @param flux convergence criteria
 * @param the maximum number of iterations
 */
void linearSolve(Matrix* A, Vector* X, Vector* B, double tol, double SOR_factor){

  double residual = 1E10;
  int iter = 0;
  int num_rows = X->getNumRows();
  Vector X_old(num_rows);
  double* x_old = X_old.getArray();
  int* IA = A->getIA();
  int* JA = A->getJA();
  double* DIAG = A->getDIAG();
  double* a = A->getA();
  double* x = X->getArray();
  double* b = B->getArray();
  int num_x = A->getNumX();
  int num_y = A->getNumY();
  int num_z = A->getNumZ();
  int num_g = A->getNumGroups();
  int row, col;

  while (iter < 1000){

    /* Pass new flux to old flux */
    X->copyTo(&X_old);

    /* Iteration over red/black cells */
    for (int color = 0; color < 2; color++){
      for (int cz = 0; cz < num_z; cz++){
        #pragma omp parallel for private(row, col)
        for (int cy = 0; cy < num_y; cy++){
          for (int cx = 0; cx < num_x; cx++){

            /* check for correct color */
            if (((cx % 2)+(cy % 2)+(cz % 2)) % 2 == color){

              for (int g = 0; g < num_g; g++){

                row = ((cz*num_y + cy)*num_x + cx)*num_g + g;
                
                /* Over-relax the xx array */
                x[row] = (1.0 - SOR_factor) * x[row];
                
                for (int i = IA[row]; i < IA[row+1]; i++){
                  
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
    }

    /* Compute the average residual */
    for (int i = 0; i < num_rows; i++){
      if (x[i] != 0.0)
        x_old[i] = pow((x[i] - x_old[i]) / x[i], 2);
      else
        x_old[i] = 0.0;
    }

    residual = sqrt(pairwise_sum(x_old, num_rows) / num_rows);

    /* Increment the interations counter */
    iter++;

    log_printf(INFO, "GS iter: %i, res: %f", iter, residual);

    if (residual < tol && iter > 10)
      break;
  }

  log_printf(DEBUG, "linear solver iterations: %i", iter);
}


/**
 * @brief Multiply matrix by vector (i.e., y = M *x).
 * @param matrix source matrix
 * @param vector_x x vector
 * @param vector_y y vector
 * @param num_blocks number of cell blocks in M matrix.
 * @param block_width number of elements in cell blocks in M matrix.
 */
void matrixMultiplication(Matrix* A, Vector* X, Vector* B){  

  int* IA = A->getIA();
  int* JA = A->getJA();
  double* a = A->getA();
  double* x = X->getArray();
  double* b = B->getArray();
  int num_rows = X->getNumRows();
  std::fill_n(b, num_rows, 0.0);
  
  #pragma omp parallel for
  for (int row = 0; row < num_rows; row++){
    for (int i = IA[row]; i < IA[row+1]; i++)
      b[row] += a[i] * x[JA[i]];
  }
}
