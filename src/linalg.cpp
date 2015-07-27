#include "linalg.h"

/**
 * @brief Solves a generalized eigenvalue problem using the Power method.
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a tolerance used
 *          for both the power method and linear solve convergence (tol), and
 *          a successive over-relaxation factor (SOR_factor) and computes the
 *          dominant eigenvalue and eigenvector using the Power method. The
 *          eigenvalue is returned along and the input X Vector is modified in
 *          place to be the corresponding eigenvector.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 * @return k_eff the dominant eigenvalue
 */
FP_PRECISION eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, FP_PRECISION tol,
                             FP_PRECISION SOR_factor) {

  log_printf(INFO, "computing eigenvalue...");

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != M->getNumX() || A->getNumX() != X->getNumX())
    log_printf(ERROR, "Cannot compute eigenvalue with different x dimensions"
               " for the A matrix, M matrix, and X vector: (%i, %i, %i)",
               A->getNumX(), M->getNumX(), X->getNumX());
  else if (A->getNumX() != M->getNumY() || A->getNumY() != X->getNumY())
    log_printf(ERROR, "Cannot compute eigenvalue with different y dimensions"
               " for the A matrix, M matrix, and X vector: (%i, %i, %i)",
               A->getNumY(), M->getNumY(), X->getNumY());
  else if (A->getNumGroups() != M->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log_printf(ERROR, "Cannot compute eigenvalue with different num groups"
               " for the A matrix, M matrix, and X vector: (%i, %i, %i)",
               A->getNumGroups(), M->getNumGroups(), X->getNumGroups());

  /* Initialize variables */
  int num_rows = X->getNumRows();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  Vector old_source(num_x, num_y, num_z, num_groups);
  Vector new_source(num_x, num_y, num_z, num_groups);
  FP_PRECISION residual, _k_eff;
  int iter;

  /* Compute and normalize the initial source */
  matrixMultiplication(M, X, &old_source);
  old_source.scaleByValue(num_rows / old_source.getSum());
  X->scaleByValue(num_rows / old_source.getSum());

  /* Power iteration diffusion solver */
  for (iter = 0; iter < 25000; iter++) {

    /* Solve phi = A^-1 * old_source */
    linearSolve(A, M, X, &old_source, tol*1e1, SOR_factor);

    /* Compute the new source */
    matrixMultiplication(M, X, &new_source);

    /* Compute and set keff */
    _k_eff = new_source.getSum() / num_rows;

    /* Scale the old source by keff */
    old_source.scaleByValue(_k_eff);

    /* Compute the residual */
    residual = computeRMSE(&new_source, &old_source, true);

    /* Normalize the new source to have an average value of 1.0 */
    new_source.scaleByValue(num_rows / new_source.getSum());
    new_source.copyTo(&old_source);

    log_printf(INFO, "CMFD iter: %i, keff: %f, error: %g",
               iter, _k_eff, residual);

    /* Check for convergence */
    if (residual < tol && iter > 10)
      break;
  }

  log_printf(INFO, "eigenvalue solve iterations: %i", iter);

  return _k_eff;
}


/**
 * @brief Solves a linear system using Red-Black Gauss Seidel with
 *        successive over-relaxation.
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a source Vector (B),
 *          a source convergence tolerance (tol) and a successive
 *          over-relaxation factor (SOR_factor) and computes the
 *          solution to the linear system. The input X Vector is modified in
 *          place to be the solution vector.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object
 * @param B the source Vector object
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 */
void linearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, FP_PRECISION tol,
                 FP_PRECISION SOR_factor) {

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != B->getNumX() || A->getNumX() != X->getNumX() ||
      A->getNumX() != M->getNumX())
    log_printf(ERROR, "Cannot perform linear solve with different x dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%i, %i, %i, %i)", A->getNumX(), M->getNumX(),
               B->getNumX(), X->getNumX());
  else if (A->getNumX() != B->getNumY() || A->getNumY() != X->getNumY() ||
           A->getNumX() != M->getNumX())
    log_printf(ERROR, "Cannot perform linear solve with different y dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%i, %i, %i, %i)", A->getNumY(), M->getNumY(),
               B->getNumY(), X->getNumY());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups() ||
           A->getNumGroups() != M->getNumGroups())
    log_printf(ERROR, "Cannot perform linear solve with different num groups"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%i, %i, %i, %i)", A->getNumGroups(), M->getNumGroups(),
               B->getNumGroups(), X->getNumGroups());

  /* Initialize variables */
  FP_PRECISION residual;
  int iter = 0;
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
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
  Vector old_source(num_x, num_y, num_z, num_groups);
  Vector new_source(num_x, num_y, num_z, num_groups);

  /* Compute initial source */
  matrixMultiplication(M, X, &old_source);

  while (iter < 1000) {

    /* Pass new flux to old flux */
    X->copyTo(&X_old);

    /* Iteration over red/black cells */
    for (int color = 0; color < 2; color++) {
      for (int cz = 0; cz < num_z; cz++){
        #pragma omp parallel for private(row, col)
        for (int cy = 0; cy < num_y; cy++) {
          for (int cx = 0; cx < num_x; cx++) {

            /* check for correct color */
            if (((cx % 2)+(cy % 2)+(cz % 2)) % 2 == color){
              
              for (int g = 0; g < num_groups; g++) {
                
                row = ((cz*num_y + cy)*num_x + cx)*num_groups + g;
                
                /* Over-relax the x array */
                x[row] = (1.0 - SOR_factor) * x[row];
                
                for (int i = IA[row]; i < IA[row+1]; i++) {
                  
                  /* Get the column index */
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

    /* Compute the new source */
    matrixMultiplication(M, X, &new_source);

    /* Compute the residual */
    residual = computeRMSE(&new_source, &old_source, true);

    /* Copy the new source to the old source */
    new_source.copyTo(&old_source);

    /* Increment the interations counter */
    iter++;

    log_printf(INFO, "SOR iter: %i, res: %g", iter, residual);

    if (residual < tol && iter > 10)
      break;
  }

  log_printf(INFO, "linear solve iterations: %i", iter);
}


/**
 * @brief Performs a matrix vector multiplication.
 * @details This function takes in a Matrix (A), a variable Vector (X),
 *          and a solution Vector (B) and computes the matrix vector product.
 *          The solution Vector is modified in place.
 * @param A a Matrix object
 * @param X the variable Vector object
 * @param B the solution Vector object
 */
void matrixMultiplication(Matrix* A, Vector* X, Vector* B) {

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != B->getNumX() || A->getNumX() != X->getNumX())
    log_printf(ERROR, "Cannot perform matrix multiplication  with different x "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%i, %i, %i)", A->getNumX(), B->getNumX(), X->getNumX());
  else if (A->getNumX() != B->getNumY() || A->getNumY() != X->getNumY())
    log_printf(ERROR, "Cannot perform matrix multiplication with different y "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%i, %i, %i)", A->getNumY(), B->getNumY(), X->getNumY());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log_printf(ERROR, "Cannot perform matrix multiplication with different "
               "num groups for the A matrix, B vector, and X vector: "
               "(%i, %i, %i)", A->getNumGroups(), B->getNumGroups(),
               X->getNumGroups());

  B->setAll(0.0);
  int* IA = A->getIA();
  int* JA = A->getJA();
  FP_PRECISION* a = A->getA();
  FP_PRECISION* x = X->getArray();
  FP_PRECISION* b = B->getArray();
  int num_rows = X->getNumRows();

  #pragma omp parallel for
  for (int row = 0; row < num_rows; row++) {
    for (int i = IA[row]; i < IA[row+1]; i++)
      b[row] += a[i] * x[JA[i]];
  }
}


/**
 * @brief Computes the Root Mean Square Error of two Vectors.
 * @details This function takes in two vectors (X and Y) and computes the
 *          Root Mean Square Error of the Vector Y with respect to Vector X.
 *          A boolean (integrated) must also be given to indicate whether the
 *          operation the vector should be group-wise integrated before
 *          performing the RMSE operation.
 * @param X a Vector object
 * @param Y a second Vector object
 * @param integrated a boolean indicating whether to group-wise integrate.
 */
FP_PRECISION computeRMSE(Vector* X, Vector* Y, bool integrated) {

  /* Check for consistency of vector dimensions */
  if (X->getNumX() != Y->getNumX() || X->getNumY() != Y->getNumY() ||
      X->getNumGroups() != Y->getNumGroups())
    log_printf(ERROR, "Cannot compute RMSE with different vector dimensions: "
               "(%i, %i, %i) and (%i, %i, %i)",
               X->getNumX(), X->getNumY(), X->getNumGroups(),
               Y->getNumX(), Y->getNumY(), Y->getNumGroups());

  FP_PRECISION rmse;
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();

  if (integrated) {

    FP_PRECISION new_source, old_source;
    Vector residual(num_x, num_y, num_z);

    /* Compute the RMSE */
    #pragma omp parallel for private(new_source, old_source)
    for (int i = 0; i < num_x*num_y*num_z; i++) {
      new_source = 0.0;
      old_source = 0.0;
      for (int g = 0; g < num_groups; g++) {
        new_source += X->getValue(i, g);
        old_source += Y->getValue(i, g);
      }

      if (new_source != 0.0)
        residual.setValue(i, 0, pow((new_source - old_source) / new_source, 2));
    }

    rmse = sqrt(residual.getSum() / (num_x * num_y * num_z));
  }
  else{

    Vector residual(num_x, num_y, num_z, num_groups);

    /* Compute the RMSE */
    #pragma omp parallel for
    for (int i = 0; i < num_x*num_y*num_z; i++) {
      for (int g = 0; g < num_groups; g++) {
        if (X->getValue(i, g) != 0.0)
          residual.setValue(i, g, pow((X->getValue(i, g) - Y->getValue(i, g)) /
                                      X->getValue(i, g), 2));
      }
    }

    rmse = sqrt(residual.getSum() / (num_x * num_y * num_z * num_groups));
  }

  return rmse;
}
