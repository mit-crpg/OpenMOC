#include "linalg.h"
#include <fstream>
#include <fenv.h>


/**
 * @brief Solves a generalized eigenvalue problem using a power iteration method
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a tolerance used
 *          for both the power method and linear solve convergence (tol), and
 *          a successive over-relaxation factor (SOR_factor) and computes the
 *          dominant eigenvalue and eigenvector using the Power method. The
 *          eigenvalue is returned and the input X Vector is modified in
 *          place to be the corresponding eigenvector.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object
 * @param k_eff initial k_effective
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 * @param convergence_data a summary of how to solver converged
 * @param comm an MPI communicator for the domain-decomposed solver
 * @return k_eff the dominant eigenvalue
 */
double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double k_eff,
                             double tol, double SOR_factor,
                             ConvergenceData* convergence_data,
                             DomainCommunicator* comm) {

  log_printf(INFO, "Computing the Matrix-Vector eigenvalue...");
  tol = std::max(MIN_LINALG_TOLERANCE, tol);

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != M->getNumX() || A->getNumX() != X->getNumX())
    log_printf(ERROR, "Cannot compute the Matrix-Vector eigenvalue with "
               "different x dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumX(), M->getNumX(), X->getNumX());
  else if (A->getNumY() != M->getNumY() || A->getNumY() != X->getNumY())
    log_printf(ERROR, "Cannot compute the Matrix-Vector eigenvalue with "
               "different y dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumY(), M->getNumY(), X->getNumY());
  else if (A->getNumZ() != M->getNumZ() || A->getNumZ() != X->getNumZ())
    log_printf(ERROR, "Cannot compute the Matrix-Vector eigenvalue with "
               "different z dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumZ(), M->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != M->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log_printf(ERROR, "Cannot compute the Matrix-Vector eigenvalue with "
               "different num groups  for the A matrix, M matrix, and X vector:"
               " (%d, %d, %d)", A->getNumGroups(), M->getNumGroups(),
               X->getNumGroups());

  /* Initialize variables */
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups);
  double residual;
  int iter;

  /* Compute and normalize the initial source */
  matrixMultiplication(M, X, &old_source);
  double old_source_sum = old_source.getSum();
  int num_rows = X->getNumRows();
#ifdef MPIx
  if (comm != NULL) {
    double temp_sum = old_source_sum;
    MPI_Allreduce(&temp_sum, &old_source_sum, 1, MPI_DOUBLE, MPI_SUM,
        comm->_MPI_cart);
    int temp_rows = num_rows;
    MPI_Allreduce(&temp_rows, &num_rows, 1, MPI_INT, MPI_SUM,
        comm->_MPI_cart);
  }
#endif
  old_source.scaleByValue(num_rows / old_source_sum);
  X->scaleByValue(num_rows * k_eff / old_source_sum);

  /* Power iteration Matrix-Vector solver */
  double initial_residual = 0;
  bool solver_failure = false;
  for (iter = 0; iter < MAX_LINALG_POWER_ITERATIONS; iter++) {

    /* Solve X = A^-1 * old_source */
    bool converged = false;
    if (!solver_failure)
      converged = linearSolve(A, M, X, &old_source, tol*1e-1, SOR_factor,
                              convergence_data, comm);

    /* If the solver failed, try the diagonally dominant solver */
    if (!converged) {
      solver_failure = true;
      converged = ddLinearSolve(A, M, X, &old_source, tol*1e-1, SOR_factor,
                                convergence_data, comm);
    }

    /* Check for divergence */
    if (!converged)
      return -1.0;

    /* Compute the new source */
    matrixMultiplication(M, X, &new_source);

    /* Compute the sum of new sources */
    double new_source_sum = new_source.getSum();
#ifdef MPIx
    if (comm != NULL) {
      double temp_sum = new_source_sum;
      MPI_Allreduce(&temp_sum, &new_source_sum, 1, MPI_DOUBLE, MPI_SUM,
          comm->_MPI_cart);
    }
#endif

    /* Compute and set keff */
    k_eff = new_source_sum / num_rows;

    /* Scale the new source by keff */
    new_source.scaleByValue(1.0 / k_eff);

    /* Compute the residual */
    residual = computeRMSE(&new_source, &old_source, true, comm);
    if (iter == 0) {
      initial_residual = residual;
      if (initial_residual < 1e-14)
        initial_residual = 1e-10;
      if (convergence_data != NULL) {
        convergence_data->cmfd_res_1 = residual;
        convergence_data->linear_iters_1 = convergence_data->linear_iters_end;
        convergence_data->linear_res_1 = convergence_data->linear_res_end;
      }
    }

    /* Copy the new source to the old source */
    new_source.copyTo(&old_source);

    log_printf(INFO_ONCE, "Matrix-Vector eigenvalue iter: %d, keff: %f, residual: "
               "%3.2e", iter, k_eff, residual);

    /* Check for convergence */
    if ((residual / initial_residual < 0.03 || residual < MIN_LINALG_TOLERANCE)
        && iter > MIN_LINALG_POWER_ITERATIONS) {
      if (convergence_data != NULL) {
        convergence_data->cmfd_res_end = residual;
        convergence_data->cmfd_iters = iter;
      }
      break;
    }
  }

  log_printf(INFO_ONCE, "Matrix-Vector eigenvalue solve iterations: %d", iter);
  if (iter == MAX_LINALG_POWER_ITERATIONS)
    log_printf(ERROR, "Eigenvalue solve failed to converge in %d iterations",
               iter);

  return k_eff;
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
 * @param convergence_data a summary of how to solver converged
 * @param comm an MPI communicator for the domain-decomposed solver
 */
bool linearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                 double SOR_factor, ConvergenceData* convergence_data,
                 DomainCommunicator* comm) {

  bool success = true;
  tol = std::max(MIN_LINALG_TOLERANCE, tol);

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != B->getNumX() || A->getNumX() != X->getNumX() ||
      A->getNumX() != M->getNumX())
    log_printf(ERROR, "Cannot perform linear solve with different x dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumX(), M->getNumX(),
               B->getNumX(), X->getNumX());
  else if (A->getNumY() != B->getNumY() || A->getNumY() != X->getNumY() ||
           A->getNumY() != M->getNumY())
    log_printf(ERROR, "Cannot perform linear solve with different y dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumY(), M->getNumY(),
               B->getNumY(), X->getNumY());
  else if (A->getNumZ() != B->getNumZ() || A->getNumZ() != X->getNumZ() ||
           A->getNumZ() != M->getNumZ())
    log_printf(ERROR, "Cannot perform linear solve with different z dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumZ(), M->getNumZ(),
               B->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups() ||
           A->getNumGroups() != M->getNumGroups())
    log_printf(ERROR, "Cannot perform linear solve with different num groups"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumGroups(), M->getNumGroups(),
               B->getNumGroups(), X->getNumGroups());

  /* Initialize variables */
  double residual;
  double min_residual = 1e6;
  int iter = 0;
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  int num_rows = X->getNumRows();
  //Vector X_old(cell_locks, num_x, num_y, num_z, num_groups);  //FIXME delete
  //CMFD_PRECISION* x_old = X_old.getArray();
  int* IA = A->getIA();
  int* JA = A->getJA();
  CMFD_PRECISION* DIAG = A->getDiag();
  CMFD_PRECISION* a = A->getA();
  CMFD_PRECISION* x = X->getArray();
  CMFD_PRECISION* b = B->getArray();
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups);

  /* Compute initial source */
  matrixMultiplication(M, X, &old_source);

  // Initialize communication buffers
  int* coupling_sizes = NULL;
  int** coupling_indexes = NULL;
  CMFD_PRECISION** coupling_coeffs = NULL;
  CMFD_PRECISION** coupling_fluxes = NULL;

  double initial_residual = 0;
  while (iter < MAX_LINEAR_SOLVE_ITERATIONS) {

    /* Pass new flux to old flux */
    //X->copyTo(&X_old);

    // Iteration over red/black cells
    for (int color = 0; color < 2; color++) {
      int offset = 0;
#ifdef MPIx
      getCouplingTerms(comm, color, coupling_sizes, coupling_indexes,
                       coupling_coeffs, coupling_fluxes, x, offset);
#endif

#pragma omp parallel for collapse(2)
      for (int iz=0; iz < num_z; iz++) {
        for (int iy=0; iy < num_y; iy++) {
          for (int ix=(iy+iz+color+offset)%2; ix < num_x; ix+=2) {

            int cell = (iz*num_y + iy)*num_x + ix;
            int row_start = cell*num_groups;

            /* Find index into communicator buffers for cells on surfaces */
            bool on_surface = (iz==0) || (iz==num_z-1) || (iy==0) || (iy==num_y-1)
                 || (ix==0) || (ix==num_x-1);
            int domain_surface_index = -1;
            if (comm != NULL && on_surface)
              domain_surface_index = comm->mapLocalToSurface[cell];

            /* Contribution of off-diagonal terms, hard to SIMD vectorize */
            for (int g=0; g < num_groups; g++) {

              int row = row_start + g;
              x[row] = (1.0 - SOR_factor) * x[row] * (DIAG[row] / SOR_factor);

              if (fabs(DIAG[row]) < FLT_EPSILON )
                  log_printf(ERROR, "A zero has been found on the diagonal of "
                             "the CMFD matrix cell [%d,%d,%d]=%d, group %d",
                             ix, iy, iz, cell, g);

              for (int i = IA[row]; i < IA[row+1]; i++) {

                // Get the column index
                int col = JA[i];
                if (row != col)
                  x[row] -= a[i] * x[col];
                else
                  x[row] += b[row];
              }
#ifdef MPIx
              // Contribution of off node fluxes
              if (comm != NULL && on_surface) {
                int row_surf = domain_surface_index * num_groups + g;
                for (int i = 0; i < coupling_sizes[row_surf]; i++) {
                  int idx = coupling_indexes[row_surf][i] * num_groups + g;
                  int domain = comm->domains[color][row_surf][i];
                  CMFD_PRECISION flux = coupling_fluxes[domain][idx];
                  x[row] -= coupling_coeffs[row_surf][i] * flux;
                }
              }
#endif
              // Perform these operations separately, for performance
              x[row] *= (SOR_factor / DIAG[row]);
            }
          }
        }
      }
    }

    // Compute the new source
    matrixMultiplication(M, X, &new_source);

    // Compute the initial residual
    feclearexcept (FE_ALL_EXCEPT);
    if (iter == 0) {
      residual = computeRMSE(&new_source, &old_source, true, comm);
      if (convergence_data != NULL)
        convergence_data->linear_res_end = residual;
      initial_residual = residual;
    }

    // Increment the interations counter
    iter++;

    double ratio_residuals = 1;
    if (initial_residual > 0)
      ratio_residuals = residual / initial_residual;
    log_printf(DEBUG, "SOR iter: %d, residual: %3.2e, initial residual: %3.2e"
               ", ratio = %3.2e, tolerance: %3.2e, end? %d", iter, residual,
               initial_residual, ratio_residuals, tol,
               (ratio_residuals < 0.1 || residual < tol) &&
               iter > MIN_LINEAR_SOLVE_ITERATIONS);

    // Compute residual only after minimum iteration number
    if (iter > MIN_LINEAR_SOLVE_ITERATIONS) {

      residual = computeRMSE(&new_source, &old_source, true, comm);

      // Record current minimum residual
      if (residual < min_residual)
        min_residual = residual;

      // Check for going off the rails
      int raised = fetestexcept (FE_INVALID);
      if ((residual > 1e3 * min_residual && min_residual > 1e-10) || raised) {
        log_printf(WARNING, "linear solve divergent : res %e min_res %e NaN? %d",
                   residual, min_residual, raised);
        if (convergence_data != NULL)
          convergence_data->linear_iters_end = iter;
        success = false;
        break;
      }

      // Check for convergence
      if (residual / initial_residual < 0.1 || residual < tol) {
        if (convergence_data != NULL)
          convergence_data->linear_iters_end = iter;
        break;
      }
    }

    // Copy the new source to the old source
    if (iter > MIN_LINEAR_SOLVE_ITERATIONS - 1 &&
        iter < MAX_LINEAR_SOLVE_ITERATIONS)
      new_source.copyTo(&old_source);
  }

  log_printf(DEBUG, "linear solve iterations: %d", iter);

  // Check if the maximum iterations were reached
  if (iter == MAX_LINEAR_SOLVE_ITERATIONS) {
    matrixMultiplication(M, X, &new_source);
    double residual = computeRMSE(&new_source, &old_source, true, comm);
    log_printf(INFO, "Ratio = %3.2e, tol = %3.2e", residual / initial_residual,
               tol);
    log_printf(NORMAL, "Linear solve failed to converge in %d iterations with "
               "initial residual %3.2e and final residual %3.2e", iter,
               initial_residual, residual);
    success = false;
  }
  return success;
}


#ifdef MPIx
/**
 * @brief Get coupling fluxes and other information from neighbors.
 *        The information are transfered by reference.
 * @param comm Structure for communication of fluxes between neighbor domains
 * @param color red or black color
 * @param coupling_sizes Number of connecting neighbors for each surface cell
 * @param coupling_indexes Surface numbers of connecting neighbors for each
 *        surface cell
 * @param coupling_coeffs Coupling coeffs between connecting neighbors and
 *        itself for each surface cell
 * @param coupling_fluxes Fluxes of connecting neighbors for each surface cell
 * @param curr_fluxes CMFD cell fluxes of current iteration
 * @param offset Sum of the starting CMFD global indexes of a domain, for
 *        calculation of the color
 */
void getCouplingTerms(DomainCommunicator* comm, int color, int*& coupling_sizes,
                      int**& coupling_indexes, CMFD_PRECISION**& coupling_coeffs,
                      CMFD_PRECISION**& coupling_fluxes,
                      CMFD_PRECISION* curr_fluxes, int& offset) {
  if (comm != NULL) {
    coupling_sizes = comm->num_connections[color];
    coupling_indexes = comm->indexes[color];
    coupling_coeffs = comm->coupling_coeffs[color];
    coupling_fluxes = comm->fluxes[color];

    offset = comm->_offset;

    MPI_Request requests[2*NUM_FACES];

    int nx = comm->_local_num_x;
    int ny = comm->_local_num_y;
    int nz = comm->_local_num_z;

    int ng = comm->num_groups;

    // Get numerical precision for communication
    MPI_Datatype flux_type;
    if (sizeof(CMFD_PRECISION) == 4)
      flux_type = MPI_FLOAT;
    else
      flux_type = MPI_DOUBLE;

    int sizes[NUM_FACES];
    for (int coord=0; coord < 3; coord++) {
      for (int d=0; d<2; d++) {

        int dir = 2*d-1;
        int surf = coord + 3*d;
        int op_surf = surf - 3*dir;
        int source, dest;

        MPI_Cart_shift(comm->_MPI_cart, coord, dir, &source, &dest);

        // Pack MPI buffer
        int size = 0;
        if (surf == SURFACE_X_MIN) {
          size = ny * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < ny; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*ny+j)+g] =
                  curr_fluxes[ng*((i*ny + j)*nx) + g];
        }
        else if (surf == SURFACE_X_MAX) {
          size = ny * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < ny; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*ny+j)+g] =
                  curr_fluxes[ng*((i*ny + j)*nx + nx-1) + g];
        }
        else if (surf == SURFACE_Y_MIN) {
          size = nx * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx*ny + j) + g];
        }
        else if (surf == SURFACE_Y_MAX) {
          size = nx * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx*ny + j + nx*(ny-1)) + g];
        }
        else if (surf == SURFACE_Z_MIN) {
          size = nx * ny * ng;
          for (int i=0; i < ny; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx + j)+g];
        }
        else if (surf == SURFACE_Z_MAX) {
          size = nx * ny * ng;
          for (int i=0; i < ny; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx + j + nx*ny*(nz-1)) + g];
        }

        sizes[surf] = size;

        // Post send
        MPI_Isend(comm->buffer[surf], size, flux_type,
                  dest, 0, comm->_MPI_cart, &requests[2*surf]);

        // Post receive
        MPI_Irecv(&comm->buffer[op_surf][size], size, flux_type,
                  source, 0, comm->_MPI_cart, &requests[2*op_surf+1]);
      }
    }

    // Block for communication round to complete
    bool round_complete = false;
    while (!round_complete) {

      round_complete = true;
      int flag;
      MPI_Status send_stat;
      MPI_Status recv_stat;

      for (int coord=0; coord < 3; coord++) {
        for (int d=0; d<2; d++) {
          int surf = coord + 3*d;

          MPI_Test(&requests[2*surf], &flag, &send_stat);
          if (flag == 0)
            round_complete = false;

          MPI_Test(&requests[2*surf+1], &flag, &recv_stat);
          if (flag == 0)
            round_complete = false;
          else
            // Copy received data into coupling_fluxes
            for (int i=0; i < sizes[surf]; i++)
              coupling_fluxes[surf][i] = comm->buffer[surf][sizes[surf]+i];
        }
      }
    }
  }
}
#endif


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
               "(%d, %d, %d)", A->getNumX(), B->getNumX(), X->getNumX());
  else if (A->getNumY() != B->getNumY() || A->getNumY() != X->getNumY())
    log_printf(ERROR, "Cannot perform matrix multiplication with different y "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%d, %d, %d)", A->getNumY(), B->getNumY(), X->getNumY());
  else if (A->getNumZ() != B->getNumZ() || A->getNumZ() != X->getNumZ())
    log_printf(ERROR, "Cannot perform matrix multiplication with different z "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%d, %d, %d)", A->getNumZ(), B->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log_printf(ERROR, "Cannot perform matrix multiplication with different "
               "num groups for the A matrix, M matrix, and X vector: "
               "(%d, %d, %d)", A->getNumGroups(), B->getNumGroups(),
               X->getNumGroups());

  B->setAll(0.0);
  int* IA = A->getIA();
  int* JA = A->getJA();
  CMFD_PRECISION* a = A->getA();
  CMFD_PRECISION* x = X->getArray();
  CMFD_PRECISION* b = B->getArray();
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
 *          The boolean integrated must also be given to indicate whether the
 *          operation on the vector should be group-wise integrated before
 *          performing the RMSE operation.
 * @param X a Vector object
 * @param Y a second Vector object
 * @param integrated a boolean indicating whether to group-wise integrate.
 * @param comm an MPI communicator to exchange residuals between domains
 */
double computeRMSE(Vector* X, Vector* Y, bool integrated,
                   DomainCommunicator* comm) {

  /* Check for consistency of vector dimensions */
  if (X->getNumX() != Y->getNumX() || X->getNumY() != Y->getNumY() ||
      X->getNumZ() != Y->getNumZ() || X->getNumGroups() != Y->getNumGroups())
    log_printf(ERROR, "Cannot compute RMSE with different vector dimensions: "
               "(%d, %d, %d, %d) and (%d, %d, %d, %d)",
               X->getNumX(), X->getNumY(), X->getNumZ(), X->getNumGroups(),
               Y->getNumX(), Y->getNumY(), Y->getNumZ(), Y->getNumGroups());

  double rmse;
  double sum_residuals = 0;
  int norm;
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  omp_lock_t* cell_locks = X->getCellLocks();

  if (integrated) {

    double new_source, old_source;
    CMFD_PRECISION residual[num_x * num_y * num_z]
         __attribute__ ((aligned(VEC_ALIGNMENT)));  //FIXME Overflow for large cases?
    memset(residual, 0, num_x * num_y * num_z * sizeof(CMFD_PRECISION));

    /* Compute the RMSE */
#pragma omp parallel for private(new_source, old_source)
    for (int i = 0; i < num_x*num_y*num_z; i++) {
      new_source = 0.0;
      old_source = 0.0;
#pragma omp simd reduction(+:new_source,old_source)
      for (int g = 0; g < num_groups; g++) {
        new_source += X->getValue(i, g);
        old_source += Y->getValue(i, g);
      }
      if (fabs(old_source) > FLUX_EPSILON)
        residual[i] = pow((new_source - old_source) / old_source, 2);
    }

    // Sum residuals
#pragma omp simd reduction(+:sum_residuals) aligned(residual)
    for (int i = 0; i < num_x*num_y*num_z; i++)
      sum_residuals += residual[i];

    norm = num_x * num_y * num_z;
  }
  else {

    //FIXME Incurs a memory allocation, uses unnecessary locks
    Vector residual(cell_locks, num_x, num_y, num_z, num_groups);

    /* Compute the RMSE */
#pragma omp parallel for
    for (int i = 0; i < num_x*num_y*num_z; i++) {
      for (int g = 0; g < num_groups; g++) {
        if (fabs(X->getValue(i, g)) > FLUX_EPSILON)
          residual.setValue(i, g, pow((X->getValue(i, g) - Y->getValue(i, g)) /
                                      X->getValue(i, g), 2));
      }
    }
    sum_residuals = residual.getSum();
    norm = num_x * num_y * num_z * num_groups;
  }

#ifdef MPIx
  if (comm != NULL) {
    double temp_residual = sum_residuals;
    MPI_Allreduce(&temp_residual, &sum_residuals, 1, MPI_DOUBLE, MPI_SUM,
        comm->_MPI_cart);
    int temp_norm = norm;
    MPI_Allreduce(&temp_norm, &norm, 1, MPI_INT, MPI_SUM, comm->_MPI_cart);
  }
#endif

  /* Error check residual components */
  if (sum_residuals < 0.0) {
    log_printf(WARNING, "CMFD Residual mean square error %6.4f less than zero",
               sum_residuals);
    sum_residuals = 0.0;
  }
  if (norm <= 0) {
    log_printf(WARNING, "CMFD residual norm %d less than one", norm);
    norm = 1;
  }

  /* Compute RMS residual error */
  rmse = sqrt(sum_residuals / norm);

#ifdef MPIx
  if (comm != NULL) {
    double temp_rmse = rmse;
    MPI_Allreduce(&temp_rmse, &rmse, 1, MPI_DOUBLE, MPI_MAX, comm->_MPI_cart);
  }
#endif

  return rmse;
}


/**
 * @brief Solves a linear system using the linear solver above, but makes the
 *        loss and streaming matrix diagonally dominant first, to increase
 *        likelihood of convergence.
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a source Vector (B),
 *          a source convergence tolerance (tol) and a successive
 *          over-relaxation factor (SOR_factor) and makes (A) diagonally
 *          dominant before calling the linear solve routine to compute the
 *          solution to the linear system. The input X Vector is modified in
 *          place to be the solution vector. The transformation to make (A)
 *          diagonally dominant is compensated by another matrix multiplication.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object
 * @param B the source Vector object
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 * @param convergence_data a summary of the convergence performance
 * @param comm a communicator for exchanging data through MPI
 */
bool ddLinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                   double SOR_factor, ConvergenceData* convergence_data,
                   DomainCommunicator* comm) {

  /* Create vector for stabilizing flux */
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  int num_rows = X->getNumRows();

  Vector dd(cell_locks, num_x, num_y, num_z, num_groups);
  dd.setAll(0.0);

  CMFD_PRECISION* dd_array = dd.getArray();
  CMFD_PRECISION* x = X->getArray();

  /* Stabilize matrix A to be diagonally dominant */
  CMFD_PRECISION* a = A->getA();
  CMFD_PRECISION* a_diag = A->getDiag();
  int* IA = A->getIA();
  int* JA = A->getJA();

  // Initialize communication buffers
  int* coupling_sizes = NULL;
  int** coupling_indexes = NULL;
  CMFD_PRECISION** coupling_coeffs = NULL;
  CMFD_PRECISION** coupling_fluxes = NULL;

  /* Loop over cells */
  for (int color = 0; color < 2; color++) {
    int offset = 0;
#ifdef MPIx
    getCouplingTerms(comm, color, coupling_sizes, coupling_indexes,
                     coupling_coeffs, coupling_fluxes, x, offset);
#endif
#pragma omp parallel for collapse(2)
    for (int iz=0; iz < num_z; iz++) {
      for (int iy=0; iy < num_y; iy++) {
        for (int ix=(iy+iz+color+offset)%2; ix < num_x; ix+=2) {

          int cell = (iz*num_y + iy)*num_x + ix;

          /* Find index into communicator buffers for cells on surfaces */
          bool on_surface = (iz==0) || (iz==num_z-1) || (iy==0) || (iy==num_y-1)
               || (ix==0) || (ix==num_x-1);
          int domain_surface_index = -1;
          if (comm != NULL && on_surface)
            domain_surface_index = comm->mapLocalToSurface[cell];

          /* Determine whether each group's row is diagonally dominant */
          for (int e = 0; e < num_groups; e++) {

            int row = cell * num_groups + e;

            /* Add local off-diagonal elements */
            double row_sum = 0.0;
            int diag_ind = -1;
            for (int idx = IA[row]; idx < IA[row+1]; idx++) {
              if (JA[idx] != row)
                row_sum += fabs(a[idx]);
              else
                diag_ind = idx;
            }

            /* Add off-node off-diagonal elements */
#ifdef MPIx
            if (comm != NULL && on_surface) {
              int row_surf = domain_surface_index * num_groups + e;
              int* coupling_sizes = comm->num_connections[color];
              CMFD_PRECISION** coupling_coeffs = comm->coupling_coeffs[color];
              for (int idx=0; idx < coupling_sizes[row_surf]; idx++)
                row_sum += fabs(coupling_coeffs[row_surf][idx]);
            }
#endif

            /* Check for diagonal dominance */
            if (row_sum > a[diag_ind])
              dd.incrementValue(cell, e, row_sum - a[diag_ind]);
          }
        }
      }
    }
  }

  /* Adjust matrix A to be diagonally dominant */
#pragma omp parallel for
  for (int i=0; i < num_x*num_y*num_z; i++) {
    for (int e=0; e < num_groups; e++) {
      A->incrementValue(i, e, i, e, dd.getValue(i,e));
    }
  }

  /* Create a vector for the remainder right hand side */
  Vector RHS(cell_locks, num_x, num_y, num_z, num_groups);
  CMFD_PRECISION* rhs_array = RHS.getArray();
  CMFD_PRECISION* b = B->getArray();

  /* Keep track of sources */
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups);

  /* Calculate source */
  matrixMultiplication(M, X, &new_source);

  /* Iterate to get the total solution */
  double initial_residual = 0.0;
  double residual = 0.0;
  double min_residual = 1e6;
  for (int iter=0; iter < MAX_LINEAR_SOLVE_ITERATIONS; iter++) {

    // Copy the new source to the old source
    new_source.copyTo(&old_source);

    for (int row=0; row < num_rows; row++)
      rhs_array[row] = b[row] + dd_array[row] * x[row];

    bool converged = linearSolve(A, M, X, &RHS, tol, SOR_factor,
                                 convergence_data, comm);
    if (!converged)
      log_printf(ERROR, "Stabilized linear solver inner iteration failed"
                 " to converge");

    // Compute the new source
    matrixMultiplication(M, X, &new_source);

    // Compute the residual
    residual = computeRMSE(&new_source, &old_source, true, comm);
    if (iter == 0){
      initial_residual = residual;
      if (initial_residual < 1e-14)
        initial_residual = 1e-10;
    }

    // Record current minimum residual
    if (residual < min_residual)
      min_residual = residual;

    // Check for going off the rails
    if (residual > 1e3 * min_residual && min_residual > 1e-10) {
      log_printf(WARNING, "Inner linear solve divergent.");
      log_printf(NORMAL, "Residual = %6.4e, Min Res = %6.4e", residual, min_residual);
      return false;
    }

    // Check for convergence
    if (residual / initial_residual < 0.1 || residual < tol)
      break;
  }

  /* Reset matrix A */
#pragma omp parallel for
  for (int i=0; i < num_x*num_y*num_z; i++) {
    for (int e=0; e < num_groups; e++) {
      A->incrementValue(i, e, i, e, -dd.getValue(i,e));
    }
  }

  return true;
}
