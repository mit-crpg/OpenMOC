#include "linalg.h"
#include <fstream>

/**
 * @brief Solves a generalized eigenvalue problem using the Power method.
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
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 * @return k_eff the dominant eigenvalue
 */
FP_PRECISION eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double k_eff,
                             FP_PRECISION tol, FP_PRECISION SOR_factor,
                             ConvergenceData* convergence_data,
                             DomainCommunicator* comm) {

  log_printf(INFO, "Computing the Matrix-Vector eigenvalue...");

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
  int num_rows = X->getNumRows();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups);
  FP_PRECISION residual;
  int iter;

  /* Compute and normalize the initial source */
  matrixMultiplication(M, X, &old_source);
  double old_source_sum = old_source.getSum();
  old_source.scaleByValue(num_rows / old_source_sum);
  X->scaleByValue(num_rows * k_eff / old_source_sum);

  if (comm != NULL) {
    int rank;
    MPI_Comm_rank(comm->_MPI_cart, &rank);
    if (rank == 0) {
      std::cout << "BEFORE" << std::endl;
      std::cout << "flux:" << std::endl;
      X->printString();
      std::cout << "source:" << std::endl;
      old_source.printString();
    }
  }

  /* Power iteration Matrix-Vector solver */
  double initial_residual = 0;
  for (iter = 0; iter < MAX_LINALG_POWER_ITERATIONS; iter++) {

    //FIXME
    Matrix* A_temp = A;
    Matrix* M_temp = M;
    Vector* X_temp = X;
    Vector* old_source_temp = &old_source;
    if (comm != NULL) {

      int nx = num_x / comm->_num_domains_x;
      int ny = num_y / comm->_num_domains_y;
      int nz = num_z / comm->_num_domains_z;

      A_temp = new Matrix(A->getCellLocks(), nx, ny, nz, A->getNumGroups());
      M_temp = new Matrix(M->getCellLocks(), nx, ny, nz, M->getNumGroups());
      X_temp = new Vector(X->getCellLocks(), nx, ny, nz, X->getNumGroups());
      old_source_temp = new Vector(old_source.getCellLocks(), nx, ny, nz,
                          old_source.getNumGroups());

      int ng = A->getNumGroups();
      int num_small_cells = nx * ny * nz;

      int x_start = comm->_domain_idx_x * nx;
      int x_end = (comm->_domain_idx_x + 1) * nx;
      int y_start = comm->_domain_idx_y * ny;
      int y_end = (comm->_domain_idx_y + 1) * ny;
      int z_start = comm->_domain_idx_z * nz;
      int z_end = (comm->_domain_idx_z + 1) * nz;

      comm->num_connections = new int*[2];
      comm->indexes = new int**[2];
      comm->domains = new int**[2];
      comm->fluxes = new FP_PRECISION**[2];
      comm->coupling_coeffs = new FP_PRECISION**[2];
      comm->buffer = new FP_PRECISION[2*num_small_cells*ng*NUM_FACES];
      comm->buffer_size = num_small_cells*ng;
      for (int rb=0; rb<2; rb++) {
        comm->num_connections[rb] = new int[num_small_cells*ng];
        comm->indexes[rb] = new int*[num_small_cells*ng];
        comm->domains[rb] = new int*[num_small_cells*ng];
        comm->fluxes[rb] = new FP_PRECISION*[num_small_cells*ng];
        comm->coupling_coeffs[rb] = new FP_PRECISION*[num_small_cells*ng];
        //printf("At %d allocating %d\n", rb, num_small_cells * ng);
        for (int nsc=0; nsc < num_small_cells * ng; nsc++) {
          comm->num_connections[rb][nsc] = 0;
          comm->indexes[rb][nsc] = new int[NUM_FACES];
          comm->domains[rb][nsc] = new int[NUM_FACES];
          comm->fluxes[rb][nsc] = new FP_PRECISION[NUM_FACES];
          comm->coupling_coeffs[rb][nsc] = new FP_PRECISION[NUM_FACES];
          for (int f=0; f < NUM_FACES; f++)
            comm->fluxes[rb][nsc][f] = X->getValue(nsc/ng, nsc%ng);
        }
      }

      comm->_offset = (x_start + y_start + z_start) % 2;

      Matrix* new_mats[2] = {A_temp, M_temp};
      Matrix* old_mats[2] = {A, M};

      Vector* new_vecs[2] = {X_temp, old_source_temp};
      Vector* old_vecs[2] = {X, &old_source};

      for (int mi=0; mi < 2; mi++) {

        Matrix* new_mat = new_mats[mi];
        Matrix* old_mat = old_mats[mi];

        Vector* new_vec = new_vecs[mi];
        Vector* old_vec = old_vecs[mi];

        int* IA = old_mat->getIA();
        int* JA = old_mat->getJA();

        for (int ix=x_start; ix < x_end; ix++) {
          for (int iy=y_start; iy < y_end; iy++) {
            for (int iz=z_start; iz < z_end; iz++) {
              int n = ((iz-z_start)*ny + (iy-y_start))*nx + ix - x_start;
              int n_global = (iz*num_y + iy)*num_x + ix;
              for (int g=0; g < ng; g++) {
                int row_global = n_global * ng + g;
                int row = n * ng + g;

                new_vec->setValue(n,g,old_vec->getValue(n_global,g));
                for (int i = IA[row_global]; i < IA[row_global+1]; i++) {

                  // Get the column index
                  int col = JA[i];
                  int j = col / ng;
                  int gp = col % ng;

                  int jx = (j % (num_x * num_y)) % num_x;
                  int jy = (j % (num_x * num_y)) / num_x;
                  int jz = j / (num_x * num_y);

                  int domain = -1;
                  //FIXME: missing corners
                  if (jx < x_start)
                    domain = SURFACE_X_MIN;
                  else if (jx >= x_end)
                    domain = SURFACE_X_MAX;
                  else if (jy < y_start)
                    domain = SURFACE_Y_MIN;
                  else if (jy >= y_end)
                    domain = SURFACE_Y_MAX;
                  else if (jz < z_start)
                    domain = SURFACE_Z_MIN;
                  else if (jz >= z_end)
                    domain = SURFACE_Z_MAX;

                  FP_PRECISION val = old_mat->getValue(j, gp, n_global, g);
                  int jp = ((jz-z_start)*ny + (jy-y_start))*nx + jx - x_start;
                  if (domain == -1) {
                    new_mat->setValue(jp, gp, n, g, val);
                  }
                  else {

                    /*
                    printf("Found off node at %d global (%d, %d, %d) local "
                           "(%d, %d, %d) group %d.\n From global index "
                           "(%d, %d, %d) group %d.\n", row, jx, jy, jz, jx-x_start,
                           jy-y_start, jz-z_start, gp, ix, iy, iz, g);
                    */

                    int color = (ix + iy + iz) % 2;
                    int idx = comm->num_connections[color][row];

                    comm->domains[color][row][idx] = domain;

                    if (domain == SURFACE_X_MIN)
                      comm->indexes[color][row][idx] = n + nx-1;
                    else if (domain == SURFACE_X_MAX)
                      comm->indexes[color][row][idx] = n - nx+1;
                    else if (domain == SURFACE_Y_MIN)
                      comm->indexes[color][row][idx] = n + nx*(ny-1);
                    else if (domain == SURFACE_Y_MAX)
                      comm->indexes[color][row][idx] = n - nx*(ny-1);
                    else if (domain == SURFACE_Z_MIN)
                      comm->indexes[color][row][idx] = n + nx * ny * (nz-1);
                    else if (domain == SURFACE_Z_MAX)
                      comm->indexes[color][row][idx] = n - nx * ny * (nz-1);

                    if (idx >= NUM_FACES) {
                      std::cout << "OUT OF BOUNDS AT " << row << std::endl;
                      std::cout << "-------> " << ix << ", " << iy << ", " <<
                        iz << " (" << g << ", " << gp << ")" << std::endl;
                    }

                    comm->coupling_coeffs[color][row][idx] = val;

                    comm->num_connections[color][row]++;
                  }
                }
              }
            }
          }
        }
      }
    }

    /*
    A->printString();
    A_temp->printString();
    for (int c=0; c<2; c++) {
      for (int r=0; r < A_temp->getNumRows(); r++) {
        for (int idx=0; idx < comm->num_connections[c][r]; idx++) {
          printf("CONNECTION (%d, %d) -> %6.4f at domain %d\n", c, r,
                 comm->coupling_coeffs[c][r][idx], comm->domains[c][r][idx]);
        }
      }
    }
    */


    /* Solve X = A^-1 * old_source */
    linearSolve(A_temp, M_temp, X_temp, old_source_temp, tol*1e-2, SOR_factor,
                convergence_data, comm);

    //FIXME
    if (comm != NULL) {
      int nx = num_x / comm->_num_domains_x;
      int ny = num_y / comm->_num_domains_y;
      int nz = num_z / comm->_num_domains_z;
      int ng = A->getNumGroups();
      int num_small_cells = nx * ny * nz;

      int x_start = comm->_domain_idx_x * nx;
      int x_end = (comm->_domain_idx_x + 1) * nx;
      int y_start = comm->_domain_idx_y * ny;
      int y_end = (comm->_domain_idx_y + 1) * ny;
      int z_start = comm->_domain_idx_z * nz;
      int z_end = (comm->_domain_idx_z + 1) * nz;


      Vector* new_vecs[2] = {X_temp, old_source_temp};
      Vector* old_vecs[2] = {X, &old_source};
      for (int vi=0; vi<2; vi++) {
        Vector* new_vec = new_vecs[vi];
        Vector* old_vec = old_vecs[vi];
        //old_vec->setAll(0.0);
#ifdef MPIx
        MPI_Datatype flux_type;
        if (sizeof(FP_PRECISION) == 4)
          flux_type = MPI_FLOAT;
        else
          flux_type = MPI_DOUBLE;


        FP_PRECISION* send_arr = new FP_PRECISION[num_x*num_y*num_z*ng];
        for (int ii=0; ii<num_x*num_y*num_z*ng; ii++)
          send_arr[ii] = 0.0;
        for (int ix=x_start; ix < x_end; ix++) {
          for (int iy=y_start; iy < y_end; iy++) {
            for (int iz=z_start; iz < z_end; iz++) {
              int n = ((iz-z_start)*ny + (iy-y_start))*nx + ix - x_start;
              int n_global = (iz*num_y + iy)*num_x + ix;
              for (int g=0; g < ng; g++) {
                old_vec->setValue(n_global,g,new_vec->getValue(n,g));
                send_arr[n_global*ng+g] = new_vec->getValue(n,g);
              }
            }
          }
        }
        if (comm->_num_domains_x * comm->_num_domains_y * comm->_num_domains_z > 1)
          MPI_Allreduce(send_arr, old_vec->getArray(), num_x * num_y * num_z * ng,
                        flux_type, MPI_SUM, comm->_MPI_cart);
        delete [] send_arr;
#endif
      }

      for (int rb=0; rb<2; rb++) {
        //printf("At %d deleting %d\n", rb, num_small_cells * ng);
        for (int nsc=0; nsc < num_small_cells * ng; nsc++) {
          delete [] comm->indexes[rb][nsc];
          delete [] comm->domains[rb][nsc];
          delete [] comm->coupling_coeffs[rb][nsc];
          delete [] comm->fluxes[rb][nsc];
        }
        delete [] comm->num_connections[rb];
        delete [] comm->indexes[rb];
        delete [] comm->domains[rb];
        delete [] comm->coupling_coeffs[rb];
        delete [] comm->fluxes[rb];
      }

      delete [] comm->num_connections;
      delete [] comm->indexes;
      delete [] comm->fluxes;
      delete [] comm->coupling_coeffs;

      delete A_temp;
      delete M_temp;
      delete X_temp;
      delete old_source_temp;
    }

    if (comm->stop) {
      X->printString();
      old_source.printString();
      MPI_Barrier(comm->_MPI_cart);
      //exit(1);
    }
    //exit(1);
    /*
    if (comm != NULL) {
      int rank;
      MPI_Comm_rank(comm->_MPI_cart, &rank);
      if (rank == 0)
        X->printString();
    }
    */


    /* Compute the new source */
    matrixMultiplication(M, X, &new_source);

    /* Compute and set keff */
    k_eff = new_source.getSum() / num_rows;

    /* Scale the new source by keff */
    new_source.scaleByValue(1.0 / k_eff);

    /* Compute the residual */
    residual = computeRMSE(&new_source, &old_source, true, iter);
    if (iter == 0) {
      initial_residual = residual;
      if (convergence_data != NULL) {
        convergence_data->cmfd_res_1 = residual;
        convergence_data->linear_iters_1 = convergence_data->linear_iters_end;
        convergence_data->linear_res_1 = convergence_data->linear_res_end;
      }
    }

    /* Copy the new source to the old source */
    new_source.copyTo(&old_source);

    log_printf(INFO, "Matrix-Vector eigenvalue iter: %d, keff: %f, residual: "
               "%3.2e", iter, k_eff, residual);

    /* Check for convergence */
    //FIXME
    if ((residual < tol || residual / initial_residual < 0.01)
           && iter > MIN_LINALG_POWER_ITERATIONS) {
      if (convergence_data != NULL) {
        convergence_data->cmfd_res_end = residual;
        convergence_data->cmfd_iters = iter;
      }
      break;
    }
  }

  log_printf(INFO, "Matrix-Vector eigenvalue solve iterations: %d", iter);
  if (iter == MAX_LINALG_POWER_ITERATIONS)
    log_printf(WARNING, "Eigenvalue solve failed to converge in %d iterations",
               iter);

  //FIXME
  /*
  if (comm != NULL) {
    int rank;
    MPI_Comm_rank(comm->_MPI_cart, &rank);
    if (rank == 0) {
      std::cout << "After" << std::endl;
      std::cout << "flux:" << std::endl;
      X->printString();
      std::cout << "source:" << std::endl;
      old_source.printString();
    }
  }
  */
  //exit(1);

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
 */
void linearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, FP_PRECISION tol,
                 FP_PRECISION SOR_factor, ConvergenceData* convergence_data,
                 DomainCommunicator* comm) {

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
  FP_PRECISION residual;
  int iter = 0;
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  int num_rows = X->getNumRows();
  Vector X_old(cell_locks, num_x, num_y, num_z, num_groups);
  FP_PRECISION* x_old = X_old.getArray();
  int* IA = A->getIA();
  int* JA = A->getJA();
  FP_PRECISION* DIAG = A->getDiag();
  FP_PRECISION* a = A->getA();
  FP_PRECISION* x = X->getArray();
  FP_PRECISION* b = B->getArray();
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups);

  /* Compute initial source */
  matrixMultiplication(M, X, &old_source);

  //
  int* coupling_sizes = NULL;
  int** coupling_indexes = NULL;
  FP_PRECISION** coupling_coeffs = NULL;
  FP_PRECISION** coupling_fluxes = NULL;

  double initial_residual = 0;
  while (iter < MAX_LINEAR_SOLVE_ITERATIONS) {

    /* Pass new flux to old flux */
    X->copyTo(&X_old);

    // Iteration over red/black cells
    for (int color = 0; color < 2; color++) {
      int offset = 0;
      getCouplingTerms(comm, color, coupling_sizes, coupling_indexes,
                       coupling_coeffs, coupling_fluxes, x, offset);
//#pragma omp parallel for collapse(2)
      /* std::cout << "NEW COLOR" << std::endl;

      for (int r=0; r < num_rows; r++) {
        for (int idx=0; idx < coupling_sizes[r]; idx++) {
          printf("In sweep connection (%d, %d) -> %6.4f at idx %d\n", color, r,
                 coupling_coeffs[r][idx], idx);
        }
      }
      */

      for (int iz=0; iz < num_z; iz++) {
        for (int iy=0; iy < num_y; iy++) {
          for (int ix=(iy+iz+color+offset)%2; ix < num_x; ix+=2) {

            int cell = (iz*num_y + iy)*num_x + ix;
            int row_start = cell*num_groups;

            for (int g=0; g < num_groups; g++) {

              int row = row_start + g;
              x[row] = (1.0 - SOR_factor) * x[row];
              for (int i = IA[row]; i < IA[row+1]; i++) {

                // Get the column index
                int col = JA[i];

                /*
                printf("At %d (%d, %d, %d) grp %d subtracting %6.4f x %6.4f / %6.4f\n",
                       row, ix, iy, iz, g, a[i], x[col], DIAG[row]);
                       */

                if (row == col)
                  x[row] += SOR_factor * b[row] / DIAG[row];
                else
                  x[row] -= SOR_factor * a[i] * x[col] / DIAG[row];
              }

              // Contribution of off node fluxes
              if (comm != NULL) {
                for (int i = 0; i < coupling_sizes[row]; i++) {
                  int idx = coupling_indexes[row][i] * num_groups + g;
                  int domain = comm->domains[color][row][i];
                  FP_PRECISION flux = coupling_fluxes[idx][domain];
                  x[row] -= SOR_factor * coupling_coeffs[row][i] * flux
                            / DIAG[row];
                  /*
                  printf("xAt (%d, %d, %d) grp %d subtracting %6.4f x %6.4f / "
                       "%6.4f\n INDEXES = (%d, %d)\n",
                       ix, iy, iz, g, coupling_coeffs[row][i], flux,
                       DIAG[row], idx, domain);
                  std::cout << "Fluxes:" << std::endl;
                  for (int ii=0; ii < M->getNumRows(); ii++)
                    std::cout << coupling_fluxes[ii][domain] << " ";
                  std::cout << std::endl;
                  */
                }
              }
            }
          }
        }
      }
    }

    // Compute the new source
    matrixMultiplication(M, X, &new_source);

    // Compute the residual
    residual = computeRMSE(&new_source, &old_source, true, 1); //FIXME
    if (iter == 0) {
      if (convergence_data != NULL)
        convergence_data->linear_res_end = residual;
      initial_residual = residual;
    }

    // Copy the new source to the old source
    new_source.copyTo(&old_source);

    // Increment the interations counter
    iter++;

    log_printf(INFO, "SOR iter: %d, residual: %f", iter, residual);

    //FIXME
    if ((residual < tol || residual / initial_residual < 0.1)
         && iter > MIN_LINEAR_SOLVE_ITERATIONS) {
      if (convergence_data != NULL)
        convergence_data->linear_iters_end = iter;
      break;
    }
  }

  log_printf(INFO, "linear solve iterations: %d", iter);

  // Check if the maximum iterations were reached
  if (iter == MAX_LINEAR_SOLVE_ITERATIONS) {
    log_printf(WARNING, "Linear solve failed to converge in %d iterations",
               iter);

    for (int i=0; i < num_x*num_y*num_z*num_groups; i++) {
      if (x[i] < 0.0)
        x[i] = 0.0;
    }
    X->scaleByValue(num_rows / X->getSum());
  }
}
/*
    for (int color = 0; color < 2; color++) {
      for (int oct = 0; oct < 8; oct++) {
///#pragma omp parallel for private(row, col) collapse(3)
        for (int cz = (oct / 4) * num_z/2; cz < (oct / 4 + 1) * num_z/2; cz++) {
          for (int cy = (oct % 4 / 2) * num_y/2;
              cy < (oct % 4 / 2 + 1) * num_y/2; cy++) {
            for (int cx = (oct % 4 % 2) * num_x/2;
                cx < (oct % 4 % 2 + 1) * num_x/2; cx++) {

              if (((cx % 2)+(cy % 2)+(cz % 2)) % 2 == color) {

                for (int g = 0; g < num_groups; g++) {

                  int row = ((cz*num_y + cy)*num_x + cx)*num_groups + g;


                  x[row] = (1.0 - SOR_factor) * x[row];
                  std::cout << row << std::endl;
                  printf("(%d, %d, %d) -> %d oct %d group %d row %d\n", cx, cy, cz,
                      color, oct, g, row);

                  for (int i = IA[row]; i < IA[row+1]; i++) {

                    int col = JA[i];

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
    }
    //exit(0);

    matrixMultiplication(M, X, &new_source);

    residual = computeRMSE(&new_source, &old_source, true, 1);
    if (iter == 0) {
      if (convergence_data != NULL)
        convergence_data->linear_res_end = residual;
      initial_residual = residual;
    }

    new_source.copyTo(&old_source);

    iter++;

    log_printf(INFO, "SOR iter: %d, residual: %f", iter, residual);

    //FIXME
    if ((residual < tol || residual / initial_residual < 0.1)
         && iter > MIN_LINEAR_SOLVE_ITERATIONS) {
      if (convergence_data != NULL)
        convergence_data->linear_iters_end = iter;
      break;
    }
  }

  log_printf(INFO, "linear solve iterations: %d", iter);

  if (iter == MAX_LINEAR_SOLVE_ITERATIONS) {
    log_printf(WARNING, "Linear solve failed to converge in %d iterations",
               iter);

    for (int i=0; i < num_x*num_y*num_z*num_groups; i++) {
      if (x[i] < 0.0)
        x[i] = 0.0;
    }
    X->scaleByValue(num_rows / X->getSum());
  }
}
*/


//FIXME
void getCouplingTerms(DomainCommunicator* comm, int color, int*& coupling_sizes,
                      int**& coupling_indexes, FP_PRECISION**& coupling_coeffs,
                      FP_PRECISION**& coupling_fluxes, FP_PRECISION* curr_fluxes,
                      int& offset) {

  if (comm != NULL) {
    coupling_sizes = comm->num_connections[color];
    coupling_indexes = comm->indexes[color];
    coupling_coeffs = comm->coupling_coeffs[color];
    coupling_fluxes = comm->fluxes[color];

    offset = comm->_offset;

    int buffer_size = comm->buffer_size;
    MPI_Request requests[2*NUM_FACES];

    for (int coord=0; coord < 3; coord++) {
      for (int d=0; d<2; d++) {

        MPI_Datatype flux_type;
        if (sizeof(FP_PRECISION) == 4)
          flux_type = MPI_FLOAT;
        else
          flux_type = MPI_DOUBLE;

        int dir = 2*d-1;
        int surf = coord + 3*d;
        int op_surf = surf - 3*dir;
        int source, dest;

        MPI_Cart_shift(comm->_MPI_cart, coord, dir, &source, &dest);

        // Pack MPI buffer
        //MPI_Barrier(comm->_MPI_cart);
        for (int i=0; i < buffer_size; i++) {
          comm->buffer[2*surf*buffer_size+i] = curr_fluxes[i];
        }

        // Post send
        /*
        std::cout << "Sending:" << std::endl;
        for (int i=0; i < buffer_size; i++) {
          std::cout << " " << coupling_fluxes[i][surf];
        }
        std::cout << std::endl;
        */
        MPI_Isend(&comm->buffer[2*surf*buffer_size], buffer_size, flux_type,
                  dest, 0, comm->_MPI_cart, &requests[2*surf]);

        // Post receive
        MPI_Irecv(&comm->buffer[(2*op_surf+1)*buffer_size], buffer_size, flux_type,
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
            for (int i=0; i < buffer_size; i++)
              coupling_fluxes[i][surf] = comm->buffer[(2*surf+1)*buffer_size+i];
        }
      }
    }
    /*
    for (int coord=0; coord < 3; coord++) {
      for (int d=0; d<2; d++) {
        int surf = coord + 3*d;
        MPI_Status stat;

         int flag;
          MPI_Test(&requests[2*surf+1], &flag, &stat);
          if (flag != 0) {
            std::cout << "Received:" << std::endl;
            for (int i=0; i < buffer_size; i++) {
              std::cout << " " << coupling_fluxes[i][surf];
            std::cout << std::endl;
            }
          }
      }
    }
    */
  }
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
 *          The boolean integrated must also be given to indicate whether the
 *          operation on the vector should be group-wise integrated before
 *          performing the RMSE operation.
 * @param X a Vector object
 * @param Y a second Vector object
 * @param integrated a boolean indicating whether to group-wise integrate.
 */
double computeRMSE(Vector* X, Vector* Y, bool integrated, int it,
                   DomainCommunicator* comm) {

  /* Check for consistency of vector dimensions */
  if (X->getNumX() != Y->getNumX() || X->getNumY() != Y->getNumY() ||
      X->getNumZ() != Y->getNumZ() || X->getNumGroups() != Y->getNumGroups())
    log_printf(ERROR, "Cannot compute RMSE with different vector dimensions: "
               "(%d, %d, %d, %d) and (%d, %d, %d, %d)",
               X->getNumX(), X->getNumY(), X->getNumZ(), X->getNumGroups(),
               Y->getNumX(), Y->getNumY(), Y->getNumZ(), Y->getNumGroups());


  double rmse;
  double sum_residuals;
  int norm;
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  omp_lock_t* cell_locks = X->getCellLocks();

  if (integrated) {

    FP_PRECISION new_source, old_source;
    Vector residual(cell_locks, num_x, num_y, num_z, 1);

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
        residual.setValue(i, 0, pow((new_source - old_source) / old_source, 2));
    }
    sum_residuals = residual.getSum();
    norm = num_x * num_y * num_z;

  }
  else {

    Vector residual(cell_locks, num_x, num_y, num_z, num_groups);

    /* Compute the RMSE */
    #pragma omp parallel for
    for (int i = 0; i < num_x*num_y*num_z; i++) {
      for (int g = 0; g < num_groups; g++) {
        if (X->getValue(i, g) != 0.0)
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
    double temp_norm = norm;
    MPI_Allreduce(&temp_norm, &norm, 1, MPI_INT, MPI_SUM, comm->_MPI_cart);
    rmse = sqrt(sum_residuals / norm);
    double temp_rmse = 0;
    MPI_Allreduce(&temp_rmse, &rmse, 1, MPI_DOUBLE, MPI_MAX, comm->_MPI_cart);
  }
#endif


  return rmse;
}
