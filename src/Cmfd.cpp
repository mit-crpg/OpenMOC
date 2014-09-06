#include "Cmfd.h"

/**
 * @brief Constructor initializes boundaries and variables that describe
 *          the Cmfd object.
 * @details The construcor initializes the many variables that describe
 *          the CMFD Mesh, the solve method, and flux type. If solve
 *          method is DIFFUSION, the fsr volumes, Materials, and fluxes
 *          are initialized.
 * @param geometry pointer to the Geometry
 * @param criteria convergence criteria on keff
 */
Cmfd::Cmfd() {

  /* Initialize Geometry and Mesh-related attribute */
  _quad = NULL;
  _SOR_factor = 1.0;

  /* Global variables used in solving CMFD problem */
  _source_convergence_threshold = 1E-7;
  _num_x = 1;
  _num_y = 1;
  _width = 0.;
  _height = 0.;
  _cell_width = 0.;
  _cell_height = 0.;
  _flux_update_on = true;
  _optically_thick = true;
  _SOR_factor = 1.0;
  _num_FSRs = 0;
  _relax_factor = 0.6;

  /* Energy group problem parameters */
  _num_moc_groups = 0;
  _num_cmfd_groups = 0;

  /* Set matrices and arrays to NULL */
  _A = NULL;
  _M = NULL;
  _flux_temp = NULL;
  _old_source = NULL;
  _new_source = NULL;
  _group_indices = NULL;
  _group_indices_map = NULL;

  /* Initialize boundaries to be reflective */
  _boundaries = new boundaryType[4];
  _boundaries[0] = REFLECTIVE;
  _boundaries[1] = REFLECTIVE;
  _boundaries[2] = REFLECTIVE;
  _boundaries[3] = REFLECTIVE;

}


/**
 * @brief Destructor deletes arrays of A and M row insertion arrays.
 */
Cmfd::~Cmfd() {

  /* Delete matrix and vector objects */

  if (_M != NULL){
    for (int i = 0; i < _num_x*_num_y; i++)
      delete [] _M[i];

    delete [] _M;
  }

  if (_A != NULL){
    for (int i = 0; i < _num_x*_num_y; i++)
      delete [] _A[i];

    delete [] _A;
  }

  if (_old_flux != NULL)
    delete [] _old_flux;

  if (_new_flux != NULL)
    delete [] _new_flux;

  if (_flux_temp != NULL)
    delete [] _flux_temp;

  if (_old_source != NULL)
    delete [] _old_source;

  if (_new_source != NULL)
    delete [] _new_source;
}


/**
 * @brief Set the number of Mesh cells in a row.
 * @param number of Mesh cells in a row
 */
void Cmfd::setNumX(int num_x){

  if (num_x < 1)
    log_printf(ERROR, "The number of lattice cells in the x direction "
               "must be > 0. Input value: %i", num_x);

  _num_x = num_x;
  if (_width != 0.)
    _cell_width = _width / _num_x;
}


/**
 * @brief Set the number of Mesh cells in a column
 * @param number of Mesh cells in a column
 */
void Cmfd::setNumY(int num_y){

  if (num_y < 1)
    log_printf(ERROR, "The number of lattice cells in the y direction "
               "must be > 0. Input value: %i", num_y);

  _num_y = num_y;
  if (_height != 0.)
    _cell_height = _height / _num_y;
}


/**
 * @brief Get the number of Mesh cells in a row.
 * @return number of Mesh cells in a row
 */
int Cmfd::getNumX(){
  return _num_x;
}


/**
 * @brief Get the number of Mesh cells in a column
 * @return number of Mesh cells in a column
 */
int Cmfd::getNumY(){
  return _num_y;
}


/**
 * @brief Set Mesh width.
 * @param width physical width of Mesh
 */
void Cmfd::setWidth(double width){
  _width = width;
  if (_num_x != 0)
    _cell_width = _width / _num_x;
}


/**
 * @brief Set Mesh height.
 * @param height physical height of Mesh
 */
void Cmfd::setHeight(double height){
  _height = height;
  if (_num_y != 0)
    _cell_height = _height / _num_y;
}


/**
 * @brief Create cross-sections and fluxes for each Cmfd cell by
 *        energy condensing and volume averaging cross sections from
 *        the MOC sweep.
 */
void Cmfd::computeXS(){

  log_printf(INFO, "Computing CMFD cross-sections...");

  /* Split corner currents to side surfaces */
  splitCorners();

  /* Initialize variables for FSR properties*/
  FP_PRECISION volume, flux, abs, tot, nu_fis, chi, dif_coef;
  FP_PRECISION* scat;

  /* Initialize tallies for each parameter */
  FP_PRECISION abs_tally, nu_fis_tally, dif_tally, rxn_tally;
  FP_PRECISION vol_tally, tot_tally, neut_prod_tally;
  FP_PRECISION scat_tally[_num_cmfd_groups];
  FP_PRECISION chi_tally[_num_cmfd_groups];

  /* Pointers to material objects */
  Material* fsr_material;
  Material* cell_material;

  /* Loop over cmfd cells */
  #pragma omp parallel for private(volume, flux, abs, tot, nu_fis, chi, \
    dif_coef, scat, abs_tally, nu_fis_tally, dif_tally, rxn_tally,  \
    vol_tally, tot_tally, scat_tally, fsr_material, cell_material, \
    neut_prod_tally, chi_tally)
  for (int i = 0; i < _num_x * _num_y; i++){

    cell_material = _materials[i];
    std::vector<int>::iterator iter;

    /* Loop over CMFD coarse energy groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Zero tallies for this group */
      abs_tally = 0;
      nu_fis_tally = 0;
      dif_tally = 0;
      rxn_tally = 0;
      vol_tally = 0;
      tot_tally = 0;
      neut_prod_tally = 0.0;

      /* Zero each group-to-group scattering tally */
      for (int g = 0; g < _num_cmfd_groups; g++){
        scat_tally[g] = 0;
        chi_tally[g] = 0.0;
      }

      /* Loop over FSRs in cmfd cell */
      for (iter = _cell_fsrs.at(i).begin();
        iter != _cell_fsrs.at(i).end(); ++iter){

        fsr_material = _FSR_materials[*iter];
        volume = _FSR_volumes[*iter];
        scat = fsr_material->getSigmaS();
        vol_tally += volume;

        /* Chi tallies */
        for (int b = 0; b < _num_cmfd_groups; b++){
          chi = 0.0;

          /* Compute the chi for group b */
          for (int h = _group_indices[b]; h < _group_indices[b+1]; h++)
            chi += fsr_material->getChi()[h];

          for (int h = 0; h < _num_moc_groups; h++){
            chi_tally[b] += chi * fsr_material->getNuSigmaF()[h] *
                _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
            neut_prod_tally += chi * fsr_material->getNuSigmaF()[h] *
                _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
          }
        }

        /* Loop over MOC energy groups within this CMFD coarse group */
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){

          /* Gets FSR volume, material, and cross sections */
          flux = _FSR_fluxes[(*iter)*_num_moc_groups+h];
          abs = fsr_material->getSigmaA()[h];
          tot = fsr_material->getSigmaT()[h];
          dif_coef = fsr_material->getDifCoef()[h];
          nu_fis = fsr_material->getNuSigmaF()[h];

          /* If Material has a diffusion coefficient, use it; otherwise
           * estimate diffusion coefficient with \f$ \frac{1}{3\Sigma_t} \f$ */
          if (fsr_material->getDifCoef()[h] > 1e-8)
            dif_tally += fsr_material->getDifCoef()[h] * flux * volume;
          else
            dif_tally += flux * volume / (3.0 * tot);

          /* Increment tallies for this group */
          abs_tally += abs * flux * volume;
          tot_tally += tot * flux * volume;
          nu_fis_tally += nu_fis * flux * volume;
          rxn_tally += flux * volume;

          /* Scattering tallies */
          for (int g = 0; g < _num_moc_groups; g++){
              scat_tally[getCmfdGroup(g)] +=
                  scat[g*_num_moc_groups+h] * flux * volume;
          }
        }
      }

      /* Set the Mesh cell properties with the tallies */
      _volumes[i] = vol_tally;
      cell_material->setSigmaAByGroup(abs_tally / rxn_tally, e+1);
      cell_material->setSigmaTByGroup(tot_tally / rxn_tally, e+1);
      cell_material->setNuSigmaFByGroup(nu_fis_tally / rxn_tally, e+1);
      cell_material->setDifCoefByGroup(dif_tally / rxn_tally, e+1);
      _old_flux[i*_num_cmfd_groups+e] = rxn_tally / vol_tally;

      /* Set chi */
      if (neut_prod_tally != 0.0)
        cell_material->setChiByGroup(chi_tally[e] / neut_prod_tally, e+1);
      else
        cell_material->setChiByGroup(0.0,e+1);

      log_printf(DEBUG, "cell: %i, group: %i, vol: %f, siga: %f, sigt: %f,"
                 " nu_sigf: %f, dif_coef: %f, flux: %f, chi: %f", i, e,
                 vol_tally, abs_tally / rxn_tally, tot_tally / rxn_tally,
                 nu_fis_tally / rxn_tally, dif_tally / rxn_tally,
                 rxn_tally / vol_tally, chi_tally[e] / (neut_prod_tally+1e-12));

      /* Set scattering xs */
      for (int g = 0; g < _num_cmfd_groups; g++){
        cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, e+1, g+1);
        log_printf(DEBUG, "scattering from %i to %i: %f", e, g,
                   scat_tally[g] / rxn_tally);
      }
    }
  }
}


/**
 * @brief Compute the diffusion coefficients:
 *          \f$ D \f$ - straight diffusion coefficient
 *          \f$ \hat{D} \f$ - surface diffusion coefficient
 *          \f$ \tilde{D} \f$ - surface diffusion coefficient correction factor
 *        for each mesh while ensuring neutron balance is achieved.
 */
void Cmfd::computeDs(){

  log_printf(INFO, "Computing CMFD diffusion coefficients...");

  FP_PRECISION d, d_next, d_hat, d_tilde;
  FP_PRECISION current, flux, flux_next, f, f_next;
  FP_PRECISION length, length_perpen, next_length_perpen;
  FP_PRECISION sense;
  int next_surface;
  int cell, cell_next;

  /* Loop over mesh cells in y direction */
  #pragma omp parallel for private(d, d_next, d_hat, d_tilde, current, flux, \
    flux_next, f, f_next, length, length_perpen, next_length_perpen, \
    sense, next_surface, cell, cell_next)
  for (int y = 0; y < _num_y; y++){

    /* Loop over Mesh cells in x direction */
    for (int x = 0; x < _num_x; x++){

      cell = y*_num_x+x;

      /* Loop over Surfaces in a cell */
      for (int surface = 0; surface < 4; surface++){

        /* Loop over groups */
        for (int e = 0; e < _num_cmfd_groups; e++){

          /* Get diffusivity and flux for Mesh cell */
          d = _materials[cell]->getDifCoef()[e];
          flux = _old_flux[cell*_num_cmfd_groups+e];
          cell_next = getCellNext(cell, surface);

          /* Set halfspace sense of the Surface */
          if (surface == 0 || surface == 3)
            sense = -1.0;
          else
            sense = 1.0;

          /* Set the length of this Surface and the perpendicular Surface */
          if (surface == 0 || surface== 2){
            length = _cell_height;
            length_perpen = _cell_width;
          }
          else if (surface == 1 || surface == 3){
            length = _cell_width;
            length_perpen = _cell_height;
          }

          /* Compute the optical thickness correction factor */
          f = computeDiffCorrect(d, length_perpen);

          /* If Surface is on a boundary, choose appropriate BCs */
          if (cell_next == -1){

            current = sense * _surface_currents[cell*_num_cmfd_groups*8 +
                                       surface*_num_cmfd_groups + e];

            /* REFLECTIVE BC */
            if (_boundaries[surface] == REFLECTIVE){

              /* Set D's */
              d_hat = 0.0;
              d_tilde = 0.0;
            }

            /* VACUUM BC */
            else if (_boundaries[surface] == VACUUM){

              /* Set D's */
              d_hat =  2 * d*f / length_perpen / (1 + 4 * d*f /
                       length_perpen);
              d_tilde = (sense * d_hat * flux - current / length) / flux;
             }
          }

          /* If Surface is an interface, use finite differencing */
          else{

            /* Set properties for cell next to Surface */
            if (surface == 0){
              next_length_perpen = _cell_width;
              next_surface = 2;
            }
            else if (surface == 1){
              next_length_perpen = _cell_height;
              next_surface = 3;
            }
            else if (surface == 2){
              next_length_perpen = _cell_width;
              next_surface = 0;
            }
            else if (surface == 3){
              next_length_perpen = _cell_height;
              next_surface = 1;
            }

            /* Set diffusion coefficient and flux for neighboring cell */
            d_next = _materials[cell_next]->getDifCoef()[e];
            flux_next = _old_flux[cell_next*_num_cmfd_groups + e];

            /* Get optical thickness correction term for meshCellNext */
            f_next = computeDiffCorrect(d_next, next_length_perpen);

            /* Compute d_hat */
            d_hat = 2.0 * d * f * d_next * f_next / (length_perpen
                    * d * f + next_length_perpen * d_next*f_next);

            /* Compute net current */
            current = sense * _surface_currents[cell*_num_cmfd_groups*8 +
                      surface*_num_cmfd_groups + e] - sense
                      * _surface_currents[cell_next*_num_cmfd_groups*8 +
                      next_surface*_num_cmfd_groups + e];

            /* Compute d_tilde */
            d_tilde = -(sense * d_hat * (flux_next - flux) +
                        current  / length) / (flux_next + flux);

            /* If the magnitude of d_tilde is greater than the magnitude of
             * d_hat, select new values d_tilde and d_hat to ensure the course
             * mesh equations are guaranteed to be diagonally dominant */
            if (fabs(d_tilde) > fabs(d_hat)){

              if (sense == -1){

                /* If d_tilde is positive */
                if (1 - fabs(d_tilde)/d_tilde < 1e-8){
                  d_hat   = - current/(2*flux*length);
                  d_tilde = - current/(2*flux*length);
                }

                /* If d_tilde is negative */
                else{
                  d_hat   = current/(2*flux_next*length);
                  d_tilde = - current/(2*flux_next*length);
                }
              }
              else{

                /* If d_tilde is positive */
                if (1 - fabs(d_tilde)/d_tilde < 1e-8){
                  d_hat   = - current/(2*flux_next*length);
                  d_tilde = - current/(2*flux_next*length);
                }

                /* If d_tilde is negative */
                else{
                  d_hat   = current/(2*flux*length);
                  d_tilde = - current/(2*flux*length);
                }
              }
            }
          }

          /* Perform underrelaxation on d_tilde */
          d_tilde =
             _materials[cell]->getDifTilde()[surface*_num_cmfd_groups + e] *
             (1 - _relax_factor) + _relax_factor * d_tilde;

          /* Set d_hat and d_tilde */
          _materials[cell]->setDifHatByGroup(d_hat, e+1, surface);
          _materials[cell]->setDifTildeByGroup(d_tilde, e+1, surface);

          log_printf(DEBUG, "cell: %i, group: %i, side: %i, flux: %f,"
                     " current: %f, d: %f, dhat: %f, dtilde: %f",
                     y*_num_x + x, e, surface, flux, current, d, d_hat,
                     d_tilde);
        }
      }
    }
  }
}



/** @brief CMFD solver that solves the diffusion problem.
 * @return k-effective the solution eigenvalue.
 */
FP_PRECISION Cmfd::computeKeff(){

  log_printf(INFO, "Running diffusion solver...");

  /* Create matrix and vector objects */
  if (_A == NULL){
    try{
    
      _M = new FP_PRECISION*[_num_x*_num_y];
      _A = new FP_PRECISION*[_num_x*_num_y];
      _old_source = new FP_PRECISION[_num_x*_num_y*_num_cmfd_groups];
      _new_source = new FP_PRECISION[_num_x*_num_y*_num_cmfd_groups];
      _volumes = new FP_PRECISION[_num_x*_num_y];

      initializeFlux();
      initializeMaterials();

      for (int i = 0; i < _num_x*_num_y; i++){
        _M[i] = new FP_PRECISION[_num_cmfd_groups*_num_cmfd_groups];
        _A[i] = new FP_PRECISION[_num_cmfd_groups*(_num_cmfd_groups+4)];
      }
    }
    catch(std::exception &e){
      log_printf(ERROR, "Could not allocate memory for the CMFD mesh objects. "
                 "Backtrace:%s", e.what());
    }
  }

  /* Initialize variables */
  FP_PRECISION sum_new, sum_old, val, residual, scale_val;
  int row;

  /* Convergence criteria on L2 norm of flux for linear solve */
  FP_PRECISION linear_solve_convergence_criteria = 1E-7;
  
  /* Compute the cross sections and surface diffusion coefficients */
  computeXS();
  computeDs();

  /* Construct matrices */
  constructMatrices();

  /* Compute and normalize the initial source */
  matrix_multiplication(_M, _old_flux, _old_source, _num_x*_num_y, 
      _num_cmfd_groups);
  sum_old = pairwise_sum(_old_source, _num_x*_num_y*_num_cmfd_groups);
  scale_val = (_num_x * _num_y * _num_cmfd_groups) / sum_old;
  vector_scale(_old_source, scale_val, _num_x*_num_y*_num_cmfd_groups);
  vector_copy(_old_flux, _new_flux, _num_x*_num_y*_num_cmfd_groups);
  vector_scale(_new_flux, scale_val, _num_x*_num_y*_num_cmfd_groups);
  sum_old = _num_x * _num_y * _num_cmfd_groups;
  
  /* Power iteration diffusion solver */
  for (int iter = 0; iter < 25000; iter++){
      
    /* Solve phi = A^-1 * old_source */
    linearSolve(_A, _new_flux, _old_source, linear_solve_convergence_criteria);
      
    /* Compute the new source */
    matrix_multiplication(_M, _new_flux, _new_source, _num_x*_num_y,
        _num_cmfd_groups);
    sum_new = pairwise_sum(_new_source, _num_x*_num_y*_num_cmfd_groups);
    
    /* Compute and set keff */
    _k_eff = sum_new / sum_old;
    
    /* Scale the old source by keff */
    vector_scale(_old_source, _k_eff, _num_x*_num_y*_num_cmfd_groups);
    
    /* Compute the L2 norm of source error */
    residual = 0.0;    
    for (int i = 0; i < _num_x*_num_y*_num_cmfd_groups; i++){
      if (_new_source[i] != 0.0)
        residual += pow((_new_source[i] - _old_source[i]) / _new_source[i], 2);
    }

    /* Compute the average value of the residual */
    residual = sqrt(residual / (_num_x*_num_y*_num_cmfd_groups));

    /* Normalize the new source to have an average value of 1.0 */
    scale_val = (_num_x * _num_y * _num_cmfd_groups) / sum_new;
    vector_scale(_new_source, scale_val, _num_x*_num_y*_num_cmfd_groups);
    vector_copy(_new_source, _old_source, _num_x*_num_y*_num_cmfd_groups);
    
    log_printf(INFO, "CMFD iter: %i, keff: %f, error: %f", 
               iter, _k_eff, residual);
    
    /* Check for convergence */
    if (residual < _source_convergence_threshold)
      break;
  }

  /* Rescale the old and new flux */
  rescaleFlux();

  /* Update the MOC flux */
  updateMOCFlux();

  return _k_eff;
}


/**
 * @brief Solve the linear system Ax=b using Gauss Seidel with SOR.
 * @param mat pointer to A matrix
 * @param vec_x pointer to x vector
 * @param vec_b pointer to b vector
 * @param conv flux convergence criteria
 * @param max_iter the maximum number of iterations
 */
void Cmfd::linearSolve(FP_PRECISION** mat, FP_PRECISION* vec_x, 
                       FP_PRECISION* vec_b, FP_PRECISION conv, int max_iter){

  FP_PRECISION residual = 1E10;
  int row, cell;
  FP_PRECISION val;
  int iter = 0;

  while (residual > conv){

    /* Pass new flux to old flux */
    vector_copy(vec_x, _flux_temp, _num_x*_num_y*_num_cmfd_groups);

    /* Iteration over red cells */
    #pragma omp parallel for private(row, val, cell)
    for (int y = 0; y < _num_y; y++){
      for (int x = y % 2; x < _num_x; x += 2){

        cell = y*_num_x+x;

        for (int g = 0; g < _num_cmfd_groups; g++){

          row = (y*_num_x+x)*_num_cmfd_groups + g;
          val = 0.0;

          /* Previous flux term */
          val += (1.0 - _SOR_factor) * vec_x[row];

          /* Source term */
          val += _SOR_factor*vec_b[row] / mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _SOR_factor * vec_x[row - _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Bottom surface */
          if (y != 0)
            val -= _SOR_factor * vec_x[row - _num_x * _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+1] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < _num_cmfd_groups; e++){
            if (e != g)
              val -= _SOR_factor * vec_x[(y*_num_x+x)*_num_cmfd_groups+e] *
                     mat[cell][g*(_num_cmfd_groups+4)+2+e] /
                     mat[cell][g*(_num_cmfd_groups+4)+g+2];
          }

          /* Right surface */
          if (x != _num_x - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+2] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Top surface */
          if (y != _num_y - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups*_num_x] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+3] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          vec_x[row] = val;
        }
      }
    }

    /* Iteration over black cells */
    #pragma omp parallel for private(row, val, cell)
    for (int y = 0; y < _num_y; y++){
      for (int x = 1 - y % 2; x < _num_x; x += 2){

        cell = y*_num_x+x;

        for (int g = 0; g < _num_cmfd_groups; g++){

          row = (y*_num_x+x)*_num_cmfd_groups + g;
          val = 0.0;

          /* Previous flux term */
          val += (1.0 - _SOR_factor) * vec_x[row];

          /* Source term */
          val += _SOR_factor*vec_b[row] / mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _SOR_factor * vec_x[row - _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Bottom surface */
          if (y != 0)
            val -= _SOR_factor * vec_x[row - _num_x * _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+1] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < _num_cmfd_groups; e++){
            if (e != g)
              val -= _SOR_factor * vec_x[(y*_num_x+x)*_num_cmfd_groups+e] *
                     mat[cell][g*(_num_cmfd_groups+4)+2+e] /
                     mat[cell][g*(_num_cmfd_groups+4)+g+2];
          }

          /* Right surface */
          if (x != _num_x - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+2] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Top surface */
          if (y != _num_y - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups*_num_x] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+3] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          vec_x[row] = val;
        }
      }
    }

    /* Compute the average residual */
    residual = 0.0;
    for (int i = 0; i < _num_x*_num_y*_num_cmfd_groups; i++){
      if (vec_x[i] != 0.0)
        residual += pow((vec_x[i] - _flux_temp[i]) / vec_x[i], 2);
    }
    residual = pow(residual, 0.5) / (_num_x*_num_y*_num_cmfd_groups);

    /* Increment the interations counter */
    iter++;

    log_printf(DEBUG, "GS iter: %i, res: %f", iter, residual);

    if (iter >= max_iter)
      break;
  }

  log_printf(DEBUG, "linear solver iterations: %i", iter);
}


/**
 * @brief Rescale the initial and converged flux arrays.
 */
void Cmfd::rescaleFlux(){

  FP_PRECISION sum_new, sum_old, scale_val;

  /* Rescale the new and old flux to have an avg source of 1.0 */
  matrix_multiplication(_M, _new_flux, _new_source, _num_x*_num_y, 
      _num_cmfd_groups);
  sum_new = pairwise_sum(_new_source, _num_x*_num_y*_num_cmfd_groups);
  scale_val = _num_x*_num_y*_num_cmfd_groups / sum_new;
  vector_scale(_new_flux, scale_val, _num_x*_num_y*_num_cmfd_groups);
  matrix_multiplication(_M, _old_flux, _old_source, _num_x*_num_y, 
          _num_cmfd_groups);
  sum_old = pairwise_sum(_old_source, _num_x*_num_y*_num_cmfd_groups);
  scale_val = _num_x*_num_y*_num_cmfd_groups / sum_old;
  vector_scale(_old_flux, scale_val, _num_x*_num_y*_num_cmfd_groups);
}


/** @brief Fill in the values in the A matrix, M matrix, and old
 *        scalar flux vector.
 */
void Cmfd::constructMatrices(){

  log_printf(INFO,"Constructing matrices...");
    
  FP_PRECISION value, volume;
  int cell, row;
  Material* material;
  
  /* Zero _A and _M matrices */
  matrix_zero(_M, _num_cmfd_groups*_num_cmfd_groups, _num_x*_num_y);
  matrix_zero(_A, _num_cmfd_groups*(_num_cmfd_groups+4), _num_x*_num_y);
  
  /* Loop over cells */
  #pragma omp parallel for private(value, volume, cell, row, material)
  for (int y = 0; y < _num_y; y++){
    for (int x = 0; x < _num_x; x++){

      cell = y*_num_x + x;
      material = _materials[cell];
      volume = _volumes[cell];

      /* Loop over groups */
      for (int e = 0; e < _num_cmfd_groups; e++){
          
        row = cell*_num_cmfd_groups + e;
    
        /* Absorption term */
        value = material->getSigmaA()[e] * volume;
        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;
        
        /* Out (1st) and in (2nd) scattering */
        for (int g = 0; g < _num_cmfd_groups; g++){
          if (e != g){
            value = material->getSigmaS()[g*_num_cmfd_groups + e] * volume;
            _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;
            value = - material->getSigmaS()[e*_num_cmfd_groups + g] * volume;
            _A[cell][e*(_num_cmfd_groups+4)+g+2] += value;
          }
        }

        /* RIGHT SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[2*_num_cmfd_groups + e]
                - material->getDifTilde()[2*_num_cmfd_groups + e])
          * _cell_height;

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (x != _num_x - 1){
          value = - (material->getDifHat()[2*_num_cmfd_groups + e]
                  + material->getDifTilde()[2*_num_cmfd_groups + e])
                  * _cell_height;

          _A[cell][e*(_num_cmfd_groups+4)+_num_cmfd_groups+2] += value;
        }

        /* LEFT SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[e]
                + material->getDifTilde()[e])
            * _cell_height;

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (x != 0){
          value = - (material->getDifHat()[e]
                     - material->getDifTilde()[e])
              * _cell_height;

          _A[cell][e*(_num_cmfd_groups+4)] += value;
        }

        /* BOTTOM SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[1*_num_cmfd_groups + e]
                - material->getDifTilde()[1*_num_cmfd_groups + e])
                * _cell_width;

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (y != 0){
          value = - (material->getDifHat()[1*_num_cmfd_groups + e]
                  + material->getDifTilde()[1*_num_cmfd_groups + e])
              * _cell_width;

          _A[cell][e*(_num_cmfd_groups+4)+1] += value;
        }

        /* TOP SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[3*_num_cmfd_groups + e]
                + material->getDifTilde()[3*_num_cmfd_groups + e])
            * _cell_width;

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (y != _num_y - 1){
          value = - (material->getDifHat()[3*_num_cmfd_groups + e]
                  - material->getDifTilde()[3*_num_cmfd_groups + e])
                  * _cell_width;

          _A[cell][e*(_num_cmfd_groups+4)+_num_cmfd_groups+3] += value;
        }

        /* Source term */
        for (int g = 0; g < _num_cmfd_groups; g++){
          value = material->getChi()[e] * material->getNuSigmaF()[g]
                  * volume;

          _M[cell][e*_num_cmfd_groups+g] += value;
        }

        log_printf(DEBUG, "cell: %i, vol; %f", cell, volume);

        for (int i = 0; i < _num_cmfd_groups+4; i++)
          log_printf(DEBUG, "i: %i, A value: %f",
                     i, _A[cell][e*(_num_cmfd_groups+4)+i]);

        for (int i = 0; i < _num_cmfd_groups; i++)
          log_printf(DEBUG, "i: %i, M value: %f",
                     i, _M[cell][e*(_num_cmfd_groups+4)+i]);
      }
    }
  }

  log_printf(INFO, "Done constructing matrices...");
}


/**
 * @brief Update the MOC flux in each FSR.
 */
void Cmfd::updateMOCFlux(){

  log_printf(INFO, "Updating MOC flux...");

  /* Initialize variables */
  FP_PRECISION old_cell_flux, new_cell_flux;

  /* Loop over mesh cells */
  #pragma omp parallel for private(old_cell_flux, new_cell_flux)
  for (int i = 0; i < _num_x*_num_y; i++){

    std::vector<int>::iterator iter;

    /* Loop over CMFD groups */
    for (int e = 0; e < _num_cmfd_groups; e++){

      /* Get the old and new Mesh cell flux */
      old_cell_flux = _old_flux[i*_num_cmfd_groups + e];
      new_cell_flux = _new_flux[i*_num_cmfd_groups + e];

      for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){

        /* Loop over FRSs in mesh cell */
        for (iter = _cell_fsrs.at(i).begin();
          iter != _cell_fsrs.at(i).end(); ++iter) {

          /* Set new flux in FSR */
          _FSR_fluxes[*iter*_num_moc_groups+h] =
            new_cell_flux/old_cell_flux * _FSR_fluxes[*iter*_num_moc_groups+h];

          log_printf(DEBUG, "Updating flux in FSR: %i, cell: %i, group: "
                     "%i, ratio: %f", *iter ,i, h,
                      new_cell_flux / old_cell_flux);
        }
      }
    }
  }
}


/**
 * @brief Compute diffusion correction factors to correct diffusion
 *        coefficients in optically thick mesh cells.
 * @param d old diffusion coefficient
 * @param h height of cell
 * @return correction factor
 */
FP_PRECISION Cmfd::computeDiffCorrect(FP_PRECISION d, FP_PRECISION h){

  if (_optically_thick){

    /* Initialize variables */
    FP_PRECISION alpha, mu, expon;
    FP_PRECISION rho, F;
    rho = 0.0;

    /* Loop over polar angles */
    for (int p = 0; p < _quad->getNumPolarAngles(); p++){
      mu = cos(asin(_quad->getSinTheta(p)));
      expon = exp(- h / (3 * d * mu));
      alpha = (1 + expon) / (1 - expon) - 2 * (3 * d * mu) / h;
      rho += mu * _quad->getWeight(p) * alpha;
    }

    /* Compute correction factor, F */
    F = 1.0 + h * rho / (2 * d);

    return F;
  }
  else
    return 1.0;
}


/**
 * @brief Set the FSR materials array pointer.
 * @param FSR_materials pointer to FSR_materials array
 */
void Cmfd::setFSRMaterials(Material** FSR_materials){
  _FSR_materials = FSR_materials;
}


/**
 * @brief Set the pointer to the array of FSR_volumes.
 * @param FSR_volumes array of FSR volumes
 */
void Cmfd::setFSRVolumes(FP_PRECISION* FSR_volumes){
  _FSR_volumes = FSR_volumes;
}


/**
 * @brief Set pointer to FSR flux array.
 * @param scalar_flux pointer to FSR flux array
 */
void Cmfd::setFSRFluxes(FP_PRECISION* scalar_flux){
  _FSR_fluxes = scalar_flux;
}


/**
 * @brief Set successive over-relaxation relaxation factor.
 * @param SOR_factor
 */
void Cmfd::setSORRelaxationFactor(FP_PRECISION SOR_factor){
    
  if (SOR_factor <= 0.0 || SOR_factor >= 2.0)
    log_printf(ERROR, "The successive over-relaxation relaxation factor "
        "must be > 0 and < 2. Input value: %i", SOR_factor);

  _SOR_factor = SOR_factor;
}


/**
 * @brief Get the number of coarse CMFD energy groups.
 * @return the number of CMFD energy groups
 */
int Cmfd::getNumCmfdGroups(){
  return _num_cmfd_groups;
}


/**
 * @brief Get the CMFD group given an MOC group.
 * @param group the MOC energy group
 * @return the CMFD energy group
 */
int Cmfd::getCmfdGroup(int group){
    return _group_indices_map[group];
}


/**
 * @brief Set the CMFD energy group structure.
 * @param an array of the CMFD group boundaries
 * @param the length of the group_indices array
 */
void Cmfd::setGroupStructure(int* group_indices, int length_group_indices){
    
  _num_cmfd_groups = length_group_indices - 1;

  /* Allocate memory */
  if (_group_indices == NULL){
    _group_indices = new int[length_group_indices];
  }
    
  if (group_indices == NULL){
    for (int i = 0; i < length_group_indices; i++){
      _group_indices[i] = i;
    }
  }
  else{
    if (group_indices[0] != 1)
      log_printf(ERROR, "The first value in group indices must be 1!");    

    /* Set first group indice to 0 */
    _group_indices[0] = 0;
        
    /* Set MOC group bounds for rest of CMFD energy groups */
    for (int i = 1; i < length_group_indices; i++){
      /* Check that the group indices are always increasing */
      if (group_indices[i] <= group_indices[i-1])
        log_printf(ERROR, "The group indices must be increasing!");
            
      _group_indices[i] = group_indices[i] - 1;
      log_printf(INFO, "group indices %i: %i", i, _group_indices[i]+1);
    }
  }
}


/**
 * @brief Initialize the flux arrays.
 */
void Cmfd::initializeFlux(){

  /* Allocate memory for fluxes and volumes */
  try{
    _new_flux = new FP_PRECISION[_num_x*_num_y*_num_cmfd_groups];
    _old_flux = new FP_PRECISION[_num_x*_num_y*_num_cmfd_groups];
    _flux_temp = new FP_PRECISION[_num_x*_num_y*_num_cmfd_groups];
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh cell fluxes, "
               "lengths, and volumes. Backtrace:%s", e.what());
  }

  /* Set initial Mesh cell flux to 1.0 and allocate memory for FSR vectors */
  for (int y = 0; y < _num_y; y++){
    for (int x = 0; x < _num_x; x++){
      for (int g = 0; g < _num_cmfd_groups; g++){
        _new_flux[(y*_num_x+x)*_num_cmfd_groups + g] = 1.0;
        _old_flux[(y*_num_x+x)*_num_cmfd_groups + g] = 1.0;
        _flux_temp[(y*_num_x+x)*_num_cmfd_groups + g] = 1.0;
      }
    }
  }
}


/**
 * @brief Initialize the CMFD materials.
 */
void Cmfd::initializeMaterials(){

  Material* material;

  try{
    _materials = new Material*[_num_x*_num_y];

    for (int y = 0; y < _num_y; y++){
      for (int x = 0; x < _num_x; x++){
        material = new Material(y*_num_x+x);
        material->setNumEnergyGroups(_num_cmfd_groups);
        _materials[y*_num_x+x] = material;
      }
    }
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh cell materials. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initializes the vector of vectors that links cmfd cells with FSRs.
 * @details This method is called by the geometry once the cmfd mesh has been
 *          initialized by the geometry. This method allocates a vector for
 *          each cmfd cell that is used to store the FSR ids contained within
 *          that cell.
 */
void Cmfd::initializeCellMap(){

  /* Allocate memory for mesh cell FSR vectors */
  for (int y = 0; y < _num_y; y++){
    for (int x = 0; x < _num_x; x++){
      std::vector<int> *fsrs = new std::vector<int>;
      _cell_fsrs.push_back(*fsrs);
    }
  }
}


/**
 * @brief Initialize and set array that links the MOC energy groups to the
 *        cmfd energy groups.
 * @details This method initializes the _group_indices_map, which is a 1D array
 *          of length _num_moc_groups that maps the MOC energy groups to cmfd 
 *          energy groups. The indices into _group_indices_map are the MOC
 *          energy groups and the values are the cmfd energy groups.
 */
void Cmfd::initializeGroupMap(){

  /* Allocate memory */
  if (_group_indices_map == NULL){
    _group_indices_map = new int[_num_moc_groups];
  }    
    
  /* Create group indices map */
  for (int e = 0; e < _num_cmfd_groups; e++){
    for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){
      _group_indices_map[h] = e;
    }
  }
}



/**
 * @brief Find the cmfd surface that a LocalCoords object lies on.
 * @details If the coords is not on a surface, -1 is returned. Otherwise,
 *        the surface ID is returned. 
 * @param The CMFD cell ID that the local coords is in.
 * @param The coords being evaluated.
 * @return The surface ID.
 */
int Cmfd::findCmfdSurface(int cell, LocalCoords* coords){
  Point* point = coords->getHighestLevel()->getPoint();
  return _lattice->getLatticeSurface(cell, point);
}


/**
 * @brief Find the cmfd cell that a LocalCoords object is in. 
 * @param The coords being evaluated.
 * @return The CMFD cell ID.
 */
int Cmfd::findCmfdCell(LocalCoords* coords){
  Point* point = coords->getHighestLevel()->getPoint();
  return _lattice->getLatticeCell(point);
}


/**
 * @brief The Lattice object used as the CMFD mesh. 
 * @param Pointer to the lattice object.
 */
void Cmfd::setLattice(Lattice* lattice){
    _lattice = lattice;
}


/**
 * @brief The structure of the Lattice to be used as the CMFD mesh.
 * @param The number of cells in the x direction.
 * @param The number of cells in the y direction.
 */
void Cmfd::setLatticeStructure(int num_x, int num_y){
  setNumX(num_x);
  setNumY(num_y);
}


/**
 * @brief Returns the Lattice object used as the CMFD mesh. 
 * @return A pointer to a Lattice object.
 */
Lattice* Cmfd::getLattice(){
  return _lattice;
}


/**
 * @brief Add an FSR ID to a vector that contains all the FSR IDs
 *        contained within a CMFD mesh cell.
 * @param The CMFD cell ID.
 * @param The FSR ID.
 */
void Cmfd::addFSRToCell(int cmfd_cell, int fsr_id){
  _cell_fsrs.at(cmfd_cell).push_back(fsr_id);
}


/**
 * @brief Set the number of MOC energy groups.
 * @param num_groups number of MOC energy groups
 */
void Cmfd::setNumMOCGroups(int num_groups){
  _num_moc_groups = num_groups;
}


/**
 * @brief Get the number of MOC energy groups.
 * @return the number of MOC energy groups
 */
int Cmfd::getNumMOCGroups(){
  return _num_moc_groups;
}


/**
 * @brief Get the number of CMFD cells.
 * @return the number of CMFD cells
 */
int Cmfd::getNumCells(){
  return _num_x*_num_y;
}


/**
 * @brief Set the pointer to the Mesh surface currents array.
 * @param surface_currents pointer to Mesh surface currents array
 */
void Cmfd::setSurfaceCurrents(FP_PRECISION* surface_currents){
  _surface_currents = surface_currents;
}


/**
 * @brief set the number of FSRs.
 * @param the number of FSRs
 */
void Cmfd::setNumFSRs(int num_fsrs){
    _num_FSRs = num_fsrs;
}


/** @brief Split the currents of the Mesh cell corners to the nearby surfaces.
 * @details left bottom corner -> bottom surface and left surface
 *          of mesh cell below; right bottom corner -> bottom surface
 *          and right surface of mesh cell below; right top corner ->
 *          right surface and top surface of mesh cell to the right;
 *          left top corner -> left surface and top surface of mesh
 *          cell to the left. The currents tallied on a corner is split
 *          equally to the adjoining surfaces.
 */
void Cmfd::splitCorners(){

  log_printf(INFO, "splitting corners...");
    
  int ncg = _num_cmfd_groups;

  for (int x = 0; x < _num_x; x++){
    for (int y = 0; y < _num_y; y++){
        
      /* Split the LEFT BOTTOM CORNER */
        
      /* If cell is not on left or bottom geometry edge
       * give to bottom surface and left surface of mesh cell below */
      if (x > 0 && y > 0){
    
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f",
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
          _surface_currents[((y-1)*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
          _surface_currents[(y*_num_x+x-1)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
        }
      }
      /* If cell is on left geometry edge
       * give to bottom surface and left surfaces */
      else if (x == 0 && y != 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 1*ncg + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
          _surface_currents[((y-1)*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
        }
      }
      /* If cell is on bottom geometry edge
       * give to bottom surface and left surfaces */
      else if (x != 0 && y == 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
          _surface_currents[(y*_num_x+x-1)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e];
        }
      }
      
      /* Split the RIGHT BOTTOM CORNER */
      
      /* If cell is not on right or bottom geometry edge
       * give to bottom surface and right surface of mesh cell below */
      if (x < _num_x - 1 && y > 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
          _surface_currents[((y-1)*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
          _surface_currents[(y*_num_x+x+1)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
        }
      }
      /* If cell is on right geometry edge
       * give to bottom surface and right surface */
      else if (x == _num_x - 1 && y > 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 1*ncg + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
          _surface_currents[((y-1)*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
        }
      }
      /* If cell is on bottom geometry edge
       * give to bottom surface and right surface */
      else if (x < _num_x - 1 && y == 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 2*ncg + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
          _surface_currents[(y*_num_x+x+1)*ncg*8 + 1*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e];
        }
      }
      
      /* Split the RIGHT TOP CORNER */
      
      /* If cell is not on right or top geometry edge
       * give to right surface and top surface of mesh cell to the right */
      if (x < _num_x - 1 && y < _num_y - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
          _surface_currents[(y*_num_x+x+1)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
          _surface_currents[((y+1)*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
        }
      }
      /* If cell is on right geometry edge
       * give to right surface and top surface */
      else if (x == _num_x - 1 && y != _num_y - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 3*ncg + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
          _surface_currents[((y+1)*_num_x+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
        }
      }
      /* If cell is on top geometry edge
       * give to right surface and top surface */
      else if (x != _num_x - 1 && y == _num_y - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + 2*ncg + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
          _surface_currents[(y*_num_x+x+1)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e];
        }
      }
      
      /* Split the LEFT TOP CORNER */
      
      /* If cell is not on left or top geometry edge
       * give to left surface and top surface of mesh cell to the left */
      if (x > 0 && y < _num_y - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
          _surface_currents[(y*_num_x+x-1)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
          _surface_currents[((y+1)*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
        }
      }
      /* If cell is on left geometry edge
       * give to top surface and left surface */
      else if (x == 0 && y != _num_y - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 3*ncg + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
          _surface_currents[((y+1)*_num_x+x)*ncg*8 + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
        }
      }
      /* If cell is on top geometry edge
       * give to top surface and left surface */
      else if (x != 0 && y == _num_y - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f", 
              y*_num_x+x,e, _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e]);
          _surface_currents[(y*_num_x+x)*ncg*8 + e] += 
              _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
          _surface_currents[(y*_num_x+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
          _surface_currents[(y*_num_x+x-1)*ncg*8 + 3*ncg + e] += 
              0.5 * _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e];
        }
      }
      
      for (int e = 0; e < ncg; e++){
        _surface_currents[(y*_num_x+x)*ncg*8 + 4*ncg + e] = 0.0;
        _surface_currents[(y*_num_x+x)*ncg*8 + 5*ncg + e] = 0.0;
        _surface_currents[(y*_num_x+x)*ncg*8 + 6*ncg + e] = 0.0;
        _surface_currents[(y*_num_x+x)*ncg*8 + 7*ncg + e] = 0.0;
      }
    }
  }
}


/**
 * @brief Get the ID of the Mesh cell next to given Mesh cell.
 * @param cell_num current Mesh cell ID
 * @param surface_id Mesh cell surface ID to look across for neighboring cell
 * @return neighboring Mesh cell ID
 */
int Cmfd::getCellNext(int cell_num, int surface_id){

  int cell_next = -1;

  if (surface_id == 0){
    if (cell_num % _num_x != 0)
      cell_next = cell_num - 1;
  }

  else if (surface_id == 1){
    if (cell_num / _num_x != 0)
      cell_next = cell_num - _num_x;
  }

  else if (surface_id == 2){
    if (cell_num % _num_x != _num_x - 1)
      cell_next = cell_num + 1;
  }

  else if (surface_id == 3){
    if (cell_num / _num_x != _num_y - 1)
      cell_next = cell_num + _num_x;
  }

  return cell_next;
}


/**
 * @brief Return whether optically thick diffusion correction factor is in use.
 * @return whether optically thick diffusion correction factor is in use.
 */
bool Cmfd::isOpticallyThick(){
  return _optically_thick;
}


/**
 * @brief Set whether optically thick diffusion correction factor is in use.
 * @param optically_thick boolean indicating whether optically thick diffusion
 *        correction factor is in use.
 */
void Cmfd::setOpticallyThick(bool optically_thick){
  _optically_thick = optically_thick;
}


/**
 * @brief Return the under-relaxation factor used in MOC flux updates.
 * @return _relax_factor the MOC flux under-relaxation factor.
 */
FP_PRECISION Cmfd::getMOCRelaxationFactor(){
  return _relax_factor;
}


/**
 * @brief Set the under-relaxation factor used in MOC flux updates.
 * @param relax_factor the MOC flux under-relaxation factor
 */
void Cmfd::setMOCRelaxationFactor(FP_PRECISION relax_factor){
  _relax_factor = relax_factor;
}


/**
 * @brief Set the cmfd boundary type for a given surface.
 * @details The CMFD boundary is assumed to be rectangular with 4
 *          boundary surfaces. The surfaces UIDs are 0 (left),
 *          1 (bottom), 2 (right), and 3 (top).
 * @param side the cmfd surface UID.
 * @param boundary the boundaryType of the surface.
 */
void Cmfd::setBoundary(int side, boundaryType boundary){
  _boundaries[side] = boundary;
}


/**
 * @brief Get the boundaryType for one side of the Mesh.
 * @param side the Mesh surface ID.
 * @return the boundaryType for the surface.
 */
boundaryType Cmfd::getBoundary(int side){
  return _boundaries[side];
}


/**
 * @brief Return the cmfd cell ID that an FSR lies in.
 * @detail Note that a cmfd cell is not an actual Cell object; rather, a cmfd
 *         cell is just a way of describing each of the rectangular regions
 *         that make up a cmfd lattice. Cmfd cells are numbered with 0 in the
 *         lower left corner and monotonically increasing from left to right.
 *         from left to right. For example, he indices for a 4 x 4 lattice are:
 *                  12  13  14  15
 *                  8    9  10  11
 *                  4    5   6   7
 *                  0    1   2   3
 * @param The FSR ID.
 * @return The cmfd cell ID. Return -1 if cell is not found.
 */
int Cmfd::convertFSRIdToCmfdCell(int fsr_id){

  std::vector<int>::iterator iter;    
  for (int cell=0; cell < _num_x*_num_y; cell++){

      for (iter = _cell_fsrs.at(cell).begin();
           iter != _cell_fsrs.at(cell).end(); ++iter){
          if (*iter  == fsr_id)
              return cell;
      }
  }

  return -1;  
}


/**
 * @brief Return a pointer to the vector of vectors that contains 
 *        the FSRs that lie in each cell.
 * @return Vector of vectors containing FSR IDs in each cell.
 */
std::vector< std::vector<int> > Cmfd::getCellFSRs(){
  return _cell_fsrs;
}

 
/**
 * @brief Set the vector of vectors that contains.
 *        the FSRs that lie in each cell.
 * @param Vector of vectors containing FSR IDs in each cell.
 */
void Cmfd::setCellFSRs(std::vector< std::vector<int> > cell_fsrs){
  _cell_fsrs = cell_fsrs;
}


/**
 * @brief Set flag indicating whether to update the MOC flux.
 * @param Flag saying whether to update MOC flux.
 */
void Cmfd::setFluxUpdateOn(bool flux_update_on){
  _flux_update_on = flux_update_on;
}


/**
 * @brief Get flag indicating whether to update the MOC flux.
 * @return Flag saying whether to update MOC flux.
 */
bool Cmfd::isFluxUpdateOn(){
 return _flux_update_on;
}


/**
 * @brief Sets the threshold for cmfd source convergence (>0)
 * @param source_thresh the threshold for source convergence
 */
void Cmfd::setSourceConvergenceThreshold(FP_PRECISION source_thresh) {

  if (source_thresh <= 0.0)
    log_printf(ERROR, "Unable to set the cmfd source convergence threshold to"
              " %f since the threshold must be positive.", source_thresh);

  _source_convergence_threshold = source_thresh;
}


/**
 * @brief Sets the type of polar angle quadrature set to use (ie, TABUCHI
 *        or LEONARD).
 * @param quadrature_type the polar angle quadrature type
 * @param num_polar the number of polar angles
 */
void Cmfd::setPolarQuadrature(quadratureType quadrature_type, int num_polar) {

  /* Deletes the old Quadrature if one existed */
  if (_quad != NULL)
    delete _quad;

  _quad = new Quadrature(quadrature_type, num_polar);
}
