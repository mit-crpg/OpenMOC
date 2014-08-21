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
Cmfd::Cmfd(double criteria) {

  /* Initialize Geometry and Mesh-related attribute */
  _quad = new Quadrature(TABUCHI);
  _SOR_factor = 1.0;

  /* Global variables used in solving CMFD problem */
  _l2_norm = 1.0;
  _conv_criteria = criteria;
  _cx = 1;
  _cy = 1;
  _width = 0.;
  _height = 0.;
  _cell_width = 0.;
  _cell_height = 0.;
  _update_flux = true;
  _overlay_mesh = true;
  _optically_thick = true;
  _SOR_factor = 1.0;
  _num_FSRs = 0;
  _relax_factor = 0.6;

  /* Energy group problem parameters */
  _num_groups = 0;
  _num_cmfd_groups = 0;

  /* matrices */
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
    for (int i = 0; i < _cx*_cy; i++)
      delete [] _M[i];

    delete [] _M;
  }

  if (_A != NULL){
    for (int i = 0; i < _cx*_cy; i++)
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
 * @param cx number of Mesh cells in a row
 */
void Cmfd::setNumX(int cx){

  if (cx < 1)
    log_printf(ERROR, "The number of lattice cells in the x direction "
               "must be > 0. Input value: %i", cx);

  _cx = cx;
  if (_width != 0.)
    _cell_width = _width / _cx;
}


/**
 * @brief Set the number of Mesh cells in a column
 * @param cy number of Mesh cells in a column
 */
void Cmfd::setNumY(int cy){

  if (cy < 1)
    log_printf(ERROR, "The number of lattice cells in the y direction "
               "must be > 0. Input value: %i", cy);

  _cy = cy;
  if (_height != 0.)
    _cell_height = _height / _cy;
}


/**
 * @brief Get the number of Mesh cells in a row.
 * @return number of Mesh cells in a row
 */
int Cmfd::getNumX(){
  return _cx;
}


/**
 * @brief Get the number of Mesh cells in a column
 * @return number of Mesh cells in a column
 */
int Cmfd::getNumY(){
  return _cy;
}


/**
 * @brief Set Mesh width.
 * @param width physical width of Mesh
 */
void Cmfd::setWidth(double width){
  _width = width;
  if (_cx != 0)
    _cell_width = _width / _cx;
}


/**
 * @brief Set Mesh height.
 * @param height physical height of Mesh
 */
void Cmfd::setHeight(double height){
  _height = height;
  if (_cy != 0)
    _cell_height = _height / _cy;
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
  double volume, flux, abs, tot, nu_fis, chi, dif_coef;
  FP_PRECISION* scat;

  /* Initialize tallies for each parameter */
  double abs_tally, nu_fis_tally, dif_tally, rxn_tally;
  double vol_tally, tot_tally, neut_prod_tally;
  double scat_tally[_num_cmfd_groups];
  double chi_groups[_num_cmfd_groups];
  double chi_tally[_num_cmfd_groups];

  /* Pointers to material objects */
  Material* fsr_material;
  Material* cell_material;

  /* Loop over Mesh cells */
  #pragma omp parallel for private(volume, flux, abs, tot, nu_fis, chi, \
    dif_coef, scat, abs_tally, nu_fis_tally, dif_tally, rxn_tally,  \
    vol_tally, tot_tally, scat_tally, fsr_material, cell_material, \
    neut_prod_tally, chi_tally, chi_groups)
  for (int i = 0; i < _cx * _cy; i++){

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

      /* Loop over FSRs in Mesh cell */
      for (iter = _cell_fsrs.at(i).begin();
        iter != _cell_fsrs.at(i).end(); ++iter){

        fsr_material = _FSR_materials[*iter];
        volume = _FSR_volumes[*iter];
        scat = fsr_material->getSigmaS();
        vol_tally += volume;

        /* Chi tallies */
        for (int b = 0; b < _num_cmfd_groups; b++){
          chi_groups[b] = 0.0;

          for (int h = _group_indices[b]; h < _group_indices[b+1]; h++)
            chi_groups[b] += fsr_material->getChi()[h];

          for (int h = 0; h < _num_groups; h++){
            chi_tally[b] += chi_groups[b] * fsr_material->getNuSigmaF()[h] *
                            _FSR_fluxes[(*iter)*_num_groups+h] * volume;
            neut_prod_tally += chi_groups[b] *
                               fsr_material->getNuSigmaF()[h] *
                               _FSR_fluxes[(*iter)*_num_groups+h] * volume;
          }
        }

        /* Loop over fine energy groups within this CMFD coarse group */
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){

          /* Gets FSR volume, material, and cross sections */
          flux = _FSR_fluxes[(*iter)*_num_groups+h];
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
          for (int g = 0; g < _num_groups; g++){
              scat_tally[getCmfdGroup(g)] +=
                  scat[g*_num_groups+h] * flux * volume;
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

      /* set chi */
      if (neut_prod_tally != 0.0)
        cell_material->setChiByGroup(chi_tally[e] / neut_prod_tally, e+1);
      else
        cell_material->setChiByGroup(0.0,e+1);

      log_printf(DEBUG, "cell: %i, group: %i, vol: %f, siga: %f, sigt: %f,"
                 " nu_sigf: %f, dif_coef: %f, flux: %f, chi: %f", i, e,
                 vol_tally, abs_tally / rxn_tally, tot_tally / rxn_tally,
                 nu_fis_tally / rxn_tally, dif_tally / rxn_tally,
                 rxn_tally / vol_tally, chi_tally[e] / (neut_prod_tally+1e-12));

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

  double d, d_next, d_hat, d_tilde;
  double current, flux, flux_next, f, f_next;
  double length, length_perpen, next_length_perpen;
  double sense;
  int next_surface;
  int cell, cell_next;

  /* Loop over mesh cells in y direction */
  #pragma omp parallel for private(d, d_next, d_hat, d_tilde, current, flux, \
    flux_next, f, f_next, length, length_perpen, next_length_perpen, \
    sense, next_surface, cell, cell_next)
  for (int y = 0; y < _cy; y++){

    /* Loop over Mesh cells in x direction */
    for (int x = 0; x < _cx; x++){

      cell = y*_cx+x;

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

            current = sense * _currents[cell*_num_cmfd_groups*8 +
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
            current = sense * _currents[cell*_num_cmfd_groups*8 +
                      surface*_num_cmfd_groups + e] - sense
                      * _currents[cell_next*_num_cmfd_groups*8 +
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
                     y*_cx + x, e, surface, flux, current, d, d_hat,
                     d_tilde);
        }
      }
    }
  }
}



/** @brief CMFD solver that solves the diffusion problem.
 * @return k-effective the solution eigenvalue.
 */
double Cmfd::computeKeff(){

  log_printf(INFO, "Running diffusion solver...");

  /* Create matrix and vector objects */
  if (_A == NULL){
    try{
    
      _M = new double*[_cx*_cy];
      _A = new double*[_cx*_cy];
      _old_source = new double[_cx*_cy*_num_cmfd_groups];
      _new_source = new double[_cx*_cy*_num_cmfd_groups];
      _volumes = new double[_cx*_cy];

      initializeFlux();
      initializeMaterials();

      for (int i = 0; i < _cx*_cy; i++){
        _M[i] = new double[_num_cmfd_groups*_num_cmfd_groups];
        _A[i] = new double[_num_cmfd_groups*(_num_cmfd_groups+4)];
      }
    }
    catch(std::exception &e){
      log_printf(ERROR, "Could not allocate memory for the CMFD mesh objects. "
                 "Backtrace:%s", e.what());
    }
  }

  /* Initialize variables */
  double sum_new, sum_old, val, norm, scale_val;
  double flux_conv = 1e-8;
  int row;

  /* Compute the cross sections and surface diffusion coefficients */
  computeXS();
  computeDs();

  /* Construct matrices */
  constructMatrices();

  /* Compute and normalize the initial source */
  matMultM(_M, _old_flux, _old_source);
  sum_old = vecSum(_old_source);
  scale_val = (_cx * _cy * _num_cmfd_groups) / sum_old;
  vecScale(_old_source, scale_val);
  vecCopy(_old_flux, _new_flux);
  vecScale(_new_flux, scale_val);
  sum_old = _cx * _cy * _num_cmfd_groups;
  
  /* Power iteration diffusion solver */
  for (int iter = 0; iter < 25000; iter++){
      
    /* Solve phi = A^-1 * old_source */
    linearSolve(_A, _new_flux, _old_source, flux_conv);
      
    /* Compute the new source */
    matMultM(_M, _new_flux, _new_source);
    sum_new = vecSum(_new_source);
    
    /* Compute and set keff */
    _k_eff = sum_new / sum_old;
    
    /* Scale the old source by keff */
    vecScale(_old_source, _k_eff);
    
    /* Compute the L2 norm of source error */
    norm = 0.0;
    
    for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++){
      if (_new_source[i] != 0.0)
        norm += pow((_new_source[i] - _old_source[i]) / _new_source[i], 2);
    }

    norm = sqrt(norm / (_cx*_cy*_num_cmfd_groups));
    
    scale_val = (_cx * _cy * _num_cmfd_groups) / sum_new;
    vecScale(_new_source, scale_val);
    vecCopy(_new_source, _old_source);
    
    log_printf(INFO, "GS POWER iter: %i, keff: %f, error: %f",
               iter, _k_eff, norm);
    
    /* Check for convergence */
    if (norm < _conv_criteria)
      break;
  }

  /* rescale the old and new flux */
  rescaleFlux();

  /* update the MOC flux */
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
void Cmfd::linearSolve(double** mat, double* vec_x, double* vec_b,
                       double conv, int max_iter){

  double norm = 1e10;
  int row, cell;
  double val;
  int iter = 0;

  while (norm > conv){

    /* Pass new flux to old flux */
    vecCopy(vec_x, _flux_temp);

    /* Iteration over red cells */
    #pragma omp parallel for private(row, val, cell)
    for (int y = 0; y < _cy; y++){
      for (int x = y % 2; x < _cx; x += 2){

        cell = y*_cx+x;

        for (int g = 0; g < _num_cmfd_groups; g++){

          row = (y*_cx+x)*_num_cmfd_groups + g;
          val = 0.0;

          /* Previous flux term */
          val += (1.0 - _SOR_factor) * vec_x[row];

          /* Source term */
          val += _SOR_factor * vec_b[row] / mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _SOR_factor * vec_x[row - _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Bottom surface */
          if (y != 0)
            val -= _SOR_factor * vec_x[row - _cx * _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+1] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < _num_cmfd_groups; e++){
            if (e != g)
              val -= _SOR_factor * vec_x[(y*_cx+x)*_num_cmfd_groups+e] *
                     mat[cell][g*(_num_cmfd_groups+4)+2+e] /
                     mat[cell][g*(_num_cmfd_groups+4)+g+2];
          }

          /* Right surface */
          if (x != _cx - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+2] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Top surface */
          if (y != _cy - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups*_cx] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+3] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          vec_x[row] = val;
        }
      }
    }

    /* Iteration over black cells */
    #pragma omp parallel for private(row, val, cell)
    for (int y = 0; y < _cy; y++){
      for (int x = 1 - y % 2; x < _cx; x += 2){

        cell = y*_cx+x;

        for (int g = 0; g < _num_cmfd_groups; g++){

          row = (y*_cx+x)*_num_cmfd_groups + g;
          val = 0.0;

          /* Previous flux term */
          val += (1.0 - _SOR_factor) * vec_x[row];

          /* Source term */
          val += _SOR_factor * vec_b[row] / mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _SOR_factor * vec_x[row - _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Bottom surface */
          if (y != 0)
            val -= _SOR_factor * vec_x[row - _cx * _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+1] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < _num_cmfd_groups; e++){
            if (e != g)
              val -= _SOR_factor * vec_x[(y*_cx+x)*_num_cmfd_groups+e] *
                     mat[cell][g*(_num_cmfd_groups+4)+2+e] /
                     mat[cell][g*(_num_cmfd_groups+4)+g+2];
          }

          /* Right surface */
          if (x != _cx - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+2] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Top surface */
          if (y != _cy - 1)
            val -= _SOR_factor * vec_x[row + _num_cmfd_groups*_cx] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+3] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          vec_x[row] = val;
        }
      }
    }

    norm = 0.0;

    for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++){
      if (vec_x[i] != 0.0)
        norm += pow((vec_x[i] - _flux_temp[i]) / vec_x[i], 2);
    }

    norm = pow(norm, 0.5) / (_cx*_cy*_num_cmfd_groups);

    iter++;

    log_printf(DEBUG, "GS iter: %i, norm: %f", iter, norm);

    if (iter >= max_iter)
      break;
  }

  log_printf(DEBUG, "linear solver iterations: %i", iter);
}


/**
 * @brief Rescale the initial and converged flux arrays.
 */
void Cmfd::rescaleFlux(){

  double sum_new, sum_old, scale_val;

  /* Rescale the new and old flux to have an avg source of 1.0 */
  matMultM(_M, _new_flux, _new_source);
  sum_new = vecSum(_new_source);
  scale_val = _cx*_cy*_num_cmfd_groups / sum_new;
  vecScale(_new_flux, scale_val);
  matMultM(_M, _old_flux, _old_source);
  sum_old = vecSum(_old_source);
  scale_val = _cx*_cy*_num_cmfd_groups / sum_old;
  vecScale(_old_flux, scale_val);
}


/**
 * @brief Scale vectgor by some scalar value.
 * @param vec vector to be scaled
 * @param scale_val value to scale vector
 */
void Cmfd::vecScale(double* vec, double scale_val){

  #pragma omp parallel for
  for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++)
    vec[i] *= scale_val;
}


/**
 * @brief Set every element in vector to some value.
 * @param vec vector to be set
 * @param val value to set each element
 */
void Cmfd::vecSet(double* vec, double val){

  #pragma omp parallel for
  for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++)
    vec[i] = val;
}


/**
 * @brief Multiply matrix by vector (i.e., y = M *x).
 * @param mat source matrix
 * @param vec_x x vector
 * @param vec_y y vector
 */
void Cmfd::matMultM(double** mat, double* vec_x, double* vec_y){

  vecSet(vec_y, 0.0);

  for (int i = 0; i < _cx*_cy; i++){
    for (int g = 0; g < _num_cmfd_groups; g++){
      for (int e = 0; e < _num_cmfd_groups; e++){
        vec_y[i*_num_cmfd_groups+g] += mat[i][g*_num_cmfd_groups+e] *
                                   vec_x[i*_num_cmfd_groups+e];
      }
    }
  }
}


/**
 * @brief Sum all elements in a vector.
 * @param vec vector to be summed
 * @return sum of vector elements
 */
double Cmfd::vecSum(double* vec){

  double sum = 0.0;

  for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++)
    sum += vec[i];

    return sum;
}


/**
 * @brief Copy a vector to another vector.
 * @param vec_from vector to be copied
 * @param vec_to vector to receive copied data
 */
void Cmfd::vecCopy(double* vec_from, double* vec_to){

  #pragma omp parallel for
  for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++)
    vec_to[i] = vec_from[i];
}


/**
 * @brief Assign all elements in a matrix to zero.
 * @param mat matrix to be zeroed
 * @param width width of matrix row
 */
void Cmfd::matZero(double** mat, int width){

  #pragma omp parallel for
  for (int i = 0; i < _cx*_cy; i++){
    for (int g = 0; g < _num_cmfd_groups*width; g++)
        mat[i][g] = 0.0;
  }
}


/** @brief Fill in the values in the A matrix, M matrix, and old
 *        scalar flux vector.
 */
void Cmfd::constructMatrices(){

  log_printf(INFO,"Constructing matrices...");
    
  double value, volume;
  int cell, row;
  Material* material;
  
  /* Zero _A and _M matrices */
  matZero(_M, _num_cmfd_groups);
  matZero(_A, _num_cmfd_groups+4);
  
  /* Loop over cells */
  #pragma omp parallel for private(value, volume, cell, row, material)
  for (int y = 0; y < _cy; y++){
    for (int x = 0; x < _cx; x++){

      cell = y*_cx + x;
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
        if (x != _cx - 1){
          value = - (material->getDifHat()[2*_num_cmfd_groups + e]
                  + material->getDifTilde()[2*_num_cmfd_groups + e])
                  * _cell_height;

          _A[cell][e*(_num_cmfd_groups+4)+_num_cmfd_groups+2] += value;
        }

        /* LEFT SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[0*_num_cmfd_groups + e]
                + material->getDifTilde()[0*_num_cmfd_groups + e])
            * _cell_height;

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (x != 0){
          value = - (material->getDifHat()[0*_num_cmfd_groups + e]
                     - material->getDifTilde()[0*_num_cmfd_groups + e])
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
        if (y != _cy - 1){
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

        log_printf(DEBUG, "cel: %i, vol; %f", cell, volume);

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
  double old_cell_flux, new_cell_flux;

  /* Loop over mesh cells */
  #pragma omp parallel for private(old_cell_flux, new_cell_flux)
  for (int i = 0; i < _cx*_cy; i++){

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
          _FSR_fluxes[*iter*_num_groups+h] =
             new_cell_flux / old_cell_flux * _FSR_fluxes[*iter*_num_groups+h];

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
double Cmfd::computeDiffCorrect(double d, double h){

  if (_optically_thick){

    /* Initialize variables */
    double alpha, mu, expon;
    double rho, F;
    rho = 0.0;

    /* Loop over polar angles */
    for (int p = 0; p < 3; p++){
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
 * @brief Return the eigenvalue \f$ k_{eff} \f$.
 * @return the eigenvalue \f$ k_{eff} \f$
 */
double Cmfd::getKeff(){
  return _k_eff;
}


/**
 * @brief Set the FSR materials array pointer.
 * @param FSR_materials pointer to FSR_materials array
 */
void Cmfd::setFSRMaterials(Material** FSR_materials){
  _FSR_materials = FSR_materials;
}


/**
 * @brief Set the fsr volumes by summing the volumes of the FSRs contained
 *        in each Mesh cell.
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
void Cmfd::setSORRelaxationFactor(double SOR_factor){
    
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


void Cmfd::setGroupStructure(int* group_indices, int ncg){
    
  _num_cmfd_groups = ncg - 1;

  /* allocate memory */
  if (_group_indices == NULL){
    _group_indices = new int[ncg];
  }
    
  if (group_indices == NULL){
    for (int i = 0; i < ncg; i++){
      _group_indices[i] = i;
    }
  }
  else{
    if (group_indices[0] != 1)
      log_printf(ERROR, "The first value in group indices must be 1!");    

    _group_indices[0] = 1;
        
    for (int i = 1; i < ncg; i++){
      /* check that the group indices are always increasing */
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
    _new_flux = new double[_cx*_cy*_num_cmfd_groups];
    _old_flux = new double[_cx*_cy*_num_cmfd_groups];
    _flux_temp = new double[_cx*_cy*_num_cmfd_groups];
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh cell fluxes, "
               "lengths, and volumes. Backtrace:%s", e.what());
  }

  /* Set initial Mesh cell flux to 1.0 and allocate memory for FSR vectors */
  for (int y = 0; y < _cy; y++){
    for (int x = 0; x < _cx; x++){
      for (int g = 0; g < _num_cmfd_groups; g++){
        _new_flux[(y*_cx+x)*_num_cmfd_groups + g] = 1.0;
        _old_flux[(y*_cx+x)*_num_cmfd_groups + g] = 1.0;
        _flux_temp[(y*_cx+x)*_num_cmfd_groups + g] = 1.0;
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
    _materials = new Material*[_cx*_cy];

    for (int y = 0; y < _cy; y++){
      for (int x = 0; x < _cx; x++){
        material = new Material(y*_cx+x);
        material->setNumEnergyGroups(_num_cmfd_groups);
        _materials[y*_cx+x] = material;
      }
    }
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh cell materials. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initializes the Cmfd by allocating memory for various arrays.
 * @details This method is called by the geometry once the width of
 *          the Mesh has been determined. This method allocates memory
 *          for the volumes, old and new flux arrays, lengths, bounds,
 *          and cell FSR vectors.
 */
void Cmfd::initializeCellMap(){

  /* Set initial Mesh cell flux to 1.0 and allocate memory for FSR vectors */
  for (int y = 0; y < _cy; y++){
    for (int x = 0; x < _cx; x++){
          
      /* Allocate memory for FSR vector */
      std::vector<int> *fsrs = new std::vector<int>;
      _cell_fsrs.push_back(*fsrs);
    }
  }
}


/**
 * @brief Initializes the Cmfd by allocating memory for various arrays.
 * @details This method is called by the geometry once the width of
 *          the Mesh has been determined. This method allocates memory
 *          for the volumes, old and new flux arrays, lengths, bounds,
 *          and cell FSR vectors.
 */
void Cmfd::initializeGroupMap(){

  /* allocate memory */
  if (_group_indices_map == NULL){
    _group_indices_map = new int[_num_groups];
  }    
    
  /* create group indices map */
  for (int e = 0; e < _num_cmfd_groups; e++){
    for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){
      _group_indices_map[h] = e;
    }
  }
}



/**
 * @brief Find the cmfd surface that a LocalCoords object likes on.
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
void Cmfd::setLatticeStructure(int cx, int cy){
  setNumX(cx);
  setNumY(cy);
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
 * @brief Set the number of energy groups.
 * @param num_groups number of energy groups
 */
void Cmfd::setNumGroups(int num_groups){
  _num_groups = num_groups;
}


/**
 * @brief Get the number of energy groups.
 * @return the number of energy groups
 */
int Cmfd::getNumGroups(){
  return _num_groups;
}


/**
 * @brief Get the number of CMFD cells.
 * @return the number of CMFD cells
 */
int Cmfd::getNumCells(){
  return _cx*_cy;
}


/**
 * @brief Set the pointer to the Mesh surface currents array.
 * @param surface_currents pointer to Mesh surface currents array
 */
void Cmfd::setSurfaceCurrents(double* surface_currents){
  _currents = surface_currents;
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

  for (int x = 0; x < _cx; x++){
    for (int y = 0; y < _cy; y++){
        
      /* split the LEFT BOTTOM CORNER */
        
      /* if cell is not on left or bottom geometry edge
       * give to bottom surface and left surface of mesh cell below */
      if (x > 0 && y > 0){
    
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f",
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 4*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
          _currents[((y-1)*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
          _currents[(y*_cx+x-1)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
        }
      }
      /* if cell is on left geometry edge
       * give to bottom surface and left surfaces */
      else if (x == 0 && y != 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 4*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 1*ncg + e] += 
              _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
          _currents[((y-1)*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
        }
      }
      /* if cell is on bottom geometry edge
       * give to bottom surface and left surfaces */
      else if (x != 0 && y == 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 4*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + e] += 
              _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
          _currents[(y*_cx+x-1)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 4*ncg + e];
        }
      }
      
      /* split the RIGHT BOTTOM CORNER */
      
      /* if cell is not on right or bottom geometry edge
       * give to bottom surface and right surface of mesh cell below */
      if (x < _cx - 1 && y > 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f", 
        y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 5*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
          _currents[((y-1)*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
          _currents[(y*_cx+x+1)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
        }
      }
      /* if cell is on right geometry edge
       * give to bottom surface and right surface */
      else if (x == _cx - 1 && y > 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 5*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 1*ncg + e] += 
              _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
          _currents[((y-1)*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
        }
      }
      /* if cell is on bottom geometry edge
       * give to bottom surface and right surface */
      else if (x < _cx - 1 && y == 0){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 5*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 2*ncg + e] += 
              _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
          _currents[(y*_cx+x+1)*ncg*8 + 1*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 5*ncg + e];
        }
      }
      
      /* split the RIGHT TOP CORNER */
      
      /* if cell is not on right or top geometry edge
       * give to right surface and top surface of mesh cell to the right */
      if (x < _cx - 1 && y < _cy - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 6*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
          _currents[(y*_cx+x+1)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
          _currents[((y+1)*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
        }
      }
      /* if cell is on right geometry edge
       * give to right surface and top surface */
      else if (x == _cx - 1 && y != _cy - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 6*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 3*ncg + e] += 
              _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
          _currents[((y+1)*_cx+x)*ncg*8 + 2*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
        }
      }
      /* if cell is on top geometry edge
       * give to right surface and top surface */
      else if (x != _cx - 1 && y == _cy - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 6*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + 2*ncg + e] += 
              _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
          _currents[(y*_cx+x+1)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 6*ncg + e];
        }
      }
      
      /* split the LEFT TOP CORNER */
      
      /* if cell is not on left or top geometry edge
       * give to left surface and top surface of mesh cell to the left */
      if (x > 0 && y < _cy - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f", 
                     y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 7*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
          _currents[(y*_cx+x-1)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
          _currents[((y+1)*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
        }
      }
      /* if cell is on left geometry edge
       * give to top surface and left surface */
      else if (x == 0 && y != _cy - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f", 
        y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 7*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 3*ncg + e] += 
              _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
          _currents[((y+1)*_cx+x)*ncg*8 + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
        }
      }
      /* if cell is on top geometry edge
       * give to top surface and left surface */
      else if (x != 0 && y == _cy - 1){
        for (int e = 0; e < ncg; e++){
          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f", 
        y*_cx+x,e, _currents[(y*_cx+x)*ncg*8 + 7*ncg + e]);
          _currents[(y*_cx+x)*ncg*8 + e] += 
              _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
          _currents[(y*_cx+x)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
          _currents[(y*_cx+x-1)*ncg*8 + 3*ncg + e] += 
              0.5 * _currents[(y*_cx+x)*ncg*8 + 7*ncg + e];
        }
      }
      
      for (int e = 0; e < ncg; e++){
        _currents[(y*_cx+x)*ncg*8 + 4*ncg + e] = 0.0;
        _currents[(y*_cx+x)*ncg*8 + 5*ncg + e] = 0.0;
        _currents[(y*_cx+x)*ncg*8 + 6*ncg + e] = 0.0;
        _currents[(y*_cx+x)*ncg*8 + 7*ncg + e] = 0.0;
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
    if (cell_num % _cx != 0)
      cell_next = cell_num - 1;
  }

  else if (surface_id == 1){
    if (cell_num / _cx != 0)
      cell_next = cell_num - _cx;
  }

  else if (surface_id == 2){
    if (cell_num % _cx != _cx - 1)
      cell_next = cell_num + 1;
  }

  else if (surface_id == 3){
    if (cell_num / _cx != _cy - 1)
      cell_next = cell_num + _cx;
  }

  return cell_next;
}


/**
 * @brief Return whether CMFD mesh is overlaid on geometry.
 * @return true if mesh overlay is on and false otherwise
 */
bool Cmfd::getOverlayMesh(){
  return _overlay_mesh;
}


/**
 * @brief Set whether CMFD mesh is overlayed on geometry.
 * @param overlay_mesh flag to toggle CMFD mesh overaid
 */
void Cmfd::setOverlayMesh(bool overlay_mesh){
  _overlay_mesh = overlay_mesh;
}


/**
 * @brief Return whether optically thick diffusion correction factor is in use.
 * @return whether optically thick diffusion correction factor is in use
 */
bool Cmfd::getOpticallyThick(){
  return _optically_thick;
}


/**
 * @brief Return whether optically thick diffusion correction factor is in use.
 * @return whether optically thick diffusion correction factor is in use
 */
void Cmfd::setOpticallyThick(bool optically_thick){
  _optically_thick = optically_thick;
}


/**
 * @brief Return the under-relaxation factor used in MOC flux updates.
 * @return _relax_factor the MOC flux under-relaxation factor
 */
double Cmfd::getMOCRelaxationFactor(){
  return _relax_factor;
}


/**
 * @brief Set the under-relaxation factor used in MOC flux updates.
 * @param relax_factor the MOC flux under-relaxation factor
 */
void Cmfd::setMOCRelaxationFactor(double relax_factor){
  _relax_factor = relax_factor;
}


/**
 * @brief Set the Mesh boundary type for left surface.
 * @param side the Mesh surface UID
 * @param boundary the boundaryType of the surface
 */
void Cmfd::setBoundary(int side, boundaryType boundary){
  _boundaries[side] = boundary;
}


/**
 * @brief Get the boundaryType for one side of the Mesh.
 * @param side the Mesh surface ID
 * @return the boundaryType for the surface
 */
boundaryType Cmfd::getBoundary(int side){
  return _boundaries[side];
}


/**
 * @brief Return the Cell ID that an FSR lies in.
 * @param The FSR ID.
 * @return The Cell ID. Return -1 if cell is not found.
 */
int Cmfd::convertFSRIdToCmfdCell(int fsr_id){

  std::vector<int>::iterator iter;    
  for (int cell=0; cell < _cx*_cy; cell++){

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
void Cmfd::setUpdateFlux(bool update_flux){
  _update_flux = update_flux;
}


/**
 * @brief Get flag indicating whether to update the MOC flux.
 * @return Flag saying whether to update MOC flux.
 */
bool Cmfd::getUpdateFlux(){
 return _update_flux;
}
