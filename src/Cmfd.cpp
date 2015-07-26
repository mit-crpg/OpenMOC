#include "Cmfd.h"

/**
 * @brief Constructor initializes boundaries and variables that describe
 *          the Cmfd object.
 * @details The construcor initializes the many variables that describe
 *          the CMFD mesh and are used to solve the nonlinear diffusion
 *          acceleration problem.
 */
Cmfd::Cmfd() {

  /* Initialize Geometry and Mesh-related attribute */
  _polar_quad = NULL;
  _geometry = NULL;
  _SOR_factor = 1.0;

  /* Global variables used in solving CMFD problem */
  _source_convergence_threshold = 1E-6;
  _num_x = 1;
  _num_y = 1;
  _width = 0.;
  _height = 0.;
  _cell_width = 0.;
  _cell_height = 0.;
  _flux_update_on = true;
  _centroid_update_on = true;
  _k_nearest = 3;
  _optically_thick = false;
  _SOR_factor = 1.0;
  _num_FSRs = 0;
  _relax_factor = 0.6;

  /* Energy group and polar angle problem parameters */
  _num_moc_groups = 0;
  _num_cmfd_groups = 0;
  _num_polar = 0;

  /* Set matrices and arrays to NULL */
  _A = NULL;
  _M = NULL;
  _old_flux = NULL;
  _new_flux = NULL;
  _old_source = NULL;
  _new_source = NULL;
  _group_indices = NULL;
  _group_indices_map = NULL;
  _surface_currents = NULL;
  _volumes = NULL;

  /* Initialize boundaries to be reflective */
  _boundaries = new boundaryType[4];
  _boundaries[SURFACE_X_MIN] = REFLECTIVE;
  _boundaries[SURFACE_X_MAX] = REFLECTIVE;
  _boundaries[SURFACE_Y_MIN] = REFLECTIVE;
  _boundaries[SURFACE_Y_MAX] = REFLECTIVE;
}


/**
 * @brief Destructor.
 */
Cmfd::~Cmfd() {

  if (_boundaries != NULL)
    delete [] _boundaries;

  if (_group_indices != NULL)
    delete [] _group_indices;

  if (_group_indices_map != NULL)
    delete [] _group_indices_map;
}


/**
 * @brief Set the number of Mesh cells in a row.
 * @param num_x number of Mesh cells in a row
 */
void Cmfd::setNumX(int num_x) {

  if (num_x < 1)
    log_printf(ERROR, "The number of lattice cells in the x direction "
               "must be > 0. Input value: %i", num_x);

  _num_x = num_x;
  if (_width != 0.)
    _cell_width = _width / _num_x;
}


/**
 * @brief Set the number of Mesh cells in a column
 * @param num_y number of Mesh cells in a column
 */
void Cmfd::setNumY(int num_y) {

  if (num_y < 1)
    log_printf(ERROR, "The number of lattice cells in the y direction "
               "must be > 0. Input value: %i", num_y);

  _num_y = num_y;
  if (_height != 0.)
    _cell_height = _height / _num_y;
}


/**
 * @brief Get the number of Mesh cells in a row.
 * @return The number of Mesh cells in a row
 */
int Cmfd::getNumX() {
  return _num_x;
}


/**
 * @brief Get the number of Mesh cells in a column
 * @return The number of Mesh cells in a column
 */
int Cmfd::getNumY() {
  return _num_y;
}


/**
 * @brief Set Mesh width.
 * @param width Physical width of Mesh
 */
void Cmfd::setWidth(double width) {
  _width = width;
  if (_num_x != 0)
    _cell_width = _width / _num_x;
}


/**
 * @brief Set Mesh height.
 * @param height Physical height of Mesh
 */
void Cmfd::setHeight(double height) {
  _height = height;
  if (_num_y != 0)
    _cell_height = _height / _num_y;
}


/**
 * @brief Create cross-sections and fluxes for each Cmfd cell by
 *        energy condensing and volume averaging cross sections from
 *        the MOC sweep.
 * @details This method performs a cell-wise energy condensation and volume
 *         average of the cross sections of the fine, unstructured FSR mesh.
 *         The cross sections are condense such that all reaction rates and
 *         the neutron produciton rate from fission are conserved. It is
 *         important to note that the volume averaging is performed before
 *         energy condensation in order to properly collapse the diffusion
 *         coefficients.
 */
void Cmfd::computeXS() {

  log_printf(INFO, "Computing CMFD cross-sections...");

  /* Split corner currents to side surfaces */
  splitCorners();

  /* Initialize variables for FSR properties*/
  FP_PRECISION volume, flux, abs, tot, nu_fis, chi, dif_coef;
  FP_PRECISION* scat;

  /* Initialize tallies for each parameter */
  FP_PRECISION abs_tally, nu_fis_tally, dif_tally, rxn_tally;
  FP_PRECISION vol_tally, tot_tally, neut_prod_tally;
  FP_PRECISION trans_tally_group, rxn_tally_group;
  FP_PRECISION scat_tally[_num_cmfd_groups];
  FP_PRECISION chi_tally[_num_cmfd_groups];

  /* Pointers to material objects */
  Material* fsr_material;
  Material* cell_material;

  /* Loop over cmfd cells */
  #pragma omp parallel for private(volume, flux, abs, tot, nu_fis, chi, \
    dif_coef, scat, abs_tally, nu_fis_tally, dif_tally, rxn_tally,  \
    vol_tally, tot_tally, scat_tally, fsr_material, cell_material, \
    neut_prod_tally, chi_tally, trans_tally_group, rxn_tally_group)
  for (int i = 0; i < _num_x * _num_y; i++) {

    cell_material = _materials[i];
    std::vector<int>::iterator iter;

    /* Loop over CMFD coarse energy groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Zero tallies for this group */
      abs_tally = 0.0;
      nu_fis_tally = 0.0;
      dif_tally = 0.0;
      rxn_tally = 0.0;
      vol_tally = 0.0;
      tot_tally = 0.0;
      neut_prod_tally = 0.0;

      /* Zero each group-to-group scattering tally */
      for (int g = 0; g < _num_cmfd_groups; g++) {
        scat_tally[g] = 0;
        chi_tally[g] = 0.0;
      }

      /* Loop over FSRs in cmfd cell to compute chi */
      for (iter = _cell_fsrs.at(i).begin();
           iter != _cell_fsrs.at(i).end(); ++iter) {

        fsr_material = _FSR_materials[*iter];
        volume = _FSR_volumes[*iter];

        /* Chi tallies */
        for (int b = 0; b < _num_cmfd_groups; b++) {
          chi = 0.0;
              
          /* Compute the chi for group b */
          for (int h = _group_indices[b]; h < _group_indices[b+1]; h++)
            chi += fsr_material->getChi()[h];

          for (int h = 0; h < _num_moc_groups; h++) {
            chi_tally[b] += chi * fsr_material->getNuSigmaF()[h] *
                _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
            neut_prod_tally += chi * fsr_material->getNuSigmaF()[h] *
                _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
          }
        }
      }

      /* Loop over MOC energy groups within this CMFD coarse group */
      for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {

        /* Reset transport, rxn, and vol tally for this MOC group */
        trans_tally_group = 0.0;
        rxn_tally_group = 0.0;
        vol_tally = 0.0;

        /* Loop over FSRs in cmfd cell */
        for (iter = _cell_fsrs.at(i).begin();
             iter != _cell_fsrs.at(i).end(); ++iter) {

          fsr_material = _FSR_materials[*iter];
          volume = _FSR_volumes[*iter];
          scat = fsr_material->getSigmaS();
          vol_tally += volume;

          /* Gets FSR volume, material, and cross sections */
          flux = _FSR_fluxes[(*iter)*_num_moc_groups+h];
          abs = fsr_material->getSigmaA()[h];
          tot = fsr_material->getSigmaT()[h];
          nu_fis = fsr_material->getNuSigmaF()[h];

          /* Increment tallies for this group */
          abs_tally += abs * flux * volume;
          tot_tally += tot * flux * volume;
          nu_fis_tally += nu_fis * flux * volume;
          rxn_tally += flux * volume;
          trans_tally_group += tot * flux * volume;
          rxn_tally_group += flux * volume;

          /* Scattering tallies */
          for (int g = 0; g < _num_moc_groups; g++) {
            scat_tally[getCmfdGroup(g)] +=
                scat[g*_num_moc_groups+h] * flux * volume;
          }
        }

        /* Energy collapse diffusion coefficient */
        dif_tally += rxn_tally_group / 
            (3.0 * (trans_tally_group / rxn_tally_group));
      }

      /* Set the Mesh cell properties with the tallies */
      _volumes->setValue(i, 0, vol_tally);
      cell_material->setSigmaAByGroup(abs_tally / rxn_tally, e+1);
      cell_material->setSigmaTByGroup(tot_tally / rxn_tally, e+1);
      cell_material->setNuSigmaFByGroup(nu_fis_tally / rxn_tally, e+1);
      cell_material->setDifCoefByGroup(dif_tally / rxn_tally, e+1);
      _old_flux->setValue(i, e, rxn_tally / vol_tally);

      /* Set chi */
      if (neut_prod_tally != 0.0)
        cell_material->setChiByGroup(chi_tally[e] / neut_prod_tally, e+1);
      else
        cell_material->setChiByGroup(0.0,e+1);

      log_printf(DEBUG, "cell: %i, group: %i, vol: %e, siga: %e, sigt: %e,"
                 " nu_sigf: %e, dif_coef: %e, flux: %e, chi: %e", i, e,
                 vol_tally, abs_tally / rxn_tally, tot_tally / rxn_tally,
                 nu_fis_tally / rxn_tally, dif_tally / rxn_tally,
                 rxn_tally / vol_tally, cell_material->getChi()[e]);

      /* Set scattering xs */
      for (int g = 0; g < _num_cmfd_groups; g++) {
        cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, e+1, g+1);
        log_printf(DEBUG, "scattering from %i to %i: %e", e, g,
                   scat_tally[g] / rxn_tally);
      }
    }
  }
}


/**
 * @brief Compute the surface diffusion coefficients for each cell surface.
 * @details This method uses finite differencing to compute the surface
 *         diffusion coefficients (\f$ \hat{D} \f$) for each cell surface.
 *         Additionally, the surface diffusion coefficent correction factors
 *         (\f$ \tilde{D} \f$) are computed in order to conserve the net
 *         leakage out of each mesh cell. This serves to preserve both local
 *         and global neutron balance.
 * @param moc_iteration MOC iteration number
 */
void Cmfd::computeDs(int moc_iteration) {

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
  for (int y = 0; y < _num_y; y++) {

    /* Loop over Mesh cells in x direction */
    for (int x = 0; x < _num_x; x++) {

      cell = y*_num_x+x;

      /* Loop over side surfaces in a cell */
      for (int surface = 0; surface < 4; surface++) {

        /* Loop over groups */
        for (int e = 0; e < _num_cmfd_groups; e++) {

          /* Get diffusivity and flux for Mesh cell */
          d = _materials[cell]->getDifCoef()[e];
          flux = _old_flux->getValue(cell, e);
          cell_next = getCellNext(cell, surface);

          /* Set halfspace sense of the Surface */
          if (surface == SURFACE_X_MIN || surface == SURFACE_Y_MIN)
            sense = -1.0;
          else
            sense = 1.0;

          /* Set the length of this Surface and the perpendicular Surface */
          if (surface == SURFACE_X_MIN || surface== SURFACE_X_MAX) {
            length = _cell_height;
            length_perpen = _cell_width;
          }
          else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX) {
            length = _cell_width;
            length_perpen = _cell_height;
          }

          /* Compute the optical thickness correction factor */
          f = computeDiffCorrect(d, length_perpen);

          /* If Surface is on a boundary, choose appropriate BCs */
          if (cell_next == -1) {

            current = sense * _surface_currents->getValue
              (cell, surface*_num_cmfd_groups + e);

            /* REFLECTIVE BC */
            if (_boundaries[surface] == REFLECTIVE) {

              /* Set D's */
              d_hat = 0.0;
              d_tilde = 0.0;
            }

            /* VACUUM BC */
            else if (_boundaries[surface] == VACUUM) {

              /* Set D's */
              d_hat =  2 * d*f / length_perpen / (1 + 4 * d*f /
                       length_perpen);
              d_tilde = (sense * d_hat * flux - current / length) / flux;
             }
          }

          /* If Surface is an interface, use finite differencing */
          else{

            /* Set properties for cell next to Surface */
            if (surface == SURFACE_X_MIN) {
              next_length_perpen = _cell_width;
              next_surface = SURFACE_X_MAX;
            }
            else if (surface == SURFACE_Y_MIN) {
              next_length_perpen = _cell_height;
              next_surface = SURFACE_Y_MAX;
            }
            else if (surface == SURFACE_X_MAX) {
              next_length_perpen = _cell_width;
              next_surface = SURFACE_X_MIN;
            }
            else if (surface == SURFACE_Y_MAX) {
              next_length_perpen = _cell_height;
              next_surface = SURFACE_Y_MIN;
            }

            /* Set diffusion coefficient and flux for neighboring cell */
            d_next = _materials[cell_next]->getDifCoef()[e];
            flux_next = _old_flux->getValue(cell_next, e);

            /* Get optical thickness correction term for meshCellNext */
            f_next = computeDiffCorrect(d_next, next_length_perpen);

            /* Compute d_hat */
            d_hat = 2.0 * d * f * d_next * f_next / (length_perpen
                    * d * f + next_length_perpen * d_next*f_next);

            /* Compute net current */
            current = sense * _surface_currents->getValue
              (cell, surface*_num_cmfd_groups + e)
              - sense * _surface_currents->getValue
              (cell_next, next_surface*_num_cmfd_groups + e);

            /* Compute d_tilde */
            d_tilde = -(sense * d_hat * (flux_next - flux) +
                        current  / length) / (flux_next + flux);

            /* If the magnitude of d_tilde is greater than the magnitude of
             * d_hat, select new values d_tilde and d_hat to ensure the course
             * mesh equations are guaranteed to be diagonally dominant */
            if (fabs(d_tilde) > fabs(d_hat) && moc_iteration != 0) {

              if (sense == -1) {

                /* If d_tilde is positive */
                if (1 - fabs(d_tilde)/d_tilde < 1e-8) {
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
                if (1 - fabs(d_tilde)/d_tilde < 1e-8) {
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

          /* Perform underrelaxation on d_tilde. If first MOC iteration, solve 
           * the diffusion problem without correcting currents */
          if (moc_iteration == 0)
            d_tilde = 0.0;
          else
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



/**
 * @brief Solve the nonlinear diffusion acceleration problem to accelerate the
 *        convergence of the MOC problem.
 * @details This method uses the information from the last MOC transport sweep
 *         and solves a simplified nonlinear diffusion problem. The diffusion
 *         problem is tightly converged and the solution is used to update the
 *         the solution of the MOC problem.
 *  @param moc_iteration MOC iteration number
 *  @return The dominant eigenvalue of the nonlinear diffusion problem
 */
FP_PRECISION Cmfd::computeKeff(int moc_iteration) {

  log_printf(INFO, "Running diffusion solver...");

  /* Create matrix and vector objects */
  if (_A == NULL) {
    try{

      /* Allocate memory for matrix and vector objects */
      _M = new Matrix(_num_x, _num_y, _num_cmfd_groups);
      _A = new Matrix(_num_x, _num_y, _num_cmfd_groups);
      _old_source = new Vector(_num_x, _num_y, _num_cmfd_groups);
      _new_source = new Vector(_num_x, _num_y, _num_cmfd_groups);
      _old_flux = new Vector(_num_x, _num_y, _num_cmfd_groups);
      _new_flux = new Vector(_num_x, _num_y, _num_cmfd_groups);
      _volumes = new Vector(_num_x, _num_y, 1);

      /* Initialize flux and materials */
      initializeMaterials();
    }
    catch(std::exception &e) {
      log_printf(ERROR, "Could not allocate memory for the CMFD mesh objects. "
                 "Backtrace:%s", e.what());
    }
  }
  
  /* Compute the cross sections and surface diffusion coefficients */
  computeXS();
  computeDs(moc_iteration);

  /* Construct matrices */
  constructMatrices();

  /* Copy old flux to new flux */
  _old_flux->copyTo(_new_flux);

  /* Solve the eigenvalue problem */
  _k_eff = eigenvalueSolve(_A, _M, _new_flux, _source_convergence_threshold,
                           _SOR_factor);

  /* Rescale the old and new flux */
  rescaleFlux();

  /* Update the MOC flux */
  updateMOCFlux();

  return _k_eff;
}


/**
 * @brief Rescale the initial and converged flux arrays.
 * @details The diffusion problem is a generalized eigenvalue problem and
 *         therefore the solution is independent of flux level. This method
 *         rescales the input flux and converged flux to both have an average
 *         fission source of 1.0 in each group in each cell.
 */
void Cmfd::rescaleFlux() {

  /* Rescale the new and old flux to have an avg source of 1.0 */
  matrixMultiplication(_M, _new_flux, _new_source);
  matrixMultiplication(_M, _old_flux, _old_source);

  _new_flux->scaleByValue(_num_x*_num_y*_num_cmfd_groups /
                          _new_source->getSum());
  _old_flux->scaleByValue(_num_x*_num_y*_num_cmfd_groups /
                          _old_source->getSum());
}


/**
 * @brief Construct the loss + streaming matrix (A) and the fission gain
 *         matrix (M) in preparation for solving the eigenvalue problem.
 * @details This method loops over all mesh cells and energy groups and
 *         accumulates the iteraction and streaming terms into their
 *         approipriate positions in the loss + streaming matrix and
 *         fission gain matrix.
 */
void Cmfd::constructMatrices() {

  log_printf(INFO,"Constructing matrices...");
    
  FP_PRECISION value, volume;
  int cell;
  Material* material;
  
  /* Zero _A and _M matrices */
  _A->clear();
  _M->clear();
  
  /* Loop over cells */
  #pragma omp parallel for private(value, volume, cell, material)
  for (int y = 0; y < _num_y; y++) {
    for (int x = 0; x < _num_x; x++) {

      cell = y*_num_x + x;
      material = _materials[cell];
      volume =_volumes->getValue(cell, 0);

      /* Loop over groups */
      for (int e = 0; e < _num_cmfd_groups; e++) {
          
        /* Absorption term */
        value = material->getSigmaA()[e] * volume;
        _A->incrementValue(cell, e, cell, e, value);
        
        /* Out (1st) and in (2nd) scattering */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          if (e != g) {
            value = material->getSigmaS()[g*_num_cmfd_groups + e] * volume;
            _A->incrementValue(cell, e, cell, e, value);
            value = - material->getSigmaS()[e*_num_cmfd_groups + g] * volume;
            _A->incrementValue(cell, g, cell, e, value);
          }
        }

        /* SURFACE_X_MIN */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[e] + material->getDifTilde()[e])
          * _cell_height;
        _A->incrementValue(cell, e, cell, e, value);        

        /* Set transport term on off diagonal */
        if (x != 0) {
          value = - (material->getDifHat()[SURFACE_X_MIN*_num_cmfd_groups + e]
                     - material->getDifTilde()
                     [SURFACE_X_MIN*_num_cmfd_groups + e]) * _cell_height;
          _A->incrementValue(cell-1, e, cell, e, value);          
        }
        else if (_boundaries[SURFACE_X_MIN] == PERIODIC) {
          value = - (material->getDifHat()[SURFACE_X_MIN*_num_cmfd_groups + e]
                     - material->getDifTilde()
                     [SURFACE_X_MIN*_num_cmfd_groups + e])
            * _cell_height;
          
          _A->incrementValue(cell + (_num_x-1), e, cell, e, value);
        }
        
        /* SURFACE_X_MAX */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[SURFACE_X_MAX*_num_cmfd_groups + e]
                - material->getDifTilde()[SURFACE_X_MAX*_num_cmfd_groups + e])
          * _cell_height;
        _A->incrementValue(cell, e, cell, e, value);

        /* Set transport term on off diagonal */
        if (x != _num_x - 1) {
          value = - (material->getDifHat()[SURFACE_X_MAX*_num_cmfd_groups + e]
                     + material->getDifTilde()
                     [SURFACE_X_MAX*_num_cmfd_groups + e]) * _cell_height;
          _A->incrementValue(cell+1, e, cell, e, value);            
        }
        else if (_boundaries[SURFACE_X_MAX] == PERIODIC) {
          value = - (material->getDifHat()[SURFACE_X_MAX*_num_cmfd_groups + e]
                     + material->getDifTilde()
                     [SURFACE_X_MAX*_num_cmfd_groups + e]) * _cell_height;
          
          _A->incrementValue(cell - (_num_x-1), e, cell, e, value);
        }
        
        /* SURFACE_Y_MIN */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[SURFACE_Y_MIN*_num_cmfd_groups + e]
                 + material->getDifTilde()[SURFACE_Y_MIN*_num_cmfd_groups + e])
          * _cell_width;
        _A->incrementValue(cell, e, cell, e, value);

        /* Set transport term on off diagonal */
        if (y != 0) {
          value = - (material->getDifHat()[SURFACE_Y_MIN*_num_cmfd_groups + e]
                     - material->getDifTilde()
                     [SURFACE_Y_MIN*_num_cmfd_groups + e]) * _cell_width;
          _A->incrementValue(cell-_num_x, e, cell, e, value);
        }
        else if (_boundaries[SURFACE_Y_MIN] == PERIODIC) {
          value = - (material->getDifHat()[SURFACE_Y_MIN*_num_cmfd_groups + e]
                     - material->getDifTilde()
                     [SURFACE_Y_MIN*_num_cmfd_groups + e]) * _cell_width;
          
          _A->incrementValue(cell + _num_x*(_num_y-1), e, cell, e, value);
        }

        /* SURFACE_Y_MAX */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[SURFACE_Y_MAX*_num_cmfd_groups + e]
                 - material->getDifTilde()[SURFACE_Y_MAX*_num_cmfd_groups + e])
          * _cell_width;
        _A->incrementValue(cell, e, cell, e, value);

        /* Set transport term on off diagonal */
        if (y != _num_y - 1) {
          value = - (material->getDifHat()[SURFACE_Y_MAX*_num_cmfd_groups + e]
                     + material->getDifTilde()
                     [SURFACE_Y_MAX*_num_cmfd_groups + e]) * _cell_width;
          _A->incrementValue(cell+_num_x, e, cell, e, value);
        }
        else if (_boundaries[SURFACE_Y_MAX] == PERIODIC) {
          value = - (material->getDifHat()[SURFACE_Y_MAX*_num_cmfd_groups + e]
                     + material->getDifTilde()
                     [SURFACE_Y_MAX*_num_cmfd_groups + e])
            * _cell_width;
          
          _A->incrementValue(cell - _num_x*(_num_y-1), e, cell, e, value);
        }

        /* Source term */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          value = material->getChi()[e] * material->getNuSigmaF()[g]
                  * volume;
          _M->incrementValue(cell, g, cell, e, value);
        }
      }
    }
  }

  log_printf(INFO, "Done constructing matrices...");
}


/**
 * @brief Update the MOC flux in each FSR.
 * @details This method uses the condensed flux from the last MOC transport
 *         sweep and the converged flux from the eigenvalue problem to
 *         update the MOC flux in each FSR.
 */
void Cmfd::updateMOCFlux() {

  log_printf(INFO, "Updating MOC flux...");

  /* Loop over mesh cells */
  #pragma omp parallel for
  for (int i = 0; i < _num_x*_num_y; i++) {

    std::vector<int>::iterator iter;

    /* Loop over CMFD groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {
  
      for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {

        /* Loop over FRSs in mesh cell */
        for (iter = _cell_fsrs.at(i).begin();
             iter != _cell_fsrs.at(i).end(); ++iter) {

          /* Set new flux in FSR */
          _FSR_fluxes[*iter*_num_moc_groups+h] = getUpdateRatio(i,h, *iter)
             * _FSR_fluxes[*iter*_num_moc_groups+h];

          log_printf(DEBUG, "Updating flux in FSR: %i, cell: %i, group: "
                     "%i, ratio: %f", *iter ,i, h,
                     getUpdateRatio(i,h, *iter));
        }
      }
    }
  }
}


/**
 * @brief Compute a diffusion coefficient correction factor for optically
 *        thick regions.
 * @param d Diffusion coefficient before applying correction factor
 * @param h Width of the cell in the direction of interest
 * @return The diffusion coefficient correction factor
 */
FP_PRECISION Cmfd::computeDiffCorrect(FP_PRECISION d, FP_PRECISION h) {

  if (_optically_thick) {

    /* Initialize variables */
    FP_PRECISION alpha, mu, expon;
    FP_PRECISION rho, F;
    rho = 0.0;

    /* Loop over polar angles */
    for (int p = 0; p < _num_polar; p++) {
      mu = cos(asin(_polar_quad->getSinTheta(p)));
      expon = exp(- h / (3 * d * mu));
      alpha = (1 + expon) / (1 - expon) - 2 * (3 * d * mu) / h;
      rho += mu * _polar_quad->getWeight(p) * alpha;
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
 * @param FSR_materials Pointer to FSR_materials array
 */
void Cmfd::setFSRMaterials(Material** FSR_materials) {
  _FSR_materials = FSR_materials;
}


/**
 * @brief Set the pointer to the array of FSR_volumes.
 * @param FSR_volumes Array of FSR volumes
 */
void Cmfd::setFSRVolumes(FP_PRECISION* FSR_volumes) {
  _FSR_volumes = FSR_volumes;
}


/**
 * @brief Set pointer to FSR flux array.
 * @param scalar_flux Pointer to FSR flux array
 */
void Cmfd::setFSRFluxes(FP_PRECISION* scalar_flux) {
  _FSR_fluxes = scalar_flux;
}


/**
 * @brief Set the successive over-relaxation relaxation factor for the
 *        linear solve within the diffusion eigenvalue solve.
 * @param SOR_factor over-relaxation factor
 */
void Cmfd::setSORRelaxationFactor(FP_PRECISION SOR_factor) {

  if (SOR_factor <= 0.0 || SOR_factor >= 2.0)
    log_printf(ERROR, "The successive over-relaxation relaxation factor "
        "must be > 0 and < 2. Input value: %i", SOR_factor);

  _SOR_factor = SOR_factor;
}


/**
 * @brief Get the number of coarse CMFD energy groups.
 * @return The number of CMFD energy groups
 */
int Cmfd::getNumCmfdGroups() {
  return _num_cmfd_groups;
}


/**
 * @brief Get the CMFD group given an MOC group.
 * @param group The MOC energy group
 * @return The CMFD energy group
 */
int Cmfd::getCmfdGroup(int group) {
  return _group_indices_map[group];
}


/**
 * @brief Set the CMFD energy group structure. 
 * @details CMFD does not necessarily need to have the same energy group 
 *          structure as the MOC problem. This function can be used to set 
 *          a sparse energy group structure to speed up the CMFD solve.
 * @param group_indices An array of the CMFD group boundaries
 * @param length_group_indices The length of the group_indices array
 */
void Cmfd::setGroupStructure(int* group_indices, int length_group_indices) {
    
  _num_cmfd_groups = length_group_indices - 1;

  /* Allocate memory */
  if (_group_indices == NULL) {
    _group_indices = new int[length_group_indices];
  }
    
  if (group_indices == NULL) {
    for (int i = 0; i < length_group_indices; i++) {
      _group_indices[i] = i;
    }
  }
  else{
    if (group_indices[0] != 1)
      log_printf(ERROR, "The first value in group indices must be 1!");    

    /* Set first group indice to 0 */
    _group_indices[0] = 0;
        
    /* Set MOC group bounds for rest of CMFD energy groups */
    for (int i = 1; i < length_group_indices; i++) {
      /* Check that the group indices are always increasing */
      if (group_indices[i] <= group_indices[i-1])
        log_printf(ERROR, "The group indices must be increasing!");
            
      _group_indices[i] = group_indices[i] - 1;
      log_printf(INFO, "group indices %i: %i", i, group_indices[i]);
    }
  }
}


/**
 * @brief Initialize the CMFD materials.
 */
void Cmfd::initializeMaterials() {

  Material* material;

  try{
    _materials = new Material*[_num_x*_num_y];

    for (int y = 0; y < _num_y; y++) {
      for (int x = 0; x < _num_x; x++) {
        material = new Material(y*_num_x+x);
        material->setNumEnergyGroups(_num_cmfd_groups);
        _materials[y*_num_x+x] = material;
      }
    }
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the Mesh cell materials. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initializes Cmfd surface currents Vector prior to first MOC iteration.
 */
void Cmfd::initializeSurfaceCurrents() {

  /* Delete old Cmfd surface currents array it it exists */
  if (_surface_currents != NULL)
    delete [] _surface_currents;

  /* Allocate memory for the Cmfd Mesh surface currents array */
  int num_mesh_cells = _num_x * _num_y;
  int size = num_mesh_cells * _num_cmfd_groups * NUM_SURFACES;
  _surface_currents = new Vector(_num_x, _num_y,
                                 _num_cmfd_groups * NUM_SURFACES);

  return;
}


/**
 * @brief Initializes the vector of vectors that links CMFD cells with FSRs.
 * @details This method is called by the geometry once the CMFD mesh has been
 *          initialized by the geometry. This method allocates a vector for
 *          each CMFD cell that is used to store the FSR ids contained within
 *          that cell.
 */
void Cmfd::initializeCellMap() {

  /* Allocate memory for mesh cell FSR vectors */
  for (int y = 0; y < _num_y; y++) {
    for (int x = 0; x < _num_x; x++) {
      std::vector<int> *fsrs = new std::vector<int>;
      _cell_fsrs.push_back(*fsrs);
    }
  }
}


/**
 * @brief Initialize and set array that links the MOC energy groups to the
 *        CMFD energy groups.
 * @details This method initializes the _group_indices_map, which is a 1D array
 *          of length _num_moc_groups that maps the MOC energy groups to CMFD 
 *          energy groups. The indices into _group_indices_map are the MOC
 *          energy groups and the values are the CMFD energy groups.
 */
void Cmfd::initializeGroupMap() {

  /* Allocate memory */
  if (_group_indices_map == NULL) {
    _group_indices_map = new int[_num_moc_groups];
  }    
    
  /* Create group indices map */
  for (int e = 0; e < _num_cmfd_groups; e++) {
    for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {
      _group_indices_map[h] = e;
    }
  }
}



/**
 * @brief Find the cmfd surface that a LocalCoords object lies on.
 * @details If the coords is not on a surface, -1 is returned. Otherwise,
 *        the surface ID is returned. 
 * @param cell The CMFD cell ID that the local coords is in.
 * @param coords The coords being evaluated.
 * @return The surface ID.
 */
int Cmfd::findCmfdSurface(int cell, LocalCoords* coords) {
  Point* point = coords->getHighestLevel()->getPoint();
  return _lattice->getLatticeSurface(cell, point);
}


/**
 * @brief Find the CMFD cell that a LocalCoords object is in. 
 * @param coords The coords being evaluated.
 * @return The CMFD cell ID.
 */
int Cmfd::findCmfdCell(LocalCoords* coords) {
  Point* point = coords->getHighestLevel()->getPoint();
  return _lattice->getLatticeCell(point);
}


/**
 * @brief The Lattice object used as the CMFD mesh. 
 * @param lattice Pointer to the lattice object.
 */
void Cmfd::setLattice(Lattice* lattice) {
    _lattice = lattice;
}


/**
 * @brief The structure of the Lattice to be used as the CMFD mesh.
 * @param num_x The number of cells in the x direction.
 * @param num_y The number of cells in the y direction.
 */
void Cmfd::setLatticeStructure(int num_x, int num_y) {
  setNumX(num_x);
  setNumY(num_y);
}


/**
 * @brief Returns the Lattice object used as the CMFD mesh. 
 * @return A pointer to a Lattice object.
 */
Lattice* Cmfd::getLattice() {
  return _lattice;
}


/**
 * @brief Add an FSR ID to a vector that contains all the FSR IDs
 *        contained within a CMFD mesh cell.
 * @param cmfd_cell The CMFD cell ID.
 * @param fsr_id The FSR ID.
 */
void Cmfd::addFSRToCell(int cmfd_cell, int fsr_id) {
  _cell_fsrs.at(cmfd_cell).push_back(fsr_id);
}


/**
 * @brief Set the number of MOC energy groups.
 * @param num_groups number of MOC energy groups
 */
void Cmfd::setNumMOCGroups(int num_groups) {
  _num_moc_groups = num_groups;
}


/**
 * @brief Get the number of MOC energy groups.
 * @return The number of MOC energy groups
 */
int Cmfd::getNumMOCGroups() {
  return _num_moc_groups;
}


/**
 * @brief Get the number of CMFD cells.
 * @return The number of CMFD cells
 */
int Cmfd::getNumCells() {
  return _num_x*_num_y;
}


/**
 * @brief Set the number of FSRs.
 * @param num_fsrs The number of FSRs
 */
void Cmfd::setNumFSRs(int num_fsrs) {
  _num_FSRs = num_fsrs;
}


/** @brief Split the currents of the Mesh cell corners to the nearby surfaces.
 * @details This metohd splits the current from each corner evenly to the two
 *          adjoining surfaces.
 */
void Cmfd::splitCorners() {

  log_printf(INFO, "splitting corners...");
    
  FP_PRECISION current;

  #pragma omp parallel for private(current)
  for (int y = 0; y < _num_y; y++) {
    for (int x = 0; x < _num_x; x++) {
      int cell = y*_num_x + x;
      for (int s = 4; s < NUM_SURFACES; s++){

        /* SURFACE_X_MIN_Y_MIN */
        if (s == SURFACE_X_MIN_Y_MIN){
          for (int e=0; e > _num_cmfd_groups; e++){
            current = 0.5*_surface_currents->getValue
              (cell, s*_num_cmfd_groups + e);
            
            _surface_currents->incrementValue
              (cell, SURFACE_X_MIN*_num_cmfd_groups + e, current);
            _surface_currents->incrementValue
              (cell, SURFACE_Y_MIN*_num_cmfd_groups + e, current);
            
            if (x != 0)
              _surface_currents->incrementValue
                (cell-1, SURFACE_Y_MIN*_num_cmfd_groups + e, current);
            if (y != 0)
              _surface_currents->incrementValue
                (cell-_num_x, SURFACE_X_MIN*_num_cmfd_groups + e, current);
          }
        }
        
        /* SURFACE_X_MAX_Y_MIN */
        else if (s == SURFACE_X_MAX_Y_MIN){
          for (int e=0; e > _num_cmfd_groups; e++){
            current = 0.5*_surface_currents->getValue
              (cell, s*_num_cmfd_groups + e);
            
            _surface_currents->incrementValue
              (cell, SURFACE_X_MAX*_num_cmfd_groups + e, current);
            _surface_currents->incrementValue
              (cell, SURFACE_Y_MIN*_num_cmfd_groups + e, current);
            
            if (x != _num_x - 1)
              _surface_currents->incrementValue
                (cell+1, SURFACE_Y_MIN*_num_cmfd_groups + e, current);
            if (y != 0)
              _surface_currents->incrementValue
                (cell-_num_x, SURFACE_X_MAX*_num_cmfd_groups + e, current);
          }
        }
        
        /* SURFACE_X_MIN_Y_MAX */
        if (s == SURFACE_X_MIN_Y_MAX){
          for (int e=0; e > _num_cmfd_groups; e++){
            current = 0.5*_surface_currents->getValue
              (cell, s*_num_cmfd_groups + e);
            
            _surface_currents->incrementValue
              (cell, SURFACE_X_MIN*_num_cmfd_groups + e, current);
            _surface_currents->incrementValue
              (cell, SURFACE_Y_MAX*_num_cmfd_groups + e, current);
            
            if (x != 0)
              _surface_currents->incrementValue
                (cell-1, SURFACE_Y_MAX*_num_cmfd_groups + e, current);
            if (y != _num_y - 1)
              _surface_currents->incrementValue
                (cell+_num_x, SURFACE_X_MIN*_num_cmfd_groups + e, current);
          }
        }
        
        /* SURFACE_X_MAX_Y_MAX */
        if (s == SURFACE_X_MAX_Y_MAX){
          for (int e=0; e > _num_cmfd_groups; e++){
            current = 0.5*_surface_currents->getValue
              (cell, s*_num_cmfd_groups + e);
            
            _surface_currents->incrementValue
              (cell, SURFACE_X_MAX*_num_cmfd_groups + e, current);
            _surface_currents->incrementValue
              (cell, SURFACE_Y_MAX*_num_cmfd_groups + e, current);
            
            if (x != _num_x - 1)
              _surface_currents->incrementValue
                (cell+1, SURFACE_Y_MAX*_num_cmfd_groups + e, current);
            if (y != _num_y - 1)
              _surface_currents->incrementValue
                (cell+_num_x, SURFACE_X_MAX*_num_cmfd_groups + e, current);
          }
        }

        /* Zero out corner currents */
        for (int e=0; e > _num_cmfd_groups; e++)
          _surface_currents->setValue(cell, s*_num_cmfd_groups + e, 0.0);
      }
    }
  }
}


/**
 * @brief Get the ID of the Mesh cell next to a given Mesh cell across a
 *        given surface.
 * @param cell_num Current Mesh cell ID
 * @param surface_id CMFD cell surface ID to look across for neighboring cell
 * @return Neighboring CMFD cell ID
 */
int Cmfd::getCellNext(int cell_num, int surface_id) {

  int cell_next = -1;

  if (surface_id == SURFACE_X_MIN) {
    if (cell_num % _num_x != 0)
      cell_next = cell_num - 1;
    else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
      cell_next = cell_num + (_num_x-1);
  }
  else if (surface_id == SURFACE_Y_MIN) {
    if (cell_num / _num_x != 0)
      cell_next = cell_num - _num_x;
    else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
      cell_next = cell_num + _num_x*(_num_y-1);
  }
  else if (surface_id == SURFACE_X_MAX) {
    if (cell_num % _num_x != _num_x - 1)
      cell_next = cell_num + 1;
    else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
      cell_next = cell_num - (_num_x-1);
  }
  else if (surface_id == SURFACE_Y_MAX) {
    if (cell_num / _num_x != _num_y - 1)
      cell_next = cell_num + _num_x;
    else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
      cell_next = cell_num - _num_x*(_num_y-1);
  }

  return cell_next;
}


/**
 * @brief Return whether optically thick diffusion correction factor is in use.
 * @return Whether optically thick diffusion correction factor is in use.
 */
bool Cmfd::isOpticallyThick() {
  return _optically_thick;
}


/**
 * @brief Set whether optically thick diffusion correction factor is in use.
 * @param optically_thick Boolean indicating whether optically thick diffusion
 *        correction factor is in use.
 */
void Cmfd::setOpticallyThick(bool optically_thick) {
  _optically_thick = optically_thick;
}


/**
 * @brief Return the under-relaxation factor used in MOC updates.
 * @return The MOC current under-relaxation factor
 */
FP_PRECISION Cmfd::getMOCRelaxationFactor() {
  return _relax_factor;
}


/**
 * @brief Set the under-relaxation factor used in MOC updates.
 * @param relax_factor The MOC current under-relaxation factor
 */
void Cmfd::setMOCRelaxationFactor(FP_PRECISION relax_factor) {
  _relax_factor = relax_factor;
}


/**
 * @brief Set the CMFD boundary type for a given surface.
 * @details The CMFD boundary is assumed to be rectangular with the
 *          surfaces identified by constants in the constants.h file.
 * @param side The CMFD surface UID.
 * @param boundary The boundaryType of the surface.
 */
void Cmfd::setBoundary(int side, boundaryType boundary) {
  _boundaries[side] = boundary;
}


/**
 * @brief Get the boundaryType for one side of the CMFD mesh.
 * @param side The CMFD mesh surface ID.
 * @return The boundaryType for the surface.
 */
int Cmfd::getBoundary(int side) {
  return _boundaries[side];
}


/**
 * @brief Return the CMFD cell ID that an FSR lies in.
 * @details Note that a CMFD cell is not an actual Cell object; rather, a CMFD
 *         cell is just a way of describing each of the rectangular regions
 *         that make up a CMFD lattice. CMFD cells are numbered with 0 in the
 *         lower left corner and monotonically increasing from left to right
 *         and from bottom to top. For example, the indices for a 4 x 4
 *         lattice are:
 *                  12  13  14  15
 *                  8    9  10  11
 *                  4    5   6   7
 *                  0    1   2   3
 * @param fsr_id The FSR ID.
 * @return The CMFD cell ID. Return -1 if cell is not found.
 */
int Cmfd::convertFSRIdToCmfdCell(int fsr_id) {

  std::vector<int>::iterator iter;    
  for (int cell=0; cell < _num_x*_num_y; cell++) {

    for (iter = _cell_fsrs.at(cell).begin();
         iter != _cell_fsrs.at(cell).end(); ++iter) {
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
std::vector< std::vector<int> >* Cmfd::getCellFSRs() {
  return &_cell_fsrs;
}

 
/**
 * @brief Set the vector of vectors that contains the FSRs that lie in
 *        each cell.
 * @param cell_fsrs Vector of vectors containing FSR IDs in each cell.
 */
void Cmfd::setCellFSRs(std::vector< std::vector<int> >* cell_fsrs) {

  if (!_cell_fsrs.empty()) {
    std::vector< std::vector<int> >::iterator iter;
    for (iter = _cell_fsrs.begin(); iter != _cell_fsrs.end(); ++iter)
      iter->clear();
    _cell_fsrs.clear();
  }

  _cell_fsrs = *cell_fsrs;
}


/**
 * @brief Set flag indicating whether to update the MOC flux.
 * @param flux_update_on Boolean saying whether to update MOC flux.
 */
void Cmfd::setFluxUpdateOn(bool flux_update_on) {
  _flux_update_on = flux_update_on;
}


/**
 * @brief Get flag indicating whether to update the MOC flux.
 * @return Boolean saying whether to update MOC flux.
 */
bool Cmfd::isFluxUpdateOn() {
 return _flux_update_on;
}


/**
 * @brief Set flag indicating whether to use FSR centroids to update
 *        the MOC flux.
 * @param centroid_update_on Flag saying whether to use centroids to
 *        update MOC flux.
 */
void Cmfd::setCentroidUpdateOn(bool centroid_update_on){
  _centroid_update_on = centroid_update_on;
}


/**
 * @brief Get flag indicating whether to use FSR centroids to update
 *        the MOC flux.
 * @return Flag saying whether to use centroids to update MOC flux.
 */
bool Cmfd::isCentroidUpdateOn(){
 return _centroid_update_on;
}


/**
 * @brief Sets the threshold for CMFD source convergence (>0)
 * @param the threshold for source convergence
 */
void Cmfd::setSourceConvergenceThreshold(FP_PRECISION source_thresh) {

  if (source_thresh <= 0.0)
    log_printf(ERROR, "Unable to set the cmfd source convergence threshold to"
              " %f since the threshold must be positive.", source_thresh);

  _source_convergence_threshold = source_thresh;
}


/**
 * @brief Sets the PolarQuad object in use by the MOC Solver.
 * @param polar_quad A PolarQuad object pointer from the Solver
 */
void Cmfd::setPolarQuadrature(PolarQuad* polar_quad) {

  /* Deletes the old Quadrature if one existed */
  if (_polar_quad != NULL)
    delete _polar_quad;

  _polar_quad = polar_quad;
  _num_polar = polar_quad->getNumPolarAngles();
}


/**
 * @brief Get the ratio of new to old CMFD cell flux in a given CMFD cell
 *        and CMFD energy group containing a given MOC energy group.
 * @param cmfd_cell The CMFD cell of interest.
 * @param moc_group The MOC energy group of interest.
 * @return The ratio of new to old CMFD cell flux.
 */
FP_PRECISION Cmfd::getFluxRatio(int cmfd_cell, int moc_group) {

  int cmfd_group = _group_indices_map[moc_group];
  FP_PRECISION old_flux = _old_flux->getValue(cmfd_cell, cmfd_group);
  FP_PRECISION new_flux = _new_flux->getValue(cmfd_cell, cmfd_group);
  return new_flux / old_flux;
}


/**
 * @brief Generate the k-nearest neighbor CMFD cell stencil for each FSR.
 * @detail This method finds the k-nearest CMFD cell stencil for each FSR
 *         and saves the stencil, ordered from the closest-to-furthest
 *         CMFD cell, in the _k_nearest_stencils map. The stencil of cells
 *         surrounding the current cell is defined as:
 *
 *                             6 7 8
 *                             3 4 5
 *                             0 1 2
 *
 *         where 4 is the given CMFD cell. If the cell is on the edge or corner
 *         of the geometry and there are less than k nearest neighbor cells,
 *         k is reduced to the number of neighbor cells for that instance.
 */
void Cmfd::generateKNearestStencils(){

  std::vector< std::pair<int, FP_PRECISION> >::iterator iter2;
  std::vector<int>::iterator iter;
  Point* centroid;

  /* Loop over mesh cells */
  for (int i = 0; i < _num_x*_num_y; i++){

    /* Loop over FRSs in mesh cell */
    for (iter = _cell_fsrs.at(i).begin();
         iter != _cell_fsrs.at(i).end(); ++iter) {

      /* Get centroid */
      centroid = _geometry->getFSRCentroid(*iter);

      /* Create new stencil */
      std::vector< std::pair<int, FP_PRECISION> > *stencil =
        new std::vector< std::pair<int, FP_PRECISION> >;
      _k_nearest_stencils[*iter] = (*stencil);

      /* Get distance to all cells that touch current cell */
      for (int j=0; j <= NUM_SURFACES; j++)
        _k_nearest_stencils[*iter]
          .push_back(std::make_pair<int, FP_PRECISION>
                     (int(j), getDistanceToCentroid(centroid, i, j)));

      /* Sort the distances */
      std::sort(_k_nearest_stencils[*iter].begin(),
                _k_nearest_stencils[*iter].end(), stencilCompare);

      /* Remove non-existent cells */
      iter2 = _k_nearest_stencils[*iter].begin();
      while (iter2 != _k_nearest_stencils[*iter].end()){
        if (iter2->second == std::numeric_limits<FP_PRECISION>::max())
          iter2 = _k_nearest_stencils[*iter].erase(iter2++);
        else
          ++iter2;
      }

      /* Resize stencil to be of size <= _k_nearest */
      _k_nearest_stencils[*iter].resize
        (std::min(_k_nearest, int(_k_nearest_stencils[*iter].size())));
    }
  }
}


/**
 * @brief Get the ratio used to update the FSR flux after converging CMFD.
 * @detail This method takes in a cmfd cell, a MOC energy group, and a FSR
 *         and returns the ratio used to update the FSR flux. There are two
 *         methods that can be used to update the flux, conventional and
 *         k-nearest centroid updating. The k-nearest centroid updating uses
 *         the k-nearest cells (with k between 1 and 9) of the current CMFD
 *         cell and the 8 neighboring CMFD cells. The stencil of cells
 *         surrounding the current cell is defined as:
 *
 *                             6 7 8
 *                             3 4 5
 *                             0 1 2
 *
 *         where 4 is the given CMFD cell. If the cell is on the edge or corner
 *         of the geometry and there are less than k nearest neighbor cells,
 *         k is reduced to the number of neighbor cells for that instance.
 * @param cmfd_cell The cmfd cell containing the FSR.
 * @param moc_group The MOC energy group being updated.
 * @param fsr The fsr being updated.
 * @return the ratio used to update the FSR flux.
 */
FP_PRECISION Cmfd::getUpdateRatio(int cmfd_cell, int moc_group, int fsr){

  FP_PRECISION ratio = 0.0;
  FP_PRECISION total_distance = 1.e-10;
  std::vector< std::pair<int, FP_PRECISION> >::iterator iter;

  if (_centroid_update_on){

    /* Compute the total distance for the stencil */
    for (iter = _k_nearest_stencils[fsr].begin();
         iter < _k_nearest_stencils[fsr].end(); ++iter)
      total_distance += iter->second;

    /* Compute the ratio */
    for (iter = _k_nearest_stencils[fsr].begin();
         iter != _k_nearest_stencils[fsr].end(); ++iter){

      /* LOWER LEFT CORNER */
      if (iter->first == 0)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell - _num_x - 1, moc_group);

      /* BOTTOM SIDE */
      else if (iter->first == 1)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell - _num_x, moc_group);

      /* LOWER RIGHT CORNER */
      else if (iter->first == 2)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell - _num_x + 1, moc_group);

      /* LEFT SIDE */
      else if (iter->first == 3)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell - 1, moc_group);

      /* RIGHT SIDE */
      else if (iter->first == 5)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell + 1, moc_group);

      /* UPPER LEFT CORNER */
      else if (iter->first == 6)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell + _num_x - 1, moc_group);

      /* TOP SIDE */
      else if (iter->first == 7)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell + _num_x, moc_group);

      /* TOP RIGHT CORNER */
      else if (iter->first == 8)
        ratio += (1.0 - iter->second/total_distance) *
          getFluxRatio(cmfd_cell + _num_x + 1, moc_group);
    }

    /* INTERNAL */
    if (_k_nearest_stencils[fsr].size() == 1)
      ratio += getFluxRatio(cmfd_cell, moc_group);
    else{
      ratio += (1.0 - _k_nearest_stencils[fsr][0].second/total_distance) *
        getFluxRatio(cmfd_cell, moc_group);
      ratio /= (_k_nearest_stencils[fsr].size() - 1);
    }
  }
  else
    ratio = getFluxRatio(cmfd_cell, moc_group);

  return ratio;
}


/**
 * @brief Get the distances from an FSR centroid to a given cmfd cell.
 * @detail This method takes in a FSR centroid, a cmfd cell, and a stencil index
 *         to a cell located in the 9-point stencil encompassing the cmfd
 *         cell an all its possible neighbors. The CMFD cell stencil is:
 *
 *                             6 7 8
 *                             3 4 5
 *                             0 1 2
 *
 *         where 4 is the given CMFD cell. If a CMFD edge or corner cells is
 *         given and the stencil indexed cell lies outside the geometry, the
 *         maximum allowable FP_PRECISION value is returned.
 * @param centroid The numerical centroid an FSR in the cell.
 * @param cell The cmfd cell containing the FSR.
 * @param stencil_index The index of the cell in the stencil that we want to
 *        get the distance from.
 * @return the distance from the CMFD cell centroid to the FSR centroid.
 */
FP_PRECISION Cmfd::getDistanceToCentroid(Point* centroid, int cell,
                                         int stencil_index){

  int x = cell % _num_x;
  int y = cell / _num_x;
  FP_PRECISION dist_x, dist_y;
  bool found = false;

  /* LOWER LEFT CORNER */
  if (x > 0 && y > 0 && stencil_index == 0){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x - 0.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y - 0.5)*_cell_height), 2.0);
    found = true;
  }
  
  /* BOTTOM SIDE */
  else if (y > 0 && stencil_index == 1){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x + 0.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y - 0.5)*_cell_height), 2.0);
    found = true;
  }

  /* LOWER RIGHT CORNER */
  else if (x < _num_x - 1 && y > 0 && stencil_index == 2){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x + 1.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y - 0.5)*_cell_height), 2.0);
    found = true;
  }

  /* LEFT SIDE */
  else if (x > 0 && stencil_index == 3){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x - 0.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y + 0.5)*_cell_height), 2.0);
    found = true;
  }

  /* CURRENT */
  else if (stencil_index == 4){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x + 0.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y + 0.5)*_cell_height), 2.0);
    found = true;
  }

  /* RIGHT SIDE */
  else if (x < _num_x - 1 && stencil_index == 5){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x + 1.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y + 0.5)*_cell_height), 2.0);
    found = true;
  }

  /* UPPER LEFT CORNER */
  else if (x > 0 && y < _num_y - 1 && stencil_index == 6){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x - 0.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y + 1.5)*_cell_height), 2.0);
    found = true;
  }

  /* TOP SIDE */
  else if (y < _num_y - 1 && stencil_index == 7){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x + 0.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y + 1.5)*_cell_height), 2.0);
    found = true;
  }

  /* UPPER RIGHT CORNER */
  else if (x < _num_x - 1 && y < _num_y - 1 && stencil_index == 8){
    dist_x = pow(centroid->getX() - (-_width/2.0+(x + 1.5)*_cell_width), 2.0);
    dist_y = pow(centroid->getY() - (-_height/2.0+(y + 1.5)*_cell_height), 2.0);
    found = true;
  }

  if (found)
    return pow(dist_x + dist_y, 0.5);
  else
    return std::numeric_limits<FP_PRECISION>::max();
}


/**
 * @brief Update the MOC boundary fluxes.
 * @details The MOC boundary fluxes are updated using the P0 approximation.
 *          With this approximation, the boundary fluxes are updated using
 *          the ratio of new to old flux for the cell that the outgoing flux
 *          from the track enters.
 * @param tracks 2D array of Tracks
 * @param boundary_flux Array of boundary fluxes
 * @return The number of Tracks
 */
void Cmfd::updateBoundaryFlux(Track** tracks, FP_PRECISION* boundary_flux, 
			      int num_tracks) {

  segment* segments;
  segment* curr_segment;
  int num_segments;
  int bc;
  FP_PRECISION* track_flux;
  FP_PRECISION ratio;
  int cmfd_cell;
  
  log_printf(INFO, "updating boundary flux");

  /* Loop over Tracks */
  for (int i=0; i < num_tracks; i++) {
      
    num_segments = tracks[i]->getNumSegments();
    segments = tracks[i]->getSegments();

    /* Update boundary flux in forward direction */
    bc = (int)tracks[i]->getBCOut();
    curr_segment = &segments[0];
    track_flux = &boundary_flux[i*2*_num_moc_groups*_num_polar];
    cmfd_cell = convertFSRIdToCmfdCell(curr_segment->_region_id);
    
    if (bc) {
      for (int e=0; e < _num_moc_groups; e++) {
        for (int p=0; p < _num_polar; p++) {
          track_flux[p*_num_moc_groups+e] *= getFluxRatio(cmfd_cell, e);
        }
      }
    }

    /* Update boundary flux in backwards direction */
    bc = (int)tracks[i]->getBCIn();
    curr_segment = &segments[num_segments-1];
    track_flux = &boundary_flux[(i*2 + 1)*_num_moc_groups*_num_polar];
    
    if (bc) {
      for (int e=0; e < _num_moc_groups; e++) {
        for (int p=0; p < _num_polar; p++) {
          track_flux[p*_num_moc_groups+e] *= getFluxRatio(cmfd_cell, e);
        }
      }
    }
  }
}


/** @brief Set a pointer to the Geometry.
 * @param goemetry A pointer to a Geometry object.
 */
void Cmfd::setGeometry(Geometry* geometry){
  _geometry = geometry;
}


/** @brief Set a number of k-nearest neighbor cells to use in updating
 *         the FSR flux.
 * @param k_nearest The number of nearest neighbor CMFD cells.
 */
void Cmfd::setKNearest(int k_nearest){

  if (_k_nearest < 1 || k_nearest > 9)
    log_printf(ERROR, "Unable to set CMFD k-nearest to %i. k-nearest "
               "must be between 1 and 9.", k_nearest);
  else
    _k_nearest = k_nearest;
}


/**
 * @brief Zero the surface currents for each mesh cell and energy group.
 */
void Cmfd::zeroSurfaceCurrents() {
  _surface_currents->clear();
}


/**
 * @brief Tallies the current contribution from this segment across the
 *        the appropriate CMFD mesh cell surface.
 * @param curr_segment the current Track segment
 * @param track_flux the outgoing angular flux for this segment
 * @param polar_weights array of polar weights for some azimuthal angle
 * @param fwd boolean indicating direction of integration along segment
 */
void Cmfd::tallySurfaceCurrent(segment* curr_segment, FP_PRECISION* track_flux, 
                               FP_PRECISION* polar_weights, bool fwd) {

  FP_PRECISION surf_current;
  int surf_id;

  if (curr_segment->_cmfd_surface_fwd != -1 && fwd) {

    surf_id = curr_segment->_cmfd_surface_fwd % NUM_SURFACES;
    int cell_id = curr_segment->_cmfd_surface_fwd / NUM_SURFACES;
    
    for (int e=0; e < _num_moc_groups; e++) {
      surf_current = 0.;

      int g = getCmfdGroup(e);

      for (int p=0; p < _num_polar; p++)
        surf_current += track_flux(p,e) * polar_weights[p] / 2.;

      /* Increment current (polar and azimuthal weighted flux, group) */
      _surface_currents->incrementValue
        (cell_id, surf_id*_num_cmfd_groups + g, surf_current);
    }
  }
  else if (curr_segment->_cmfd_surface_bwd != -1 && !fwd) {

    surf_id = curr_segment->_cmfd_surface_bwd % NUM_SURFACES;
    int cell_id = curr_segment->_cmfd_surface_bwd / NUM_SURFACES;
    
    for (int e=0; e < _num_moc_groups; e++) {
      surf_current = 0.;

      int g = getCmfdGroup(e);
      
      for (int p=0; p < _num_polar; p++)
        surf_current += track_flux(p,e) * polar_weights[p] / 2.;

      /* Increment current (polar and azimuthal weighted flux, group) */
      _surface_currents->incrementValue
        (cell_id, surf_id*_num_cmfd_groups + g, surf_current);
    }
  }
}
