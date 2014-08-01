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
Cmfd::Cmfd(Geometry* geometry, double criteria) {

  /* Initialize Geometry and Mesh-related attribute */
  _geometry = geometry;
  _mesh = geometry->getMesh();
  _num_FSRs = _mesh->getNumFSRs();
  _cx = _mesh->getCellsX();
  _cy = _mesh->getCellsY();
  _quad = new Quadrature(TABUCHI);
  _timer = new Timer();
  _omega = 1.0;

  /* Boolean and Enum flags to toggle features */
  _solve_method = _mesh->getSolveType();
  _flux_type = PRIMAL;
  _eigen_method = POWER;

  /* Global variables used in solving CMFD problem */
  _l2_norm = 1.0;
  _conv_criteria = criteria;

  /* Energy group problem parameters */
  _num_groups = _geometry->getNumEnergyGroups();
  _num_cmfd_groups = 0;

  /* matrices */
  _A = NULL;
  _M = NULL;
  _AM = NULL;
  _phi_temp = NULL;
  _old_source = NULL;
  _new_source = NULL;
  _group_indices = NULL;
  _group_indices_map = NULL;


  /* If solving diffusion problem, create arrays for FSR parameters */
  if (_solve_method == DIFFUSION){
    try{
      _FSR_volumes = (FP_PRECISION*)calloc(_num_FSRs, sizeof(FP_PRECISION));
      _FSR_materials = new Material*[_num_FSRs];
      _FSR_fluxes =
           (FP_PRECISION*)calloc(_num_FSRs*_num_groups, sizeof(FP_PRECISION));
    }
    catch(std::exception &e){
      log_printf(ERROR, "Could not allocate memory for the CMFD diffusion "
                 "Mesh objects. Backtrace:%s", e.what());
    }

    _mesh->initializeSurfaceCurrents();
    initializeFSRs();
  }
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

  if (_AM != NULL){
    for (int i = 0; i < _cx*_cy; i++)
      delete [] _AM[i];

    delete [] _AM;
  }

  if (_phi_temp != NULL)
    delete [] _phi_temp;

  if (_old_source != NULL)
    delete [] _old_source;

  if (_new_source != NULL)
    delete [] _new_source;
}


/**
 * @brief Create cross-sections and fluxes for each Cmfd cell by
 *        energy condensing and volume averaging cross sections from
 *        the MOC sweep.
 */
void Cmfd::computeXS(){

  log_printf(INFO, "Computing CMFD cross-sections...");

  /* Split corner currents to side surfaces */
  _mesh->splitCorners();

  /* Initialize variables for FSR properties*/
  double volume, flux, abs, tot, nu_fis, chi, dif_coef;
  FP_PRECISION* scat;
  Material** materials = _mesh->getMaterials();
  double* fluxes = _mesh->getFluxes(PRIMAL);

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

    cell_material = materials[i];
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
      for (iter = _mesh->getCellFSRs()->at(i).begin();
        iter != _mesh->getCellFSRs()->at(i).end(); ++iter){

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
      _mesh->setVolume(vol_tally, i);
      cell_material->setSigmaAByGroup(abs_tally / rxn_tally, e+1);
      cell_material->setSigmaTByGroup(tot_tally / rxn_tally, e+1);
      cell_material->setNuSigmaFByGroup(nu_fis_tally / rxn_tally, e+1);
      cell_material->setDifCoefByGroup(dif_tally / rxn_tally, e+1);
      fluxes[i*_num_cmfd_groups+e] = rxn_tally / vol_tally;

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
        cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, g+1, e+1);
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

  Material** materials = _mesh->getMaterials();
  double* cell_flux = _mesh->getFluxes(PRIMAL);
  double* lengths_y = _mesh->getLengthsY();
  double* lengths_x = _mesh->getLengthsX();
  double* currents = _mesh->getCurrents();

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
          d = materials[cell]->getDifCoef()[e];
          flux = cell_flux[cell*_num_cmfd_groups+e];
          cell_next = _mesh->getCellNext(cell, surface);

          /* Set halfspace sense of the Surface */
          if (surface == 0 || surface == 3)
            sense = -1.0;
          else
            sense = 1.0;

          /* Set the length of this Surface and the perpendicular Surface */
          if (surface == 0 || surface== 2){
            length = lengths_y[y];
            length_perpen = lengths_x[x];
          }
          else if (surface == 1 || surface == 3){
            length = lengths_x[x];
            length_perpen = lengths_y[y];
          }

          /* Compute the optical thickness correction factor */
          f = computeDiffCorrect(d, length_perpen);

          /* If Surface is on a boundary, choose appropriate BCs */
          if (cell_next == -1){

            current = sense * currents[cell*_num_cmfd_groups*8 +
                                       surface*_num_cmfd_groups + e];

            /* REFLECTIVE BC */
            if (_mesh->getBoundary(surface) == REFLECTIVE){

              /* Set D's */
              d_hat = 0.0;
              d_tilde = 0.0;
            }

            /* VACUUM BC */
            else if (_mesh->getBoundary(surface) == VACUUM){

              /* Set D's */
              d_hat =  2 * d*f / length_perpen / (1 + 4 * d*f /
                       length_perpen);

              if (_solve_method == MOC)
                d_tilde = (sense * d_hat * flux - current / length) / flux;
              else
                d_tilde = 0.0;
             }

            /* ZERO_FLUX BC */
            else if (_mesh->getBoundary(surface) == ZERO_FLUX){

              /* Set D's */
              d_hat = 2 * d*f / length_perpen;
              d_tilde = 0.0;
            }
          }

          /* If Surface is an interface, use finite differencing */
          else{

            /* Set properties for cell next to Surface */
            if (surface == 0){
              next_length_perpen = lengths_x[cell_next % _cx];
              next_surface = 2;
            }
            else if (surface == 1){
              next_length_perpen = lengths_y[cell_next / _cx];
              next_surface = 3;
            }
            else if (surface == 2){
              next_length_perpen = lengths_x[cell_next % _cx];
              next_surface = 0;
            }
            else if (surface == 3){
              next_length_perpen = lengths_y[cell_next / _cx];
              next_surface = 1;
            }

            /* Set diffusion coefficient and flux for neighboring cell */
            d_next = materials[cell_next]->getDifCoef()[e];
            flux_next = cell_flux[cell_next*_num_cmfd_groups + e];

            /* Get optical thickness correction term for meshCellNext */
            f_next = computeDiffCorrect(d_next, next_length_perpen);

            /* Compute d_hat */
            d_hat = 2.0 * d * f * d_next * f_next / (length_perpen
                    * d * f + next_length_perpen * d_next*f_next);

            /* Compute net current */
            current = sense * currents[cell*_num_cmfd_groups*8 +
                      surface*_num_cmfd_groups + e] - sense
                      * currents[cell_next*_num_cmfd_groups*8 +
                      next_surface*_num_cmfd_groups + e];

            /* Compute d_tilde */
            if (_solve_method == MOC)
              d_tilde = -(sense * d_hat * (flux_next - flux) +
                        current  / length) / (flux_next + flux);
            else
              d_tilde = 0.0;

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
             materials[cell]->getDifTilde()[surface*_num_cmfd_groups + e] *
             (1 - _mesh->getRelaxFactor()) + _mesh->getRelaxFactor() * d_tilde;

          /* Set d_hat and d_tilde */
          materials[cell]->setDifHatByGroup(d_hat, e+1, surface);
          materials[cell]->setDifTildeByGroup(d_tilde, e+1, surface);

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
 * @return k-effective
 */
double Cmfd::computeKeff(){

  log_printf(INFO, "Running diffusion solver...");

  /* Create matrix and vector objects */
  if (_A == NULL){
    try{
    
      if (_num_cmfd_groups == 0)
        createGroupStructure(NULL, _num_groups+1);

      _AM = NULL;
      _M = new double*[_cx*_cy];
      _A = new double*[_cx*_cy];
      _phi_temp = new double[_cx*_cy*_num_cmfd_groups];
      _old_source = new double[_cx*_cy*_num_cmfd_groups];
      _new_source = new double[_cx*_cy*_num_cmfd_groups];

      _mesh->setNumGroups(_num_cmfd_groups);

      _mesh->initializeFlux();

      if (_solve_method == MOC)
        _mesh->initializeMaterialsMOC();

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

  /* If solving diffusion problem, initialize timer */
  if (_solve_method == DIFFUSION)
    _timer->startTimer();

  /* Initialize variables */
  double sum_new, sum_old, val, norm, scale_val;
  double flux_conv = 1e-8;
  int row;
  double* phi_old = _mesh->getFluxes(PRIMAL);
  double* phi_new = _mesh->getFluxes(PRIMAL_UPDATE);

  /* Compute the cross sections and surface diffusion coefficients */
  if (_solve_method == MOC)
    computeXS();

  computeDs();

  /* Construct matrices */
  constructMatrices();

  if (_eigen_method == POWER){

    /* Compute and normalize the initial source */
    matMultM(_M, phi_old, _old_source);
    sum_old = vecSum(_old_source);
    scale_val = (_cx * _cy * _num_cmfd_groups) / sum_old;
    vecScale(_old_source, scale_val);
    vecCopy(phi_old, phi_new);
    vecScale(phi_new, scale_val);
    sum_old = _cx * _cy * _num_cmfd_groups;

    /* Power iteration diffusion solver */
    for (int iter = 0; iter < 25000; iter++){

      /* Solve phi = A^-1 * old_source */
      linearSolve(_A, phi_new, _old_source, flux_conv);

      /* Compute the new source */
      matMultM(_M, phi_new, _new_source);
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
  }
  else{

    /* Allocate memory for AM matrix */
    if (_AM == NULL){
      log_printf(INFO, "Allocating memory for AM");

      try{
        _AM = new double*[_cx*_cy];

        for (int i = 0; i < _cx*_cy; i++){
          _AM[i] = new double[_num_cmfd_groups*(_num_cmfd_groups+4)];
        }
      }
      catch(std::exception &e){
        log_printf(ERROR, "Could not allocate memory for the _AM matrix."
                   " Backtrace:%s", e.what());
      }
      log_printf(INFO, "Done allocating memory for AM");
    }

    int max_iter = 100;
    double shift, offset;
    int iter = 0;
    norm = 1.0;

    /* Copy old flux to new */
    vecCopy(phi_old, phi_new);

    /* Compute the initial k_eff */
    _k_eff = rayleighQuotient(phi_new, _new_source, _phi_temp);

    /* Compute and normalize the initial source */
    matMultM(_M, phi_new, _old_source);
    sum_old = vecSum(_old_source);
    scale_val = (_cx * _cy * _num_cmfd_groups) / sum_old;
    vecScale(_old_source, scale_val);
    vecScale(phi_new, scale_val);
    sum_old = _cx * _cy * _num_cmfd_groups;

    shift = 1.5;

    while (norm > 1e-5){

      /* Reconstruct _AM */
      matSubtract(_AM, _A, 1.0/shift, _M);

      /* Solve inverse system */
      linearSolve(_AM, phi_new, _old_source, flux_conv);

      /* Compute new flux */
      vecScale(phi_new, vecMax(phi_new));
      matMultM(_M, phi_new, _new_source);
      sum_new = vecSum(_new_source);
      vecScale(_new_source, (_cx*_cy*_num_cmfd_groups) / sum_new);
      vecScale(phi_new, (_cx*_cy*_num_cmfd_groups) / sum_new);

      /* Compute new eigenvalue */
      _k_eff = rayleighQuotient(phi_new, _new_source, _phi_temp);

      /* Compute the L2 norm of source error */
      norm = 0.0;
      for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++){
        if (_new_source[i] != 0.0)
          norm += pow((_new_source[i] - _old_source[i]) / _new_source[i], 2);
      }

      norm = pow(norm, 0.5);
      norm = norm / (_cx*_cy*_num_cmfd_groups);
      vecCopy(_new_source, _old_source);

      iter++;

      log_printf(INFO, "iter: %i, k_eff: %f, norm: %f", iter, _k_eff, norm);
    }

    offset = 0.05;
    log_printf(INFO, "offset set to: %f", offset);

    while (norm > _conv_criteria){

      /* Reconstruct _AM */
      matSubtract(_AM, _A, 1.0/(_k_eff + offset), _M);

      /* Solve inverse system */
      linearSolve(_AM, phi_new, _old_source, flux_conv);

      /* compute new flux */
      vecScale(phi_new, vecMax(phi_new));
      matMultM(_M, phi_new, _new_source);
      sum_new = vecSum(_new_source);
      vecScale(_new_source, (_cx*_cy*_num_cmfd_groups) / sum_new);
      vecScale(phi_new, (_cx*_cy*_num_cmfd_groups) / sum_new);

      /* compute new eigenvalue */
      _k_eff = rayleighQuotient(phi_new, _new_source, _phi_temp);

      /* compute the L2 norm of source error */
      norm = 0.0;
      for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++){
        if (_new_source[i] != 0.0)
          norm += pow((_new_source[i] - _old_source[i]) / _new_source[i], 2);
      }

      norm = pow(norm, 0.5);
      norm = norm / (_cx*_cy*_num_cmfd_groups);
      vecCopy(_new_source, _old_source);

      iter++;

      log_printf(INFO, "iter: %i, k_eff: %f, norm: %f", iter, _k_eff, norm);
    }
  }

  /* rescale the old and new flux */
  rescaleFlux();

  /* update the MOC flux */
  if (_solve_method == MOC)
    updateMOCFlux();

  if (_flux_type == ADJOINT)
    vecCopy(phi_new, _mesh->getFluxes(ADJOINT));

  /* If solving diffusion problem, print timing results */
  if (_solve_method == DIFFUSION){
    std::string msg_string;
    log_printf(TITLE, "TIMING REPORT");
    _timer->stopTimer();
    _timer->recordSplit("Total time to solve diffusion eigenvalue problem");

    double tot_time = _timer->getSplit("Total time to solve diffusion "
                                       "eigenvalue problem");
    msg_string = "Total time to solve diffusion eigenvalue problem";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);
  }

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
    vecCopy(vec_x, _phi_temp);

    /* Iteration over red cells */
    #pragma omp parallel for private(row, val, cell)
    for (int y = 0; y < _cy; y++){
      for (int x = y % 2; x < _cx; x += 2){

        cell = y*_cx+x;

        for (int g = 0; g < _num_cmfd_groups; g++){

          row = (y*_cx+x)*_num_cmfd_groups + g;
          val = 0.0;

          /* Previous flux term */
          val += (1.0 - _omega) * vec_x[row];

          /* Source term */
          val += _omega * vec_b[row] / mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _omega * vec_x[row - _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Bottom surface */
          if (y != _cy - 1)
            val -= _omega * vec_x[row + _cx * _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+1] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < _num_cmfd_groups; e++){
            if (e != g)
              val -= _omega * vec_x[(y*_cx+x)*_num_cmfd_groups+e] *
                     mat[cell][g*(_num_cmfd_groups+4)+2+e] /
                     mat[cell][g*(_num_cmfd_groups+4)+g+2];
          }

          /* Right surface */
          if (x != _cx - 1)
            val -= _omega * vec_x[row + _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+2] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Top surface */
          if (y != 0)
            val -= _omega * vec_x[row - _num_cmfd_groups*_cx] *
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
          val += (1.0 - _omega) * vec_x[row];

          /* Source term */
          val += _omega * vec_b[row] / mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Left surface */
          if (x != 0)
            val -= _omega * vec_x[row - _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Bottom surface */
          if (y != _cy - 1)
            val -= _omega * vec_x[row + _cx * _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+1] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Group-to-group */
          for (int e = 0; e < _num_cmfd_groups; e++){
            if (e != g)
              val -= _omega * vec_x[(y*_cx+x)*_num_cmfd_groups+e] *
                     mat[cell][g*(_num_cmfd_groups+4)+2+e] /
                     mat[cell][g*(_num_cmfd_groups+4)+g+2];
          }

          /* Right surface */
          if (x != _cx - 1)
            val -= _omega * vec_x[row + _num_cmfd_groups] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+2] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          /* Top surface */
          if (y != 0)
            val -= _omega * vec_x[row - _num_cmfd_groups*_cx] *
                   mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+3] /
                   mat[cell][g*(_num_cmfd_groups+4)+g+2];

          vec_x[row] = val;
        }
      }
    }

    norm = 0.0;

    for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++){
      if (vec_x[i] != 0.0)
        norm += pow((vec_x[i] - _phi_temp[i]) / vec_x[i], 2);
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
  double* phi_old = _mesh->getFluxes(PRIMAL);
  double* phi_new = _mesh->getFluxes(PRIMAL_UPDATE);

  /* Rescale the new and old flux to have an avg source of 1.0 */
  matMultM(_M, phi_new, _new_source);
  sum_new = vecSum(_new_source);
  scale_val = _cx*_cy*_num_cmfd_groups / sum_new;
  vecScale(phi_new, scale_val);
  matMultM(_M, phi_old, _old_source);
  sum_old = vecSum(_old_source);
  scale_val = _cx*_cy*_num_cmfd_groups / sum_old;
  vecScale(phi_old, scale_val);
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
 * @brief Normalize vector to have avg source of 1.0.
 * @param mat source matrix
 * @param vec vector to be normalized
 */
void Cmfd::vecNormal(double** mat, double* vec){\
  double source, scale_val;
  matMultM(mat, vec, _phi_temp);
  source = vecSum(_phi_temp);
  scale_val = (_cx*_cy*_num_cmfd_groups) / source;
  vecScale(vec, scale_val);
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

  double* heights = _mesh->getLengthsY();
  double* widths = _mesh->getLengthsX();

  /* Zero _A and _M matrices */
  matZero(_M, _num_cmfd_groups);
  matZero(_A, _num_cmfd_groups+4);

  /* Loop over cells */
  #pragma omp parallel for private(value, volume, cell, row, material)
  for (int y = 0; y < _cy; y++){
    for (int x = 0; x < _cx; x++){

      cell = y*_cx + x;
      material = _mesh->getMaterials()[cell];
      volume = _mesh->getVolumes()[cell];

      /* Loop over groups */
      for (int e = 0; e < _num_cmfd_groups; e++){

        row = cell*_num_cmfd_groups + e;

        /* Absorption term */
        value = material->getSigmaA()[e] * volume;
        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Out (1st) and in (2nd) scattering */
        if (_flux_type == PRIMAL){
          for (int g = 0; g < _num_cmfd_groups; g++){
            if (e != g){
              value = material->getSigmaS()[g*_num_cmfd_groups + e] * volume;
              _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;
              value = - material->getSigmaS()[e*_num_cmfd_groups + g] * volume;
              _A[cell][e*(_num_cmfd_groups+4)+g+2] += value;
            }
          }
        }
        else{
          for (int g = 0; g < _num_cmfd_groups; g++){
            if (e != g){
              value = material->getSigmaS()[e*_num_cmfd_groups + g] * volume;
              _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;
              value = - material->getSigmaS()[g*_num_cmfd_groups + e] * volume;
              _A[cell][e*(_num_cmfd_groups+4)+g+2] += value;
            }
          }
        }

        /* RIGHT SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[2*_num_cmfd_groups + e]
                - material->getDifTilde()[2*_num_cmfd_groups + e])
          * heights[cell / _cx];

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (x != _cx - 1){
          value = - (material->getDifHat()[2*_num_cmfd_groups + e]
                  + material->getDifTilde()[2*_num_cmfd_groups + e])
                  * heights[cell / _cx];

          _A[cell][e*(_num_cmfd_groups+4)+_num_cmfd_groups+2] += value;
        }

        /* LEFT SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[0*_num_cmfd_groups + e]
                + material->getDifTilde()[0*_num_cmfd_groups + e])
                * heights[cell / _cx];

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (x != 0){
          value = - (material->getDifHat()[0*_num_cmfd_groups + e]
                  - material->getDifTilde()[0*_num_cmfd_groups + e])
                  * heights[cell / _cx];

          _A[cell][e*(_num_cmfd_groups+4)] += value;
        }

        /* BOTTOM SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[1*_num_cmfd_groups + e]
                - material->getDifTilde()[1*_num_cmfd_groups + e])
                * widths[cell % _cx];

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (y != _cy - 1){
          value = - (material->getDifHat()[1*_num_cmfd_groups + e]
                  + material->getDifTilde()[1*_num_cmfd_groups + e])
                  * widths[cell % _cx];

          _A[cell][e*(_num_cmfd_groups+4)+1] += value;
        }

        /* TOP SURFACE */

        /* Set transport term on diagonal */
        value = (material->getDifHat()[3*_num_cmfd_groups + e]
                + material->getDifTilde()[3*_num_cmfd_groups + e])
                * widths[cell % _cx];

        _A[cell][e*(_num_cmfd_groups+4)+e+2] += value;

        /* Set transport term on off diagonal */
        if (y != 0){
          value = - (material->getDifHat()[3*_num_cmfd_groups + e]
                  - material->getDifTilde()[3*_num_cmfd_groups + e])
                  * widths[cell % _cx];

          _A[cell][e*(_num_cmfd_groups+4)+_num_cmfd_groups+3] += value;
        }

        /* Source term */
        for (int g = 0; g < _num_cmfd_groups; g++){
          value = material->getChi()[e] * material->getNuSigmaF()[g]
                  * volume;

          if (_flux_type == PRIMAL)
            _M[cell][e*_num_cmfd_groups+g] += value;
          else
            _M[cell][g*_num_cmfd_groups+e] += value;
        }

        log_printf(DEBUG, "cel: %i, vol; %f", cell, _mesh->getVolumes()[cell]);

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
  double* old_flux = _mesh->getFluxes(PRIMAL);
  double* new_flux = _mesh->getFluxes(PRIMAL_UPDATE);
  double old_cell_flux, new_cell_flux;

  /* Loop over mesh cells */
  #pragma omp parallel for private(old_cell_flux, new_cell_flux)
  for (int i = 0; i < _cx*_cy; i++){

    std::vector<int>::iterator iter;

    /* Loop over CMFD groups */
    for (int e = 0; e < _num_cmfd_groups; e++){

      /* Get the old and new Mesh cell flux */
      old_cell_flux = old_flux[i*_num_cmfd_groups + e];
      new_cell_flux = new_flux[i*_num_cmfd_groups + e];

      for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){

        /* Loop over FRSs in mesh cell */
        for (iter = _mesh->getCellFSRs()->at(i).begin();
          iter != _mesh->getCellFSRs()->at(i).end(); ++iter) {

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

  if (_mesh->getOpticallyThick() && _solve_method == MOC){

    /* Initialize variables */
    double alpha, mu, expon;
    double rho, F;
    rho = 0.0;

    /* Loop over polar angles */
    for (int p = 0; p < 3; p++){
      mu = cos(asin(_quad->getSinTheta(p)));
      expon = exp(- h / (3 * d * mu));
      alpha = (1 + expon) / (1 - expon) - 2 * mu / h;
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
 * @brief Initialize the FSRs.
 */
void Cmfd::initializeFSRs(){

  log_printf(INFO, "Initialize FSRs...");

  /* Initialize variables */
  int fsr_id;
  CellBasic* cell;
  Material* material;
  Universe* univ_zero = _geometry->getUniverse(0);
  double* heights = _mesh->getLengthsY();
  double* widths = _mesh->getLengthsX();

  for (int i = 0; i < _cx * _cy; i++){

    /* Get mesh cell and fsr volume */
    fsr_id = _mesh->getCellFSRs()->at(i).front();
    _FSR_volumes[fsr_id] = heights[i / _cx] * widths[i % _cx];

    /* Initialize the fsr fluxes to 1.0 */
    for (int e = 0; e < _num_groups; e++)
      _FSR_fluxes[fsr_id*_num_groups+e] = 1.0;

    /* Get the cell corresponding to this FSR from the geometry */
    cell = _geometry->findCellContainingFSR(fsr_id);

    /* Get the cell's material and assign it to the FSR */
    material = _geometry->getMaterial(cell->getMaterial());
    _FSR_materials[fsr_id] = material;

    log_printf(DEBUG, "cell %i with FSR id = %d has cell id = %d and "
               "material id = %d and volume = %f", i, fsr_id, cell->getId(),
               _FSR_materials[fsr_id]->getUid(), _FSR_volumes[fsr_id]);

  }

  log_printf(INFO, "Done initializing FSRs");
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

  std::vector<int>::iterator iter;
  double volume;

  /* Set volume of Mesh cells */
  for (int y = 0; y < _cy; y++){
    for (int x = 0; x < _cx; x++){
      volume = 0.0;

        for (iter = _mesh->getCellFSRs()->at(y*_cx+x).begin();
             iter != _mesh->getCellFSRs()->at(y*_cx+x).end(); ++iter)
          volume += _FSR_volumes[*iter];

        _mesh->setVolume(volume, y*_cx+x);

        log_printf(DEBUG, "set cell %i volume to: %f", y*_cx+x, volume);
    }
  }
}


/**
 * @brief Set pointer to FSR flux array.
 * @param scalar_flux pointer to FSR flux array
 */
void Cmfd::setFSRFluxes(FP_PRECISION* scalar_flux){
  _FSR_fluxes = scalar_flux;
}


/**
 * @brief Get pointer to the Mesh object.
 * @return pointer to mesh
 */
Mesh* Cmfd::getMesh(){
  return _mesh;
}


/**
 * @brief Set the flux type (PRIMAL or ADJOINT).
 * @param flux_type char string representing enum for flux type
 */
void Cmfd::setFluxType(const char* flux_type){

  if (strcmp("PRIMAL", flux_type) == 0)
    _flux_type = PRIMAL;
  else if (strcmp("ADJOINT", flux_type) == 0)
    _flux_type = ADJOINT;
  else
    log_printf(ERROR, "Could not recognize flux type: "
               " the options are PRIMAL and ADJOINT");
}


/**
 * @brief Set the eigenvalue solution method (POWER or WIELANDT).
 * @param eigen_method char string representing enum for eigen method
 */
void Cmfd::setEigenMethod(const char* eigen_method){

  if (strcmp("POWER", eigen_method) == 0)
    _eigen_method = POWER;
  else if (strcmp("WIELANDT", eigen_method) == 0)
    _eigen_method = WIELANDT;
  else
    log_printf(ERROR, "Could not recognize eigen method: "
               " the options are POWER and WIELANDT");
}


/**
 * @brief Dump a vector to the console.
 * @param vec vector to be dumped
 * @param length length of vector
 */
void Cmfd::dumpVec(double* vec, int length){

  log_printf(NORMAL, "dumping vector...");

  for (int i = 0; i < length; i++)
    log_printf(NORMAL, "cell: %i, value: %f", i, vec[i]);

  log_printf(NORMAL, "done dumping vector...");
}


/**
 * @brief Set the SOR factor.
 * @param omega
 */
void Cmfd::setOmega(double omega){
  _omega = omega;
}


/**
 * @brief Computes the Rayleigh quotient.
 * @param x
 * @param snew
 * @param sold
 * @return
 */
double Cmfd::rayleighQuotient(double* x, double* snew, double* sold){

  double numer = 0.0;
  double denom = 0.0;

  matMultA(_A, x, sold);
  matMultM(_M, x, snew);

  for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++){
    numer += x[i]*snew[i];
    denom += x[i]*sold[i];
  }

  return numer/denom;
}


/**
 * @brief Multiply matrix by vector (i.e., y = M *x).
 * @param mat source matrix
 * @param vec_x x vector
 * @param vec_y y vector
 */
void Cmfd::matMultA(double** mat, double* vec_x, double* vec_y){

  vecSet(vec_y, 0.0);
  int row, cell;

  for (int y = 0; y < _cy; y++){
    for (int x = 0; x < _cx; x++){

      cell = y*_cx+x;

      for (int g = 0; g < _num_cmfd_groups; g++){

        row = cell*_num_cmfd_groups + g;

        if (x != 0)
          vec_y[row] += mat[cell][g*(_num_cmfd_groups+4)] *
                        vec_x[(cell-1)*_num_cmfd_groups+g];

        if (y != _cy - 1)
          vec_y[row] += mat[cell][g*(_num_cmfd_groups+4)+1] *
                        vec_x[(cell+_cx)*_num_cmfd_groups+g];

        if (x != _cx - 1)
          vec_y[row] += mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+2] *
                        vec_x[(cell+1)*_num_cmfd_groups+g];

        if (y != 0)
          vec_y[row] += mat[cell][g*(_num_cmfd_groups+4)+_num_cmfd_groups+3] *
                        vec_x[(cell-_cx)*_num_cmfd_groups+g];

        for (int e = 0; e < _num_cmfd_groups; e++)
          vec_y[row] += mat[cell][g*(_num_cmfd_groups+4)+2+e] *
                        vec_x[cell*_num_cmfd_groups+e];
      }
    }
  }
}


/**
 * @brief
 * @param AM
 * @param A
 * @param omega
 * @param M
 */
void Cmfd::matSubtract(double** AM, double** A, double omega, double** M){

  /* Copy A to AM */
  for (int i = 0; i < _cx*_cy; i++){
    for (int g = 0; g < _num_cmfd_groups*(_num_cmfd_groups+4); g++)
      AM[i][g] = A[i][g];
  }

  for (int i = 0; i < _cx*_cy; i++){
    for (int e = 0; e < _num_cmfd_groups; e++){
      for (int g = 0; g < _num_cmfd_groups; g++)
        AM[i][g*(_num_cmfd_groups+4)+e+2] -= omega*M[i][g*_num_cmfd_groups+e];
    }
  }
}


/**
 * @brief Finds and returns the maximum element in a vector.
 * @param vec the vector of interest
 * @return the maximum element in the vector
 */
double Cmfd::vecMax(double* vec){

  double max = vec[0];

  for (int i = 0; i < _cx*_cy*_num_cmfd_groups; i++)
    max = std::max(max, vec[i]);

  return max;
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
  * @brief Create the CMFD coarse energy group structure.
  */
void Cmfd::createGroupStructure(int* group_indices, int ncg){

    _num_cmfd_groups = ncg - 1;

    /* allocate memory */
    if (_group_indices == NULL){
        _group_indices = new int[ncg];
        _group_indices_map = new int[_num_groups];
    }

    if (group_indices == NULL){
        for (int i = 0; i < ncg; i++){
            _group_indices[i] = i;
        }
    }
    else{
        /* check that the group indices span the group space */
        if (group_indices[0] != 0 || group_indices[ncg-1] != _num_groups)
            log_printf(ERROR, "The first and last indicies of group structure "
                       " must be 0 and the number of MOC energy groups");
        
        _group_indices[0] = 0;
        
        for (int i = 1; i < ncg; i++){
            /* check that the group indices are always increasing */
            if (group_indices[i] <= group_indices[i-1])
                log_printf(ERROR, "The group indices must be increasing!");
            
            _group_indices[i] = group_indices[i];
            log_printf(INFO, "group indices %i: %i", i, _group_indices[i]);
        }
    }
    
    /* create group indices map */
    for (int e = 0; e < _num_cmfd_groups; e++){
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){
            _group_indices_map[h] = e;
        }
    }
}
