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
  _quadrature = NULL;
  _geometry = NULL;
  _materials = NULL;

  /* Global variables used in solving CMFD problem */
  _num_x = 1;
  _num_y = 1;
  _width_x = 0.;
  _width_y = 0.;
  _x_min = 0.;
  _y_min = 0.;
  _cell_width_x = 0.;
  _cell_width_y = 0.;
  _flux_update_on = true;
  _centroid_update_on = true;
  _k_nearest = 3;
  _SOR_factor = 1.0;
  _num_FSRs = 0;

  /* Energy group and polar angle problem parameters */
  _num_moc_groups = 0;
  _num_cmfd_groups = 0;
  _num_polar_2 = 0;

  /* Set matrices and arrays to NULL */
  _A = NULL;
  _M = NULL;
  _old_flux = NULL;
  _new_flux = NULL;
  _flux_ratio = NULL;
  _old_source = NULL;
  _new_source = NULL;
  _group_indices = NULL;
  _group_indices_map = NULL;
  _user_group_indices = false;
  _surface_currents = NULL;
  _cell_locks = NULL;
  _volumes = NULL;
  _lattice = NULL;

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

  if (_cell_locks != NULL)
    delete [] _cell_locks;

  if (_boundaries != NULL)
    delete [] _boundaries;

  if (_group_indices != NULL)
    delete [] _group_indices;

  if (_group_indices_map != NULL)
    delete [] _group_indices_map;

  /* Delete the Matrix and Vector objects */
  if (_M != NULL)
    delete _M;

  if (_A != NULL)
    delete _A;

  if (_old_source != NULL)
    delete _old_source;

  if (_new_source != NULL)
    delete _new_source;

  if (_old_flux != NULL)
    delete _old_flux;

  if (_new_flux != NULL)
    delete _new_flux;

  if (_flux_ratio != NULL)
    delete _flux_ratio;

  if (_surface_currents != NULL)
    delete _surface_currents;

  if (_volumes != NULL)
    delete _volumes;

  /* Delete Cmfd materials array */
  if (_materials != NULL) {
    for (int i=0; i < _num_x * _num_y; i++)
      delete _materials[i];
  }
  delete [] _materials;

  /* Delete the Cmfd lattice */
  if (_lattice != NULL)
    delete _lattice;

  /* Clear the _cell_fsrs vector of vectors */
  std::vector< std::vector<int> >::iterator iter1;
  for (iter1 = _cell_fsrs.begin(); iter1 != _cell_fsrs.end(); ++iter1)
    iter1->clear();
  _cell_fsrs.clear();

  /* Clear the _k_nearest_stencils map of vectors */
  std::map<int, std::vector< std::pair<int, FP_PRECISION> > >::iterator iter2;
  for (iter2 = _k_nearest_stencils.begin(); iter2 != _k_nearest_stencils.end();
       ++iter2)
    iter2->second.clear();
  _k_nearest_stencils.clear();
}


/**
 * @brief Set the number of Mesh cells in a row.
 * @param num_x Number of Mesh cells in a row
 */
void Cmfd::setNumX(int num_x) {

  if (num_x < 1)
    log_printf(ERROR, "The number of lattice cells in the x direction "
               "must be > 0. Input value: %d", num_x);

  _num_x = num_x;
  if (_width_x != 0.)
    _cell_width_x = _width_x / _num_x;
}


/**
 * @brief Set the number of Mesh cells in a column
 * @param num_y Number of Mesh cells in a column
 */
void Cmfd::setNumY(int num_y) {

  if (num_y < 1)
    log_printf(ERROR, "The number of lattice cells in the y direction "
               "must be > 0. Input value: %d", num_y);

  _num_y = num_y;
  if (_width_y != 0.)
    _cell_width_y = _width_y / _num_y;
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
 * @brief Set Mesh width in the x-direction
 * @param width Physical width of Mesh in the x-direction
 */
void Cmfd::setWidthX(double width) {
  _width_x = width;
  if (_num_x != 0)
    _cell_width_x = _width_x / _num_x;
}


/**
 * @brief Set Mesh width in the y-direction
 * @param width Physical width of Mesh in the y-direction
 */
void Cmfd::setWidthY(double width) {
  _width_y = width;
  if (_num_y != 0)
    _cell_width_y = _width_y / _num_y;
}


/**
 * @brief Collapse cross-sections and fluxes for each Cmfd cell by
 *        energy condensing and volume averaging cross sections from
 *        the MOC sweep.
 * @details This method performs a cell-wise energy condensation and volume
 *          average of the cross sections of the fine, unstructured FSR mesh.
 *          The cross sections are condensed such that all reaction rates and
 *          the neutron production rate from fission are conserved. It is
 *          important to note that the volume averaging is performed before
 *          energy condensation in order to properly collapse the diffusion
 *          coefficients.
 */
void Cmfd::collapseXS() {

  log_printf(DEBUG, "Collapsing cross-sections onto CMFD mesh...");

  /* Split edge currents to side surfaces */
  splitEdgeCurrents();

#pragma omp parallel
  {

    /* Initialize variables for FSR properties*/
    FP_PRECISION volume, flux, tot, nu_fis, chi;
    FP_PRECISION* scat;

    /* Initialize tallies for each parameter */
    FP_PRECISION nu_fis_tally, rxn_tally;
    FP_PRECISION vol_tally, tot_tally, neut_prod_tally;
    FP_PRECISION scat_tally[_num_cmfd_groups];
    FP_PRECISION chi_tally[_num_cmfd_groups];

    /* Pointers to material objects */
    Material* fsr_material;
    Material* cell_material;

    /* Loop over cmfd cells */
#pragma omp for
    for (int i = 0; i < _num_x * _num_y; i++) {

      cell_material = _materials[i];
      std::vector<int>::iterator iter;

      /* Loop over CMFD coarse energy groups */
      for (int e = 0; e < _num_cmfd_groups; e++) {

        /* Zero tallies for this group */
        nu_fis_tally = 0.0;
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
            for (int h = _group_indices[b]; h < _group_indices[b + 1]; h++)
              chi += fsr_material->getChiByGroup(h+1);

            for (int h = 0; h < _num_moc_groups; h++) {
              chi_tally[b] += chi * fsr_material->getNuSigmaFByGroup(h+1) *
                  _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
              neut_prod_tally += chi * fsr_material->getNuSigmaFByGroup(h+1) *
                  _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
            }
          }
        }

        /* Loop over MOC energy groups within this CMFD coarse group */
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {

          /* Reset volume tally for this MOC group */
          vol_tally = 0.0;

          /* Loop over FSRs in cmfd cell */
          for (iter = _cell_fsrs.at(i).begin();
               iter != _cell_fsrs.at(i).end(); ++iter) {

            /* Gets FSR volume, material, and cross sections */
            fsr_material = _FSR_materials[*iter];
            volume = _FSR_volumes[*iter];
            scat = fsr_material->getSigmaS();
            flux = _FSR_fluxes[(*iter)*_num_moc_groups+h];
            tot = fsr_material->getSigmaTByGroup(h+1);
            nu_fis = fsr_material->getNuSigmaFByGroup(h+1);

            /* Increment tallies for this group */
            tot_tally += tot * flux * volume;
            nu_fis_tally += nu_fis * flux * volume;
            rxn_tally += flux * volume;
            vol_tally += volume;

            /* Scattering tallies */
            for (int g = 0; g < _num_moc_groups; g++) {
              scat_tally[getCmfdGroup(g)] +=
                  scat[g*_num_moc_groups+h] * flux * volume;
            }
          }
        }

        /* Set the Mesh cell properties with the tallies */
        _volumes->setValue(i, 0, vol_tally);
        cell_material->setSigmaTByGroup(tot_tally / rxn_tally, e + 1);
        cell_material->setNuSigmaFByGroup(nu_fis_tally / rxn_tally, e + 1);
        _old_flux->setValue(i, e, rxn_tally / vol_tally);

        /* Set chi */
        if (neut_prod_tally != 0.0)
          cell_material->setChiByGroup(chi_tally[e] / neut_prod_tally, e + 1);
        else
          cell_material->setChiByGroup(0.0, e + 1);

        /* Set scattering xs */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, e + 1,
                                          g + 1);
        }
      }
    }
  }
}


/**
 * @brief Computes the diffusion coefficient for a given cmfd cell and cmfd
 *        energy group.
 * @details This method computes the diffusion coefficient for a cmfd cell and
 *          cmfd energy group by spatially collapsing the total/transport xs
 *          in each FSR contained within the cmfd cell and then energy
 *          collapsing the diffusion coefficient (\f$1 / (3 * \Sigma_t)\f$) for
 *          all MOC groups in the given cmfd energy group.
 * @param cmfd_cell A Cmfd cell
 * @param group A Cmfd energy group
 * @return The diffusion coefficient
 */
FP_PRECISION Cmfd::getDiffusionCoefficient(int cmfd_cell, int group) {

  /* Pointers to material objects */
  Material* fsr_material;
  Material* cell_material = _materials[cmfd_cell];
  std::vector<int>::iterator iter;

  /* Zero tallies for this group */
  FP_PRECISION dif_tally = 0.0;
  FP_PRECISION rxn_tally = 0.0;
  FP_PRECISION trans_tally_group, rxn_tally_group;
  FP_PRECISION tot, flux, volume;

  /* Loop over MOC energy groups within this CMFD coarse group */
  for (int h = _group_indices[group]; h < _group_indices[group+1]; h++) {

    /* Reset transport and rxn tally for this MOC group */
    trans_tally_group = 0.0;
    rxn_tally_group = 0.0;

    /* Loop over FSRs in cmfd cell */
    for (iter = _cell_fsrs.at(cmfd_cell).begin();
         iter != _cell_fsrs.at(cmfd_cell).end(); ++iter) {

      fsr_material = _FSR_materials[*iter];
      volume = _FSR_volumes[*iter];
      flux = _FSR_fluxes[(*iter)*_num_moc_groups+h];
      tot = fsr_material->getSigmaTByGroup(h+1);

      /* Increment tallies for this group */
      rxn_tally += flux * volume;
      trans_tally_group += tot * flux * volume;
      rxn_tally_group += flux * volume;
    }

    /* Energy collapse diffusion coefficient */
    dif_tally += rxn_tally_group /
        (3.0 * (trans_tally_group / rxn_tally_group));
  }

  return dif_tally / rxn_tally;
}


/**
 * @brief Compute the surface diffusion coefficient for a given cmfd cell,
 *        cell surface, and group.
 * @details This method uses finite differencing to compute the surface
 *          diffusion coefficient (\f$ \hat{D} \f$) or surface diffusion
 *          coefficient correction (\f$ \tilde{D} \f$) for a given cmfd cell,
 *          cell surface, and cmfd energy group. If the MOC iteration is zero,
 *          (\f$ \tilde{D} \f$) is returned as zero. Since (\f$ \hat{D} \f$) and
 *          (\f$ \tilde{D} \f$) are dependent on each other, they must be
 *          computed together; therefore, the boolean correction is used to
 *          indicate which value is to be returned.
 * @param cmfd_cell A Cmfd cell
 * @param surface A surface of the Cmfd cell
 * @param group A Cmfd energy group
 * @param moc_iteration MOC iteration number
 * @param correction Boolean indicating whether (\f$ \hat{D} \f$) or
 *                   (\f$ \tilde{D} \f$) is to be returned
 * @return The surface diffusion coefficient, (\f$ \hat{D} \f$) or
 *         (\f$ \tilde{D} \f$)
 */
FP_PRECISION Cmfd::getSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                                  int group, int moc_iteration,
                                                  bool correction) {

  FP_PRECISION dif_surf, dif_surf_corr;
  FP_PRECISION current, current_out, current_in;
  FP_PRECISION flux_next;

  /* Get diffusivity and flux for Mesh cell */
  FP_PRECISION dif_coef = getDiffusionCoefficient(cmfd_cell, group);
  FP_PRECISION flux = _old_flux->getValue(cmfd_cell, group);
  int cmfd_cell_next = getCellNext(cmfd_cell, surface);
  FP_PRECISION delta_interface = getSurfaceWidth(surface);
  FP_PRECISION delta = getPerpendicularSurfaceWidth(surface);
  int sense = getSense(surface);

  /* Correct the diffusion coefficient with Larsen's effective diffusion
   * coefficient correction factor */
  dif_coef *= computeLarsensEDCFactor(dif_coef, delta);

  /* If surface is on a boundary, choose appropriate BC */
  if (cmfd_cell_next == -1) {

    /* REFLECTIVE BC */
    if (_boundaries[surface] == REFLECTIVE) {
      dif_surf = 0.0;
      dif_surf_corr = 0.0;
    }

    /* VACUUM BC */
    else if (_boundaries[surface] == VACUUM) {

      /* Compute the surface-averaged current leaving the cell */
      current_out = sense * _surface_currents->getValue
          (cmfd_cell, surface*_num_cmfd_groups + group) / delta_interface;

      /* Set the surface diffusion coefficient and MOC correction */
      dif_surf =  2 * dif_coef / delta / (1 + 4 * dif_coef / delta);
      dif_surf_corr = (sense * dif_surf * flux - current_out) / flux;
    }
  }

  /* If surface is an interface, use finite differencing */
  else{

    /* Get the surface index for the surface in the neighboring cell */
    int surface_next = (surface + NUM_FACES / 2) % NUM_FACES;

    /* Set diffusion coefficient and flux for the neighboring cell */
    FP_PRECISION dif_coef_next = getDiffusionCoefficient(cmfd_cell_next, group);
    flux_next = _old_flux->getValue(cmfd_cell_next, group);

    /* Correct the diffusion coefficient with Larsen's effective diffusion
     * coefficient correction factor */
    dif_coef_next *= computeLarsensEDCFactor(dif_coef_next, delta);

    /* Compute the surface diffusion coefficient */
    dif_surf = 2.0 * dif_coef * dif_coef_next
        / (delta * dif_coef + delta * dif_coef_next);

    /* Get the outward current on surface */
    current_out = _surface_currents->getValue
        (cmfd_cell, surface*_num_cmfd_groups + group);

    /* Get the inward current on the surface */
    current_in = _surface_currents->getValue
        (cmfd_cell_next, surface_next*_num_cmfd_groups + group);

    /* Compute the surface-averaged net current across the surface */
    current = sense * (current_out - current_in) / delta_interface;

    /* Compute the surface diffusion coefficient correction */
    dif_surf_corr = -(sense * dif_surf * (flux_next - flux) + current)
        / (flux_next + flux);

    /* If the magnitude of dif_surf_corr is greater than the magnitude of
     * dif_surf, select new values dif_surf_corr and dif_surf to ensure the
     * coarse mesh equations are guaranteed to be diagonally dominant */
    if (fabs(dif_surf_corr) > dif_surf && moc_iteration != 0) {

      /* Get the sign of dif_surf_corr */
      int sign = (int) round(fabs(dif_surf_corr) / dif_surf_corr);

      /* Compute the surface diffusion coefficient while preserving the
       * the surface-averaged net current across the surface and the signs
       * of the surface diffusion coefficients. */
      if (sense == sign)
        dif_surf = fabs(current / (2 * flux_next));
      else
        dif_surf = fabs(current / (2 * flux));

      /* Set dif_surf_corr to have the same magnitude as dif_surf,
       *  but the same sign as originially computed. */
      dif_surf_corr = sign * dif_surf;
    }
  }

  /* If it is the first MOC iteration, solve the straight diffusion problem
   * with no MOC correction */
  if (moc_iteration == 0)
    dif_surf_corr = 0.0;

  /* Determine which surface diffusion coefficient is corrected */
  if (correction)
    return dif_surf_corr;
  else
    return dif_surf;
}


/**
 * @brief Solve the nonlinear diffusion acceleration problem to accelerate the
 *        convergence of the MOC problem.
 * @details This method uses the information from the last MOC transport sweep
 *          and solves a simplified nonlinear diffusion problem. The diffusion
 *          problem is tightly converged and the solution is used to update the
 *          the solution of the MOC problem.
 *  @param moc_iteration MOC iteration number
 *  @return The dominant eigenvalue of the nonlinear diffusion problem
 */
FP_PRECISION Cmfd::computeKeff(int moc_iteration) {

  log_printf(DEBUG, "Running diffusion solver...");

  /* Create matrix and vector objects */
  if (_A == NULL) {
    log_printf(ERROR, "Unable to compute k-eff in Cmfd since the Cmfd "
               "linear algebra matrices and arrays have not been created.");
  }

  /* Collapse the cross sections onto the CMFD mesh */
  collapseXS();

  /* Construct matrices */
  constructMatrices(moc_iteration);

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
 *          therefore the solution is independent of flux level. This method
 *          rescales the input flux and converged flux to both have an average
 *          fission source of 1.0 in each group in each cell.
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
 *          accumulates the iteraction and streaming terms into their
 *          approipriate positions in the loss + streaming matrix and
 *          fission gain matrix.
 */
void Cmfd::constructMatrices(int moc_iteration) {

  log_printf(DEBUG,"Constructing matrices...");

  /* Zero _A and _M matrices */
  _A->clear();
  _M->clear();

#pragma omp parallel
  {

    FP_PRECISION value, volume, delta;
    FP_PRECISION dif_surf, dif_surf_corr;
    int sense;
    Material* material;

    /* Loop over cells */
#pragma omp for
    for (int i = 0; i < _num_x*_num_y; i++) {

      material = _materials[i];
      volume = _volumes->getValue(i, 0);

      /* Loop over groups */
      for (int e = 0; e < _num_cmfd_groups; e++) {

        /* Net removal term */
        value = material->getSigmaTByGroup(e+1) * volume;
        _A->incrementValue(i, e, i, e, value);

        /* Scattering gain from all groups */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          value = - material->getSigmaSByGroup(g+1, e+1) * volume;
          _A->incrementValue(i, g, i, e, value);
        }

        /* Streaming to neighboring cells */
        for (int s = 0; s < NUM_FACES; s++) {

          sense = getSense(s);
          delta = getSurfaceWidth(s);

          /* Set transport term on diagonal */
          dif_surf = getSurfaceDiffusionCoefficient(
              i, s, e, moc_iteration, false);
          dif_surf_corr = getSurfaceDiffusionCoefficient(
              i, s, e, moc_iteration, true);

          /* Set the diagonal term */
          value = (dif_surf - sense * dif_surf_corr) * delta;
          _A->incrementValue(i, e, i, e, value);

          /* Set the off diagonal term */
          if (getCellNext(i, s) != -1) {
            value = - (dif_surf + sense * dif_surf_corr) * delta;
            _A->incrementValue(getCellNext(i, s), e, i, e, value);
          }
        }

        /* Fission source term */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          value = material->getChiByGroup(e+1)
              * material->getNuSigmaFByGroup(g+1) * volume;
          _M->incrementValue(i, g, i, e, value);
        }
      }
    }
  }
}


/**
 * @brief Update the MOC flux in each FSR.
 * @details This method uses the condensed flux from the last MOC transport
 *          sweep and the converged flux from the eigenvalue problem to
 *          update the MOC flux in each FSR.
 */
void Cmfd::updateMOCFlux() {

  log_printf(DEBUG, "Updating MOC flux...");

  /* Precompute the CMFD flux ratios */
#pragma omp parallel for
  for (int i = 0; i < _num_x * _num_y; i++) {
    for (int e = 0; e < _num_cmfd_groups; e++)
      _flux_ratio->setValue(i, e, _new_flux->getValue(i, e)
                            / _old_flux->getValue(i, e));
  }

  /* Loop over mesh cells */
#pragma omp parallel for
  for (int i = 0; i < _num_x * _num_y; i++) {

    std::vector<int>::iterator iter;

    /* Loop over CMFD groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Loop over FRSs in mesh cell */
      for (iter = _cell_fsrs.at(i).begin();
           iter != _cell_fsrs.at(i).end(); ++iter) {

        FP_PRECISION update_ratio = getUpdateRatio(i, e, *iter);

        for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++) {

          /* Update FSR flux using ratio of old and new CMFD flux */
          _FSR_fluxes[*iter*_num_moc_groups + h] = update_ratio
            * _FSR_fluxes[*iter*_num_moc_groups + h];

          log_printf(DEBUG, "Updating flux in FSR: %d, cell: %d, MOC group: "
            "%d, CMFD group: %d, ratio: %f", *iter ,i, h, e, update_ratio);
        }
      }
    }
  }
}


/**
 * @brief Compute Larsen's effective diffusion coefficient correction factor.
 * @details By conserving reaction and leakage rates within cells, CMFD
 *          guarantees preservation of area-averaged scalar fluxes and net
 *          surface currents from the MOC fixed source iteration if the CMFD
 *          equations can be converged. However, when the MOC mesh cell size
 *          becomes significantly larger than the neutron mean free path in that
 *          cell, the step characteristics no longer preserve the linear
 *          infinite medium solution to the transport equation. While the
 *          surface diffusion coefficient correction term in CMFD is guaranteed
 *          to preserve reaction rates and surface net currents for any choice
 *          of diffusion coefficient, convergence (and convergence rate) of the
 *          nonlinear iteration acceleration of CMFD is affected by the choice
 *          of diffusion coefficient. All flat source methods, when applied for
 *          thick optical meshes, artificially distribute neutrons in space.
 *          This is the reason that Larsen’s effective diffusion coefficient is
 *          useful in assuring that the CMFD acceleration equations have a
 *          diffusion coefficient (on the flux gradient term) that is
 *          consistent, not with the physical transport problem, but with the
 *          transport problem that is being accelerated by the CMFD equations.
 *          Larsen’s effective diffusion coefficient is precisely this term in
 *          the one-dimensional limit. The following publications provide
 *          further background on how this term is derived and used:
 *
 *            [1] E. Larsen, "Infinite Medium Solutions to the transport
 *                equation, Sn discretization schemes, and the diffusion
 *                approximation", M&C 2001.
 *            [2] S. Shaner, "Transient Method of Characteristics via the
 *                Adiabatic, Theta, and Multigrid Amplitude Function Methods",
 *                Masters Thesis, MIT 2014.
 * @param dif_coef Diffusion coefficient before applying correction factor
 * @param delta Width of the cell in the direction of interest
 * @return The diffusion coefficient correction factor
 */
FP_PRECISION Cmfd::computeLarsensEDCFactor(FP_PRECISION dif_coef,
                                           FP_PRECISION delta) {

  /* Initialize variables */
  FP_PRECISION alpha, mu, expon;
  FP_PRECISION rho = 0.0;

  /* Loop over polar angles */
  for (int p = 0; p < _num_polar_2; p++) {
    mu = cos(asin(_quadrature->getSinTheta(0, p)));
    expon = exp(-delta / (3 * dif_coef * mu));
    alpha = (1 + expon) / (1 - expon) - 2 * (3 * dif_coef * mu) / delta;
    rho += 2.0 * mu * _quadrature->getPolarWeight(0, p) * alpha;
  }

  /* Compute the correction factor */
  FP_PRECISION correction = 1.0 + delta * rho / (2 * dif_coef);

  return correction;
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
 * @brief Set the successive over-relaxation factor for the
 *        linear solve within the diffusion eigenvalue solve.
 * @param SOR_factor Over-relaxation factor
 */
void Cmfd::setSORRelaxationFactor(FP_PRECISION SOR_factor) {

  if (SOR_factor <= 0.0 || SOR_factor >= 2.0)
    log_printf(ERROR, "The successive over-relaxation factor "
        "must be > 0 and < 2. Input value: %d", SOR_factor);

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
 * @brief Set a coarse energy group structure for CMFD.
 * @details CMFD does not necessarily need to have the same energy group
 *          structure as the MOC problem. This function can be used to set
 *          a sparse energy group structure to speed up the CMFD solve. An
 *          example of how this may be called from Python to use a coarse
 *          2-group CMFD structure atop a fine 7-group MOC structure is
 *          illustrated below:
 *
 * @code
 *          cmfd.setGroupStructure([[1,2,3], [4,5,6,7]])
 * @endcode
 *
 * @param group_indices A nested vector of MOC-to-CMFD group mapping
 */
void Cmfd::setGroupStructure(std::vector< std::vector<int> > group_indices) {

  _user_group_indices = true;

  /* Delete old group indices array if it exists */
  if (_group_indices != NULL)
    delete [] _group_indices;

  /* Allocate memory for new group indices */
  _num_cmfd_groups = group_indices.size();
  _group_indices = new int[_num_cmfd_groups+1];

  /* Initialize first group index to 0 */
  int last_moc_group = 0;

  /* Set MOC group bounds for rest of CMFD energy groups */
  for (int i=0; i < _num_cmfd_groups; i++) {
    for (int j=0; j < group_indices[i].size(); j++) {
      if (group_indices[i][j] <= last_moc_group)
	log_printf(ERROR, "The CMFD coarse group indices are not "
		   "monotonically increasing");
      last_moc_group = group_indices[i][j];
    }
    _group_indices[i] = group_indices[i][0] - 1;
    log_printf(DEBUG, "CMFD group indices %d: %d", i, _group_indices[i]);
  }

  /* Set the last group index */
  _group_indices[_num_cmfd_groups] =
    group_indices[_num_cmfd_groups-1].back();
  log_printf(DEBUG, "CMFD group indices %d: %d",
	     _num_cmfd_groups, _group_indices[_num_cmfd_groups]);
}


/**
 * @brief Initialize the CMFD materials.
 */
void Cmfd::initializeMaterials() {

  Material* material;

  /* Delete old Cmfd surface currents vector if it exists */
  if (_materials != NULL)
    delete [] _materials;

  try{
    _materials = new Material*[_num_x*_num_y];

    for (int y = 0; y < _num_y; y++) {
      for (int x = 0; x < _num_x; x++) {
        material = new Material(y*_num_x + x);
        material->setNumEnergyGroups(_num_cmfd_groups);
        _materials[y*_num_x + x] = material;
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
void Cmfd::initializeCurrents() {

  /* Delete old Cmfd surface currents vector if it exists */
  if (_surface_currents != NULL)
    delete _surface_currents;

  /* Allocate memory for the Cmfd Mesh surface currents Vectors */
  _surface_currents = new Vector(_cell_locks, _num_x, _num_y,
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
    for (int x = 0; x < _num_x; x++)
      _cell_fsrs.push_back(std::vector<int>());
  }
}


/**
 * @brief Initialize and set array that links the MOC energy groups to the
 *        CMFD energy groups.
 * @details This method initializes the _group_indices_map, which is a 1D array
 *           of length _num_moc_groups that maps the MOC energy groups to CMFD
 *           energy groups. The indices into _group_indices_map are the MOC
 *           energy groups and the values are the CMFD energy groups.
 */
void Cmfd::initializeGroupMap() {

  /* Setup one-to-one fine-to-coarse group map if not specified by user */
  if (!_user_group_indices) {
    _num_cmfd_groups = _num_moc_groups;

    /* Delete old group indices array if it exists */
    if (_group_indices != NULL)
      delete [] _group_indices;

    /* Allocate memory for new group indices */
    _group_indices = new int[_num_cmfd_groups+1];

    /* Populate a 1-to-1 mapping from MOC to CMFD groups */
    for (int i = 0; i <= _num_cmfd_groups; i++) {
      _group_indices[i] = i;
    }
  }
  else {
    if (_num_moc_groups != _group_indices[_num_cmfd_groups])
      log_printf(ERROR, "The CMFD coarse group mapping is specified for "
		 "%d groups, but the MOC problem contains %d groups",
		 _group_indices[_num_cmfd_groups], _num_moc_groups);
  }

  /* Delete old group indices map if it exists */
  if (_group_indices_map != NULL)
    delete [] _group_indices_map;

  /* Allocate memory for new group indices map */
  _group_indices_map = new int[_num_moc_groups];

  /* Create group indices map */
  for (int e = 0; e < _num_cmfd_groups; e++) {
    for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++) {
      _group_indices_map[h] = e;
    }
  }
}



/**
 * @brief Find the cmfd surface that a LocalCoords object lies on.
 * @details If the coords is not on a surface, -1 is returned. Otherwise,
 *          the surface ID is returned.
 * @param cell_id The CMFD cell ID that the local coords is in.
 * @param coords The coords being evaluated.
 * @return The surface ID.
 */
int Cmfd::findCmfdSurface(int cell_id, LocalCoords* coords) {
  Point* point = coords->getHighestLevel()->getPoint();
  return _lattice->getLatticeSurface(cell_id, point);
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
 * @param cell_id The CMFD cell ID.
 * @param fsr_id The FSR ID.
 */
void Cmfd::addFSRToCell(int cell_id, int fsr_id) {
  _cell_fsrs.at(cell_id).push_back(fsr_id);
}


/**
 * @brief Set the number of MOC energy groups.
 * @param num_groups Number of MOC energy groups
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
  return _num_x * _num_y;
}


/**
 * @brief Set the number of FSRs.
 * @param num_fsrs The number of FSRs
 */
void Cmfd::setNumFSRs(int num_fsrs) {
  _num_FSRs = num_fsrs;
}


/**
 * @brief Split the currents of the Mesh cell edges to the adjacent faces.
 * @details This method takes the currents tallied across the edges (or corners)
 *          of a cmfd cell and splits them evenly across the adjacent faces
 *          (locations 1 and 2). In order to transport the current through to
 *          the diagonal cell, the current is also tallied on the surfaces of
 *          the adjacent cells (locations 3 and 4). Essentially, the tracks that
 *          cross through edges are split into two half-weight tracks as shown
 *          in the illustration below:
 *
 *                                       |    /
 *                                       | __/_
 *                                       |/   /
 *                                     3 /   /
 *                                      /|  / 4
 *                   ------------------/-+-/------------------
 *                                  1 /  |/
 *                                   /   / 2
 *                                  /___/|
 *                                   /   |
 *                                  /    |
 *
 */
void Cmfd::splitEdgeCurrents() {

  log_printf(DEBUG, "Splitting CMFD edge currents...");

  int ncg = _num_cmfd_groups;
  int nf = NUM_FACES;
  int ne = NUM_EDGES;
  int ns = NUM_SURFACES;

#pragma omp parallel
  {

    FP_PRECISION current;
    std::vector<int> surfaces;
    std::vector<int>::iterator iter;
    int cell, surface;

#pragma omp for
    for (int i = 0; i < _num_x * _num_y; i++) {
      for (int e = nf - 1; e < nf + ne; e++) {

        /* Get the surfaces to split this edge's current to */
        getEdgeSplitSurfaces(i, e, &surfaces);

        for (int g=0; g > ncg; g++) {
          /* Divide edge current by 2 since we will split to 2 surfaces,
           * which propagate through 2 surfaces */
          current = _surface_currents->getValue(i, e * ncg + g) / 2;

          /* Increment current for faces adjacent to this edge */
          for (iter = surfaces.begin(); iter != surfaces.end(); ++iter) {
            cell = (*iter) / ns;
            surface = (*iter) % ns;
            _surface_currents->incrementValue(cell, surface * ncg + g, current);
          }
        }
      }
    }
  }
}


/**
 * @brief Get the faces to split the currents of the Mesh cell edges.
 * @details The process by which the current of tracks passing through edges
 *          is split is described in the comment for Cmfd::splitEdgeCurrents().
 *          This method takes in the cell and edge that is being split as well
 *          as a std::vector used to store the IDs of surfaces that are crossed
 *          by the partial-weight tracks. This method properly accounts for
 *          crossings on the geometry boundaries by applying the corresponding
 *          boundary conditions to split the currents.
 * @param cell The CMFD cell ID that the edge is in.
 * @param edge The edge that the track crosses through.
 * @param surfaces A std::vector that is populated with the IDs of surfaces that
 *        are crossed by partial-weight tracks.
 */
void Cmfd::getEdgeSplitSurfaces(int cell, int edge,
                                std::vector<int>* surfaces) {

  surfaces->clear();
  int x = cell % _num_x;
  int y = cell / _num_x;
  int ns = NUM_SURFACES;

  if (edge == SURFACE_X_MIN_Y_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1))* ns + SURFACE_Y_MIN);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Y_MIN);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MIN);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_X_MIN);
  }
  else if (edge == SURFACE_X_MAX_Y_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1))* ns + SURFACE_Y_MIN);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Y_MIN);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MAX);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_X_MAX);
  }
  else if (edge == SURFACE_X_MIN_Y_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1))* ns + SURFACE_Y_MAX);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Y_MAX);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MIN);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_X_MIN);
  }
  else if (edge == SURFACE_X_MAX_Y_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1))* ns + SURFACE_Y_MAX);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Y_MAX);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MAX);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_X_MAX);
  }
}


/**
 * @brief Get the ID of the Mesh cell next to a given Mesh cell across a
 *        given surface.
 * @param cell_id Current Mesh cell ID
 * @param surface_id CMFD cell surface ID to look across for neighboring cell
 * @return Neighboring CMFD cell ID
 */
int Cmfd::getCellNext(int cell_id, int surface_id) {

  int cell_next_id = -1;
  int x = cell_id % _num_x;
  int y = cell_id / _num_x;

  if (surface_id == SURFACE_X_MIN) {
    if (x != 0)
      cell_next_id = cell_id - 1;
    else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
      cell_next_id = cell_id + (_num_x - 1);
  }
  else if (surface_id == SURFACE_Y_MIN) {
    if (y != 0)
      cell_next_id = cell_id - _num_x;
    else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
      cell_next_id = cell_id + _num_x * (_num_y - 1);
  }
  else if (surface_id == SURFACE_X_MAX) {
    if (x != _num_x - 1)
      cell_next_id = cell_id + 1;
    else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
      cell_next_id = cell_id - (_num_x - 1);
  }
  else if (surface_id == SURFACE_Y_MAX) {
    if (y != _num_y - 1)
      cell_next_id = cell_id + _num_x;
    else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
      cell_next_id = cell_id - _num_x * (_num_y - 1);
  }

  return cell_next_id;
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
 *          cell is just a way of describing each of the rectangular regions
 *          that make up a CMFD lattice. CMFD cells are numbered with 0 in the
 *          lower left corner and monotonically increasing from left to right
 *          and from bottom to top. For example, the indices for a 4 x 4
 *          lattice are:
 *                  12  13  14  15
 *                  8    9  10  11
 *                  4    5   6   7
 *                  0    1   2   3
 * @param fsr_id The FSR ID.
 * @return The CMFD cell ID. Return -1 if cell is not found.
 */
int Cmfd::convertFSRIdToCmfdCell(int fsr_id) {

  std::vector<int>::iterator iter;
  for (int cell_id=0; cell_id < _num_x * _num_y; cell_id++) {

    for (iter = _cell_fsrs.at(cell_id).begin();
         iter != _cell_fsrs.at(cell_id).end(); ++iter) {
      if (*iter  == fsr_id)
        return cell_id;
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
void Cmfd::setCentroidUpdateOn(bool centroid_update_on) {
  _centroid_update_on = centroid_update_on;
}


/**
 * @brief Get flag indicating whether to use FSR centroids to update
 *        the MOC flux.
 * @return Flag saying whether to use centroids to update MOC flux.
 */
bool Cmfd::isCentroidUpdateOn() {
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
 * @brief Sets the Quadrature object in use by the MOC Solver.
 * @param quadrature A Quadrature object pointer from the Solver
 */
void Cmfd::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
  _num_polar_2 = quadrature->getNumPolarAngles() / 2;
}


/**
 * @brief Generate the k-nearest neighbor CMFD cell stencil for each FSR.
 * @details This method finds the k-nearest CMFD cell stencil for each FSR
 *          and saves the stencil, ordered from the closest-to-furthest
 *          CMFD cell, in the _k_nearest_stencils map. The stencil of cells
 *          surrounding the current cell is defined as:
 *
 *                             6 7 8
 *                             3 4 5
 *                             0 1 2
 *
 *          where 4 is the given CMFD cell. If the cell is on the edge or corner
 *          of the geometry and there are less than k nearest neighbor cells,
 *          k is reduced to the number of neighbor cells for that instance.
 */
void Cmfd::generateKNearestStencils() {

  if (!_centroid_update_on)
    return;

  std::vector< std::pair<int, FP_PRECISION> >::iterator stencil_iter;
  std::vector<int>::iterator fsr_iter;
  Point* centroid;
  int fsr_id;

  /* Loop over mesh cells */
  for (int i = 0; i < _num_x*_num_y; i++) {

    /* Loop over FRSs in mesh cell */
    for (fsr_iter = _cell_fsrs.at(i).begin();
         fsr_iter != _cell_fsrs.at(i).end(); ++fsr_iter) {

      fsr_id = *fsr_iter;

      /* Get centroid */
      centroid = _geometry->getFSRCentroid(fsr_id);

      /* Create new stencil */
      _k_nearest_stencils[fsr_id] =
        std::vector< std::pair<int, FP_PRECISION> >();

      /* Get distance to all cells that touch current cell */
      for (int j=0; j <= NUM_SURFACES; j++)
        _k_nearest_stencils[fsr_id]
          .push_back(std::make_pair<int, FP_PRECISION>
                     (int(j), getDistanceToCentroid(centroid, i, j)));

      /* Sort the distances */
      std::sort(_k_nearest_stencils[fsr_id].begin(),
                _k_nearest_stencils[fsr_id].end(), stencilCompare);

      /* Remove ghost cells that are outside the geometry boundaries */
      stencil_iter = _k_nearest_stencils[fsr_id].begin();
      while (stencil_iter != _k_nearest_stencils[fsr_id].end()) {
        if (stencil_iter->second == std::numeric_limits<FP_PRECISION>::max())
          stencil_iter = _k_nearest_stencils[fsr_id].erase(stencil_iter++);
        else
          ++stencil_iter;
      }

      /* Resize stencil to be of size <= _k_nearest */
      _k_nearest_stencils[fsr_id].resize
        (std::min(_k_nearest, int(_k_nearest_stencils[fsr_id].size())));
    }
  }

  /* Precompute (1.0 - cell distance / total distance) of each FSR centroid to
   * its k-nearest CMFD cells */
  FP_PRECISION total_distance;
  for (int i=0; i < _num_FSRs; i++) {
    total_distance = 1.e-10;

    /* Compute the total distance of each FSR centroid to its k-nearest CMFD
     * cells */
    for (stencil_iter = _k_nearest_stencils[i].begin();
         stencil_iter < _k_nearest_stencils[i].end(); ++stencil_iter)
      total_distance += stencil_iter->second;

    /* Reset the second stencil value to
     * (1.0 - cell_distance / total_distance) */
    for (stencil_iter = _k_nearest_stencils[i].begin();
         stencil_iter < _k_nearest_stencils[i].end(); ++stencil_iter)
      stencil_iter->second = 1.0 - stencil_iter->second / total_distance;
  }
}


/**
 * @brief Get the ID of the Mesh cell given a stencil ID and Mesh cell ID.
 * @details The stencil of cells surrounding the current cell is defined as:
 *
 *                             6 7 8
 *                             3 4 5
 *                             0 1 2
 *
 * @param cell_id Current Mesh cell ID
 * @param stencil_id CMFD cell stencil ID
 * @return Neighboring CMFD cell ID
 */
int Cmfd::getCellByStencil(int cell_id, int stencil_id) {

  int cell_next_id = -1;
  int x = cell_id % _num_x;
  int y = cell_id / _num_x;

  if (stencil_id == 0) {
    if (x != 0 && y != 0)
      cell_next_id = cell_id - _num_x - 1;
  }
  else if (stencil_id == 1) {
    if (y != 0)
      cell_next_id = cell_id - _num_x;
    else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
      cell_next_id = cell_id + _num_x * (_num_y - 1);
  }
  else if (stencil_id == 2) {
    if (x != _num_x - 1 && y != 0)
      cell_next_id = cell_id - _num_x + 1;
  }
  else if (stencil_id == 3) {
    if (x != 0)
      cell_next_id = cell_id - 1;
    else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
      cell_next_id = cell_id + (_num_x - 1);
  }
  else if (stencil_id == 4) {
    cell_next_id = cell_id;
  }
  else if (stencil_id == 5) {
    if (x != _num_x - 1)
      cell_next_id = cell_id + 1;
    else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
      cell_next_id = cell_id - (_num_x - 1);
  }
  else if (stencil_id == 6) {
    if (x != 0 && y != _num_y - 1)
      cell_next_id = cell_id + _num_x - 1;
  }
  else if (stencil_id == 7) {
    if (y != _num_y - 1)
      cell_next_id = cell_id + _num_x;
    else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
      cell_next_id = cell_id - _num_x * (_num_y - 1);
  }
  else if (stencil_id == 8) {
    if (x != _num_x - 1 && y != _num_y - 1)
      cell_next_id = cell_id + _num_x + 1;
  }

  return cell_next_id;
}


/**
 * @brief Get the ratio used to update the FSR flux after converging CMFD.
 * @details This method takes in a cmfd cell, a MOC energy group, and a FSR
 *          and returns the ratio used to update the FSR flux. There are two
 *          methods that can be used to update the flux, conventional and
 *          k-nearest centroid updating. The k-nearest centroid updating uses
 *          the k-nearest cells (with k between 1 and 9) of the current CMFD
 *          cell and the 8 neighboring CMFD cells. The stencil of cells
 *          surrounding the current cell is defined as:
 *
 *                             6 7 8
 *                             3 4 5
 *                             0 1 2
 *
 *          where 4 is the given CMFD cell. If the cell is on the edge or corner
 *          of the geometry and there are less than k nearest neighbor cells,
 *          k is reduced to the number of neighbor cells for that instance.
 * @param cell_id The cmfd cell ID containing the FSR.
 * @param group The CMFD energy group being updated.
 * @param fsr The fsr being updated.
 * @return the ratio used to update the FSR flux.
 */
FP_PRECISION Cmfd::getUpdateRatio(int cell_id, int group, int fsr) {

  FP_PRECISION ratio = 0.0;
  std::vector< std::pair<int, FP_PRECISION> >::iterator iter;
  int cell_next_id;

  if (_centroid_update_on) {

    /* Compute the ratio for all the surrounding cells */
    for (iter = _k_nearest_stencils[fsr].begin();
         iter != _k_nearest_stencils[fsr].end(); ++iter) {

      if (iter->first != 4) {
        cell_next_id = getCellByStencil(cell_id, iter->first);
        ratio += iter->second * _flux_ratio->getValue(cell_next_id, group);
      }
    }

    /* INTERNAL */
    if (_k_nearest_stencils[fsr].size() == 1)
      ratio += _flux_ratio->getValue(cell_id, group);
    else {
      ratio += _k_nearest_stencils[fsr][0].second *
        _flux_ratio->getValue(cell_id, group);
      ratio /= (_k_nearest_stencils[fsr].size() - 1);
    }
  }
  else
    ratio = _flux_ratio->getValue(cell_id, group);

  return ratio;
}


/**
 * @brief Get the distances from an FSR centroid to a given cmfd cell.
 * @details This method takes in a FSR centroid, a cmfd cell, and a stencil index
 *          to a cell located in the 9-point stencil encompassing the cmfd
 *          cell an all its possible neighbors. The CMFD cell stencil is:
 *
 *                             6 7 8
 *                             3 4 5
 *                             0 1 2
 *
 *          where 4 is the given CMFD cell. If a CMFD edge or corner cells is
 *          given and the stencil indexed cell lies outside the geometry, the
 *          maximum allowable FP_PRECISION value is returned.
 * @param centroid The numerical centroid an FSR in the cell.
 * @param cell_id The cmfd cell containing the FSR.
 * @param stencil_index The index of the cell in the stencil that we want to
 *        get the distance from.
 * @return the distance from the CMFD cell centroid to the FSR centroid.
 */
FP_PRECISION Cmfd::getDistanceToCentroid(Point* centroid, int cell_id,
                                         int stencil_index) {

  int x = cell_id % _num_x;
  int y = cell_id / _num_x;
  FP_PRECISION dist_x, dist_y;
  bool found = false;
  FP_PRECISION centroid_x = centroid->getX();
  FP_PRECISION centroid_y = centroid->getY();

  /* LOWER LEFT CORNER */
  if (x > 0 && y > 0 && stencil_index == 0) {
    dist_x = pow(centroid_x - (_x_min + (x - 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y - 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* BOTTOM SIDE */
  else if (y > 0 && stencil_index == 1) {
    dist_x = pow(centroid_x - (_x_min + (x + 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y - 0.5)*_cell_width_y),2.0);
    found = true;
  }

  /* LOWER RIGHT CORNER */
  else if (x < _num_x - 1 && y > 0 && stencil_index == 2) {
    dist_x = pow(centroid_x - (_x_min + (x + 1.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y - 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* LEFT SIDE */
  else if (x > 0 && stencil_index == 3) {
    dist_x = pow(centroid_x - (_x_min + (x - 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y + 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* CURRENT */
  else if (stencil_index == 4) {
    dist_x = pow(centroid_x - (_x_min + (x + 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y + 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* RIGHT SIDE */
  else if (x < _num_x - 1 && stencil_index == 5) {
    dist_x = pow(centroid_x - (_x_min + (x + 1.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y + 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* UPPER LEFT CORNER */
  else if (x > 0 && y < _num_y - 1 && stencil_index == 6) {
    dist_x = pow(centroid_x - (_x_min + (x - 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y + 1.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* TOP SIDE */
  else if (y < _num_y - 1 && stencil_index == 7) {
    dist_x = pow(centroid_x - (_x_min + (x + 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y + 1.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* UPPER RIGHT CORNER */
  else if (x < _num_x - 1 && y < _num_y - 1 && stencil_index == 8) {
    dist_x = pow(centroid_x - (_x_min + (x + 1.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (_y_min + (y + 1.5)*_cell_width_y), 2.0);
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
  int cell_id;

  log_printf(DEBUG, "Updating boundary flux...");

  /* Loop over Tracks */
  for (int i=0; i < num_tracks; i++) {

    num_segments = tracks[i]->getNumSegments();
    segments = tracks[i]->getSegments();

    /* Update boundary flux in forward direction */
    bc = (int)tracks[i]->getBCIn();
    curr_segment = &segments[0];
    track_flux = &boundary_flux[i*2*_num_moc_groups*_num_polar_2];
    cell_id = convertFSRIdToCmfdCell(curr_segment->_region_id);

    if (bc) {
      for (int e=0; e < _num_moc_groups; e++) {
        for (int p=0; p < _num_polar_2; p++) {
          track_flux[p*_num_moc_groups + e] *= _flux_ratio->getValue
            (cell_id, e);
        }
      }
    }

    /* Update boundary flux in backwards direction */
    bc = (int)tracks[i]->getBCOut();
    curr_segment = &segments[num_segments - 1];
    track_flux = &boundary_flux[(i*2 + 1)*_num_moc_groups*_num_polar_2];
    cell_id = convertFSRIdToCmfdCell(curr_segment->_region_id);

    if (bc) {
      for (int e=0; e < _num_moc_groups; e++) {
        for (int p=0; p < _num_polar_2; p++) {
          track_flux[p*_num_moc_groups + e] *= _flux_ratio->getValue
            (cell_id, e);
        }
      }
    }
  }
}


/** @brief Set a pointer to the Geometry.
 * @param goemetry A pointer to a Geometry object.
 */
void Cmfd::setGeometry(Geometry* geometry) {
  _geometry = geometry;
}


/** @brief Set a number of k-nearest neighbor cells to use in updating
 *         the FSR flux.
 * @param k_nearest The number of nearest neighbor CMFD cells.
 */
void Cmfd::setKNearest(int k_nearest) {

  if (_k_nearest < 1 || k_nearest > 9)
    log_printf(ERROR, "Unable to set CMFD k-nearest to %d. k-nearest "
               "must be between 1 and 9.", k_nearest);
  else
    _k_nearest = k_nearest;
}


/**
 * @brief Zero the surface currents for each mesh cell and energy
 *        group.
 */
void Cmfd::zeroCurrents() {
  _surface_currents->clear();
}


/**
 * @brief Tallies the current contribution from this segment across the
 *        the appropriate CMFD mesh cell surface.
 * @param curr_segment The current Track segment
 * @param track_flux The outgoing angular flux for this segment
 * @param azim_index Azimuthal angle index of the current Track
 * @param fwd Boolean indicating direction of integration along segment
 */
void Cmfd::tallyCurrent(segment* curr_segment, FP_PRECISION* track_flux,
                        int azim_index, bool fwd) {

  int surf_id, cell_id;
  int ncg = _num_cmfd_groups;
  FP_PRECISION currents[_num_cmfd_groups];
  memset(currents, 0.0, sizeof(FP_PRECISION) * _num_cmfd_groups);

  if (fwd) {
    if (curr_segment->_cmfd_surface_fwd != -1) {

      surf_id = curr_segment->_cmfd_surface_fwd % NUM_SURFACES;
      cell_id = curr_segment->_cmfd_surface_fwd / NUM_SURFACES;

      for (int e=0; e < _num_moc_groups; e++) {
        int g = getCmfdGroup(e);

        for (int p=0; p < _num_polar_2; p++)
          currents[g] += track_flux(p, e) *
                         _quadrature->getWeightInline(azim_index, p);
      }

      /* Increment currents */
      _surface_currents->incrementValues
          (cell_id, surf_id*ncg, (surf_id+1)*ncg - 1, currents);
    }
  }
  else {
    if (curr_segment->_cmfd_surface_bwd != -1) {

      surf_id = curr_segment->_cmfd_surface_bwd % NUM_SURFACES;
      cell_id = curr_segment->_cmfd_surface_bwd / NUM_SURFACES;

      for (int e=0; e < _num_moc_groups; e++) {

        int g = getCmfdGroup(e);

        for (int p=0; p < _num_polar_2; p++)
          currents[g] += track_flux(p, e) *
                         _quadrature->getWeightInline(azim_index, p);
      }

      /* Increment currents */
      _surface_currents->incrementValues
          (cell_id, surf_id*ncg, (surf_id+1)*ncg - 1, currents);
    }
  }
}


/**
 * @brief Initialize the Matrix and Vector objects, k-nearest stencils, the
 *        CMFD cell currents and MOC materials.
 */
void Cmfd::initialize() {

  int num_cells = getNumCells();

  /* Delete old Matrix and Vector objects if they exist */
  if (_M != NULL)
    delete _M;
  if (_A != NULL)
    delete _A;
  if (_old_source != NULL)
    delete _old_source;
  if (_new_source != NULL)
    delete _new_source;
  if (_old_flux != NULL)
    delete _old_flux;
  if (_new_flux != NULL)
    delete _new_flux;
  if (_flux_ratio != NULL)
    delete _flux_ratio;
  if (_volumes != NULL)
    delete _volumes;
  if (_cell_locks != NULL)
    delete [] _cell_locks;

  try {

    /* Allocate array of OpenMP locks for each CMFD cell */
    _cell_locks = new omp_lock_t[num_cells];

    /* Loop over all cells to initialize OpenMP locks */
#pragma omp parallel for schedule(guided)
    for (int r=0; r < num_cells; r++)
      omp_init_lock(&_cell_locks[r]);

    /* Allocate memory for matrix and vector objects */
    _M = new Matrix(_cell_locks, _num_x, _num_y, _num_cmfd_groups);
    _A = new Matrix(_cell_locks, _num_x, _num_y, _num_cmfd_groups);
    _old_source = new Vector(_cell_locks, _num_x, _num_y, _num_cmfd_groups);
    _new_source = new Vector(_cell_locks, _num_x, _num_y, _num_cmfd_groups);
    _old_flux = new Vector(_cell_locks, _num_x, _num_y, _num_cmfd_groups);
    _new_flux = new Vector(_cell_locks, _num_x, _num_y, _num_cmfd_groups);
    _flux_ratio = new Vector(_cell_locks, _num_x, _num_y, _num_cmfd_groups);
    _volumes = new Vector(_cell_locks, _num_x, _num_y, 1);

    /* Set the minimum x and y values of the geometry */
    _x_min = _geometry->getMinX();
    _y_min = _geometry->getMinY();

    /* Initialize k-nearest stencils, currents, flux, and materials */
    generateKNearestStencils();
    initializeCurrents();
    initializeMaterials();
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the CMFD mesh objects. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initialize the CMFD lattice.
 */
void Cmfd::initializeLattice(Point* offset) {

  /* Delete old lattice if it exists */
  if (_lattice != NULL)
    delete _lattice;

  /* Initialize the lattice */
  _lattice = new Lattice();
  _lattice->setNumX(_num_x);
  _lattice->setNumY(_num_y);
  _lattice->setWidth(_cell_width_x, _cell_width_y);
  _lattice->setOffset(offset->getX(), offset->getY(), 0.0);
}


/**
 * @brief Returns the width of a given surface
 * @param surface A surface index, from 0 to NUM_SURFACES - 1
 * @return The surface width
 */
FP_PRECISION Cmfd::getSurfaceWidth(int surface) {

  FP_PRECISION width;

  if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
    width = _cell_width_y;
  else
    width = _cell_width_x;

  return width;
}


/**
 * @brief Returns the width of the surface perpendicular to a given surface
 * @param surface A surface index, from 0 to NUM_SURFACES - 1
 * @return The perpendicular surface width
 */
FP_PRECISION Cmfd::getPerpendicularSurfaceWidth(int surface) {

  if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
    return _cell_width_x;
  else
    return _cell_width_y;
}


/**
 * @brief Returns the sense of a given surface
 * @details The sense of minimum surfaces (e.g. SURFACE_X_MIN) is defined to be
 *          -1 while maximum surfaces (e.g. SURFACE_X_MAX) are defined to have a
 *          sense of +1. This is based on the current exiting a cell from a
 *          minimum surface being in the direction of negative net current and
 *          the current leaving a cell from a maximum surface being in the
 *          direction of positive net current.
 * @param surface A surface index, from 0 to NUM_SURFACES - 1
 * @return The sense of the surface
 */
int Cmfd::getSense(int surface) {

  if (surface == SURFACE_X_MIN || surface == SURFACE_Y_MIN)
    return -1;
  else
    return 1;
}
