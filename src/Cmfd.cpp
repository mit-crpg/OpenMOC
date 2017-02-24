#include "Cmfd.h"

/**
 * @brief Constructor initializes boundaries and variables that describe
 *          the CMFD object.
 * @details The construcor initializes the many variables that describe
 *          the CMFD mesh and are used to solve the nonlinear diffusion
 *          acceleration problem.
 */
Cmfd::Cmfd() {

  /* Initialize Geometry and Mesh-related attribute */
  _quadrature = NULL;
  _geometry = NULL;
  _materials = NULL;

  //FIXME
  _convergence_data = NULL;
  _domain_communicator = NULL;

  /* Global variables used in solving CMFD problem */
  _source_convergence_threshold = 1E-5;
  _num_x = 1;
  _num_y = 1;
  _num_z = 1;
  _local_num_x = 1;
  _local_num_y = 1;
  _local_num_z = 1;
  _width_x = 0.;
  _width_y = 0.;
  _width_z = 0.;
  _cell_width_x = 0.;
  _cell_width_y = 0.;
  _cell_width_z = 0.;
  _flux_update_on = true;
  _centroid_update_on = true;
  _use_axial_interpolation = false;
  _k_nearest = 1;
  _SOR_factor = 1.0;
  _num_FSRs = 0;
  _solve_3D = false;
  _total_tally_size = 0;
  _tallies_allocated = false;
  _linear_source = false;

  /* Energy group and polar angle problem parameters */
  _num_moc_groups = 0;
  _num_cmfd_groups = 0;
  _num_polar = 0;

  /* Set matrices and arrays to NULL */
  _A = NULL;
  _M = NULL;
  _k_eff = 1.0;
  _relaxation_factor = 1.0;
  _old_flux = NULL;
  _new_flux = NULL;
  _old_dif_surf_corr = NULL;
  _old_source = NULL;
  _new_source = NULL;
  _flux_moments = NULL;
  _group_indices = NULL;
  _group_indices_map = NULL;
  _user_group_indices = false;
  _surface_currents = NULL;
  _cell_locks = NULL;
  _volumes = NULL;
  _lattice = NULL;
  _azim_spacings = NULL;
  _polar_spacings = NULL;

  /* Initialize boundaries to be reflective */
  _boundaries = new boundaryType[6];
  _boundaries[SURFACE_X_MIN] = REFLECTIVE;
  _boundaries[SURFACE_X_MAX] = REFLECTIVE;
  _boundaries[SURFACE_Y_MIN] = REFLECTIVE;
  _boundaries[SURFACE_Y_MAX] = REFLECTIVE;
  _boundaries[SURFACE_Z_MIN] = REFLECTIVE;
  _boundaries[SURFACE_Z_MAX] = REFLECTIVE;

  /* Initialize CMFD timer */
  _timer = new Timer();
}


/**
 * @brief Destructor deletes arrays of A and M row insertion arrays.
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

  if (_surface_currents != NULL)
    delete _surface_currents;

  if (_volumes != NULL)
    delete _volumes;

  if (_azim_spacings != NULL)
    delete [] _azim_spacings;

  if (_polar_spacings != NULL) {
    for (int a=0; a < _quadrature->getNumAzimAngles()/4; a++)
      delete [] _polar_spacings[a];
    delete [] _polar_spacings;
  }

  /* Delete CMFD materials array */
  if (_materials != NULL) {
    for (int i=0; i < _num_x * _num_y * _num_z; i++)
      delete _materials[i];
  }
  delete [] _materials;

  /* Delete the CMFD lattice */
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

  /* Delete tally information */
  if (_tallies_allocated) {

    delete [] _tally_memory;
    delete [] _nu_fission_tally;
    delete [] _reaction_tally;
    delete [] _volume_tally;
    delete [] _total_tally;
    delete [] _neutron_production_tally;
    delete [] _diffusion_tally;

    for (int i=0; i < _num_x * _num_y * _num_z; i++) {
      delete [] _scattering_tally[i];
      delete [] _chi_tally[i];
    }
    delete [] _scattering_tally;
    delete [] _chi_tally;
#ifdef MPIx
    if (_geometry->isDomainDecomposed()) {
      delete [] _tally_buffer;
      delete [] _currents_buffer;
    }
#endif
  }

  /* TODO: clean, document */
  int num_cells_local = _local_num_x * _local_num_y * _local_num_z;
  if (_domain_communicator != NULL) {
    for (int rb=0; rb<2; rb++) {
      for (int nsc=0; nsc < num_cells_local * _num_cmfd_groups; nsc++) {
        delete [] _domain_communicator->indexes[rb][nsc];
        delete [] _domain_communicator->domains[rb][nsc];
        delete [] _domain_communicator->coupling_coeffs[rb][nsc];
        delete [] _domain_communicator->fluxes[rb][nsc];
      }
      delete [] _domain_communicator->num_connections[rb];
      delete [] _domain_communicator->indexes[rb];
      delete [] _domain_communicator->domains[rb];
      delete [] _domain_communicator->coupling_coeffs[rb];
      delete [] _domain_communicator->fluxes[rb];
    }

    delete [] _domain_communicator->num_connections;
    delete [] _domain_communicator->indexes;
    delete [] _domain_communicator->fluxes;
    delete [] _domain_communicator->coupling_coeffs;
    delete _domain_communicator;
  }

  for (int r=0; r < _num_FSRs; r++)
    delete [] _axial_interpolants.at(r);
}


/**
 * @brief Set the number of Mesh cells in a row.
 * @param number of Mesh cells in a row
 */
void Cmfd::setNumX(int num_x) {

  if (num_x < 1)
    log_printf(ERROR, "The number of lattice cells in the x direction "
               "must be > 0. Input value: %i", num_x);

  _num_x = num_x;
  _local_num_x = _num_x;
  if (_domain_communicator != NULL)
    _local_num_x = _num_x / _domain_communicator->_num_domains_x;
  if (_width_x != 0.)
    _cell_width_x = _width_x / _num_x;
}


/**
 * @brief Set the number of Mesh cells in a column
 * @param number of Mesh cells in a column
 */
void Cmfd::setNumY(int num_y) {

  if (num_y < 1)
    log_printf(ERROR, "The number of lattice cells in the y direction "
               "must be > 0. Input value: %i", num_y);

  _num_y = num_y;
  _local_num_y = _num_y;
  if (_domain_communicator != NULL)
    _local_num_y = _num_y / _domain_communicator->_num_domains_y;
  if (_width_y != 0.)
    _cell_width_y = _width_y / _num_y;
}


/**
 * @brief Set the number of Mesh cells in a column
 * @param number of Mesh cells in a column
 */
void Cmfd::setNumZ(int num_z) {

  if (num_z < 1)
    log_printf(ERROR, "The number of lattice cells in the z direction "
               "must be > 0. Input value: %i", num_z);

  _num_z = num_z;
  _local_num_z = _num_z;
  if (_domain_communicator != NULL)
    _local_num_z = _num_z / _domain_communicator->_num_domains_z;
  if (_width_z != 0.)
    _cell_width_z = _width_z / _num_z;
  if (_width_z == std::numeric_limits<double>::infinity())
    _cell_width_z = 1.0;
}


/**
 * @brief Get the number of Mesh cells in a row.
 * @return number of Mesh cells in a row
 */
int Cmfd::getNumX() {
  return _num_x;
}


/**
 * @brief Get the number of Mesh cells in a column
 * @return number of Mesh cells in a column
 */
int Cmfd::getNumY() {
  return _num_y;
}


int Cmfd::getNumZ() {
  return _num_z;
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
 * @brief Set Mesh width in the z-direction
 * @param width Physical width of Mesh in the z-direction
 */
void Cmfd::setWidthZ(double width) {
  _width_z = width;
  if (_num_z != 0)
    _cell_width_z = _width_z / _num_z;
}


#ifdef MPIx
//TODO: document
void Cmfd::setNumDomains(int num_x, int num_y, int num_z) {

  if (_domain_communicator == NULL) {
    _domain_communicator = new DomainCommunicator;
    _geometry->getMPICart();
    _domain_communicator->_MPI_cart = _geometry->getMPICart();
  }

  _domain_communicator->_num_domains_x = num_x;
  _domain_communicator->_num_domains_y = num_y;
  _domain_communicator->_num_domains_z = num_z;

  _local_num_x = _num_x / num_x;
  _local_num_y = _num_y / num_y;
  _local_num_z = _num_z / num_z;
}


//TODO: document
void Cmfd::setDomainIndexes(int idx_x, int idx_y, int idx_z) {

  if (_domain_communicator == NULL) {
    _domain_communicator = new DomainCommunicator;
    _geometry->getMPICart();
    _domain_communicator->_MPI_cart = _geometry->getMPICart();
  }

  _domain_communicator->stop = false;
  _domain_communicator->_domain_idx_x = idx_x;
  _domain_communicator->_domain_idx_y = idx_y;
  _domain_communicator->_domain_idx_z = idx_z;
}
#endif


/**
 * @brief Collapse cross-sections and fluxes for each CMFD cell by
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

  log_printf(INFO, "Collapsing cross-sections onto CMFD mesh...");

  /* Check to see that CMFD tallies have been allocated */
  if (!_tallies_allocated)
    log_printf(ERROR, "Tallies need to be allocated before collapsing "
               "cross-sections");

  /* Reduce currents across domains if necessary */
#ifdef MPIx
  MPI_Comm comm;
  MPI_Datatype precision;
  if (_geometry->isDomainDecomposed()) {

    /* Start recording MPI communication time */
    _timer->startTimer();

    /* Retrieve MPI information */
    comm = _geometry->getMPICart();
    if (sizeof(FP_PRECISION) == 4)
      precision = MPI_FLOAT;
    else
      precision = MPI_DOUBLE;

    /* Reduce currents form all domains */
    FP_PRECISION* currents_array = _surface_currents->getArray();
    long num_currents = _num_x * _num_y * _num_z * _num_cmfd_groups
        * NUM_FACES;

    long num_messages = num_currents / CMFD_BUFFER_SIZE + 1;
    for (int i=0; i < num_messages; i++) {
      long min_idx = i * CMFD_BUFFER_SIZE;
      long max_idx = (i+1) * CMFD_BUFFER_SIZE;
      max_idx = std::min(max_idx, num_currents);
      int num_elements = max_idx - min_idx;
      FP_PRECISION* send_array = &currents_array[i*CMFD_BUFFER_SIZE];
      if (num_elements > 0)
        MPI_Allreduce(send_array, _currents_buffer, num_elements,
                      precision, MPI_SUM, comm);
      for (int j=0; j < num_elements; j++)
        currents_array[i*CMFD_BUFFER_SIZE+j] = _currents_buffer[j];
    }

    /* Tally MPI communication time */
    _timer->stopTimer();
    _timer->recordSplit("Total MPI communication time");
  }
#endif

  /* Split vertex and edge currents to side surfaces */
  splitVertexCurrents();
  splitEdgeCurrents();

  //TODO: clean
  bool neg_fluxes = false;

#pragma omp parallel
  {

    /* Initialize variables for FSR properties*/
    FP_PRECISION volume, flux, tot, nu_fis, chi;
    FP_PRECISION* scat;

    /* Pointers to material objects */
    Material* fsr_material;
    Material* cell_material;

    /* Loop over CMFD cells */
#pragma omp for
    for (int i = 0; i < _num_x * _num_y * _num_z; i++) {

      cell_material = _materials[i];
      std::vector<int>::iterator iter;

      /* Loop over CMFD coarse energy groups */
      for (int e = 0; e < _num_cmfd_groups; e++) {

        /* Zero tallies for this group */
        _nu_fission_tally[i][e] = 0.0;
        _reaction_tally[i][e] = 0.0;
        _volume_tally[i][e] = 0.0;
        _total_tally[i][e] = 0.0;
        _neutron_production_tally[i][e] = 0.0;
        _diffusion_tally[i][e] = 0.0;

        /* Zero each group-to-group scattering tally */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          _scattering_tally[i][e][g] = 0.0;
          _chi_tally[i][e][g] = 0.0;
        }

        /* Loop over FSRs in CMFD cell to compute chi */
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
              _chi_tally[i][e][b] += chi *
                  fsr_material->getNuSigmaFByGroup(h+1) *
                  _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
              _neutron_production_tally[i][e] += chi *
                  fsr_material->getNuSigmaFByGroup(h+1) *
                  _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;
            }
          }
        }

        /* Loop over MOC energy groups within this CMFD coarse group */
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {

          /* Reset volume tally for this MOC group */
          _volume_tally[i][e] = 0.0;
          FP_PRECISION rxn_tally_group = 0.0;
          FP_PRECISION trans_tally_group = 0.0;

          /* Loop over FSRs in CMFD cell */
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
            _total_tally[i][e] += tot * flux * volume;
            _nu_fission_tally[i][e] += nu_fis * flux * volume;
            _reaction_tally[i][e] += flux * volume;
            _volume_tally[i][e] += volume;

            /* Increment diffusion MOC group-wise tallies */
            rxn_tally_group += flux * volume;
            trans_tally_group += tot * flux * volume;

            /* Scattering tallies */
            for (int g = 0; g < _num_moc_groups; g++) {
              _scattering_tally[i][e][getCmfdGroup(g)] +=
                  scat[g*_num_moc_groups+h] * flux * volume;
            }
          }
          if (rxn_tally_group != 0 && trans_tally_group != 0) {
            FP_PRECISION flux_avg_sigma_t = trans_tally_group /
                rxn_tally_group;
            _diffusion_tally[i][e] += rxn_tally_group /
                (3.0 * flux_avg_sigma_t);
          }
        }
      }
    }
  }

  /* Reduce tallies across domains if necessary */
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {

    /* Start recording MPI communication time */
    _timer->startTimer();

    /* Communicate tallies for XS condensation */
    long num_messages = _total_tally_size / CMFD_BUFFER_SIZE + 1;
    for (int i=0; i < num_messages; i++) {
      long min_idx = i * CMFD_BUFFER_SIZE;
      long max_idx = (i+1) * CMFD_BUFFER_SIZE;
      max_idx = std::min(max_idx, _total_tally_size);
      int num_elements = max_idx - min_idx;
      FP_PRECISION* send_array = &_tally_memory[i*CMFD_BUFFER_SIZE];
      if (num_elements > 0)
        MPI_Allreduce(send_array, _tally_buffer, num_elements,
                      precision, MPI_SUM, comm);
      for (int j=0; j < num_elements; j++)
        _tally_memory[i*CMFD_BUFFER_SIZE+j] = _tally_buffer[j];
    }

    /* Tally MPI communication time */
    _timer->stopTimer();
    _timer->recordSplit("Total MPI communication time");
  }
#endif

  //TODO: clean
  if (neg_fluxes)
      log_printf(WARNING, "Negative FSR fluxes found in group collapse");

  /* Loop over CMFD cells and set cross sections */
#pragma parallel omp for
  for (int i = 0; i < _num_x * _num_y * _num_z; i++) {

    int ind = getLocalCMFDCell(i);

    Material* cell_material = _materials[i];

    /* Loop over CMFD coarse energy groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Load tallies at this cell and energy group */
      FP_PRECISION vol_tally = _volume_tally[i][e];
      FP_PRECISION tot_tally = _total_tally[i][e];
      FP_PRECISION nu_fis_tally = _nu_fission_tally[i][e];
      FP_PRECISION rxn_tally = _reaction_tally[i][e];
      FP_PRECISION neut_prod_tally = _neutron_production_tally[i][e];
      FP_PRECISION* scat_tally = _scattering_tally[i][e];
      FP_PRECISION* chi_tally = _chi_tally[i][e];

      /* Set the Mesh cell properties with the tallies */
      _volumes->setValue(i, 0, vol_tally);
      cell_material->setSigmaTByGroup(tot_tally / rxn_tally, e + 1);
      cell_material->setNuSigmaFByGroup(nu_fis_tally / rxn_tally, e + 1);
      if (ind != -1)
        _old_flux->setValue(ind, e, rxn_tally / vol_tally);
      _old_flux_full->setValue(i, e, rxn_tally / vol_tally);

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



/**
 * @brief Computes the diffusion coefficient for a given CMFD cell and CMFD
 *        energy group.
 * @details This method computes the diffusion coefficient for a CMFD cell and
 *          CMFD energy group by spatially collapsing the total/transport xs
 *          in each FSR contained within the CMFD cell and then energy
 *          collapsing the diffusion coefficient (\f$1 / (3 * \Sigma_t)\f$) for
 *          all MOC groups in the given CMFD energy group.
 * @param cmfd_cell A CMFD cell
 * @param group A CMFD energy group
 * @return The diffusion coefficient
 */
FP_PRECISION Cmfd::getDiffusionCoefficient(int cmfd_cell, int group) {
  return _diffusion_tally[cmfd_cell][group] /
    _reaction_tally[cmfd_cell][group];
}


/**
 * @brief Compute the surface diffusion coefficient for a given CMFD cell,
 *        cell surface, and group.
 * @details This method uses finite differencing to compute the surface
 *          diffusion coefficient (\f$ \hat{D} \f$) or surface diffusion
 *          coefficient correction (\f$ \tilde{D} \f$) for a given CMFD cell,
 *          cell surface, and CMFD energy group. If the MOC iteration is zero,
 *          (\f$ \tilde{D} \f$) is returned as zero. Since (\f$ \hat{D} \f$) and
 *          (\f$ \tilde{D} \f$) are dependent on each other, they must be
 *          computed together; therefore, the boolean correction is used to
 *          indicate which value is to be returned.
 * @param cmfd_cell A CMFD cell
 * @param surface A surface of the CMFD cell
 * @param group A CMFD energy group
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
  FP_PRECISION flux = _old_flux_full->getValue(cmfd_cell, group); //FIXME
  int cmfd_cell_next = getCellNext(cmfd_cell, surface); //FIXME
  FP_PRECISION delta_interface = getSurfaceWidth(surface);
  FP_PRECISION delta = getPerpendicularSurfaceWidth(surface);
  int sense = getSense(surface);

  /* Correct the diffusion coefficient with Larsen's effective diffusion
   * coefficient correction factor */
  if (!_linear_source)
    dif_coef *= computeLarsensEDCFactor(dif_coef, delta);

  /* If surface is on a boundary with REFLECTIVE or VACUUM BCs, choose
   * approipriate BC */
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

      /* Weight the old and new corrected diffusion coefficients by the
         relaxation factor */
      FP_PRECISION old_dif_surf_corr = _old_dif_surf_corr->getValue
          (cmfd_cell, surface*_num_cmfd_groups+group);
      dif_surf_corr = _relaxation_factor * dif_surf_corr +
          (1.0 - _relaxation_factor) * old_dif_surf_corr;
    }
  }

  /* If surface is an interface or PERIODIC BC, use finite differencing */
  else {

    /* Get the surface index for the surface in the neighboring cell */
    int surface_next = (surface + NUM_FACES / 2) % NUM_FACES;

    /* Set diffusion coefficient and flux for the neighboring cell */
    FP_PRECISION dif_coef_next = getDiffusionCoefficient(cmfd_cell_next, group);

    //FIXME
    flux_next = _old_flux_full->getValue(cmfd_cell_next, group);

    /* Correct the diffusion coefficient with Larsen's effective diffusion
     * coefficient correction factor */
    if (!_linear_source)
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

    /* Weight the old and new corrected diffusion coefficients by the
       relaxation factor */
    FP_PRECISION old_dif_surf_corr = _old_dif_surf_corr->getValue
        (cmfd_cell, surface*_num_cmfd_groups+group);
    dif_surf_corr = _relaxation_factor * dif_surf_corr +
        (1.0 - _relaxation_factor) * old_dif_surf_corr;

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

  log_printf(INFO, "Running diffusion solver...");

  /* Start recording total CMFD time */
  _timer->startTimer();

  /* Create matrix and vector objects */
  if (_A == NULL) {
    log_printf(ERROR, "Unable to compute k-eff in CMFD since the CMFD "
               "linear algebra matrices and arrays have not been created.");
  }

  /* Start recording XS collapse time */
  _timer->startTimer();

  /* Collapse the cross sections onto the CMFD mesh */
  collapseXS();

  /* Tally the XS collpase time */
  _timer->stopTimer();
  _timer->recordSplit("Total collapse time");

  /* Construct matrices */
  constructMatrices(moc_iteration);

  /* Copy old flux to new flux */
  _old_flux->copyTo(_new_flux);

  /* Start recording CMFD solve time */
  _timer->startTimer();

  /* Solve the eigenvalue problem */
  /*
  std::cout << "BEFORE" << std::endl;
  _new_flux->printString();
  std::cout << "A....." << std::endl;
  _A->printString();
  if (_domain_communicator != NULL) {
    for (int c = 0; c < 2; c++) {
      int size = _local_num_x * _local_num_y * _local_num_z * _num_cmfd_groups;
      for (int cg=0; cg < size; cg++) {
        int ns = _domain_communicator->num_connections[c][cg];
        for (int s=0; s < ns; s++) {
          std::cout << "Connection with domain " <<
            _domain_communicator->domains[c][cg][s] << " at idx " <<
            _domain_communicator->indexes[c][cg][s] << " for color " << c << " and "
            << "row " << cg << " = " <<
            _domain_communicator->coupling_coeffs[c][cg][s] << std::endl;
        }
      }
    }
  }
  std::cout << "M....." << std::endl;
  _M->printString();
  */
  _k_eff = eigenvalueSolve(_A, _M, _new_flux, _k_eff,
                           _source_convergence_threshold, _SOR_factor,
                           _convergence_data, _domain_communicator);
  /*
  std::cout << "AFTER" << std::endl;
  _new_flux->printString();
  if (moc_iteration >= 0)
    exit(0);
    */

  /* Tally the CMFD solver time */
  _timer->stopTimer();
  _timer->recordSplit("Total solver time");

  /* Rescale the old and new flux */
  rescaleFlux();

  /* Update the MOC flux */
  updateMOCFlux();

  /* Tally the total CMFD time */
  _timer->stopTimer();
  _timer->recordSplit("Total CMFD time");

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

  double new_source_sum = _new_source->getSum();
  double old_source_sum = _old_source->getSum();
#ifdef MPIx
  if (_domain_communicator != NULL) {
    //FIXME
    if (_domain_communicator->_num_domains_x *
        _domain_communicator->_num_domains_y *
        _domain_communicator->_num_domains_z > 1) {
    double temp_sum_new = new_source_sum;
    MPI_Allreduce(&temp_sum_new, &new_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                  _domain_communicator->_MPI_cart);
    double temp_sum_old = old_source_sum;
    MPI_Allreduce(&temp_sum_old, &old_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                  _domain_communicator->_MPI_cart);
    }
  }
#endif
  _new_flux->scaleByValue(1.0 / new_source_sum);
  _old_flux->scaleByValue(1.0 / old_source_sum);
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

  log_printf(INFO,"Constructing matrices...");

  /* Zero _A and _M matrices */
  _A->clear();
  _M->clear();

  //TODO: more elegant
  if (_domain_communicator != NULL) {
    int num_local_cells = _local_num_x * _local_num_y * _local_num_z;
    for (int c=0; c<2; c++) {
      for (int ncg=0; ncg < num_local_cells * _num_cmfd_groups; ncg++) {
        _domain_communicator->num_connections[c][ncg] = 0;
      }
    }
  }
#pragma omp parallel
  {

    FP_PRECISION value, volume, delta;
    FP_PRECISION dif_surf, dif_surf_corr;
    int sense;
    Material* material;

    int num_cells = _num_x * _num_y * _num_z;
    int x_start = 0;
    int y_start = 0;
    int z_start = 0;
    int x_end = _num_x;
    int y_end = _num_y;
    int z_end = _num_z;
    if (_geometry->isDomainDecomposed()) {
      if (_domain_communicator != NULL) {
        x_start = _domain_communicator->_domain_idx_x * _local_num_x;
        x_end = x_start + _local_num_x;
        y_start = _domain_communicator->_domain_idx_y * _local_num_y;
        y_end = y_start + _local_num_y;
        z_start = _domain_communicator->_domain_idx_z * _local_num_z;
        z_end = z_start + _local_num_z;
      }
    }

    /* Loop over cells */
#pragma omp for
    for (int i = 0; i < _num_x*_num_y*_num_z; i++) {

      int ind = getLocalCMFDCell(i);

      if (ind == -1)
        continue;

      //FIXME maybe???
      int color = getCellColor(i);

      material = _materials[i];
      volume = _volumes->getValue(i, 0);

      /* Loop over groups */
      for (int e = 0; e < _num_cmfd_groups; e++) {

        /* Net removal term */
        value = material->getSigmaTByGroup(e+1) * volume;
        _A->incrementValue(ind, e, ind, e, value);

        /* Scattering gain from all groups */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          value = - material->getSigmaSByGroup(g+1, e+1) * volume;
          _A->incrementValue(ind, g, ind, e, value);
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

          /* Record the corrected diffusion coefficient */
          _old_dif_surf_corr->setValue(i, s*_num_cmfd_groups+e, dif_surf_corr);

          /* Set the diagonal term */
          value = (dif_surf - sense * dif_surf_corr) * delta;
          _A->incrementValue(ind, e, ind, e, value);

          /* Set the off diagonal term */
          if (getCellNext(ind, s, false, false) != -1) {
            value = - (dif_surf + sense * dif_surf_corr) * delta;
            _A->incrementValue(getCellNext(ind, s, false, false), e, ind, e, value);
          }
          /* Check for cell in neighboring domain if applicable */
          else if (_geometry->isDomainDecomposed()) {
            if (_domain_communicator != NULL) {
              if (getCellNext(ind, s, false, true) != -1) {
                int neighbor_cell = getCellNext(ind, s, false, true);
                int row = ind * _num_cmfd_groups + e;
                int idx = _domain_communicator->num_connections[color][row];
                value = - (dif_surf + sense * dif_surf_corr) * delta;
                _domain_communicator->indexes[color][row][idx] = neighbor_cell;
                _domain_communicator->domains[color][row][idx] = s;
                _domain_communicator->coupling_coeffs[color][row][idx] = value;
                _domain_communicator->num_connections[color][row]++;
              }
            }
          }
        }

        /* Fission source term */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          value = material->getChiByGroup(e+1)
              * material->getNuSigmaFByGroup(g+1) * volume;
          _M->incrementValue(ind, g, ind, e, value);
        }
      }
    }
  }

  log_printf(INFO, "Done constructing matrices...");
}


/**
 * @brief Update the MOC flux in each FSR.
 * @details This method uses the condensed flux from the last MOC transport
 *          sweep and the converged flux from the eigenvalue problem to
 *          update the MOC flux in each FSR.
 */
void Cmfd::updateMOCFlux() {

  log_printf(INFO, "Updating MOC flux...");

  /* Set max prolongation factor */
  if (_convergence_data != NULL)
    _convergence_data->pf = 1.0;

  /* Loop over mesh cells */
#pragma omp parallel for
  for (int i = 0; i < _num_x * _num_y * _num_z; i++) {

    //FIXME
    int ind = getLocalCMFDCell(i);
    if (ind == -1)
      continue;

    std::vector<int>::iterator iter;

    /* Loop over CMFD groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Loop over FRSs in mesh cell */
      for (iter = _cell_fsrs.at(i).begin();
           iter != _cell_fsrs.at(i).end(); ++iter) {

        /* Get the update ratio */
        /*
        FP_PRECISION update_ratio = getUpdateRatio(i, e, *iter);
              */
        FP_PRECISION update_ratio = _new_flux->getValue(ind, e) /
              _old_flux->getValue(ind, e);
        if (_convergence_data != NULL)
          if (std::abs(log(update_ratio)) > std::abs(log(_convergence_data->pf)))
            _convergence_data->pf = update_ratio;

        for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++) {

          /* Update FSR flux using ratio of old and new CMFD flux */
          _FSR_fluxes[*iter*_num_moc_groups + h] *= update_ratio;

          /* Update flux moments if they were set */
          if (_linear_source) {
            _flux_moments[(*iter)*3*_num_moc_groups + h*3] *= update_ratio;
            _flux_moments[(*iter)*3*_num_moc_groups + h*3 + 1] *= update_ratio;
            _flux_moments[(*iter)*3*_num_moc_groups + h*3 + 2] *= update_ratio;
          }

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
 *          useful in ensuring that the CMFD acceleration equations have a
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
  for (int p = 0; p < _num_polar/2; p++) {
    mu = cos(asin(_quadrature->getSinTheta(0,p)));
    expon = exp(-delta / (3 * dif_coef * mu));
    alpha = (1 + expon) / (1 - expon) - 2 * (3 * dif_coef * mu) / delta;
    rho += 2.0 * mu * _quadrature->getPolarWeight(0,p) * alpha;
  }

  /* Compute the correction factor */
  FP_PRECISION correction = 1.0 + delta * rho / (2 * dif_coef);

  return correction;
}


/**
 * @brief Set the FSR materials array pointer.
 * @param pointer to FSR_materials array
 */
void Cmfd::setFSRMaterials(Material** FSR_materials) {
  _FSR_materials = FSR_materials;
}


/**
 * @brief Set the pointer to the array of FSR_volumes.
 * @param array of FSR volumes
 */
void Cmfd::setFSRVolumes(FP_PRECISION* FSR_volumes) {
  _FSR_volumes = FSR_volumes;
}


/**
 * @brief Set pointer to FSR flux array.
 * @param pointer to FSR flux array
 */
void Cmfd::setFSRFluxes(FP_PRECISION* scalar_flux) {
  _FSR_fluxes = scalar_flux;
}


/**
 * @brief Set pointer to source region flux moments array
 * @param pointer to source region flux moments array
 */
void Cmfd::setFluxMoments(FP_PRECISION* flux_moments) {
  _flux_moments = flux_moments;
  _linear_source = true;
}


/**
 * @brief Set successive over-relaxation relaxation factor.
 * @param over-relaxation factor
 */
void Cmfd::setSORRelaxationFactor(FP_PRECISION SOR_factor) {

  if (SOR_factor <= 0.0 || SOR_factor >= 2.0)
    log_printf(ERROR, "The successive over-relaxation relaxation factor "
                      "must be > 0 and < 2. Input value: %i", SOR_factor);

  _SOR_factor = SOR_factor;
}


/**
 * @brief Set the CMFD relaxation factor applied to diffusion coefficients
 * @param CMFD relaxation factor
 */
void Cmfd::setCMFDRelaxationFactor(FP_PRECISION relaxation_factor) {

  if (relaxation_factor <= 0.0 || relaxation_factor > 1.0)
    log_printf(ERROR, "The successive over-relaxation relaxation factor "
                      "must be greater than 0 and less than or equal to 1. "
                      "Input value: %i", relaxation_factor);

  _relaxation_factor = relaxation_factor;
  if (relaxation_factor != 1.0)
    log_printf(NORMAL, "CMFD relaxation factor: %6.4f", _relaxation_factor);
}


/**
 * @brief Get the number of coarse CMFD energy groups.
 * @return the number of CMFD energy groups
 */
int Cmfd::getNumCmfdGroups() {
  return _num_cmfd_groups;
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

  /* Delete old CMFD surface currents vector if it exists */
  if (_materials != NULL)
    delete [] _materials;

  try{
    _materials = new Material*[_num_x*_num_y*_num_z];

    for (int z = 0; z < _num_z; z++) {
      for (int y = 0; y < _num_y; y++) {
        for (int x = 0; x < _num_x; x++) {
          material = new Material(z*_num_x*_num_y + y*_num_x + x);
          material->setNumEnergyGroups(_num_cmfd_groups);
          _materials[z*_num_x*_num_y + y*_num_x + x] = material;
        }
      }
    }
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the Mesh cell materials. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initializes CMFD surface currents Vector prior to first MOC iteration.
 */
void Cmfd::initializeCurrents() {

  /* Delete old CMFD surface currents vector if it exists */
  if (_surface_currents != NULL)
    delete _surface_currents;

  /* Allocate memory for the CMFD Mesh surface and corner currents Vectors */
  _surface_currents = new Vector(_cell_locks, _num_x, _num_y, _num_z,
                                 _num_cmfd_groups * NUM_FACES);
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
  //FIXME
  for (int z = 0; z < _num_z; z++) {
    for (int y = 0; y < _num_y; y++) {
      for (int x = 0; x < _num_x; x++)
        _cell_fsrs.push_back(std::vector<int>());
    }
  }
}


//TODO DOCUMENT XXX
void Cmfd::allocateTallies() {

  if (_num_x*_num_y*_num_z == 0)
    log_printf(ERROR, "Zero cells in CMFD mesh. Please set CMFD mesh before "
               "initializing CMFD tallies.");

  if (_num_cmfd_groups == 0)
    log_printf(ERROR, "Zero CMFD gropus. Please set CMFD group structure "
               "before initializing CMFD tallies.");

  /* Determine tally sizes */
  int num_cells = _num_x * _num_y * _num_z;
  int integrated_tally_size = num_cells * _num_cmfd_groups;
  int total_integrated_tally_size = 6 * integrated_tally_size;
  int groupwise_tally_size = integrated_tally_size * _num_cmfd_groups;
  int total_groupwise_tally_size = 2 * groupwise_tally_size;
  _total_tally_size = total_integrated_tally_size +
      total_groupwise_tally_size;

  /* Allocate memory for tallies */
  _tally_memory = new FP_PRECISION[_total_tally_size];
  FP_PRECISION** integrated_tallies[6];
  for (int t=0; t<6; t++) {
    integrated_tallies[t] = new FP_PRECISION*[num_cells];
    for (int i=0; i < num_cells; i++) {
      int idx = i*_num_cmfd_groups + t*integrated_tally_size;
      integrated_tallies[t][i] = &_tally_memory[idx];
    }
  }
  FP_PRECISION*** groupwise_tallies[2];
  for (int t=0; t<2; t++) {
    groupwise_tallies[t] = new FP_PRECISION**[num_cells];
    for (int i=0; i < num_cells; i++) {
      groupwise_tallies[t][i] = new FP_PRECISION*[_num_cmfd_groups];
      for (int g=0; g < _num_cmfd_groups; g++) {
        int idx = i*_num_cmfd_groups*_num_cmfd_groups + g*_num_cmfd_groups
              + t*groupwise_tally_size + total_integrated_tally_size;
        groupwise_tallies[t][i][g] = &_tally_memory[idx];
      }
    }
  }

  /* Assign tallies to allocated data */
  _nu_fission_tally = integrated_tallies[0];
  _reaction_tally = integrated_tallies[1];
  _volume_tally = integrated_tallies[2];
  _total_tally = integrated_tallies[3];
  _neutron_production_tally = integrated_tallies[4];
  _diffusion_tally = integrated_tallies[5];
  _scattering_tally =  groupwise_tallies[0];
  _chi_tally =  groupwise_tallies[1];
  _tallies_allocated = true;

  /* Create copy buffer of currents and tallies for domain decomposition */
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    _tally_buffer = new FP_PRECISION[CMFD_BUFFER_SIZE];
    _currents_buffer = new FP_PRECISION[CMFD_BUFFER_SIZE];
  }
#endif
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
 * @brief Find the CMFD surface that a LocalCoords object lies on.
 * @details If the coords is not on a surface, -1 is returned. Otherwise,
 *          the surface ID is returned.
 * @param The CMFD cell ID that the local coords is in.
 * @param The coords being evaluated.
 * @return The surface ID.
 */
int Cmfd::findCmfdSurface(int cell, LocalCoords* coords) {
  Point* point = coords->getHighestLevel()->getPoint();
  return _lattice->getLatticeSurface(cell, point);
}


/**
 * @brief Find the CMFD cell that a LocalCoords object is in.
 * @param The coords being evaluated.
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
 * @param num_z The number of cells in the z direction.
 */
void Cmfd::setLatticeStructure(int num_x, int num_y, int num_z) {
  setNumX(num_x);
  setNumY(num_y);
  setNumZ(num_z);
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
 * @param The CMFD cell ID.
 * @param The FSR ID.
 */
//FIXME
void Cmfd::addFSRToCell(int cmfd_cell, int fsr_id) {
  _cell_fsrs.at(cmfd_cell).push_back(fsr_id);
}


/**
 * @brief Set the number of MOC energy groups.
 * @param number of MOC energy groups
 */
void Cmfd::setNumMOCGroups(int num_groups) {
  _num_moc_groups = num_groups;
}


/**
 * @brief Get the number of MOC energy groups.
 * @return the number of MOC energy groups
 */
int Cmfd::getNumMOCGroups() {
  return _num_moc_groups;
}


/**
 * @brief Get the number of CMFD cells.
 * @return the number of CMFD cells
 */
int Cmfd::getNumCells() {
  return _num_x * _num_y * _num_z;
}


/**
 * @brief set the number of FSRs.
 * @param the number of FSRs
 */
void Cmfd::setNumFSRs(int num_fsrs) {
  _num_FSRs = num_fsrs;
}


/**
 * @brief Split the currents of the Mesh cell vertices to the adjacent faces and
 *        edges.
 * @details This method takes the currents tallied across the vertices of a CMFD
 *          cell and splits them evenly across the adjacent faces and edges. In
 *          order to transport the current through to the diagonal cell, the
 *          current is also tallied on the edges of the adjacent cells.
 *          Essentially, the tracks that cross through vertices are split into
 *          three one-third-weight tracks as shown in the illustration below.
 *          Face crossings are denoted as "o" and edge crossings are denoted as
 *          "e". As shown, each partial-weight track first crosses a face and
 *          then crosses through an edge to get into an adjacent cell. After all
 *          partial-weight tracks reach the diagonal cell, they are recombined
 *          into one full-weight track. Note tracks are not physically split
 *          into partial-weight tracks for ray tracing; rather tracks cross
 *          through the vertices and the current through each vertex is tallied.
 *
 *                                                       . .      .
 *                                                    .   \|  .
 *                                        |    /   .   .   .
 *                                        |   / .   .       \
 *                                        |  e   .           .
 *                                        | / .           .
 *                                     .  e/           .
 *                 x -------------------.-+---------e-----------
 *                               o   o   /|       .
 *                            .   .     / |    .
 *                         .   .       /  | .
 *                          \ |       /  o|
 *                           .       /.   |
 *                        .   \    ./     |
 *                     .        .   y     z
 *                  .
 *
 */
void Cmfd::splitVertexCurrents() {

  log_printf(INFO, "Splitting CMFD vertex currents...");

  int ncg = _num_cmfd_groups;
  int nf = NUM_FACES;
  int ne = NUM_EDGES;
  int ns = NUM_SURFACES;

#pragma omp parallel
  {

    FP_PRECISION current;
    std::vector<int> surfaces;
    std::vector<int>::iterator iter;
    std::map<int, FP_PRECISION>::iterator it;
    int cell, surface;

#pragma omp for
    for (int i=0; i < _num_x * _num_y * _num_z; i++) {
      for (int v = NUM_FACES + NUM_EDGES; v < NUM_SURFACES; v++) {

        /* Check if this surface is contained in the map */
        int ind = i * NUM_SURFACES * ncg + v * ncg;
        it = _edge_corner_currents.find(ind);
        if (it == _edge_corner_currents.end())
          continue;

        getVertexSplitSurfaces(i, v, &surfaces);

        for (int g=0; g > ncg; g++) {
          /* Divide vertex current by 3 since we will split to 3 surfaces,
           * which propagate through 3 edges */
          current = _edge_corner_currents[ind+g] / 3;

          /* Increment current for faces and edges adjacent to this vertex */
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
 * @brief Split the currents of the Mesh cell edges to the adjacent faces.
 * @details This method takes the currents tallied across the edges (or corners)
 *          of a CMFD cell and splits them evenly across the adjacent surfaces
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

  log_printf(INFO, "Splitting CMFD edge currents...");

  int ncg = _num_cmfd_groups;
  int nf = NUM_FACES;
  int ne = NUM_EDGES;
  int ns = NUM_SURFACES;

#pragma omp parallel
  {

    FP_PRECISION current;
    std::vector<int> surfaces;
    std::vector<int>::iterator iter;
    std::map<int, FP_PRECISION>::iterator it;
    int cell, surface;

#pragma omp for
    for (int i=0; i < _num_x * _num_y * _num_z; i++) {
      for (int e = NUM_FACES; e < NUM_SURFACES + NUM_EDGES; e++) {

        /* Check if this surface is contained in the map */
        int ind = i * NUM_SURFACES * ncg + e * ncg;
        it = _edge_corner_currents.find(ind);
        if (it == _edge_corner_currents.end())
          continue;

        int i = ind / (_num_x * _num_y * _num_z);
        getEdgeSplitSurfaces(i, e, &surfaces);

        for (int g=0; g > ncg; g++) {
          /* Divide edge current by 2 since we will split to 2 surfaces,
           * which propagate through 2 surfaces */
          current = _edge_corner_currents[ind+g] / 2;

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
 * @brief Get the faces and edges to split the currents of the Mesh cell
 *        vertices.
 * @details The process by which the current of tracks passing through vertices
 *          is split is described in the comment for
 *          Cmfd::splitVertexCurrents(). This method takes in the cell and
 *          vertex that is being split as well as a std::vector used to store
 *          the IDs of surfaces that are crossed by the partial-weight tracks.
 *          This method properly accounts for crossings on the geometry
 *          boundaries by applying the corresponding boundary conditions to
 *          split the currents.
 * @param cell The CMFD cell ID that the vertex is in.
 * @param vertex The vertex that the track crosses through.
 * @param surfaces A std::vector that is populated with the IDs of surfaces that
 *        are crossed by partial-weight tracks.
 */
void Cmfd::getVertexSplitSurfaces(int cell, int vertex,
                                  std::vector<int>* surfaces) {

  surfaces->clear();
  int x = (cell % (_num_x * _num_y)) % _num_x;
  int y = (cell % (_num_x * _num_y)) / _num_x;
  int z = cell / (_num_x * _num_y);
  int ns = NUM_SURFACES;

  if (vertex == SURFACE_X_MIN_Y_MIN_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN_Z_MIN);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1)) * ns + SURFACE_Y_MIN_Z_MIN);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Y_MIN_Z_MIN);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Z_MIN);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MIN_Z_MIN);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_X_MIN_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Y_MIN);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MIN_Y_MIN);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_X_MIN_Y_MIN);
  }
  else if (vertex == SURFACE_X_MIN_Y_MIN_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN_Z_MAX);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1)) * ns + SURFACE_Y_MIN_Z_MAX);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Y_MIN_Z_MAX);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Z_MAX);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MIN_Z_MAX);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_X_MIN_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Y_MIN);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MIN_Y_MIN);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_X_MIN_Y_MIN);
  }
  else if (vertex == SURFACE_X_MIN_Y_MAX_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX_Z_MIN);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1)) * ns + SURFACE_Y_MAX_Z_MIN);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Y_MAX_Z_MIN);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Z_MIN);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MIN_Z_MIN);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_X_MIN_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Y_MAX);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MIN_Y_MAX);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_X_MIN_Y_MAX);
  }
  else if (vertex == SURFACE_X_MIN_Y_MAX_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX_Z_MAX);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1)) * ns + SURFACE_Y_MAX_Z_MAX);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Y_MAX_Z_MAX);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Z_MAX);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MIN_Z_MAX);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_X_MIN_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN_Y_MAX);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MIN_Y_MAX);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_X_MIN_Y_MAX);
  }
  else if (vertex == SURFACE_X_MAX_Y_MIN_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN_Z_MIN);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1)) * ns + SURFACE_Y_MIN_Z_MIN);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Y_MIN_Z_MIN);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Z_MIN);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MAX_Z_MIN);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_X_MAX_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Y_MIN);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MAX_Y_MIN);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_X_MAX_Y_MIN);
  }
  else if (vertex == SURFACE_X_MAX_Y_MIN_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN_Z_MAX);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1)) * ns + SURFACE_Y_MIN_Z_MAX);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Y_MIN_Z_MAX);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Z_MAX);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MAX_Z_MAX);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_X_MAX_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Y_MIN);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MAX_Y_MIN);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_X_MAX_Y_MIN);
  }
  else if (vertex == SURFACE_X_MAX_Y_MAX_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX_Z_MIN);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1)) * ns + SURFACE_Y_MAX_Z_MIN);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Y_MAX_Z_MIN);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Z_MIN);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MAX_Z_MIN);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_X_MAX_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Y_MAX);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MAX_Y_MAX);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_X_MAX_Y_MAX);
  }
  else if (vertex == SURFACE_X_MAX_Y_MAX_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX_Z_MAX);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1)) * ns + SURFACE_Y_MAX_Z_MAX);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Y_MAX_Z_MAX);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Z_MAX);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1)) * ns
                            + SURFACE_X_MAX_Z_MAX);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_X_MAX_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX_Y_MAX);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MAX_Y_MAX);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_X_MAX_Y_MAX);
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
  int x = (cell % (_num_x * _num_y)) % _num_x;
  int y = (cell % (_num_x * _num_y)) / _num_x;
  int z = cell / (_num_x * _num_y);
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
  else if (edge == SURFACE_X_MIN_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MIN);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1))* ns + SURFACE_Z_MIN);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MIN);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_Z_MIN);
  }
  else if (edge == SURFACE_X_MAX_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MIN);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1))* ns + SURFACE_Z_MIN);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MAX);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_X_MAX);
  }
  else if (edge == SURFACE_X_MIN_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (x == 0) {
      if (_boundaries[SURFACE_X_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MAX);
      else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
        surfaces->push_back((cell + (_num_x - 1))* ns + SURFACE_Z_MAX);
    }
    else
      surfaces->push_back((cell - 1) * ns + SURFACE_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MIN);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MIN);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_X_MIN);
  }
  else if (edge == SURFACE_X_MAX_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_X_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (x == _num_x - 1) {
      if (_boundaries[SURFACE_X_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MAX);
      else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
        surfaces->push_back((cell - (_num_x - 1))* ns + SURFACE_Z_MAX);
    }
    else
      surfaces->push_back((cell + 1) * ns + SURFACE_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_X_MAX);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_X_MAX);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_X_MAX);
  }
  else if (edge == SURFACE_Y_MIN_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MIN);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1))* ns + SURFACE_Z_MIN);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_Y_MIN);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_Y_MIN);
  }
  else if (edge == SURFACE_Y_MAX_Z_MIN) {
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MIN);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MIN);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1))* ns + SURFACE_Z_MIN);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_Z_MIN);

    if (z == 0) {
      if (_boundaries[SURFACE_Z_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX);
      else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_Y_MAX);
    }
    else
      surfaces->push_back((cell - _num_x * _num_y) * ns + SURFACE_Y_MAX);
  }
  else if (edge == SURFACE_Y_MIN_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_Y_MIN);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (y == 0) {
      if (_boundaries[SURFACE_Y_MIN] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MAX);
      else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
        surfaces->push_back((cell + _num_x * (_num_y - 1))* ns + SURFACE_Z_MAX);
    }
    else
      surfaces->push_back((cell - _num_x) * ns + SURFACE_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MIN);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_Y_MIN);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_Y_MIN);
  }
  else if (edge == SURFACE_Y_MAX_Z_MAX) {
    surfaces->push_back(cell * ns + SURFACE_Y_MAX);
    surfaces->push_back(cell * ns + SURFACE_Z_MAX);

    if (y == _num_y - 1) {
      if (_boundaries[SURFACE_Y_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Z_MAX);
      else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * (_num_y - 1))* ns + SURFACE_Z_MAX);
    }
    else
      surfaces->push_back((cell + _num_x) * ns + SURFACE_Z_MAX);

    if (z == _num_z - 1) {
      if (_boundaries[SURFACE_Z_MAX] == REFLECTIVE)
        surfaces->push_back(cell * ns + SURFACE_Y_MAX);
      else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
        surfaces->push_back((cell - _num_x * _num_y * (_num_z - 1)) * ns
                            + SURFACE_Y_MAX);
    }
    else
      surfaces->push_back((cell + _num_x * _num_y) * ns + SURFACE_Y_MAX);
  }
}


/**
 * @brief Get the ID of the Mesh cell next to given Mesh cell.
 * @param current Mesh cell ID
 * @param CMFD cell surface ID to look across for neighboring cell
 * @return neighboring CMFD cell ID
 */
//FIXME: HEADER
int Cmfd::getCellNext(int cell, int surface_id, bool global, bool neighbor) {

  int cell_next = -1;

  int x, y, z;
  int nx, ny, nz;
  int x_global, y_global, z_global;
  if (global || _domain_communicator == NULL) {
    x_global = (cell % (_num_x * _num_y)) % _num_x;
    y_global = (cell % (_num_x * _num_y)) / _num_x;
    z_global = cell / (_num_x * _num_y);
    x = x_global;
    y = y_global;
    z = z_global;
    nx = _num_x;
    ny = _num_y;
    nz = _num_z;
  }
  else {
    x = (cell % (_local_num_x * _local_num_y)) % _local_num_x;
    y = (cell % (_local_num_x * _local_num_y)) / _local_num_x;
    z = cell / (_local_num_x * _local_num_y);
    x_global = x + _domain_communicator->_domain_idx_x * _local_num_x;
    y_global = y + _domain_communicator->_domain_idx_y * _local_num_y;
    z_global = z + _domain_communicator->_domain_idx_z * _local_num_z;
    nx = _local_num_x;
    ny = _local_num_y;
    nz = _local_num_z;
  }

  if (surface_id == SURFACE_X_MIN) {
    if (x != 0)
      cell_next = cell - 1;
    else if (neighbor && !global && x_global != 0)
      cell_next = z * _local_num_y + y;
    else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
      cell_next = cell + (_num_x-1);
  }

  else if (surface_id == SURFACE_Y_MIN) {
    if (y != 0)
      cell_next = cell - nx;
    else if (neighbor && !global && y_global != 0)
      cell_next = z * _local_num_x + x;
    else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
      cell_next = cell + _num_x*(_num_y-1);
  }

  else if (surface_id == SURFACE_Z_MIN) {
    if (z != 0)
      cell_next = cell - nx*ny;
    else if (neighbor && !global && z_global != 0)
      cell_next = y * _local_num_x + x;
    else if (_boundaries[SURFACE_Z_MIN] == PERIODIC)
      cell_next = cell + _num_x*_num_y*(_num_z-1);
  }

  else if (surface_id == SURFACE_X_MAX) {
    if (x != nx - 1)
      cell_next = cell + 1;
    else if (neighbor && !global && x_global != _num_x - 1)
      cell_next = z * _local_num_y + y;
    else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
      cell_next = cell - (_num_x-1);
  }

  else if (surface_id == SURFACE_Y_MAX) {
    if (y != ny - 1)
      cell_next = cell + nx;
    else if (neighbor && !global && y_global != _num_y - 1)
      cell_next = z * _local_num_x + x;
    else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
      cell_next = cell - _num_x*(_num_y-1);
  }

  else if (surface_id == SURFACE_Z_MAX) {
    if (z != nz - 1)
      cell_next = cell + nx*ny;
    else if (neighbor && !global && z_global != _num_z - 1)
      cell_next = y * _local_num_x + x;
    else if (_boundaries[SURFACE_Z_MAX] == PERIODIC)
      cell_next = cell - _num_x*_num_y*(_num_z-1);
  }

  return cell_next;
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
 * @param the CMFD mesh surface ID.
 * @return the boundaryType for the surface.
 */
int Cmfd::getBoundary(int side) {
  return _boundaries[side];
}


/**
 * @brief Return the CMFD cell ID that an FSR lies in.
 * @details Note that a CMFD cell is not an actual Cell object; rather, a CMFD
 *          cell is just a way of describing each of the rectangular regions
 *          that make up a CMFD lattice. CMFD cells are numbered with 0 in the
 *          lower left corner and monotonically increasing from left to right.
 *          from left to right. For example, he indices for a 4 x 4 lattice are:
 *                  12  13  14  15
 *                  8    9  10  11
 *                  4    5   6   7
 *                  0    1   2   3
 * @param The FSR ID.
 * @return The CMFD cell ID. Return -1 if cell is not found.
 */
//FIXME
int Cmfd::convertFSRIdToCmfdCell(int fsr_id) {

  std::vector<int>::iterator iter;
  for (int cell_id=0; cell_id < _num_x*_num_y*_num_z;
      cell_id++) {
    for (iter = _cell_fsrs.at(cell_id).begin();
         iter != _cell_fsrs.at(cell_id).end(); ++iter) {
      if (*iter  == fsr_id)
        return cell_id;
    }
  }

  return -1;
}


//TODO document
int Cmfd::convertGlobalFSRIdToCmfdCell(int global_fsr_id) {

  /* Determine the domain and local FSR ID */
  int cmfd_cell = -1;
  if (!_geometry->isDomainDecomposed()) {
    cmfd_cell = convertFSRIdToCmfdCell(global_fsr_id);
  }
#ifdef MPIx
  else {

    int fsr_id;
    int domain;
    _geometry->getLocalFSRId(global_fsr_id, fsr_id, domain);

    /* Get the FSR centroid in the correct domain */
    int rank;
    MPI_Comm comm = _geometry->getMPICart();
    MPI_Comm_rank(comm, &rank);
    int temp_cmfd_cell = 0;
    if (rank == domain)
      temp_cmfd_cell = convertFSRIdToCmfdCell(fsr_id);

    /* Broadcast the centroid */
    MPI_Allreduce(&temp_cmfd_cell, &cmfd_cell, 1, MPI_INT, MPI_SUM, comm);
  }
#endif
  return cmfd_cell;
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
 * @brief Set the vector of vectors that contains.
 *        the FSRs that lie in each cell.
 * @param Vector of vectors containing FSR IDs in each cell.
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
 * @param Flag saying whether to update MOC flux.
 */
void Cmfd::setFluxUpdateOn(bool flux_update_on) {
  _flux_update_on = flux_update_on;
}


//TODO document
void Cmfd::setConvergenceData(ConvergenceData* convergence_data) {
  _convergence_data = convergence_data;
}


/**
 * @brief Set flag indicating whether to use axial interpolation for update
 *        ratios
 * @param Flag saying whether to use axial interpolation.
 */
void Cmfd::useAxialInterpolation(bool interpolate) {
  _use_axial_interpolation = interpolate;
}


/**
 * @brief Get flag indicating whether to update the MOC flux.
 * @return Flag saying whether to update MOC flux.
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
    log_printf(ERROR, "Unable to set the CMFD source convergence threshold to"
              " %f since the threshold must be positive.", source_thresh);

  _source_convergence_threshold = source_thresh;
}


/**
 * @brief Sets the PolarQuad object in use by the MOC Solver.
 * @param polar_quad a PolarQuad object pointer from the Solver
 */
void Cmfd::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
  _num_polar = quadrature->getNumPolarAngles();
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

  /* Number of cells in stencil */
  int num_cells_in_stencil = 9;

  /* Loop over mesh cells */
  //FIXME
  for (int i = 0; i < _num_x*_num_y*_num_z; i++) {

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
      for (int j=0; j < num_cells_in_stencil; j++)
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

  /* Initialize axial quadratic interpolant values */
  _axial_interpolants.resize(_num_FSRs);
  for (int r=0; r < _num_FSRs; r++) {
    _axial_interpolants.at(r) = new double[3];
    for (int j=0; j < 3; j++)
      _axial_interpolants.at(r)[j] = 0.0;
  }

  /* Compute axial quadratic interpolation values if requested */
  if (_use_axial_interpolation) {

    /* Calculate common factors */
    double dz = _cell_width_z;
    double dz_2 = _cell_width_z * _cell_width_z;

    /* Loop over mesh cells */
    for (int i = 0; i < _local_num_x*_local_num_y*_local_num_z; i++) {

      /* Calculate the CMFD cell z-coordinate */
      int z_ind = i / (_local_num_x * _local_num_y);
      double z_cmfd = (z_ind + 0.5) * _cell_width_z + _lattice->getMinZ();
      if (_domain_communicator != NULL)
        z_cmfd += _domain_communicator->_domain_idx_z * _cell_width_z * _local_num_z;

      /* Loop over FRSs in mesh cell */
      int i_global = getGlobalCMFDCell(i);
      for (fsr_iter = _cell_fsrs.at(i_global).begin();
           fsr_iter != _cell_fsrs.at(i_global).end(); ++fsr_iter) {

        /* Get centroid and calculate relative z-coordinate */
        fsr_id = *fsr_iter;
        Point* centroid = _geometry->getFSRCentroid(fsr_id);
        double zc = (centroid->getZ() - z_cmfd) / dz;
        if (std::abs(zc) > 0.5)
          log_printf(ERROR, "Found FSR %d with z-centroid offset in z "
                     "from CMFD cell %d by %6.4f, whereas the CMFD z-spacing"
                     " is %6.4f. Coordinates: (%6.4f, %6.4f, %6.4f), cmfd z: "
                     "%6.4f", fsr_id, i, zc*dz, dz, centroid->getX(),
                     centroid->getY(), centroid->getZ(), z_cmfd);

        /* Check that the CMFD cell is not an end cell */
        if (z_ind == 0)
          zc -= 1.0;
        else if (z_ind == _local_num_z - 1)
          zc += 1.0;

        /* Calculate components for quadratic interpolation */
        _axial_interpolants.at(fsr_id)[0] = zc * zc/2.0 - zc/2.0 - 1.0/24.0;
        _axial_interpolants.at(fsr_id)[1] = -zc * zc + 26.0/24.0;
        _axial_interpolants.at(fsr_id)[2] = zc * zc/2.0 + zc/2.0 - 1.0/24.0;
      }
    }
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
  int x = (cell_id % (_num_x * _num_y)) % _num_x;
  int y = (cell_id % (_num_x * _num_y)) / _num_x;

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
 * @details This method takes in a CMFD cell, a MOC energy group, and a FSR
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
 * @param cell_id The CMFD cell ID containing the FSR.
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

        ratio += iter->second * getFluxRatio(cell_next_id, group, fsr);
      }
    }

    /* INTERNAL */
    if (_k_nearest_stencils[fsr].size() == 1)
      ratio += getFluxRatio(cell_id, group, fsr);
    else {
      ratio += _k_nearest_stencils[fsr][0].second *
            getFluxRatio(cell_id, group, fsr);
      ratio /= (_k_nearest_stencils[fsr].size() - 1);
    }
  }
  else
    ratio = getFluxRatio(cell_id, group, fsr);

  return ratio;
}


/**
 * @brief Retreives the ratio of pre- and post- CMFD solve fluxes
 * @details The CMFD flux ratio is returned for the given FSR. A quadratic
 *          axial interpolant is used to estimate the value at the FSR.
 * @param cell_id The CMFD cell ID containing the FSR.
 * @param group The CMFD energy group being updated.
 * @param fsr The fsr being updated.
 * @return the ratio of CMFD fluxes
 */
FP_PRECISION Cmfd::getFluxRatio(int cell_id, int group, int fsr) {

  cell_id = getLocalCMFDCell(cell_id);

  double* interpolants = _axial_interpolants.at(fsr);
  double ratio = 1.0;
  if (interpolants[0] != 0 || interpolants[2] != 0) {

    int z_ind = cell_id / (_local_num_x * _local_num_y);
    int cell_mid = cell_id;
    if (z_ind == 0)
      cell_mid += _local_num_x * _local_num_y;
    else if (z_ind == _num_z - 1)
      cell_mid -= _local_num_x * _local_num_y;

    int cell_prev = cell_mid - _local_num_x * _local_num_y;
    int cell_next = cell_mid + _local_num_x * _local_num_y;

    double old_flux_prev = _old_flux->getValue(cell_prev, group);
    double new_flux_prev = _new_flux->getValue(cell_prev, group);

    double old_flux_next = _old_flux->getValue(cell_next, group);
    double new_flux_next = _new_flux->getValue(cell_next, group);

    double old_flux_mid = _old_flux->getValue(cell_mid, group);
    double new_flux_mid = _new_flux->getValue(cell_mid, group);

    double old_flux = interpolants[0] * old_flux_prev +
           interpolants[1] * old_flux_mid +
           interpolants[2] * old_flux_next;
    double new_flux = interpolants[0] * new_flux_prev +
           interpolants[1] * new_flux_mid +
           interpolants[2] * new_flux_next;

    if (old_flux != 0)
      ratio = new_flux / old_flux;

    if (ratio < 0) {

      /* Try a linear interpolation */
      double zc = sqrt(26.0/24.0 - interpolants[1]);
      if (z_ind < _num_z / 2) {
        old_flux = zc * (old_flux_mid - old_flux_prev) + old_flux_mid;
        new_flux = zc * (new_flux_mid - new_flux_prev) + new_flux_mid;
        if (old_flux != 0)
          ratio = new_flux / old_flux;
        else
          ratio = 0;
      }
      else {
        old_flux = zc * (old_flux_next - old_flux_mid) + old_flux_mid;
        new_flux = zc * (new_flux_next - new_flux_mid) + new_flux_mid;
        if (old_flux != 0)
          ratio = new_flux / old_flux;
        else
          ratio = 0;
      }

      /* Fallback: using the cell average flux ratio */
      if (ratio < 0) {
        if (_old_flux->getValue(cell_id, group) != 0.0)
          ratio = _new_flux->getValue(cell_id, group) /
                  _old_flux->getValue(cell_id, group);
        else
          ratio = 0.0;
      }
    }

    return ratio;
  }
  else {
    if (_old_flux->getValue(cell_id, group) != 0.0)
      return _new_flux->getValue(cell_id, group) /
              _old_flux->getValue(cell_id, group);
    else
      return 0.0;
  }
}


/**
 * @brief Get the distances from an FSR centroid to a given CMFD cell.
 * @details This method takes in a FSR centroid, a CMFD cell, and a stencil index
 *          to a cell located in the 9-point stencil encompassing the CMFD
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
 * @param cell_id The CMFD cell containing the FSR.
 * @param stencil_index The index of the cell in the stencil that we want to
 *        get the distance from.
 * @return the distance from the CMFD cell centroid to the FSR centroid.
 */
FP_PRECISION Cmfd::getDistanceToCentroid(Point* centroid, int cell_id,
                                         int stencil_index) {

  int x = (cell_id % (_num_x * _num_y)) % _num_x;
  int y = (cell_id % (_num_x * _num_y)) / _num_x;
  FP_PRECISION dist_x, dist_y;
  bool found = false;
  FP_PRECISION centroid_x = centroid->getX();
  FP_PRECISION centroid_y = centroid->getY();

  /* LOWER LEFT CORNER */
  if (x > 0 && y > 0 && stencil_index == 0) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x - 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y - 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* BOTTOM SIDE */
  else if (y > 0 && stencil_index == 1) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x + 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y - 0.5)*_cell_width_y),2.0);
    found = true;
  }

  /* LOWER RIGHT CORNER */
  else if (x < _num_x - 1 && y > 0 && stencil_index == 2) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x + 1.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y - 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* LEFT SIDE */
  else if (x > 0 && stencil_index == 3) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x - 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y + 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* CURRENT */
  else if (stencil_index == 4) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x + 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y + 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* RIGHT SIDE */
  else if (x < _num_x - 1 && stencil_index == 5) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x + 1.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y + 0.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* UPPER LEFT CORNER */
  else if (x > 0 && y < _num_y - 1 && stencil_index == 6) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x - 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y + 1.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* TOP SIDE */
  else if (y < _num_y - 1 && stencil_index == 7) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x + 0.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y + 1.5)*_cell_width_y), 2.0);
    found = true;
  }

  /* UPPER RIGHT CORNER */
  else if (x < _num_x - 1 && y < _num_y - 1 && stencil_index == 8) {
    dist_x = pow(centroid_x - (-_width_x/2.0 + (x + 1.5)*_cell_width_x), 2.0);
    dist_y = pow(centroid_y - (-_width_y/2.0 + (y + 1.5)*_cell_width_y), 2.0);
    found = true;
  }

  if (found)
    return pow(dist_x + dist_y, 0.5);
  else
    return std::numeric_limits<FP_PRECISION>::max();
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
    log_printf(ERROR, "Unable to set CMFD k-nearest to %i. k-nearest "
               "must be between 1 and 9.", k_nearest);
  else
    _k_nearest = k_nearest;
}


/**
 * @brief Zero the surface currents for each mesh cell and energy group.
 */
void Cmfd::zeroCurrents() {
  _surface_currents->clear();
}



/**
 * @brief Initialize the Matrix and Vector objects, k-nearest stencils, the
 *        CMFD cell currents and MOC materials.
 */
void Cmfd::initialize() {

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
  if (_old_dif_surf_corr != NULL)
    delete _old_dif_surf_corr;
  if (_volumes != NULL)
    delete _volumes;
  if (_cell_locks != NULL)
    delete [] _cell_locks;

  /* Calculate the number of elements */
  int num_cells = _local_num_x * _local_num_y * _local_num_z;
  int ncg = _num_cmfd_groups;

  try {

    /* Allocate array of OpenMP locks for each CMFD cell */
    _cell_locks = new omp_lock_t[_num_x * _num_y * _num_z]; //FIXME

    /* Loop over all cells to initialize OpenMP locks */
#pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_x*_num_y*_num_z; r++)
      omp_init_lock(&_cell_locks[r]);

    /* Allocate memory for matrix and vector objects */
    _M = new Matrix(_cell_locks, _local_num_x, _local_num_y, _local_num_z,
                    ncg);
    _A = new Matrix(_cell_locks, _local_num_x, _local_num_y, _local_num_z,
                    ncg);
    _old_source = new Vector(_cell_locks, _local_num_x, _local_num_y,
                             _local_num_z, ncg);
    _new_source = new Vector(_cell_locks, _local_num_x, _local_num_y,
                             _local_num_z, ncg);
    _old_flux = new Vector(_cell_locks, _local_num_x, _local_num_y,
                           _local_num_z, ncg);
    _old_flux_full = new Vector(_cell_locks, _num_x, _num_y, _num_z, ncg);
    _new_flux = new Vector(_cell_locks, _local_num_x, _local_num_y,
                           _local_num_z, ncg);
    _old_dif_surf_corr = new Vector(_cell_locks, _num_x, _num_y, _num_z,
                                    NUM_FACES * ncg);
    _old_dif_surf_corr->setAll(0.0);
    _volumes = new Vector(_cell_locks, _num_x, _num_y, _num_z, 1);

    /* Initialize k-nearest stencils, currents, flux, and materials */
    generateKNearestStencils();
    initializeCurrents();
    initializeMaterials();
    allocateTallies();


    /* TODO: document, clean */
    if (_domain_communicator != NULL) {
      _local_num_x = _num_x / _domain_communicator->_num_domains_x;
      _local_num_y = _num_y / _domain_communicator->_num_domains_y;
      _local_num_z = _num_z / _domain_communicator->_num_domains_z;
      int offset = _domain_communicator->_domain_idx_x * _local_num_x +
                    _domain_communicator->_domain_idx_y * _local_num_y +
                    _domain_communicator->_domain_idx_z * _local_num_z;
      _domain_communicator->_offset = offset;
      _domain_communicator->_local_num_x = _local_num_x;
      _domain_communicator->_local_num_y = _local_num_y;
      _domain_communicator->_local_num_z = _local_num_z;
      _domain_communicator->num_groups = ncg;

      int dir_sizes[3] = {num_cells / _local_num_x,  num_cells / _local_num_y,
                          num_cells / _local_num_z};

      _domain_communicator->num_connections = new int*[2];
      _domain_communicator->indexes = new int**[2];
      _domain_communicator->domains = new int**[2];
      _domain_communicator->fluxes = new FP_PRECISION**[2];
      _domain_communicator->coupling_coeffs = new FP_PRECISION**[2];
      _domain_communicator->buffer = new FP_PRECISION*[NUM_FACES];
      for (int rb=0; rb<2; rb++) {
        _domain_communicator->num_connections[rb] = new int[num_cells*ncg];
        _domain_communicator->indexes[rb] = new int*[num_cells*ncg];
        _domain_communicator->domains[rb] = new int*[num_cells*ncg];
        _domain_communicator->fluxes[rb] = new FP_PRECISION*[NUM_FACES];
        _domain_communicator->coupling_coeffs[rb] =
                            new FP_PRECISION*[num_cells*ncg];

        for (int coord=0; coord < 3; coord++) {
          for (int d=0; d < 2; d++) {
            int surf = coord + 3 * d;
            _domain_communicator->fluxes[rb][surf] =
                                new FP_PRECISION[dir_sizes[coord]*ncg];
            _domain_communicator->buffer[surf] =
                                new FP_PRECISION[2*dir_sizes[coord]*ncg];
          }
        }
        for (int nsc=0; nsc < num_cells * ncg; nsc++) {
          _domain_communicator->num_connections[rb][nsc] = 0;
          _domain_communicator->indexes[rb][nsc] = new int[NUM_FACES];
          _domain_communicator->domains[rb][nsc] = new int[NUM_FACES];
          _domain_communicator->coupling_coeffs[rb][nsc] =
                              new FP_PRECISION[NUM_FACES];
        }
      }
    }
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
  _lattice->setNumZ(_num_z);
  _lattice->setWidth(_cell_width_x, _cell_width_y, _cell_width_z);
  _lattice->setOffset(offset->getX(), offset->getY(), offset->getZ());
}


/**
 * @brief Returns the width of a given surface
 * @param surface A surface index, from 0 to NUM_FACES - 1
 * @return The surface width
 */
FP_PRECISION Cmfd::getSurfaceWidth(int surface) {

  FP_PRECISION width;

  if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
    width = _cell_width_y * _cell_width_z;
  else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX)
    width = _cell_width_x * _cell_width_z;
  else
    width = _cell_width_x * _cell_width_y;

  return width;
}


/**
 * @brief Returns the width of the surface perpendicular to a given surface
 * @param surface A surface index, from 0 to NUM_FACES - 1
 * @return The perpendicular surface width
 */
FP_PRECISION Cmfd::getPerpendicularSurfaceWidth(int surface) {

  if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
    return _cell_width_x;
  else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX)
    return _cell_width_y;
  else
    return _cell_width_z;
}


/**
 * @brief Returns the sense of a given surface
 * @details The sense of minimum surfaces (e.g. SURFACE_X_MIN) is defined to be
 *          -1 while maximum surfaces (e.g. SURFACE_X_MAX) are defined to have a
 *          sense of +1. This is based on the current exiting a cell from a
 *          minimum surface being in the direction of negative net current and
 *          the current leaving a cell from a maximum surface being in the
 *          direction of positive net current.
 * @param surface A surface index, from 0 to NUM_FACES - 1
 * @return The sense of the surface
 */
int Cmfd::getSense(int surface) {

  if (surface == SURFACE_X_MIN || surface == SURFACE_Y_MIN ||
      surface == SURFACE_Z_MIN)
    return -1;
  else
    return 1;
}


/**
 * @brief Sets a flag to indicate whether a 2D or 3D problem is being solved.
 * @param solve_3D A boolean indicate whether a 2D or 3D problem is being
 *        solved.
 */
void Cmfd::setSolve3D(bool solve_3D) {
  _solve_3D = solve_3D;
}


/**
 * @brief Sets the azimuthal spacings.
 * @param azim_spacings An array of azimuthal spacings for each azimuthal angle.
 * @param num_azim the number of azimuthal angles.
 */
void Cmfd::setAzimSpacings(FP_PRECISION* azim_spacings, int num_azim) {

  if (_azim_spacings != NULL)
    delete [] _azim_spacings;

  _azim_spacings = new FP_PRECISION[num_azim/4];

  for (int a=0; a < num_azim/4; a++)
    _azim_spacings[a] = FP_PRECISION(azim_spacings[a]);
}


/**
 * @brief Sets the polar spacings.
 * @param polar_spacings A 2D array of polar spacings for each azimuthal and
 *        polar angle combination.
 * @param num_azim the number of azimuthal angles.
 * @param num_polar the number of polar angles.
 */
void Cmfd::setPolarSpacings(FP_PRECISION** polar_spacings, int num_azim,
                            int num_polar) {

  if (_polar_spacings != NULL) {
    for (int a=0; a < num_azim/4; a++)
      delete [] _polar_spacings[a];
    delete [] _polar_spacings;
  }

  _polar_spacings = new FP_PRECISION*[num_azim/4];
  for (int a=0; a < num_azim/4; a++)
    _polar_spacings[a] = new FP_PRECISION[num_polar/2];

  for (int a=0; a < num_azim/4; a++) {
    for (int p=0; p < num_polar/2; p++)
      _polar_spacings[a][p] = FP_PRECISION(polar_spacings[a][p]);
  }
}


//TODO: document
void Cmfd::printTimerReport() {

  std::string msg_string;

  /* Get the total CMFD time */
  double tot_time = _timer->getSplit("Total CMFD time");
  msg_string = "Total CMFD computation time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);

  /* Get the total XS collapse time */
  double xs_collapse_time = _timer->getSplit("Total collapse time");
  msg_string = "Total XS collapse time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), xs_collapse_time);

  /* Get the total communication time */
  double comm_time = _timer->getSplit("Total MPI communication time");
  msg_string = "Total MPI communication time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), comm_time);

  /* Get the total solver time */
  double solver_time = _timer->getSplit("Total solver time");
  msg_string = "Total CMFD solver time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), solver_time);
}


//TODO: Document + FIXME
void Cmfd::checkNeutronBalance() {

  /* Initialize variables */
  omp_lock_t* cell_locks = _old_flux->getCellLocks();
  int num_rows = _old_flux->getNumRows();
  int num_x = _old_flux->getNumX();
  int num_y = _old_flux->getNumY();
  int num_z = _old_flux->getNumZ();
  int num_groups = _old_flux->getNumGroups();
  Vector m_phi(cell_locks, num_x, num_y, num_z, num_groups);
  Vector a_phi(cell_locks, num_x, num_y, num_z, num_groups);

  /* Compute CMFD balance */
  matrixMultiplication(_A, _old_flux, &a_phi);
  matrixMultiplication(_M, _old_flux, &m_phi);

  double max_imbalance = 0.0;
  int max_imbalance_cell = -1;
  int max_imbalance_grp = -1;

  /* Loop over CMFD cells */
  for (int i = 0; i < _num_x * _num_y * _num_z; i++) {

    int x = (i % (_num_x * _num_y)) % _num_x;
    int y = (i % (_num_x * _num_y)) / _num_x;
    int z = i / (_num_x * _num_y);

    Material* cell_material = _materials[i];

    /* Loop over CMFD coarse energy groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Initialize tallies */
      double total = 0.0;
      double in_scattering = 0.0;
      double fission = 0.0;

      /* Loop over FSRs in CMFD cell */
      for (int j = 0; j < _cell_fsrs.at(i).size(); j++) {

        int fsr_id = _cell_fsrs.at(i).at(j);
        Material* fsr_material = _FSR_materials[fsr_id];
        double volume = _FSR_volumes[fsr_id];
        FP_PRECISION* scat = fsr_material->getSigmaS();
        FP_PRECISION* flux = &_FSR_fluxes[fsr_id*_num_moc_groups];

        /* Loop over MOC energy groups within this CMFD coarse group */
        double chi = 0.0;
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++)
          chi += fsr_material->getChiByGroup(h+1);

        /* Calculate total fission and in-scattering in the FSR */
        double tot_fission = 0.0;
        for (int g = 0; g < _num_moc_groups; g++) {

          /* Tally total fission */
          double nu_fis = fsr_material->getNuSigmaFByGroup(g+1);
          tot_fission += nu_fis * flux[g] * volume;

          /* Loop over MOC energy groups within this CMFD coarse group */
          for (int h = _group_indices[e]; h < _group_indices[e+1]; h++)
            in_scattering += scat[h*_num_moc_groups + g] * flux[g] * volume;
        }

        /* Calculate fission contribution to this CMFD coarse group */
        fission += chi * tot_fission / _k_eff;

        /* Calcualte total reaction rate in this CMFD coarse group */
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {
          double tot = fsr_material->getSigmaTByGroup(h+1);
          total += tot * flux[h] * volume;
        }
      }

      /* Calculate net current out of the cell */
      double net_current = 0.0;
      for (int s = 0; s < NUM_FACES; s++) {
        int idx = s * _num_cmfd_groups + e;
        int cmfd_cell_next = getCellNext(i, s); //FIXME
        int surface_next = (s + NUM_FACES / 2) % NUM_FACES;
        int idx_next = surface_next * _num_cmfd_groups + e;
        if (cmfd_cell_next == -1) {
          if (_boundaries[s] == VACUUM)
            net_current += _surface_currents->getValue(i, idx);
        }
        else {
          net_current += _surface_currents->getValue(i, idx) -
              _surface_currents->getValue(cmfd_cell_next,idx_next);
        }
      }

      double moc_balance = in_scattering + fission - total - net_current;

      double cmfd_balance = m_phi.getValue(i, e) / _k_eff -
            a_phi.getValue(i, e);

      double tmp_imbalance = std::max(std::abs(moc_balance),
                                      std::abs(cmfd_balance));
      if (tmp_imbalance > max_imbalance){
        max_imbalance = tmp_imbalance;
        max_imbalance_cell = i;
        max_imbalance_grp = e;
      }

      if (std::abs(moc_balance - cmfd_balance) > 1e-14) {
        log_printf(NORMAL, "MOC neutron balance in cell (%d, %d, %d) for CMFD "
                   "group %d = %g", x, y, z, e, moc_balance);

        log_printf(NORMAL, "CMFD neutron balance in cell (%d, %d, %d) for CMFD "
                   "group %d = %g", x, y, z, e, cmfd_balance);

        log_printf(NORMAL, "Net neutron balance in cell (%d, %d, %d) for CMFD "
                   "group %d = %g", x, y, z, e, moc_balance - cmfd_balance);
      }
    }
  }
  log_printf(NORMAL, "Maximum neutron imbalance of %g at cell %i and group "
             "%d.", max_imbalance, max_imbalance_cell, max_imbalance_grp);
}



//TODO: REMOVE
int Cmfd::getLocalCMFDCell(int cmfd_cell) {

  int x_start = 0;
  int y_start = 0;
  int z_start = 0;
  int x_end = _num_x;
  int y_end = _num_y;
  int z_end = _num_z;
  if (_geometry->isDomainDecomposed()) {
    if (_domain_communicator != NULL) {
      x_start = _domain_communicator->_domain_idx_x * _local_num_x;
      x_end = x_start + _local_num_x;
      y_start = _domain_communicator->_domain_idx_y * _local_num_y;
      y_end = y_start + _local_num_y;
      z_start = _domain_communicator->_domain_idx_z * _local_num_z;
      z_end = z_start + _local_num_z;
    }
  }

  int ix = (cmfd_cell % (_num_x * _num_y)) % _num_x;
  int iy = (cmfd_cell % (_num_x * _num_y)) / _num_x;
  int iz = cmfd_cell / (_num_x * _num_y);

  int local_cmfd_cell;
  if (ix < x_start || ix >= x_end || iy < y_start || iy >= y_end ||
      iz < z_start || iz >= z_end)
    local_cmfd_cell = -1;
  else
    local_cmfd_cell = ((iz - z_start) * _local_num_y + iy - y_start) * _local_num_x
                      + ix - x_start;
  return local_cmfd_cell;
}

//TODO: REMOVE
int Cmfd::getGlobalCMFDCell(int cmfd_cell) {

  int x_start = 0;
  int y_start = 0;
  int z_start = 0;
  int x_end = _num_x;
  int y_end = _num_y;
  int z_end = _num_z;
  if (_geometry->isDomainDecomposed()) {
    if (_domain_communicator != NULL) {
      x_start = _domain_communicator->_domain_idx_x * _local_num_x;
      x_end = x_start + _local_num_x;
      y_start = _domain_communicator->_domain_idx_y * _local_num_y;
      y_end = y_start + _local_num_y;
      z_start = _domain_communicator->_domain_idx_z * _local_num_z;
      z_end = z_start + _local_num_z;
    }
  }

  int ix = (cmfd_cell % (_local_num_x * _local_num_y)) % _local_num_x;
  int iy = (cmfd_cell % (_local_num_x * _local_num_y)) / _local_num_x;
  int iz = cmfd_cell / (_local_num_x * _local_num_y);

  return ((iz + z_start) * _num_y + iy + y_start) * _num_x
                + ix + x_start;
}


//TODO: REMOVE
int Cmfd::getCellColor(int cmfd_cell) {
  int ix = (cmfd_cell % (_num_x * _num_y)) % _num_x;
  int iy = (cmfd_cell % (_num_x * _num_y)) / _num_x;
  int iz = cmfd_cell / (_num_x * _num_y);
  int color = (ix + iy + iz) % 2;
  return color;
}
