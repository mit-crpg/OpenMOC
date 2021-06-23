#include "Cmfd.h"

/**
 * @brief Constructor initializes boundaries and variables that describe
 *          the CMFD object.
 * @details The constructor initializes the many variables that describe
 *          the CMFD mesh and are used to solve the nonlinear diffusion
 *          acceleration problem.
 */
Cmfd::Cmfd() {

  /* Initialize Geometry and Mesh-related attribute */
  _quadrature = NULL;
  _geometry = NULL;
  _materials = NULL;

  /* Communication structs */
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
  _num_unbounded_iterations = 0;
  _flux_update_on = true;
  _centroid_update_on = false;
  _use_axial_interpolation = 0;
  _flux_limiting = true;
  _balance_sigma_t = false;
  _k_nearest = 1;
  _SOR_factor = 1.5;
  _num_FSRs = 0;
#ifndef THREED
  _SOLVE_3D = false;
#endif
  _total_tally_size = 0;
  _tallies_allocated = false;
  _domain_communicator_allocated = false;
  _linear_source = false;
  _old_dif_surf_valid = false;
  _non_uniform = false;
  _widths_adjusted_for_domains = false;

  /* Additional debug output */
  _check_neutron_balance = false;
  _print_cmfd_prolongation_ratios = false;

  /* Energy group and polar angle problem parameters */
  _num_moc_groups = 0;
  _num_cmfd_groups = 0;
  _num_backup_groups = 1;
  _num_polar = 0;
  _num_azim = 0;

  /* Set matrices and arrays to NULL */
  _A = NULL;
  _M = NULL;
  _moc_iteration = 0;
  _k_eff = 1.0;
  _relaxation_factor = 0.7;
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
  _starting_currents = NULL;
  _net_currents = NULL;
  _full_surface_currents = NULL;
  _cell_locks = NULL;
  _volumes = NULL;
  _lattice = NULL;
  _azim_spacings = NULL;
  _polar_spacings = NULL;
  _backup_cmfd = NULL;
  _cmfd_group_to_backup_group = NULL;
  _backup_group_structure.resize(0);

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

  if (_starting_currents != NULL)
    delete _starting_currents;

  if (_net_currents != NULL)
    delete _net_currents;

  if (_full_surface_currents != NULL)
    delete _full_surface_currents;

  if (_old_dif_surf_corr != NULL)
    delete _old_dif_surf_corr;

  if (_volumes != NULL)
    delete _volumes;

  if (_azim_spacings != NULL)
    delete [] _azim_spacings;

  if (_polar_spacings != NULL) {
    for (int a=0; a < _num_azim / 4; a++)
      delete [] _polar_spacings[a];
    delete [] _polar_spacings;
  }

  /* Delete CMFD materials array */
  if (_materials != NULL) {
    for (int i=0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      delete _materials[i];
    delete [] _materials;
  }

  /* Delete the CMFD lattice */
  if (_lattice != NULL)
    delete _lattice;

  /* Clear the _cell_fsrs vector of vectors */
  std::vector< std::vector<long> >::iterator iter1;
  for (iter1 = _cell_fsrs.begin(); iter1 != _cell_fsrs.end(); ++iter1)
    iter1->clear();
  _cell_fsrs.clear();

  /* Clear the _k_nearest_stencils map of vectors */
  std::map<long, std::vector< std::pair<int, double> > >::iterator iter2;
  for (iter2 = _k_nearest_stencils.begin(); iter2 != _k_nearest_stencils.end();
       ++iter2)
    iter2->second.clear();
  _k_nearest_stencils.clear();

  /* Delete tally information */
  if (_tallies_allocated) {

    delete [] _tally_memory;

    delete [] _reaction_tally;
    delete [] _volume_tally;
    delete [] _diffusion_tally;
  }

#ifdef MPIx
  /* De-allocate domain communicator */
  if (_domain_communicator != NULL) {
    if (_domain_communicator_allocated) {
      for (int rb=0; rb<2; rb++) {
        for (int f=0; f < NUM_FACES; f++) {
          delete [] _domain_communicator->indexes[rb][f];
          delete [] _domain_communicator->domains[rb][f];
          delete [] _domain_communicator->coupling_coeffs[rb][f];
          delete [] _domain_communicator->fluxes[rb][f];
        }
        delete [] _domain_communicator->num_connections[rb];
        delete [] _domain_communicator->indexes[rb];
        delete [] _domain_communicator->domains[rb];
        delete [] _domain_communicator->coupling_coeffs[rb];
        delete [] _domain_communicator->fluxes[rb];
      }

      delete [] _domain_communicator->num_connections;
      delete [] _domain_communicator->indexes;
      delete [] _domain_communicator->domains;
      delete [] _domain_communicator->fluxes;
      delete [] _domain_communicator->coupling_coeffs;
      delete _domain_communicator;

      delete [] _inter_domain_data;
      for (int s=0; s < NUM_FACES; s++) {
        delete [] _boundary_volumes[s];
        delete [] _boundary_reaction[s];
        delete [] _boundary_diffusion[s];
        delete [] _old_boundary_flux[s];
        delete [] _boundary_surface_currents[s];
      }

      delete [] _boundary_volumes;
      delete [] _boundary_reaction;
      delete [] _boundary_diffusion;
      delete [] _old_boundary_flux;
      delete [] _boundary_surface_currents;
    }
  }
#endif

  for (long r=0; r < _axial_interpolants.size(); r++)
    delete [] _axial_interpolants.at(r);

  if (_backup_cmfd != NULL)
    delete _backup_cmfd;

  delete _timer;
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
  if (!_widths_adjusted_for_domains) {
    _local_num_x = _num_x;
    if (_domain_communicator != NULL)
      _local_num_x = _num_x / _domain_communicator->_num_domains_x;
  }
}


/**
 * @brief Set the number of Mesh cells in a column.
 * @param num_y number of Mesh cells in a column
 */
void Cmfd::setNumY(int num_y) {

  if (num_y < 1)
    log_printf(ERROR, "The number of lattice cells in the y direction "
               "must be > 0. Input value: %i", num_y);

  _num_y = num_y;
  if (!_widths_adjusted_for_domains) {
    _local_num_y = _num_y;
    if (_domain_communicator != NULL)
      _local_num_y = _num_y / _domain_communicator->_num_domains_y;
  }
}


/**
 * @brief Set the number of Mesh cells in the z-direction.
 * @param num_z number of Mesh cells in the z direction
 */
void Cmfd::setNumZ(int num_z) {

  if (num_z < 1)
    log_printf(ERROR, "The number of lattice cells in the z direction "
               "must be > 0. Input value: %i", num_z);

  _num_z = num_z;
  if (!_widths_adjusted_for_domains) {
    _local_num_z = _num_z;
    if (_domain_communicator != NULL)
      _local_num_z = _num_z / _domain_communicator->_num_domains_z;
  }
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


/**
 * @brief Get the number of Mesh cells in the z-direction
 * @return number of Mesh cells in the z-direction
 */
int Cmfd::getNumZ() {
  return _num_z;
}


/**
 * @brief Get the number of Mesh cells in the z-direction in the local domain
 * @return number of Mesh cells in the z-direction in the domain
 */
int Cmfd::getLocalNumZ() {
  return _local_num_z;
}


/**
 * @brief Get the Vector of surface currents.
 * @return pointer to a vector containing the surface currents
 */
Vector* Cmfd::getLocalCurrents() {
  return _surface_currents;
}


/**
 * @brief Get the array of surface currents on the boundaries.
 * @return 3D array containing the boundary surface currents
 */
CMFD_PRECISION*** Cmfd::getBoundarySurfaceCurrents() {
  return _boundary_surface_currents;
}


/**
 * @brief Set Mesh width in the x-direction
 * @param width Physical width of Mesh in the x-direction
 */
void Cmfd::setWidthX(double width) {
  _width_x = width;
}


/**
 * @brief Set Mesh width in the y-direction
 * @param width Physical width of Mesh in the y-direction
 */
void Cmfd::setWidthY(double width) {
  _width_y = width;
}


/**
 * @brief Set Mesh width in the z-direction
 * @param width Physical width of Mesh in the z-direction
 */
void Cmfd::setWidthZ(double width) {
  _width_z = width;
}


#ifdef MPIx
/**
 * @brief Set the number of domains in each direction.
 * @param num_x number of domains in the X direction
 * @param num_y number of domains in the Y direction
 * @param num_z number of domains in the Z direction
 */
void Cmfd::setNumDomains(int num_x, int num_y, int num_z) {

  if (_domain_communicator == NULL) {
    _domain_communicator = new DomainCommunicator;
    _domain_communicator->_MPI_cart = _geometry->getMPICart();
  }

  _domain_communicator->_num_domains_x = num_x;
  _domain_communicator->_num_domains_y = num_y;
  _domain_communicator->_num_domains_z = num_z;

  _accumulate_lmx.resize(num_x + 1, 0);
  _accumulate_lmy.resize(num_y + 1, 0);
  _accumulate_lmz.resize(num_z + 1, 0);

  std::vector<std::pair<int, double> > divisions_missing_x;
  std::vector<std::pair<int, double> > divisions_missing_y;
  std::vector<std::pair<int, double> > divisions_missing_z;

  /* Find the position of domain decomposition interfaces among the non-uniform
     CMFD mesh cell boundaries, in the X direction */
  int j, j_prev;
  for (int i=0; i<num_x; i++) {
    double coord = (i + 1) * _width_x / num_x;
    for (j=1; j<_num_x+1; j++) {

      /* Keep track of index in mesh before domain boundary */
      if (_accumulate_x[j] < coord)
        j_prev = j;

      /* Exit loop if division is found */
      if (fabs(coord - _accumulate_x[j]) < FLT_EPSILON) {
        _accumulate_lmx[i+1] = j;
        break;
      }
    }
    if (j == _num_x+1)
      divisions_missing_x.push_back(std::make_pair(j_prev, coord));
  }

  /* Find the position of domain decomposition interfaces among the non-uniform
     CMFD mesh cell boundaries, in the Y direction */
  for (int i=0; i<num_y; i++) {
    double coord = (i + 1) * _width_y / num_y;
    for (j=1; j<_num_y+1; j++) {

      /* Keep track of index in mesh before domain boundary */
      if (_accumulate_y[j] < coord)
        j_prev = j;

      /* Exit loop if division is found */
      if (fabs(coord - _accumulate_y[j]) < FLT_EPSILON) {
        _accumulate_lmy[i+1] = j;
        break;
      }
    }
    if (j == _num_y+1)
      divisions_missing_y.push_back(std::make_pair(j_prev, coord));
  }

  /* Find the position of domain decomposition interfaces among the non-uniform
     CMFD mesh cell boundaries, in the Z direction */
  for (int i=0; i<num_z; i++) {
    double coord = (i + 1) * _width_z / num_z;
    for (j=1; j<_num_z+1; j++) {

      /* Keep track of index in mesh before domain boundary */
      if (_accumulate_z[j] < coord)
        j_prev = j;

      /* Exit loop if division is found */
      if (fabs(coord - _accumulate_z[j]) < FLT_EPSILON) {
        _accumulate_lmz[i+1] = j;
        break;
      }
    }
    if (j == _num_z+1)
      divisions_missing_z.push_back(std::make_pair(j_prev, coord));
  }

  /* Output the missing subdivisions */
  std::string div_miss_x = "";
  std::string div_miss_y = "";
  std::string div_miss_z = "";
  std::vector<std::pair<int, double> >::iterator iter;
  for (iter = divisions_missing_x.begin(); iter != divisions_missing_x.end();
       ++iter) {
    div_miss_x += std::to_string(iter->second);
    if (iter != std::prev(divisions_missing_x.end()))
      div_miss_x += ", ";
  }
  for (iter = divisions_missing_y.begin(); iter != divisions_missing_y.end();
       ++iter) {
    div_miss_y += std::to_string(iter->second);
    if (iter != std::prev(divisions_missing_y.end()))
      div_miss_y += ", ";
  }
  for (iter = divisions_missing_z.begin(); iter != divisions_missing_z.end();
       ++iter) {
    div_miss_z += std::to_string(iter->second);
    if (iter != std::prev(divisions_missing_z.end()))
      div_miss_z += ", ";
  }

  if (divisions_missing_x.size() > 0)
    log_printf(WARNING_ONCE, "Domain boundaries [%s] are not in CMFD X-mesh",
               div_miss_x.c_str());
  if (divisions_missing_y.size() > 0)
    log_printf(WARNING_ONCE, "Domain boundaries [%s] are not in CMFD Y-mesh",
               div_miss_y.c_str());
  if (divisions_missing_z.size() > 0)
    log_printf(WARNING_ONCE, "Domain boundaries [%s] are not in CMFD Z-mesh",
               div_miss_z.c_str());

  /* Modify the CMFD mesh, in order to have divisions at domain boundaries,
     but also to avoid small CMFD cells */
  std::string added_cells_x = "";
  std::string modified_cells_x = "";
  int cells_added = 0;
  for (iter = divisions_missing_x.begin(); iter != divisions_missing_x.end();
       ++iter) {

    j_prev = iter->first;
    double cmfd_cell_width = _accumulate_x[j_prev+1] - _accumulate_x[j_prev];
    double coord = iter->second;

    /* If the new CMFD cells are about to be very small,
       move the CMFD cell boundaries instead, to the closest one */
    double delta_below = coord - _accumulate_x[j_prev];
    double delta_above = _accumulate_x[j_prev+1] - coord;

    if (delta_below < 0.1 * cmfd_cell_width ||
        delta_above < 0.1 * cmfd_cell_width) {

      modified_cells_x += (std::to_string(coord) + " ");
      if (delta_below < delta_above) {
        _cell_widths_x[j_prev + cells_added - 1] += delta_below;
        _cell_widths_x[j_prev + cells_added] -= delta_below;
      }
      else {
        _cell_widths_x[j_prev + cells_added] -= delta_above;
        _cell_widths_x[j_prev + cells_added + 1] += delta_above;
        //FIXME Capture j_prev+1 overflow
      }
    }

    /* Else, add a new subdivision to the CMFD mesh */
    else {
      _cell_widths_x[j_prev + cells_added] = delta_below;
      _cell_widths_x.insert(_cell_widths_x.begin() + j_prev + cells_added + 1,
                            delta_above);
      added_cells_x += (std::to_string(coord) + " ");
      cells_added++;
    }
  }
  if (added_cells_x.compare("") != 0)
    log_printf(WARNING_ONCE, "New CMFD cells created for domain decomposition "
               "boundaries in the X direction at [%s] cm.",
               added_cells_x.c_str());
  if (modified_cells_x.compare("") != 0)
    log_printf(WARNING_ONCE, "CMFD mesh cell widths adjusted to fit "
               "domain decomposition boundaries in the X direction "
               "at [%s] cm.", modified_cells_x.c_str());

  /* Y-direction */
  std::string added_cells_y = "";
  std::string modified_cells_y = "";
  cells_added = 0;
  for (iter = divisions_missing_y.begin(); iter != divisions_missing_y.end();
       ++iter) {

    j_prev = iter->first;
    double cmfd_cell_width = _accumulate_y[j_prev+1] - _accumulate_y[j_prev];
    double coord = iter->second;

    /* If the new CMFD cells are about to be very small,
       move the CMFD cell boundaries instead, to the closest one */
    double delta_below = coord - _accumulate_y[j_prev];
    double delta_above = _accumulate_y[j_prev+1] - coord;

    if (delta_below < 0.1 * cmfd_cell_width ||
        delta_above < 0.1 * cmfd_cell_width) {

      modified_cells_y += (std::to_string(coord) + " ");
      if (delta_below < delta_above) {
        _cell_widths_y[j_prev + cells_added - 1] += delta_below;
        _cell_widths_y[j_prev + cells_added] -= delta_below;
      }
      else {
        _cell_widths_y[j_prev + cells_added] -= delta_above;
        _cell_widths_y[j_prev + cells_added + 1] += delta_above;
        //FIXME Capture j_prev+1 overflow
      }
    }

    /* Else, add a new subdivision to the CMFD mesh */
    else {
      _cell_widths_y[j_prev + cells_added] = delta_below;
      _cell_widths_y.insert(_cell_widths_y.begin() + j_prev + cells_added + 1,
                            delta_above);
      added_cells_y += (std::to_string(coord) + " ");
      cells_added++;
    }
  }
  if (added_cells_y.compare("") != 0)
    log_printf(WARNING_ONCE, "New CMFD cells created for domain decomposition "
               "boundaries in the Y direction at [%s] cm.",
               added_cells_y.c_str());
  if (modified_cells_y.compare("") != 0)
    log_printf(WARNING_ONCE, "CMFD mesh cell widths adjusted to fit "
               "domain decomposition boundaries in the Y direction "
               "at [%s] cm.", modified_cells_y.c_str());

  /* Z-direction */
  std::string added_cells_z = "";
  std::string modified_cells_z = "";
  cells_added = 0;
  for (iter = divisions_missing_z.begin(); iter != divisions_missing_z.end();
       ++iter) {

    j_prev = iter->first;
    double cmfd_cell_width = _accumulate_z[j_prev+1] - _accumulate_z[j_prev];
    double coord = iter->second;

    /* If the new CMFD cells are about to be very small,
       move the CMFD cell boundaries instead, to the closest one */
    double delta_below = coord - _accumulate_z[j_prev];
    double delta_above = _accumulate_z[j_prev+1] - coord;

    if (delta_below < 0.1 * cmfd_cell_width ||
        delta_above < 0.1 * cmfd_cell_width) {

      modified_cells_z += (std::to_string(coord) + " ");
      if (delta_below < delta_above) {
        _cell_widths_z[j_prev + cells_added - 1] += delta_below;
        _cell_widths_z[j_prev + cells_added] -= delta_below;
      }
      else {
        _cell_widths_z[j_prev + cells_added] -= delta_above;
        _cell_widths_z[j_prev + cells_added + 1] += delta_above;
        //FIXME Capture j_prev+1 overflow
      }
    }

    /* Else, add a new subdivision to the CMFD mesh */
    else {
      _cell_widths_z[j_prev + cells_added] = delta_below;
      _cell_widths_z.insert(_cell_widths_z.begin() + j_prev + cells_added + 1,
                            delta_above);
      added_cells_z += (std::to_string(coord) + " ");
      cells_added++;
    }
  }
  if (added_cells_z.compare("") != 0)
    log_printf(WARNING_ONCE, "New CMFD cells created for domain decomposition "
               "boundaries in the Z direction at [%s] cm.",
               added_cells_z.c_str());
  if (modified_cells_z.compare("") != 0)
    log_printf(WARNING_ONCE, "CMFD mesh cell widths adjusted to fit "
               "domain decomposition boundaries in the Z direction "
               "at [%s] cm.", modified_cells_z.c_str());

  /* Re-initialize lattice since widths have been modified */
  if (divisions_missing_x.size() > 0 || divisions_missing_y.size() > 0 ||
      divisions_missing_z.size() > 0) {
    _non_uniform = true;
    Point offset;
    offset.setXYZ(_lattice->getOffset()->getXYZ());
    initializeLattice(&offset);

    /* Re-run routine to obtain _accumulate_lmxyz arrays */
    if (!_widths_adjusted_for_domains) {
      _widths_adjusted_for_domains = true;
      setNumDomains(num_x, num_y, num_z);
    }
    else
      log_printf(ERROR, "Mesh adjustment to domain decomposition boundaries "
                 "failed. Manually add the domain boundaries in a non-"
                 "uniform CMFD mesh.");

  }
}


/**
 * @brief Set the indexes of the domain among the global lattice of domains.
 */
void Cmfd::setDomainIndexes(int idx_x, int idx_y, int idx_z) {

  if (_domain_communicator == NULL) {
    _domain_communicator = new DomainCommunicator;
    _domain_communicator->_MPI_cart = _geometry->getMPICart();
  }

  _domain_communicator->_domain_idx_x = idx_x;
  _domain_communicator->_domain_idx_y = idx_y;
  _domain_communicator->_domain_idx_z = idx_z;

  _local_num_x = _accumulate_lmx[idx_x + 1] - _accumulate_lmx[idx_x];
  _local_num_y = _accumulate_lmy[idx_y + 1] - _accumulate_lmy[idx_y];
  _local_num_z = _accumulate_lmz[idx_z + 1] - _accumulate_lmz[idx_z];
}
#endif


/**
 * @brief Collapse cross-sections and fluxes for each CMFD cell by
 *        energy condensing and volume averaging cross sections from
 *        the MOC sweep.
 * @details This method performs a cell-wise energy condensation and flux-volume
 *          average of the cross sections of the fine, unstructured FSR mesh.
 *          The cross sections are condensed such that all reaction rates and
 *          the neutron production rate from fission are conserved. It is
 *          important to note that the volume averaging is performed before
 *          energy condensation in order to properly collapse the diffusion
 *          coefficients.
 */
void Cmfd::collapseXS() {

  log_printf(INFO, "Collapsing cross-sections onto CMFD mesh...");

  /* Record net currents over cells if neutron balance of sigma-t requested */
  if (_balance_sigma_t)
    recordNetCurrents();

  /* Check to see that CMFD tallies have been allocated */
  if (!_tallies_allocated)
    log_printf(ERROR, "Tallies need to be allocated before collapsing "
               "cross-sections");

  /* Split vertex and edge currents to side surfaces */
  splitVertexCurrents();
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    communicateSplits(false);
#endif
  splitEdgeCurrents();
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    communicateSplits(true);
#endif

  /* Report number of negative currents */
  long num_negative_currents = _surface_currents->getNumNegativeValues();
  int total_negative_CMFD_current_domains = (num_negative_currents > 0);
#ifdef MPIx
  if (_domain_communicator != NULL) {
    long temp_sum_neg = num_negative_currents;
    MPI_Allreduce(&temp_sum_neg, &num_negative_currents, 1, MPI_LONG, MPI_SUM,
                  _domain_communicator->_MPI_cart);
    int temp_sum_dom = total_negative_CMFD_current_domains;
    MPI_Allreduce(&temp_sum_dom, &total_negative_CMFD_current_domains, 1,
                  MPI_INT, MPI_SUM, _domain_communicator->_MPI_cart);
  }
#endif

  if (_SOLVE_3D && num_negative_currents > 0)
    log_printf(WARNING_ONCE, "Negative CMFD currents in %ld surfaces-groups in"
               " %d domains.", num_negative_currents,
               total_negative_CMFD_current_domains);

#pragma omp parallel
  {

    /* Initialize variables for FSR properties */
    FP_PRECISION volume, flux;
    FP_PRECISION tot, nu_fis, chi;
    FP_PRECISION* scat;

    /* Allocate arrays for tallies on the stack */
    double scat_tally[_num_cmfd_groups];
    double chi_tally[_num_cmfd_groups];

    /* Pointers to material objects */
    Material* fsr_material;
    Material* cell_material;

    /* Loop over CMFD cells */
#pragma omp for
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

      std::vector<long>::iterator iter;
      cell_material = _materials[i];

      /* Zero group-wise fission terms */
      double neutron_production_tally = 0.0;
      for (int e = 0; e < _num_cmfd_groups; e++)
        chi_tally[e] = 0.0;

      /* Loop over FSRs in CMFD cell */
      for (iter = _cell_fsrs.at(i).begin();
           iter != _cell_fsrs.at(i).end(); ++iter) {

        fsr_material = _FSR_materials[*iter];
        volume = _FSR_volumes[*iter];

        /* Calculate total neutron production in the FSR */
        double neutron_production = 0.0;
        for (int h = 0; h < _num_moc_groups; h++)
          neutron_production += fsr_material->getNuSigmaFByGroup(h+1) *
              _FSR_fluxes[(*iter)*_num_moc_groups+h] * volume;

        /* Calculate contribution to all CMFD groups */
        for (int e=0; e < _num_cmfd_groups; e++) {
          chi = 0;
          for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++)
            chi += fsr_material->getChiByGroup(h+1);

          chi_tally[e] += chi * neutron_production;
        }

        /* Add to total neutron production within the CMFD cell */
        neutron_production_tally += neutron_production;
      }

      /* Set chi */
      if (fabs(neutron_production_tally) > 0) {

        /* Calculate group-wise fission contributions */
        for (int e=0; e < _num_cmfd_groups; e++)
          cell_material->setChiByGroup(chi_tally[e] / neutron_production_tally,
                                       e + 1);
      }
      else {
        /* Calculate group-wise chi to zero */
        for (int e=0; e < _num_cmfd_groups; e++)
          cell_material->setChiByGroup(0.0, e + 1);
      }

      /* Loop over CMFD coarse energy groups */
      for (int e = 0; e < _num_cmfd_groups; e++) {

        /* Zero tallies for this group */
        double nu_fission_tally = 0.0;
        double total_tally = 0.0;

        _diffusion_tally[i][e] = 0.0;
        _reaction_tally[i][e] = 0.0;
        _volume_tally[i][e] = 0.0;

        /* Zero each group-to-group scattering tally */
        for (int g = 0; g < _num_cmfd_groups; g++)
          scat_tally[g] = 0.0;

        /* Loop over MOC energy groups within this CMFD coarse group */
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {

          /* Reset volume tally for this MOC group */
          _volume_tally[i][e] = 0.0;
          double rxn_tally_group = 0.0;
          double trans_tally_group = 0.0;

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
            total_tally += tot * flux * volume;
            nu_fission_tally += nu_fis * flux * volume;
            _reaction_tally[i][e] += flux * volume;
            _volume_tally[i][e] += volume;

            /* Increment diffusion MOC group-wise tallies */
            rxn_tally_group += flux * volume;
            trans_tally_group += tot * flux * volume;

            /* Scattering tallies */
            for (int g = 0; g < _num_moc_groups; g++) {
              scat_tally[getCmfdGroup(g)] +=
                  scat[g*_num_moc_groups+h] * flux * volume;
            }
          }

          /* Condense diffusion coefficient (with homogenized transport XS) */
          if (fabs(trans_tally_group) > fabs(rxn_tally_group) * FLT_EPSILON) {
            CMFD_PRECISION flux_avg_sigma_t = trans_tally_group /
                rxn_tally_group;
            _diffusion_tally[i][e] += rxn_tally_group /
                (3.0 * flux_avg_sigma_t);
          }
        }

        /* Save cross-sections to material */
        double rxn_tally = _reaction_tally[i][e];

        if (rxn_tally <= 0) {
          int cell = getGlobalCMFDCell(i);
          int x = (cell % (_num_x * _num_y)) % _num_x;
          int y = (cell % (_num_x * _num_y)) / _num_x;
          int z = cell / (_num_x * _num_y);
          log_printf(WARNING, "Negative or zero reaction tally calculated in "
                     "CMFD cell %d [%d %d %d] in CMFD group %d : %e", cell, x,
                     y, z, e + 1, rxn_tally);

          /* Set all cross sections to be 1 */
          rxn_tally = ZERO_SIGMA_T;
          _reaction_tally[i][e] = ZERO_SIGMA_T;
          _diffusion_tally[i][e] = ZERO_SIGMA_T;
          total_tally = ZERO_SIGMA_T;
          if (nu_fission_tally != 0)
            nu_fission_tally = ZERO_SIGMA_T;

          /* Avoid excessive downscatter */
          for (int g = 0; g < _num_cmfd_groups; g++)
            scat_tally[g] = 0;
        }

        cell_material->setSigmaTByGroup(total_tally / rxn_tally, e + 1);
        cell_material->setNuSigmaFByGroup(nu_fission_tally / rxn_tally, e + 1);

        /* Set scattering xs */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, e + 1,
                                          g + 1);
        }
      }
    }
  }

#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    if (_domain_communicator != NULL) {

      /* Start recording MPI communication time */
      _timer->startTimer();

      /* Do the Ghost cell exchange */
      ghostCellExchange();

      /* Tally the MPI communication time */
      _timer->stopTimer();
      _timer->recordSplit("CMFD MPI communication time");
    }
  }
#endif

  /* Calculate (local) old fluxes and set volumes */
#pragma omp parallel for
  for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

    /* Loop over CMFD coarse energy groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Load tallies at this cell and energy group */
      double vol_tally = _volume_tally[i][e];
      double rxn_tally = _reaction_tally[i][e];
      _old_flux->setValue(i, e, rxn_tally / vol_tally);

      /* Set the Mesh cell properties with the tallies */
      _volumes->setValue(i, 0, vol_tally);
    }
  }

  /* Loop over boundary CMFD cells and set cross sections */
  if (_geometry->isDomainDecomposed()) {
#pragma omp parallel for
    for (int s=0; s < NUM_FACES; s++) {

      /* Loop over all CMFD cells on the current surface */
      std::map<int, int>::iterator it;
      for (it=_boundary_index_map.at(s).begin();
          it != _boundary_index_map.at(s).end(); ++it) {

        int idx = it->second;

        /* Loop over CMFD coarse energy groups */
        for (int e = 0; e < _num_cmfd_groups; e++) {

          /* Load tallies at this cell and energy group */
          double vol_tally = _boundary_volumes[s][idx][0];
          double rxn_tally = _boundary_reaction[s][idx][e];
          _old_boundary_flux[s][idx][e] = rxn_tally / vol_tally;
        }
      }
    }
  }

  /* Output max CMFD cell optical thickness for stability concerns */
  float max_tau = -1;
  for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++) {
    int global_cmfd_cell = getGlobalCMFDCell(i);
    CMFD_PRECISION delta_x = getPerpendicularSurfaceWidth(0, global_cmfd_cell);
    CMFD_PRECISION delta_y = getPerpendicularSurfaceWidth(1, global_cmfd_cell);
    CMFD_PRECISION delta_z = -1;
    if (_SOLVE_3D)
      delta_z = getPerpendicularSurfaceWidth(2, global_cmfd_cell);
    max_tau = std::max(max_tau, float(_materials[i]->getMaxSigmaT() *
              std::max(delta_x, std::max(delta_y, delta_z))));
  }

  float global_max_tau = max_tau;
#ifdef MPIx
  if (_domain_communicator != NULL)
    MPI_Allreduce(&global_max_tau, &max_tau, 1, MPI_FLOAT, MPI_MAX,
                  _domain_communicator->_MPI_cart);
#endif
  if (_moc_iteration == 0)
    log_printf(INFO_ONCE, "Max CMFD optical thickness in all domains %.2e",
               global_max_tau);
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
CMFD_PRECISION Cmfd::getDiffusionCoefficient(int cmfd_cell, int group) {
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
 *          computed together.
 * @param cmfd_cell A CMFD cell
 * @param surface A surface of the CMFD cell
 * @param group A CMFD energy group
 * @param dif_surf the surface diffusion coefficient \f$ \hat{D} \f$
 * @param dif_surf_corr the correction diffusion coefficient \f$ \tilde{D} \f$
 */
void Cmfd::getSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                          int group, CMFD_PRECISION& dif_surf,
                                          CMFD_PRECISION& dif_surf_corr) {

  FP_PRECISION current, current_out, current_in;
  CMFD_PRECISION flux_next;

  /* Get diffusivity and flux for Mesh cell */
  CMFD_PRECISION dif_coef = getDiffusionCoefficient(cmfd_cell, group);
  int global_cmfd_cell = getGlobalCMFDCell(cmfd_cell);
  int global_cmfd_cell_next = getCellNext(global_cmfd_cell, surface);
  CMFD_PRECISION flux = _old_flux->getValue(cmfd_cell, group);
  CMFD_PRECISION delta_interface = getSurfaceWidth(surface, global_cmfd_cell);
  CMFD_PRECISION delta = getPerpendicularSurfaceWidth(surface, global_cmfd_cell);

  CMFD_PRECISION delta_next = 0.0;
  if (global_cmfd_cell_next != -1)
    delta_next = getPerpendicularSurfaceWidth(surface, global_cmfd_cell_next);

  /* Correct the diffusion coefficient with Larsen's effective diffusion
   * coefficient correction factor */
  if (!_linear_source)
    dif_coef *= computeLarsensEDCFactor(dif_coef, delta);

  /* If surface is on a boundary with REFLECTIVE or VACUUM BCs, choose
   * appropriate BC */
  if (global_cmfd_cell_next == -1) {

    /* REFLECTIVE BC */
    if (_boundaries[surface] == REFLECTIVE) {
      dif_surf = 0.0;
      dif_surf_corr = 0.0;
    }

    /* VACUUM BC */
    else if (_boundaries[surface] == VACUUM) {

      /* Compute the surface-averaged current leaving the cell */
      current_out = _surface_currents->getValue
          (cmfd_cell, surface*_num_cmfd_groups + group) / delta_interface;

      /* Set the surface diffusion coefficient and MOC correction */
      dif_surf =  2 * dif_coef / delta / (1 + 4 * dif_coef / delta);
      dif_surf_corr = (dif_surf * flux - current_out) / flux;
    }
  }

  /* If surface is an interface or PERIODIC BC, use finite differencing */
  else {

    /* Get the surface index for the surface in the neighboring cell */
    int surface_next = (surface + NUM_FACES / 2) % NUM_FACES;

    /* Get the outward current on surface */
    current_out = _surface_currents->getValue
        (cmfd_cell, surface*_num_cmfd_groups + group);

    /* Set diffusion coefficient and flux for the neighboring cell */
    int cmfd_cell_next = getLocalCMFDCell(global_cmfd_cell_next);
    CMFD_PRECISION dif_coef_next;
    if (cmfd_cell_next == -1) {

      /* Get the currents in cells touching this boundary */
      CMFD_PRECISION** boundary_currents = _boundary_surface_currents[surface];

      int idx = _boundary_index_map.at(surface)[global_cmfd_cell_next];
      dif_coef_next = _boundary_diffusion[surface][idx][group] /
                _boundary_reaction[surface][idx][group];
      flux_next = _old_boundary_flux[surface][idx][group];

      /* Get the inward current on the surface */
      current_in = boundary_currents[idx][surface_next*_num_cmfd_groups+group];
    }
    else {

      dif_coef_next = getDiffusionCoefficient(cmfd_cell_next, group);
      flux_next = _old_flux->getValue(cmfd_cell_next, group);

      /* Get the inward current on the surface */
      current_in = _surface_currents->getValue
          (cmfd_cell_next, surface_next*_num_cmfd_groups + group);
    }

    /* Correct the diffusion coefficient with Larsen's effective diffusion
     * coefficient correction factor */
    if (!_linear_source)
      dif_coef_next *= computeLarsensEDCFactor(dif_coef_next, delta_next);

    /* Compute the surface diffusion coefficient */
    dif_surf = 2.0 * dif_coef * dif_coef_next
               / (delta_next * dif_coef + delta * dif_coef_next);

    /* Compute the surface-averaged net current across the surface */
    current = (current_out - current_in) / delta_interface;

    /* Compute the surface diffusion coefficient correction */
    dif_surf_corr = -(dif_surf * (flux_next - flux) + current)
        / (flux_next + flux);

    /* Flux limiting condition */
    if (_flux_limiting && _moc_iteration > 0) {
      double ratio = dif_surf_corr / dif_surf;
      if (std::abs(ratio) > 1.0) {

        if (current > 0.0)
          dif_surf = std::abs(current / (2.0*flux));
        else
          dif_surf = std::abs(current / (2.0*flux_next));

        dif_surf_corr = -(dif_surf * (flux_next - flux) + current)
                        / (flux_next + flux);

        /* Make sure diffusion coefficient is larger than the corrected one,
           to floating point precision */
        dif_surf = std::max(dif_surf, std::abs(dif_surf_corr));
      }
    }
  }
  /* Weight the old and new corrected diffusion coefficients by the
     relaxation factor */
  if (_old_dif_surf_valid) {
    CMFD_PRECISION old_dif_surf_corr = _old_dif_surf_corr->getValue
        (cmfd_cell, surface*_num_cmfd_groups+group);
    dif_surf_corr = _relaxation_factor * dif_surf_corr +
        (1.0 - _relaxation_factor) * old_dif_surf_corr;
  }

  /* If it is the first MOC iteration, solve the straight diffusion problem
   * with no MOC correction */
  if (_moc_iteration == 0 && !_check_neutron_balance)
    dif_surf_corr = 0.0;
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
double Cmfd::computeKeff(int moc_iteration) {

  log_printf(INFO, "Running diffusion solver...");

  /* Start recording total CMFD time */
  _timer->startTimer();

  /* Save MOC iteration number */
  _moc_iteration = moc_iteration;

  /* Create matrix and vector objects */
  if (_A == NULL) {
    log_printf(ERROR, "Unable to compute k-eff in CMFD since the CMFD "
               "linear algebra matrices and arrays have not been created.");
  }

  /* Start recording XS collapse time */
  _timer->startTimer();

  /* Copy surface currents if neutron balance check requested */
  if (_check_neutron_balance)
    copyFullSurfaceCurrents();

  /* Collapse the cross sections onto the CMFD mesh */
  collapseXS();

  /* Tally the XS collapse time */
  _timer->stopTimer();
  _timer->recordSplit("Total collapse time");

  /* Construct matrices and record time */
  _timer->startTimer();
  constructMatrices();
  _timer->stopTimer();
  _timer->recordSplit("Matrix construction time");

  /* Check neutron balance if requested */
  if (_check_neutron_balance)
    checkNeutronBalance(false, false);

  /* Copy old flux to new flux to use collapsed flux as a starting guess */
  _old_flux->copyTo(_new_flux);

  /* Start recording CMFD solve time */
  _timer->startTimer();

  /* Solve the eigenvalue problem */
  double k_eff = eigenvalueSolve(_A, _M, _new_flux, _k_eff,
                                 _source_convergence_threshold, _SOR_factor,
                                 _convergence_data, _domain_communicator);

  /* Try to use a few-group solver to remedy convergence issues */
  bool reduced_group_solution = false;
  if (fabs(k_eff + 1) < FLT_EPSILON && _num_cmfd_groups > _num_backup_groups) {

    log_printf(NORMAL, "Switching to a %d group CMFD solver on this iteration",
               _num_backup_groups);

    if (_backup_cmfd == NULL)
      initializeBackupCmfdSolver();

    copyCurrentsToBackup();
    k_eff = _backup_cmfd->computeKeff(moc_iteration);
    reduced_group_solution = true;
  }

  /* Tally the CMFD solver time */
  _timer->stopTimer();
  _timer->recordSplit("Total solver time");

  /* Check for a legitimate solve */
  if (fabs(k_eff + 1) > FLT_EPSILON)
    _k_eff = k_eff;
  else
    return _k_eff;

  /* Do not prolong if the few-group solution was used */
  if (reduced_group_solution) {
    log_printf(NORMAL, "The %d group CMFD solver was successful",
               _num_backup_groups);
    return _k_eff;
  }

  /* Rescale the old and new flux */
  rescaleFlux();

  /* Update the MOC flux */
  _timer->startTimer();
  if (_flux_update_on)
    updateMOCFlux();

  /* Tally the update and total CMFD time */
  _timer->stopTimer();
  _timer->recordSplit("Total MOC flux update time");
  _timer->stopTimer();
  _timer->recordSplit("Total CMFD time");

  /* If debugging, print CMFD prolongation factors */
  if (get_log_level() == DEBUG || _print_cmfd_prolongation_ratios)
    printProlongationFactors();

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
    double temp_sum_new = new_source_sum;
    MPI_Allreduce(&temp_sum_new, &new_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                  _domain_communicator->_MPI_cart);
    double temp_sum_old = old_source_sum;
    MPI_Allreduce(&temp_sum_old, &old_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                  _domain_communicator->_MPI_cart);
  }
#endif
  _new_flux->scaleByValue(1.0 / new_source_sum);
  _old_flux->scaleByValue(1.0 / old_source_sum);

  /* Check for negative fluxes in CMFD flux */
  long num_negative_fluxes = _new_flux->getNumNegativeValues();
  int total_negative_CMFD_flux_domains = (num_negative_fluxes > 0);
#ifdef MPIx
  if (_domain_communicator != NULL) {
    long temp_sum_neg = num_negative_fluxes;
    MPI_Allreduce(&temp_sum_neg, &num_negative_fluxes, 1, MPI_LONG, MPI_SUM,
                  _domain_communicator->_MPI_cart);
    int temp_sum_dom = total_negative_CMFD_flux_domains;
    MPI_Allreduce(&temp_sum_dom, &total_negative_CMFD_flux_domains, 1,
                  MPI_INT, MPI_SUM, _domain_communicator->_MPI_cart);
  }
#endif

  if (num_negative_fluxes > 0)
    log_printf(WARNING_ONCE, "Negative CMFD fluxes in %ld cells on %d domains.",
               num_negative_fluxes, total_negative_CMFD_flux_domains);
}


/**
 * @brief Construct the loss + streaming matrix (A) and the fission gain
 *         matrix (M) in preparation for solving the eigenvalue problem.
 * @details This method loops over all mesh cells and energy groups and
 *          accumulates the iteraction and streaming terms into their
 *          appropriate positions in the loss + streaming matrix and
 *          fission gain matrix.
 */
void Cmfd::constructMatrices() {

  log_printf(INFO, "Constructing matrices...");

  /* Zero _A and _M matrices */
  _A->clear();
  _M->clear();

  /* Zero the number of connections */
  if (_domain_communicator != NULL) {
    int num_boundary_cells = 2 * (_local_num_x * _local_num_z + _local_num_y *
                             _local_num_z + _local_num_x * _local_num_y);
    for (int c=0; c<2; c++) {
      for (int ncg=0; ncg < num_boundary_cells * _num_cmfd_groups; ncg++) {
        _domain_communicator->num_connections[c][ncg] = 0;
      }
    }
  }
#pragma omp parallel
  {

    FP_PRECISION value, volume, delta;
    CMFD_PRECISION dif_surf, dif_surf_corr;
    int sense;
    Material* material;

    int x_start = 0;
    int y_start = 0;
    int z_start = 0;
    int x_end = _num_x;
    int y_end = _num_y;
    int z_end = _num_z;
    if (_geometry->isDomainDecomposed()) {
      if (_domain_communicator != NULL) {
        x_start = _accumulate_lmx[_domain_communicator->_domain_idx_x];
        x_end = x_start + _local_num_x;
        y_start = _accumulate_lmy[_domain_communicator->_domain_idx_y];
        y_end = y_start + _local_num_y;
        z_start = _accumulate_lmz[_domain_communicator->_domain_idx_z];
        z_end = z_start + _local_num_z;
      }
    }

    /* Loop over cells */
#pragma omp for
    for (int i = 0; i < _local_num_x*_local_num_y*_local_num_z; i++) {

      int global_ind = getGlobalCMFDCell(i);
      int color = getCellColor(global_ind);

      material = _materials[i];
      volume = _volumes->getValue(i, 0);

      /* Find neighboring cells on other domains */
      int domain_surface_index = -1;
      int neighbor_cells[NUM_FACES];
      if (_geometry->isDomainDecomposed() && _domain_communicator != NULL) {

        bool on_surface = false;

        for (int s = 0; s < NUM_FACES; s++) {
          neighbor_cells[s] = -1;
          if (getCellNext(i, s, false, false) == -1) {
            neighbor_cells[s] = getCellNext(i, s, false, true);
            if (neighbor_cells[s] != -1)
              on_surface = true;
          }
        }

        /* If the cell is on a surface, find the index into the comm buffers */
        if (on_surface)
          domain_surface_index = _domain_communicator->mapLocalToSurface[i];
      }

      /* Loop over groups */
      for (int e = 0; e < _num_cmfd_groups; e++) {

        /* Net removal term */
        value = material->getSigmaTByGroup(e+1) * volume;
        _A->incrementValue(i, e, i, e, value);

        /* Re-compute diagonal if neutron re-balance requested */
        if (_balance_sigma_t) {
          enforceBalanceOnDiagonal(i, e);
        }

        /* Scattering gain from all groups */
        for (int g = 0; g < _num_cmfd_groups; g++) {
          value = - material->getSigmaSByGroup(g+1, e+1) * volume;
          if (std::abs(value) > FLT_EPSILON)
            _A->incrementValue(i, g, i, e, value);
        }

        /* Streaming to neighboring cells */
        for (int s = 0; s < NUM_FACES; s++) {

          delta = getSurfaceWidth(s, global_ind);

          /* Set transport term on diagonal */
          getSurfaceDiffusionCoefficient(i, s, e, dif_surf, dif_surf_corr);

          /* Record the corrected diffusion coefficient */
          _old_dif_surf_corr->setValue(i, s*_num_cmfd_groups+e, dif_surf_corr);

          /* Set the diagonal term */
          value = (dif_surf - dif_surf_corr) * delta;
          _A->incrementValue(i, e, i, e, value);

          /* Set the off diagonal term */
          if (getCellNext(i, s, false, false) != -1) {
            value = - (dif_surf + dif_surf_corr) * delta;
            _A->incrementValue(getCellNext(i, s, false, false), e, i, e, value);
          }

          /* Check for cell in neighboring domain if applicable */
          else if (_geometry->isDomainDecomposed()) {
            if (_domain_communicator != NULL) {
              if (neighbor_cells[s] != -1) {
                int neighbor_cell = neighbor_cells[s];
                int row = domain_surface_index * _num_cmfd_groups + e;
                //FIXME Make num_connections, indexes and domains array not
                // group dependent
                int idx = _domain_communicator->num_connections[color][row];
                value = - (dif_surf + dif_surf_corr) * delta;
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
          if (std::abs(value) > FLT_EPSILON)
            _M->incrementValue(i, g, i, e, value);
        }
      }
    }
  }

  /* Mark correction diffusion coefficient as valid for relaxation purposes */
  _old_dif_surf_valid = true;
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
  for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

    CMFD_PRECISION thread_max_update_ratio = 1;
    std::vector<long>::iterator iter;

    /* Loop over CMFD groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Loop over FRSs in mesh cell */
      for (iter = _cell_fsrs.at(i).begin();
           iter != _cell_fsrs.at(i).end(); ++iter) {

        /* Get the update ratio */
        CMFD_PRECISION update_ratio = getUpdateRatio(i, e, *iter);

        /* Limit the update ratio for stability purposes. For very low flux
           regions, update ratio may be left unrestricted a few iterations*/
        if (_moc_iteration > _num_unbounded_iterations) {
          if (update_ratio > 20.0)
            update_ratio = 20.0;
          else if (update_ratio < 0.05)
            update_ratio = 0.05;
        }

        /* Save max update ratio among fsrs and groups in a cell */
        if (_convergence_data != NULL)
            if (std::abs(log(update_ratio)) >
                std::abs(log(thread_max_update_ratio)))
              thread_max_update_ratio = update_ratio;

        for (int h = _group_indices[e]; h < _group_indices[e + 1]; h++) {

          /* Update FSR flux using ratio of old and new CMFD flux */
          _FSR_fluxes[*iter*_num_moc_groups + h] *= update_ratio;

          /* Update flux moments if they were set */
          if (_linear_source) {
            _flux_moments[(*iter)*3*_num_moc_groups + h] *= update_ratio;
            _flux_moments[(*iter)*3*_num_moc_groups + _num_moc_groups + h]
                 *= update_ratio;
            _flux_moments[(*iter)*3*_num_moc_groups + 2*_num_moc_groups + h]
                 *= update_ratio;
          }

          log_printf(DEBUG, "Updating flux in FSR: %d, cell: %d, MOC group: "
            "%d, CMFD group: %d, ratio: %f", *iter ,i, h, e, update_ratio);
        }
      }
    }

    /* Save maximum update ratio among CMFD cells, for output purposes */
    if (_convergence_data != NULL) {
#pragma omp critical
      {
        if (std::abs(log(thread_max_update_ratio)) >
            std::abs(log(_convergence_data->pf)))
              _convergence_data->pf = thread_max_update_ratio;
      }
    }
  }

  /* Compare maximum update ratios between domains and keep the maximum */
#ifdef MPIx
  if (_domain_communicator != NULL && _convergence_data != NULL) {
    double max_pf = _convergence_data->pf;
    MPI_Allreduce(&max_pf, &_convergence_data->pf, 1, MPI_DOUBLE, MPI_MAX,
                  _domain_communicator->_MPI_cart);
  }
#endif
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
 *          This is the reason that Larsen's effective diffusion coefficient is
 *          useful in ensuring that the CMFD acceleration equations have a
 *          diffusion coefficient (on the flux gradient term) that is
 *          consistent, not with the physical transport problem, but with the
 *          transport problem that is being accelerated by the CMFD equations.
 *          Larsen's effective diffusion coefficient is precisely this term in
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
CMFD_PRECISION Cmfd::computeLarsensEDCFactor(CMFD_PRECISION dif_coef,
                                             CMFD_PRECISION delta) {

  /* Initialize variables */
  CMFD_PRECISION alpha, mu, expon;
  double rho = 0.0;

  /* Loop over azimuthal angles */
  for (int a=0; a < _num_azim/2; a++) {

    CMFD_PRECISION wa = _quadrature->getAzimWeight(a);

    /* Loop over polar angles */
    for (int p = 0; p < _num_polar/2; p++) {
      mu = sqrt(1.0 - pow(_quadrature->getSinTheta(a,p), 2));
      expon = exp(-delta / (3 * dif_coef * mu));
      alpha = (1 + expon) / (1 - expon) - 2 * (3 * dif_coef * mu) / delta;
      rho += 2.0 * mu * _quadrature->getPolarWeight(a,p) * wa * alpha;
    }
  }

  /* Compute the correction factor */
  CMFD_PRECISION correction = 1.0 + delta * rho / (2 * dif_coef);

  return correction;
}


/**
 * @brief Set the FSR materials array pointer.
 * @param FSR_materials pointer to FSR_materials array
 */
void Cmfd::setFSRMaterials(Material** FSR_materials) {
  _FSR_materials = FSR_materials;
}


/**
 * @brief Set the pointer to the array of FSR_volumes.
 * @param FSR_volumes array of FSR volumes
 */
void Cmfd::setFSRVolumes(FP_PRECISION* FSR_volumes) {
  _FSR_volumes = FSR_volumes;
}


/**
 * @brief Set pointer to FSR flux array.
 * @param scalar_flux pointer to FSR flux array
 */
void Cmfd::setFSRFluxes(FP_PRECISION* scalar_flux) {
  _FSR_fluxes = scalar_flux;
}


/**
 * @brief Set pointer to FSR source array.
 * @param sources pointer to FSR source array
 */
void Cmfd::setFSRSources(FP_PRECISION* sources) {
  _FSR_sources = sources;
}


/**
 * @brief Set pointer to source region flux moments array.
 * @param flux_moments pointer to source region flux moments array
 */
void Cmfd::setFluxMoments(FP_PRECISION* flux_moments) {
  _flux_moments = flux_moments;
  _linear_source = true;
}


/**
 * @brief Set successive over-relaxation relaxation factor.
 * @param SOR_factor over-relaxation factor
 */
void Cmfd::setSORRelaxationFactor(double SOR_factor) {

  if (SOR_factor <= 0.0 || SOR_factor >= 2.0)
    log_printf(ERROR, "The successive over-relaxation relaxation factor "
                      "must be > 0 and < 2. Input value: %i", SOR_factor);

  _SOR_factor = SOR_factor;
}


/**
 * @brief Set the CMFD relaxation factor applied to diffusion coefficients.
 * @param relaxation_factor CMFD relaxation factor
 */
void Cmfd::setCMFDRelaxationFactor(double relaxation_factor) {

  if (relaxation_factor <= 0.0 || relaxation_factor > 1.0)
    log_printf(ERROR, "The CMFD relaxation factor must be greater than 0 and "
                      "less than or equal to 1. Input value: %6.4f",
                      relaxation_factor);

  _relaxation_factor = relaxation_factor;
}


/**
 * @brief Forces CMFD to check neutron balance on every solve.
 */
void Cmfd::checkBalance() {
  _check_neutron_balance = true;
}


/**
 * @brief Print the CMFD prolongation factors at every iteration.
 */
void Cmfd::printProlongation() {
  _print_cmfd_prolongation_ratios = true;
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
      if (group_indices[i][j] != last_moc_group + 1)
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
  if (_materials != NULL){
    for (int i=0; i < _local_num_x * _local_num_y * _local_num_z; i++)
      delete _materials[i];
    delete [] _materials;
  }

  /* Compute and log size in memory of Material array */
  double size = (double) (_num_cmfd_groups + 4) * _num_cmfd_groups *
              _local_num_x * _local_num_y * _local_num_z *
              sizeof(FP_PRECISION) / (double) 1e6;
  log_printf(NORMAL, "CMFD material storage per domain = %6.2f MB", size);

  try {
    _materials = new Material*[_local_num_x*_local_num_y*_local_num_z];
    for (int z = 0; z < _local_num_z; z++) {
      for (int y = 0; y < _local_num_y; y++) {
        for (int x = 0; x < _local_num_x; x++) {
          int ind = z*_local_num_x*_local_num_y + y*_local_num_x + x;
          material = new Material(ind);
          material->setNumEnergyGroups(_num_cmfd_groups);
          _materials[ind] = material;
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

  float size = _num_cmfd_groups * (NUM_FACES + 2 * _balance_sigma_t) *
              _local_num_x * _local_num_y * _local_num_z *
              sizeof(CMFD_PRECISION) / float(1e6);
  log_printf(INFO_ONCE, "CMFD surface current storage per domain = %6.2f MB",
             size);

  /* Allocate memory for the CMFD Mesh surface and corner currents Vectors */
  _surface_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                 _local_num_z, _num_cmfd_groups * NUM_FACES);

  if (_balance_sigma_t) {
    /* Allocate memory for the actual starting currents on boundary CMFD cells */
    _starting_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                    _local_num_z, _num_cmfd_groups);

    /* Allocate memory for the net currents of all CMFD cells */
    _net_currents = new Vector(_cell_locks, _local_num_x, _local_num_y,
                               _local_num_z, _num_cmfd_groups);
  }
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
  for (int z = 0; z < _local_num_z; z++) {
    for (int y = 0; y < _local_num_y; y++) {
      for (int x = 0; x < _local_num_x; x++)
        _cell_fsrs.push_back(std::vector<long>());
    }
  }
}


/**
 * @brief Allocates memory for the CMFD tallies.
 * @details This method is called by the CMFD initialization routine, and
 *          allocates memory for the diffusion, reaction and volume tallies for
 *          every CMFD cells.
 */
void Cmfd::allocateTallies() {

  if (_num_x*_num_y*_num_z == 0)
    log_printf(ERROR, "Zero cells in CMFD mesh. Please set CMFD mesh before "
               "initializing CMFD tallies.");

  if (_num_cmfd_groups == 0)
    log_printf(ERROR, "Zero CMFD groups. Please set CMFD group structure "
               "before initializing CMFD tallies.");

  /* Determine tally sizes */
  int num_cells = _num_x * _num_y * _num_z;
  int local_num_cells = _local_num_x * _local_num_y * _local_num_z;
  int tally_size = local_num_cells * _num_cmfd_groups;
  _total_tally_size = 3 * tally_size;
  _tally_memory = new CMFD_PRECISION[_total_tally_size];
  CMFD_PRECISION** all_tallies[3];
  for (int t=0; t < 3; t++) {
    all_tallies[t] = new CMFD_PRECISION*[local_num_cells];
    for (int i=0; i < local_num_cells; i++) {
      int idx = i * _num_cmfd_groups + t * tally_size;
      all_tallies[t][i] = &_tally_memory[idx];
    }
  }
  log_printf(INFO_ONCE, "CMFD tally storage per domain = %6.2f MB",
             (_total_tally_size * sizeof(CMFD_PRECISION) + 3 * local_num_cells
             * sizeof(CMFD_PRECISION*)) / float(1e6));

  /* Assign tallies to allocated data */
  _diffusion_tally = all_tallies[0];
  _reaction_tally = all_tallies[1];
  _volume_tally = all_tallies[2];
  _tallies_allocated = true;
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
 * @param cell the CMFD cell ID that the local coords is in.
 * @param coords the coords being evaluated.
 * @return The surface ID.
 */
int Cmfd::findCmfdSurface(int cell, LocalCoords* coords) {
  Point* point = coords->getHighestLevel()->getPoint();

  /* If domain decomposition, compute the global CMFD cell ID */
  if (_geometry->isDomainDecomposed())
    cell = getGlobalCMFDCell(cell);
  return _lattice->getLatticeSurface(cell, point);
}


/**
 * @brief Find the CMFD cell that a LocalCoords object is in.
 * @param coords the coords being evaluated.
 * @return The CMFD cell ID.
 */
int Cmfd::findCmfdCell(LocalCoords* coords) {
  Point* point = coords->getHighestLevel()->getPoint();
  int global_cmfd_cell = _lattice->getLatticeCell(point);

  /* If domain decomposition, compute the local CMFD cell ID*/
  if (_geometry->isDomainDecomposed()) {
    int local_cmfd_cell = getLocalCMFDCell(global_cmfd_cell);
    return local_cmfd_cell;
  }
  else
    /* Without decomposition, global and local CMFD cell ID are equal.*/
    return global_cmfd_cell;
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
 * @param cmfd_cell the CMFD cell ID.
 * @param fsr_id the FSR ID.
 */
void Cmfd::addFSRToCell(int cmfd_cell, long fsr_id) {
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
 * @param num_fsrs the number of FSRs
 */
void Cmfd::setNumFSRs(long num_fsrs) {
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
    std::map<int, CMFD_PRECISION>::iterator it;
    int cell, surface;


#pragma omp for
    for (int i=0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

      int global_id = getGlobalCMFDCell(i);

      for (int v = NUM_FACES + NUM_EDGES; v < NUM_SURFACES; v++) {

        /* Check if this surface is contained in the map */
        int ind = i * NUM_SURFACES * ncg + v * ncg;
        it = _edge_corner_currents.find(ind);
        if (it == _edge_corner_currents.end())
          continue;

        getVertexSplitSurfaces(global_id, v, &surfaces);

        //TODO Optimize since dont want to look for cells at every group
        for (int g=0; g < ncg; g++) {
          /* Divide vertex current by 3 since we will split to 3 surfaces,
           * which propagate through 3 edges */
          current = _edge_corner_currents[ind+g] / 3;

          /* Increment current for faces and edges adjacent to this vertex */
          for (iter = surfaces.begin(); iter != surfaces.end(); ++iter) {
            cell = (*iter) / ns;
            surface = (*iter) % ns;

            /* Look for the CMFD cell on-domain */
            int local_cell = getLocalCMFDCell(cell);

            if (local_cell != -1) {
              /* Add face contributions */
              if (surface < NUM_FACES) {
                _surface_currents->incrementValue(local_cell,
                                                  surface * ncg + g, current);
              }
              /* Add edges contributions */
              else {

                /* Map is accessed directly at new index, if the key doesn't
                exist, the default constructor is called, initializing at 0. */
                int new_ind = (local_cell * NUM_SURFACES + surface) * ncg + g;

                /* Add the contribution */
                omp_set_lock(&_edge_corner_lock);
                _edge_corner_currents[new_ind] += current;
                omp_unset_lock(&_edge_corner_lock);
              }
            }

            /* Look for the CMFD cell off-domain */
            else {

              /* Look for the boundary containing the cell */
              for (int s=0; s < NUM_FACES; s++) {

                std::map<int, int>::iterator it =
                  _boundary_index_map.at(s).find(cell);

                if (it != _boundary_index_map.at(s).end()) {

                  int idx = it->second;

                  /* Add the current to the off-domain split currents cell */
                  // Comment out this line to check for MOC neutron balance
                  omp_set_lock(&_edge_corner_lock);
                  _off_domain_split_currents[s][idx][surface * ncg + g] +=
                    current;
                  omp_unset_lock(&_edge_corner_lock);
                  break;
                }
              }
            }
          }
          _edge_corner_currents[ind+g] = 0.0;
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
    std::map<int, CMFD_PRECISION>::iterator it;
    int cell, surface;

#pragma omp for
    for (int i=0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

      int global_id = getGlobalCMFDCell(i);

      for (int e = NUM_FACES; e < NUM_FACES + NUM_EDGES; e++) {

        /* Check if this surface is contained in the map */
        int ind = i * NUM_SURFACES * ncg + e * ncg;
        it = _edge_corner_currents.find(ind);
        if (it == _edge_corner_currents.end())
          continue;

        getEdgeSplitSurfaces(global_id, e, &surfaces);

        //TODO Optimize since dont want to look for cells at every group
        for (int g=0; g < ncg; g++) {
          /* Divide edge current by 2 since we will split to 2 surfaces,
           * which propagate through 2 surfaces */
          current = _edge_corner_currents[ind+g] / 2;

          /* Increment current for faces and edges adjacent to this vertex */
          for (iter = surfaces.begin(); iter != surfaces.end(); ++iter) {
            cell = (*iter) / ns;
            surface = (*iter) % ns;

            /* Look for the CMFD cell on-domain */
            int local_cell = getLocalCMFDCell(cell);
            if (local_cell != -1) {
              _surface_currents->incrementValue(local_cell, surface * ncg + g,
                                                current);
            }

            /* Look for the CMFD cell off-domain */
            else {

              /* Look for the boundary containing the cell */
              for (int s=0; s < NUM_FACES; s++) {

                std::map<int, int>::iterator it =
                  _boundary_index_map.at(s).find(cell);

                if (it != _boundary_index_map.at(s).end()) {

                  int idx = it->second;

                  /* Add the current to the off-domain split currents cell */
                  // Comment out this line to check for MOC neutron balance
                  omp_set_lock(&_edge_corner_lock);
                  _off_domain_split_currents[s][idx][surface * ncg + g] +=
                    current;
                  omp_unset_lock(&_edge_corner_lock);
                  break;
                }
              }
            }
          }
          _edge_corner_currents[ind+g] = 0.0;
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

  int cell_indexes[3] = {x,y,z};
  int cell_limits[3] = {_num_x, _num_y, _num_z};

  int direction[3];
  convertSurfaceToDirection(vertex, direction);

  /* Get the partial surfaces composing the edge split */
  int remainder_surfaces[3];
  int partial_surfaces[3];
  for (int i=0; i < 3; i++) {
    int remainder_direction[3];
    int partial_direction[3];
    for (int j=0; j < 3; j++) {
      if (i == j) {
        remainder_direction[j] = 0;
        partial_direction[j] = direction[j];
      }
      else {
        remainder_direction[j] = direction[j];
        partial_direction[j] = 0;
      }
    }
    remainder_surfaces[i] = convertDirectionToSurface(remainder_direction);
    partial_surfaces[i] = convertDirectionToSurface(partial_direction);
  }

  /* Treat all partial surfaces */
  for (int i=0; i < 3; i++) {

    int remainder_surface = remainder_surfaces[i];
    int partial_surface = partial_surfaces[i];

    surfaces->push_back(cell * ns + partial_surface);

    /* Tally current on neighboring cell or appropriate boundary */
    int cell_next = getCellNext(cell, partial_surface);
    if ((cell_indexes[i] == 0 && direction[i] == -1) ||
        (cell_indexes[i] == cell_limits[i] - 1 && direction[i] == +1)) {

      if (_boundaries[partial_surface] == REFLECTIVE)
        surfaces->push_back(cell * ns + remainder_surface);
        //NOTE Comment out this line to check MOC neutron balance

      else if (_boundaries[partial_surface] == PERIODIC)
        surfaces->push_back(cell_next * ns + remainder_surface);

    }
    else
      surfaces->push_back(cell_next * ns + remainder_surface);

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

  int cell_indexes[3] = {x,y,z};
  int cell_limits[3] = {_num_x, _num_y, _num_z};

  int direction[3];
  convertSurfaceToDirection(edge, direction);

  /* Get the partial surfaces composing the edge split */
  int partial_surfaces[2];
  int opposite_surfaces[2];
  int ind = 0;
  for (int i=0; i < 3; i++) {
    if (direction[i] != 0) {
      int partial_direction[3] = {0,0,0};
      partial_direction[i] = direction[i];
      partial_surfaces[ind] = convertDirectionToSurface(partial_direction);
      ind++;
    }
  }

  /* Treat all partial surfaces */
  ind = 0;
  for (int i=0; i < 3; i++) {
    if (direction[i] != 0) {

      int partial_surface = partial_surfaces[ind];
      int other_surface = partial_surfaces[1-ind];

      surfaces->push_back(cell * ns + partial_surface);

      /* Tally current on neighboring cell or appropriate boundary */
      int cell_next = getCellNext(cell, partial_surface);
      if ((cell_indexes[i] == 0 && direction[i] == -1) ||
          (cell_indexes[i] == cell_limits[i] - 1 && direction[i] == +1)) {

        if (_boundaries[partial_surface] == REFLECTIVE)
          surfaces->push_back(cell * ns + other_surface);
          //NOTE Comment out this line to check MOC neutron balance

        else if (_boundaries[partial_surface] == PERIODIC)
          surfaces->push_back(cell_next * ns + other_surface);
      }
      else
        surfaces->push_back(cell_next * ns + other_surface);

      ind++;
    }
  }
}


/**
 * @brief Get the ID of the Mesh cell next to given Mesh cell.
 * @param cell index of the current CMFD cell
 * @param surface_id id of the surface between the current cell and the next
 * @param global work at the global (all domains together) level
 * @param neighbor give cell in neighboring domain
 * @return neighboring CMFD cell ID
 */
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
    x_global = x + _accumulate_lmx[_domain_communicator->_domain_idx_x];
    y_global = y + _accumulate_lmy[_domain_communicator->_domain_idx_y];
    z_global = z + _accumulate_lmz[_domain_communicator->_domain_idx_z];
    nx = _local_num_x;
    ny = _local_num_y;
    nz = _local_num_z;
  }

  /* Find the cell on the other side of the surface */
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
 * @param side the CMFD mesh surface ID.
 * @return the boundaryType for the surface.
 */
int Cmfd::getBoundary(int side) {
  return _boundaries[side];
}


/**
 * @brief Return the CMFD cell ID that a FSR lies in.
 * @details Note that a CMFD cell is not an actual Cell object; rather, a CMFD
 *          cell is just a way of describing each of the rectangular regions
 *          that make up a CMFD lattice. CMFD cells are numbered with 0 in the
 *          lower left corner and monotonically increasing from left to right,
 *          from bottom to top. For example, the indices for a 4 x 4 lattice
 *          are:
 *                  12  13  14  15
 *                  8    9  10  11
 *                  4    5   6   7
 *                  0    1   2   3
 * @param fsr_id the FSR ID.
 * @return The CMFD cell ID. Return -1 if cell is not found.
 */
int Cmfd::convertFSRIdToCmfdCell(long fsr_id) {

  std::vector<long>::iterator iter;
  for (int cell_id=0; cell_id < _local_num_x*_local_num_y*_local_num_z;
      cell_id++) {
    for (iter = _cell_fsrs.at(cell_id).begin();
         iter != _cell_fsrs.at(cell_id).end(); ++iter) {
      if (*iter  == fsr_id)
        return cell_id;
    }
  }

  return -1;
}


/**
 * @brief Return the CMFD cell ID that a FSR lies in.
 * @param global_fsr_id The global FSR ID.
 * @return The CMFD cell ID.
 */
int Cmfd::convertGlobalFSRIdToCmfdCell(long global_fsr_id) {

  /* Determine the domain and local FSR ID */
  int cmfd_cell = -1;
  if (!_geometry->isDomainDecomposed()) {
    cmfd_cell = convertFSRIdToCmfdCell(global_fsr_id);
  }
#ifdef MPIx
  else {

    long fsr_id;
    int domain;
    _geometry->getLocalFSRId(global_fsr_id, fsr_id, domain);

    /* Get the FSR centroid in the correct domain */
    int rank;
    MPI_Comm comm = _geometry->getMPICart();
    MPI_Comm_rank(comm, &rank);
    int temp_cmfd_cell = 0;
    if (rank == domain)
      temp_cmfd_cell = convertFSRIdToCmfdCell(fsr_id);

    /* Broadcast the temp_cmfd_cell */
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
std::vector< std::vector<long> >* Cmfd::getCellFSRs() {
  return &_cell_fsrs;
}


/**
 * @brief Set the vector of vectors that contains the FSRs that lie in each
 *        cell.
 * @param cell_fsrs vector of vectors containing FSR IDs in each cell.
 */
void Cmfd::setCellFSRs(std::vector< std::vector<long> >* cell_fsrs) {

  if (!_cell_fsrs.empty()) {
    std::vector< std::vector<long> >::iterator iter;
    for (iter = _cell_fsrs.begin(); iter != _cell_fsrs.end(); ++iter)
      iter->clear();
    _cell_fsrs.clear();
  }

  _cell_fsrs = *cell_fsrs;
}


/**
 * @brief Set flag indicating whether to update the MOC flux.
 * @param flux_update_on Flag saying whether to update MOC flux.
 */
void Cmfd::setFluxUpdateOn(bool flux_update_on) {
  _flux_update_on = flux_update_on;
}


/**
 * @brief Sets the a ConvergenceData object to record diagnostics.
 * @details The ConvergenceData object records the number of fission source
 *          and flux iterations for the CMFD solver as well as the maximum
 *          magnitude prolongation factor
 * @param convergence_data The convergence data object
 */
void Cmfd::setConvergenceData(ConvergenceData* convergence_data) {
  _convergence_data = convergence_data;
}


/**
 * @brief Set the flag indicating whether to use quadratic axial interpolation
 *        for update ratios.
 * @param interpolate flag meaning No interpolation(0), FSR axially averaged
          value(1) or centroid z-coordinate evaluated value(2)
 */
void Cmfd::useAxialInterpolation(int interpolate) {

  if (interpolate<0 || interpolate>2)
    log_printf(ERROR, "interpolate can only has value 0, 1, or 2, respectively"
               " meaning No interpolation, FSR axially averaged value or"
               " centroid z-coordinate evaluated value");
  if (interpolate==1 || interpolate==2)
    log_printf(WARNING_ONCE, "Axial interpolation CMFD prolongation may only"
               " be effective when all the FSRs are axially homogeneous");
  _use_axial_interpolation = interpolate;
}


/**
 * @brief Turns on the flux limiting condition.
 * @details If the CMFD correction diffusion coefficient is larger than the
 *          diffusion coefficient, recompute the diffusion coefficient as the
 *          ratio of current to twice the flux, and re-compute a correction
 *          diffusion coefficient.
 * @param flux_limiting whether to turn on the flux limiting condition
 */
void Cmfd::useFluxLimiting(bool flux_limiting) {
  _flux_limiting = flux_limiting;
}


/**
 * @brief Modifies the diagonal element to be consistent with the MOC solve
 * @details This function re-computes a new total cross-section x volume that
 *          maintains consistency with the MOC solution. Generally, this will
 *          not change the diagonal element at all since CMFD should be
 *          consistent with MOC. However, if negative fluxes are corrected to
 *          zero after the MOC transport sweep, there will be an inconsistency.
 *          This function modifies sigma-t so that there is consistency with
 *          the altered solution.
 * @param cmfd_cell The cmfd cell of the element to adjust
 * @param group The cmfd group of the element to adjust
 */
void Cmfd::enforceBalanceOnDiagonal(int cmfd_cell, int group) {

  /* Initialize tallies */
  Material* material = _materials[cmfd_cell];
  double cmfd_volume = _volumes->getValue(cmfd_cell, 0);

  /* Loop over FSRs in CMFD cell to tally the total neutron source */
  double moc_source = 0.0;
  for (int j = 0; j < _cell_fsrs.at(cmfd_cell).size(); j++) {

    long fsr_id = _cell_fsrs.at(cmfd_cell).at(j);
    FP_PRECISION volume = _FSR_volumes[fsr_id];

    /* Loop over MOC energy groups within this CMFD coarse group */
    for (int h = _group_indices[group]; h < _group_indices[group+1]; h++)
      moc_source += 4 * M_PI * volume *
        _FSR_sources[fsr_id * _num_moc_groups + h];
  }

  if (fabs(moc_source) < FLUX_EPSILON)
    moc_source = FLUX_EPSILON;

  /* Compute updated value */
  double flux = _old_flux->getValue(cmfd_cell, group);
  CMFD_PRECISION net_current = _net_currents->getValue(cmfd_cell, group);
  CMFD_PRECISION updated_value = (moc_source - net_current) / flux;

  if (updated_value < 0.0)
    log_printf(ERROR, "Negative Total XS of %6.4f computed in CMFD rebalance",
               updated_value);

  /* Update the diagonal element */
  _A->setValue(cmfd_cell, group, cmfd_cell, group, updated_value);
}


/**
 * @brief Rebalances the total cross section to be consistent with the MOC
 *        solution on every sweep.
 * @param balance_sigma_t Wheter to compute the rebalanced total cross-section
 */
void Cmfd::rebalanceSigmaT(bool balance_sigma_t) {
  _balance_sigma_t = balance_sigma_t;
}


/**
 * @brief Returns a flag indicating whether the sigma-t rebalance is on.
 * @return A flag indicating whether the rebalance is on
 */
bool Cmfd::isSigmaTRebalanceOn() {
  return _balance_sigma_t;
}


/**
 * @brief Get flag indicating whether to update the MOC flux.
 * @return Flag saying whether to update MOC flux
 */
bool Cmfd::isFluxUpdateOn() {
 return _flux_update_on;
}


/**
 * @brief Set flag indicating whether to use FSR centroids to update
 *        the MOC flux.
 * @param centroid_update_on Flag saying whether to use centroids to
 *        update MOC flux
 */
void Cmfd::setCentroidUpdateOn(bool centroid_update_on) {
  _centroid_update_on = centroid_update_on;
}


/**
 * @brief Get flag indicating whether to use FSR centroids to update
 *        the MOC flux.
 * @return Flag saying whether to use centroids to update MOC flux
 */
bool Cmfd::isCentroidUpdateOn() {
  return _centroid_update_on;
}


/**
 * @brief Sets the threshold for CMFD source convergence (>0).
 * @param source_thresh the threshold for source convergence
 */
void Cmfd::setSourceConvergenceThreshold(double source_thresh) {

  if (source_thresh <= 0.0)
    log_printf(ERROR, "Unable to set the CMFD source convergence threshold to"
              " %f since the threshold must be positive.", source_thresh);

  _source_convergence_threshold = source_thresh;
}


/**
 * @brief Sets the Quadrature object in use by the MOC Solver.
 * @param quadrature a Quadrature object pointer from the Solver
 */
void Cmfd::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
  _num_polar = quadrature->getNumPolarAngles();
  _num_azim = quadrature->getNumAzimAngles();
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
  std::vector< std::pair<int, double> >::iterator stencil_iter;
  std::vector<long>::iterator fsr_iter;
  Point* centroid;
  long fsr_id;

  if (_centroid_update_on){
    /* Number of cells in stencil */
    int num_cells_in_stencil = 9;

    /* Loop over mesh cells */
    for (int i = 0; i < _local_num_x*_local_num_y*_local_num_z; i++) {

      int global_ind = getGlobalCMFDCell(i);

      /* Loop over FRSs in mesh cell */
      for (fsr_iter = _cell_fsrs.at(i).begin();
          fsr_iter != _cell_fsrs.at(i).end(); ++fsr_iter) {

        fsr_id = *fsr_iter;

        /* Get centroid */
        centroid = _geometry->getFSRCentroid(fsr_id);

        /* Create new stencil */
        _k_nearest_stencils[fsr_id] =
          std::vector< std::pair<int, double> >();

        /* Get distance to all cells that touch current cell */
        for (int j=0; j < num_cells_in_stencil; j++)
          _k_nearest_stencils[fsr_id]
            .push_back(std::make_pair<int, double>
                      (int(j), getDistanceToCentroid(centroid, global_ind, i,
                                                     j)));

        /* Sort the distances */
        std::sort(_k_nearest_stencils[fsr_id].begin(),
                  _k_nearest_stencils[fsr_id].end(), stencilCompare);

        /* Remove ghost cells that are outside the geometry boundaries */
        stencil_iter = _k_nearest_stencils[fsr_id].begin();
        while (stencil_iter != _k_nearest_stencils[fsr_id].end()) {
          if (stencil_iter->second > FLT_INFINITY)
            stencil_iter = _k_nearest_stencils[fsr_id].erase(stencil_iter);
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
    double total_distance;
    for (long i=0; i < _num_FSRs; i++) {
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



  /* Compute axial quadratic interpolation values if requested */
  if (_use_axial_interpolation && _local_num_z >= 3) {

    /* Initialize axial quadratic interpolant values */
    _axial_interpolants.resize(_num_FSRs);
    for (long r=0; r < _num_FSRs; r++) {
      _axial_interpolants.at(r) = new double[3]();
    }
    log_printf(INFO_ONCE, "CMFD axial interpolation storage per domain = %6.2f"
               " MB", _num_FSRs * 3 * sizeof(double) / float(1e6));

    /* Loop over mesh cells */
    for (int i = 0; i < _local_num_x*_local_num_y*_local_num_z; i++) {

      /* Starting z number of CMFD mesh in this domain */
      int z_start = 0;
      if (_domain_communicator != NULL)
        z_start = _accumulate_lmz[_domain_communicator->_domain_idx_z];

      /* Calculate the CMFD cell z-coordinate */
      int z_ind = i / (_local_num_x * _local_num_y);

      /* The heights of neighboring three CMFD meshes for quadratic fit */
      double h0, h1, h2;

      /* The z coordinate of the mesh center of the middle-CMFD cell  */
      double z_cmfd;

      if (z_ind == 0) {
        h0 = _cell_widths_z[z_start + z_ind];
        h1 = _cell_widths_z[z_start + z_ind + 1];
        h2 = _cell_widths_z[z_start + z_ind + 2];
        z_cmfd = _accumulate_z[z_start + z_ind+1] + h1/2. + _lattice->getMinZ();
      }
      else if (z_ind == _local_num_z - 1) {
        h0 = _cell_widths_z[z_start + z_ind - 2];
        h1 = _cell_widths_z[z_start + z_ind - 1];
        h2 = _cell_widths_z[z_start + z_ind];
        z_cmfd = _accumulate_z[z_start + z_ind-1] + h1/2. + _lattice->getMinZ();
      }
      else {
        h0 = _cell_widths_z[z_start + z_ind - 1];
        h1 = _cell_widths_z[z_start + z_ind];
        h2 = _cell_widths_z[z_start + z_ind + 1];
        z_cmfd = _accumulate_z[z_start + z_ind] + h1/2. + _lattice->getMinZ();
      }

      /* Start and end relative z-coordinate of an FSR */
      double zs, ze;

      /* Loop over FRSs in mesh cell */
      int num_fissionable_FSRs = 0;
      for (fsr_iter = _cell_fsrs.at(i).begin();
           fsr_iter != _cell_fsrs.at(i).end(); ++fsr_iter) {

        /* Get centroid and calculate relative z-coordinate */
        fsr_id = *fsr_iter;
        Point* centroid = _geometry->getFSRCentroid(fsr_id);
        Point* feature_point = _geometry->getFSRPoint(fsr_id);

        double zc = (centroid->getZ() - z_cmfd) / h1;
        zs = (feature_point->getZ() - z_cmfd) / h1;
        ze = 2*zc - zs;

        /* Calculate components for quadratic interpolation of the FSR axially
           averaged value */
        if (_use_axial_interpolation ==1) {
          _axial_interpolants.at(fsr_id)[0] = (h1*(h1+h1*zc*4.0+h2*zc*8.0-
            h1*(zc*zc)*1.6E1-h1*(zs*zs)*4.0+h1*zc*zs*8.0)*(-1.0/4.0))
                                              /((h0+h1)*(h0+h1+h2));


          _axial_interpolants.at(fsr_id)[1] = (h1+h2-h1*zc*2.0)/(h1+h2)+(h1*(h1+
          h1*zc*4.0+h2*zc*8.0-h1*(zc*zc)*1.6E1-h1*(zs*zs)*4.0+h1*zc*zs*8.0)*
          (1.0/4.0))/(h2*(h0+h1))-((h1*h1)*(h1+h1*zc*4.0+h2*zc*8.0-h1*(zc*zc)*
          1.6E1-h1*(zs*zs)*4.0+h1*zc*zs*8.0)*(1.0/4.0))/(h2*(h1+h2)*(h0+h1+h2));


          _axial_interpolants.at(fsr_id)[2] = (h1*(-h1+h0*zc*8.0+h1*zc*4.0+
            h1*(zc*zc)*1.6E1+h1*(zs*zs)*4.0-h1*zc*zs*8.0)*(1.0/4.0))
                                              /((h1+h2)*(h0+h1+h2));


          log_printf(DEBUG, "CMFD-ID: %d, FSR-ID: %ld, c0= %10.6f, c1= %10.6f,"
                     " c2= %10.6f", i, fsr_id, _axial_interpolants.at(fsr_id)[0],
                     _axial_interpolants.at(fsr_id)[1],
                     _axial_interpolants.at(fsr_id)[2]);
        }

      /* Calculate components for quadratic interpolation of the centroid
         z-coordinate evaluated value. */
        else if (_use_axial_interpolation == 2) {
          _axial_interpolants.at(fsr_id)[0] = -(h1*(h1+h1*zc*4.0+h2*zc*8.0-h1*
          (zc*zc)*1.2E1))/((h0*4.0+h1*4.0)*(h0+h1+h2));


          _axial_interpolants.at(fsr_id)[1] = (-zc*(h0*(h1*h1)*1.2E1+(h0*h0)*h1*
          8.0-h1*(h2*h2)*8.0-(h1*h1)*h2*1.2E1)+h0*(h1*h1)*9.0+(h0*h0)*h1*4.0+h0*
          (h2*h2)*4.0+(h0*h0)*h2*4.0+h1*(h2*h2)*4.0+(h1*h1)*h2*9.0+(h1*h1*h1)*6.0
          -(zc*zc)*(h0*(h1*h1)*1.2E1+(h1*h1)*h2*1.2E1+(h1*h1*h1)*2.4E1)+h0*h1*h2
          *1.2E1)/((h1+h2)*(h0*4.0+h1*4.0)*(h0+h1+h2));


         _axial_interpolants.at(fsr_id)[2] = (h1*(-h1+h0*zc*8.0+h1*zc*4.0+h1*
         (zc*zc)*1.2E1))/((h1*4.0+h2*4.0)*(h0+h1+h2));


         log_printf(DEBUG, "CMFD-ID: %d, FSR-ID: %ld, c0= %10.6f, c1= %10.6f,"
                    " c2= %10.6f", i, fsr_id, _axial_interpolants.at(fsr_id)[0],
                    _axial_interpolants.at(fsr_id)[1],
                    _axial_interpolants.at(fsr_id)[2]);
        }
        /* Calculate components for quadratic interpolation of the centroid
           z-coordinate evaluted value. For uniform CMFD */
        /*_axial_interpolants.at(fsr_id)[0] = zc * zc/2.0 - zc/2.0 - 1.0/24.0;
        _axial_interpolants.at(fsr_id)[1] = -zc * zc + 26.0/24.0;
        _axial_interpolants.at(fsr_id)[2] = zc * zc/2.0 + zc/2.0 - 1.0/24.0;*/

        /* Set zero axial prolongation for cells with no fissionable material */
        if (_FSR_materials[fsr_id]->isFissionable())
          num_fissionable_FSRs++;
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
  int x = (cell_id % (_local_num_x * _local_num_y)) % _local_num_x;
  int y = (cell_id % (_local_num_x * _local_num_y)) / _local_num_x;

  if (stencil_id == 0) {
    if (x != 0 && y != 0)
      cell_next_id = cell_id - _local_num_x - 1;
  }
  else if (stencil_id == 1) {
    if (y != 0)
      cell_next_id = cell_id - _local_num_x;
    else if (_boundaries[SURFACE_Y_MIN] == PERIODIC)
      cell_next_id = cell_id + _local_num_x * (_local_num_y - 1);
  }
  else if (stencil_id == 2) {
    if (x != _local_num_x - 1 && y != 0)
      cell_next_id = cell_id - _local_num_x + 1;
  }
  else if (stencil_id == 3) {
    if (x != 0)
      cell_next_id = cell_id - 1;
    else if (_boundaries[SURFACE_X_MIN] == PERIODIC)
      cell_next_id = cell_id + (_local_num_x - 1);
  }
  else if (stencil_id == 4) {
    cell_next_id = cell_id;
  }
  else if (stencil_id == 5) {
    if (x != _local_num_x - 1)
      cell_next_id = cell_id + 1;
    else if (_boundaries[SURFACE_X_MAX] == PERIODIC)
      cell_next_id = cell_id - (_local_num_x - 1);
  }
  else if (stencil_id == 6) {
    if (x != 0 && y != _local_num_y - 1)
      cell_next_id = cell_id + _local_num_x - 1;
  }
  else if (stencil_id == 7) {
    if (y != _local_num_y - 1)
      cell_next_id = cell_id + _local_num_x;
    else if (_boundaries[SURFACE_Y_MAX] == PERIODIC)
      cell_next_id = cell_id - _local_num_x * (_local_num_y - 1);
  }
  else if (stencil_id == 8) {
    if (x != _local_num_x - 1 && y != _local_num_y - 1)
      cell_next_id = cell_id + _local_num_x + 1;
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
CMFD_PRECISION Cmfd::getUpdateRatio(int cell_id, int group, long fsr) {

  CMFD_PRECISION ratio = 0.0;
  std::vector< std::pair<int, double> >::iterator iter;
  int cell_next_id;

  if (_centroid_update_on) {

    /* Compute the ratio for all the surrounding cells */
    for (iter = _k_nearest_stencils[fsr].begin();
         iter != _k_nearest_stencils[fsr].end(); ++iter) {

      if (iter->first != 4) {
        cell_next_id = getCellByStencil(cell_id, iter->first);//cell_id is Local itself here.

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
 * @brief Retrieves the ratio of pre- and post- CMFD solve fluxes
 * @details The CMFD flux ratio is returned for the given FSR. A quadratic
 *          axial interpolant is used to estimate the value at the FSR.
 * @param cell_id The CMFD cell ID containing the FSR.
 * @param group The CMFD energy group being updated.
 * @param fsr The fsr being updated.
 * @return the ratio of CMFD fluxes
 */
CMFD_PRECISION Cmfd::getFluxRatio(int cell_id, int group, long fsr) {

  double ratio = 1.0;

  if (_use_axial_interpolation && _local_num_z >= 3) {

    /* Get pre-computed interpolation ratios */
    double* interpolants = _axial_interpolants.at(fsr);

    int z_ind = cell_id / (_local_num_x * _local_num_y);
    int cell_mid = cell_id;

    /* Shift up or down one cell if at top/bottom, interpolants corrects
       for the shift in cells */
    if (z_ind == 0)
      cell_mid += _local_num_x * _local_num_y;
    else if (z_ind == _local_num_z - 1)
      cell_mid -= _local_num_x * _local_num_y;

    /* Get cell index above and below current CMFD cell */
    int cell_prev = cell_mid - _local_num_x * _local_num_y;
    int cell_next = cell_mid + _local_num_x * _local_num_y;

    /* Get new and old fluxes in bottom/mid/top cells */
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

    if (fabs(old_flux) > FLUX_EPSILON)
      ratio = new_flux / old_flux;

    /* Fallback: using the cell average flux ratio */
    if (ratio < 0) {
      if (fabs(_old_flux->getValue(cell_id, group)) > FLUX_EPSILON)
        ratio = _new_flux->getValue(cell_id, group) /
                _old_flux->getValue(cell_id, group);
      else
        ratio = 0.0;
    }

    return ratio;
  }
  else {
    if (fabs(_old_flux->getValue(cell_id, group)) > FLUX_EPSILON)
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
 *          maximum allowable double value is returned.
 * @param centroid The numerical centroid an FSR in the cell.
 * @param cell_id The global CMFD cell containing the FSR.
 * @param local_cell_id The local CMFD id of the cell containing the fsr
 * @param stencil_index The index of the cell in the stencil that we want to
 *        get the distance from.
 * @return the distance from the CMFD cell centroid to the FSR centroid.
 */
double Cmfd::getDistanceToCentroid(Point* centroid, int cell_id,
                                   int local_cell_id, int stencil_index) {

  int x = (cell_id % (_num_x * _num_y)) % _num_x;
  int y = (cell_id % (_num_x * _num_y)) / _num_x;
  int xl = (local_cell_id % (_local_num_x * _local_num_y)) % _local_num_x;
  int yl = (local_cell_id % (_local_num_x * _local_num_y)) / _local_num_x;
  double dist_x, dist_y;
  bool found = false;
  double centroid_x = centroid->getX();
  double centroid_y = centroid->getY();

  /* The center of geometry is not always at (0,0,0), then relative coordinates
     should be used. */
  double dx = centroid_x - _lattice->getMinX();
  double dy = centroid_y - _lattice->getMinY();

  /* LOWER LEFT CORNER */
  if (xl > 0 && yl > 0 && stencil_index == 0) {
    dist_x = pow(dx - (_accumulate_x[x-1]+_cell_widths_x[x-1]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y-1]+_cell_widths_y[y-1]/2), 2.0);
    found = true;
  }

  /* BOTTOM SIDE */
  else if (yl > 0 && stencil_index == 1) {
    dist_x = pow(dx - (_accumulate_x[x  ]+_cell_widths_x[x  ]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y-1]+_cell_widths_y[y-1]/2), 2.0);
    found = true;
  }

  /* LOWER RIGHT CORNER */
  else if (xl < _local_num_x - 1 && yl > 0 && stencil_index == 2) {
    dist_x = pow(dx - (_accumulate_x[x+1]+_cell_widths_x[x+1]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y-1]+_cell_widths_y[y-1]/2), 2.0);
    found = true;
  }

  /* LEFT SIDE */
  else if (xl > 0 && stencil_index == 3) {
    dist_x = pow(dx - (_accumulate_x[x-1]+_cell_widths_x[x-1]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y  ]+_cell_widths_y[y  ]/2), 2.0);
    found = true;
  }

  /* CURRENT */
  else if (stencil_index == 4) {
    dist_x = pow(dx - (_accumulate_x[x  ]+_cell_widths_x[x  ]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y  ]+_cell_widths_y[y  ]/2), 2.0);
    found = true;
  }

  /* RIGHT SIDE */
  else if (xl < _local_num_x - 1 && stencil_index == 5) {
    dist_x = pow(dx - (_accumulate_x[x+1]+_cell_widths_x[x+1]/2 ), 2.0);
    dist_y = pow(dy - (_accumulate_y[y  ]+_cell_widths_y[y  ]/2 ), 2.0);
    found = true;
  }

  /* UPPER LEFT CORNER */
  else if (xl > 0 && yl < _local_num_y - 1 && stencil_index == 6) {
    dist_x = pow(dx - (_accumulate_x[x-1]+_cell_widths_x[x-1]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y+1]+_cell_widths_y[y+1]/2), 2.0);
    found = true;
  }

  /* TOP SIDE */
  else if (yl < _local_num_y - 1 && stencil_index == 7) {
    dist_x = pow(dx - (_accumulate_x[x  ]+_cell_widths_x[x  ]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y+1]+_cell_widths_y[y+1]/2), 2.0);
    found = true;
  }

  /* UPPER RIGHT CORNER */
  else if (xl < _local_num_x - 1 && yl < _local_num_y - 1 && stencil_index == 8) {
    dist_x = pow(dx - (_accumulate_x[x+1]+_cell_widths_x[x+1]/2), 2.0);
    dist_y = pow(dy - (_accumulate_y[y+1]+_cell_widths_y[y+1]/2), 2.0);
    found = true;
  }

  if (found)
    return pow(dist_x + dist_y, 0.5);
  else
    return std::numeric_limits<double>::max();
}


/**
 * @brief Set a pointer to the Geometry.
 * @param geometry A pointer to a Geometry object.
 */
void Cmfd::setGeometry(Geometry* geometry) {
  _geometry = geometry;
}


/**
 * @brief Set the number of iterations where the CMFD update ratios are not
 *        bounded.
 * @param unbounded number of iterations without bounds on CMFD update ratios
 */
void Cmfd::setNumUnboundedIterations(int unbounded_iterations) {
  _num_unbounded_iterations = unbounded_iterations;
}


/**
 * @brief Set a number of k-nearest neighbor cells to use in updating
 *         the FSR flux.
 * @param k_nearest The number of nearest neighbor CMFD cells.
 */
void Cmfd::setKNearest(int k_nearest) {

  if (k_nearest < 1 || k_nearest > 9)
    log_printf(ERROR, "Unable to set CMFD k-nearest to %i. k-nearest "
               "must be between 1 and 9.", k_nearest);
  else
    _k_nearest = k_nearest;

  /* Enables use of centroids for K-nearest */
  _centroid_update_on = true;
}


/**
 * @brief Zero the surface currents for each mesh cell and energy group.
 */
void Cmfd::zeroCurrents() {

  _surface_currents->clear();

  if (_balance_sigma_t) {
    _starting_currents->clear();
    _net_currents->clear();
  }

  /* Clear boundary currents */
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
#pragma omp parallel for
    for (int s=0; s < NUM_FACES; s++) {

      /* Loop over all CMFD cells on the current surface */
      std::map<int, int>::iterator it;
      for (it=_boundary_index_map.at(s).begin();
          it != _boundary_index_map.at(s).end(); ++it) {

        int idx = it->second;

        /* Loop over CMFD coarse energy groups */
        for (int e = 0; e < _num_cmfd_groups; e++) {

          /* Loop over cell faces */
          for (int f=0; f < NUM_FACES; f++)
            _boundary_surface_currents[s][idx][f*_num_cmfd_groups+e] = 0.0;

          /* Loop over all cell faces and edges */
          for (int fe=0; fe < NUM_FACES + NUM_EDGES; fe++) {
            _off_domain_split_currents[s][idx][fe*_num_cmfd_groups+e] = 0.0;
            _received_split_currents[s][idx][fe*_num_cmfd_groups+e] = 0.0;
          }
        }
      }
    }
  }
#endif
}


/**
 * @brief Initialize the Matrix and Vector objects, k-nearest stencils, the
 *        CMFD cell currents and MOC materials.
 */
void Cmfd::initialize() {

  /* Delete old Matrix and Vector objects if they exist */
  if (_A != NULL)
    delete _A;
  if (_M != NULL)
    delete _M;
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
    _cell_locks = new omp_lock_t[num_cells];

    /* Loop over all cells to initialize OpenMP locks */
#pragma omp parallel for schedule(guided)
    for (int r=0; r < num_cells; r++)
      omp_init_lock(&_cell_locks[r]);
    omp_init_lock(&_edge_corner_lock);

    /* Compute and log size in memory of CMFD matrices */
    int num_rows = _num_cmfd_groups * _local_num_x * _local_num_y *
                   _local_num_z * 2;
    // A matrix is duplicated as list of list and CSR form (allocated before
    // each transport iteration's CMFD solve)

    int num_non_zero_coeffs = 6 + 20 * _num_cmfd_groups / 70;
    // 20 is estimated number of non-zero scatters per group for 70g structure
    double size = (double) (num_rows) * num_non_zero_coeffs *
                  sizeof(CMFD_PRECISION) / (double) 1e6;
    log_printf(NORMAL, "CMFD A matrix est. storage per domain = %6.2f MB",
               size);

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
    _new_flux = new Vector(_cell_locks, _local_num_x, _local_num_y,
                           _local_num_z, ncg);
    _old_dif_surf_corr = new Vector(_cell_locks, _local_num_x, _local_num_y,
                                    _local_num_z, NUM_FACES * ncg);
    _old_dif_surf_corr->setAll(0.0);
    _volumes = new Vector(_cell_locks, _local_num_x, _local_num_y,
                          _local_num_z, 1);
    log_printf(INFO_ONCE, "CMFD flux, source and diffusion coefficient storage"
               " = %6.2f MB", num_cells * ((4 + NUM_FACES) * ncg *
               sizeof(CMFD_PRECISION) + sizeof(omp_lock_t)) / float(1e6));

    /* Initialize k-nearest stencils, currents, flux, materials and tallies */
    generateKNearestStencils();
    initializeCurrents();
    initializeMaterials();
    allocateTallies();

#ifdef MPIx
    /* Initialize domain communicator */
    if (_domain_communicator != NULL) {
      _domain_communicator->stop = false;
      int offset = _accumulate_lmx[_domain_communicator->_domain_idx_x] +
                   _accumulate_lmy[_domain_communicator->_domain_idx_y] +
                   _accumulate_lmz[_domain_communicator->_domain_idx_z];
      _domain_communicator->_offset = offset;
      _domain_communicator->_local_num_x = _local_num_x;
      _domain_communicator->_local_num_y = _local_num_y;
      _domain_communicator->_local_num_z = _local_num_z;
      _domain_communicator->num_groups = ncg;

      int dir_sizes[3] = {num_cells / _local_num_x,  num_cells / _local_num_y,
                          num_cells / _local_num_z};

      /* Count total number of cells at all faces of the domain */
      int num_per_side[3] = {_local_num_y * _local_num_z,
                          _local_num_x * _local_num_z,
                          _local_num_x * _local_num_y};
      int num_boundary_cells = 0;
      for (int s=0; s < NUM_FACES; s++)
        num_boundary_cells += num_per_side[s % 3];

      float size = 2 * num_boundary_cells * ncg * (1 + 2 * NUM_FACES) *
                    sizeof(int) + (2 * NUM_FACES * 2 * ncg * num_boundary_cells
                    + 2 * num_boundary_cells * ncg * NUM_FACES + NUM_FACES *
                    4 * ncg * num_boundary_cells) * sizeof(CMFD_PRECISION);
      log_printf(INFO_ONCE, "CMFD domain communicator size per domain = %6.2f "
                 "MB", size / 1e6);

      /* Allocate arrays to contain information about the domain's neighbors */
      _domain_communicator->num_connections = new int*[2];
      _domain_communicator->indexes = new int**[2];
      _domain_communicator->domains = new int**[2];

      /* Arrays to contain data to communicate to/receive from other domains */
      _domain_communicator->fluxes = new CMFD_PRECISION**[2];
      _domain_communicator->coupling_coeffs = new CMFD_PRECISION**[2];
      _domain_communicator->buffer = new CMFD_PRECISION*[NUM_FACES];
      for (int rb=0; rb<2; rb++) {
        _domain_communicator->num_connections[rb] = new
              int[num_boundary_cells*ncg];
        _domain_communicator->indexes[rb] = new int*[num_boundary_cells*ncg];
        _domain_communicator->domains[rb] = new int*[num_boundary_cells*ncg];
        _domain_communicator->fluxes[rb] = new CMFD_PRECISION*[NUM_FACES];
        _domain_communicator->coupling_coeffs[rb] =
                            new CMFD_PRECISION*[num_boundary_cells*ncg];

        for (int coord=0; coord < 3; coord++) {
          for (int d=0; d < 2; d++) {
            int surf = coord + 3 * d;
            _domain_communicator->fluxes[rb][surf] =
                                new CMFD_PRECISION[dir_sizes[coord]*ncg];
            _domain_communicator->buffer[surf] =
                                new CMFD_PRECISION[2*dir_sizes[coord]*ncg];
          }
        }
        for (int nsc=0; nsc < num_boundary_cells * ncg; nsc++) {
          _domain_communicator->num_connections[rb][nsc] = 0;
          _domain_communicator->indexes[rb][nsc] = new int[NUM_FACES];
          _domain_communicator->domains[rb][nsc] = new int[NUM_FACES];
          _domain_communicator->coupling_coeffs[rb][nsc] =
                              new CMFD_PRECISION[NUM_FACES];
        }
        _domain_communicator_allocated = true;
      }

      /* Create map of continuous indexes for the domain's 6 faces : same
         order as the surfaces (XMIN=0, YMIN=1...) */
      int count = 0;
      for (int iz=0; iz < _local_num_z; iz++) {
        for (int iy=0; iy < _local_num_y; iy++) {
          int cell = (iz*_local_num_y + iy)*_local_num_x;
          _domain_communicator->mapLocalToSurface[cell] = count++;
        }
      }
      for (int iz=0; iz < _local_num_z; iz++) {
        for (int ix=1; ix < _local_num_x; ix++) {
          int cell = (iz*_local_num_y)*_local_num_x + ix;
          _domain_communicator->mapLocalToSurface[cell] = count++;
        }
      }
      for (int iy=1; iy < _local_num_y; iy++) {
        for (int ix=1; ix < _local_num_x; ix++) {
          int cell = (iy)*_local_num_x + ix;
          _domain_communicator->mapLocalToSurface[cell] = count++;
        }
      }
      for (int iz=1; iz < _local_num_z; iz++) {
        for (int iy=1; iy < _local_num_y; iy++) {
          int cell = (iz*_local_num_y + iy)*_local_num_x + _local_num_x - 1;
          _domain_communicator->mapLocalToSurface[cell] = count++;
        }
      }
      for (int iz=1; iz < _local_num_z; iz++) {
        for (int ix=1; ix < _local_num_x-1; ix++) {
          int cell = (iz*_local_num_y + _local_num_y - 1)*_local_num_x + ix;
          _domain_communicator->mapLocalToSurface[cell] = count++;
        }
      }
      for (int iy=1; iy < _local_num_y-1; iy++) {
        for (int ix=1; ix < _local_num_x-1; ix++) {
          int cell = ((_local_num_z - 1)*_local_num_y + iy)*_local_num_x + ix;
          _domain_communicator->mapLocalToSurface[cell] = count++;
        }
      }

      /* Allocate communication buffers for CMFD matrix construction */
      int storage_per_cell = ((2 + NUM_FACES) * ncg + 1);
      int internal = ncg * num_boundary_cells;
      int comm_data_size = storage_per_cell * num_boundary_cells;

      //NOTE Rank 0 is at a corner
      log_printf(INFO_ONCE, "CMFD communication buffers size per domain = "
                 "%6.2f MB", (4 * comm_data_size + internal) *
                 sizeof(CMFD_PRECISION) / float(1e6));

      _inter_domain_data = new CMFD_PRECISION[comm_data_size + internal];
      _send_domain_data = new CMFD_PRECISION[comm_data_size];

      /* Allocate memory for communication of off-domain quantities */
      _domain_data_by_surface = new CMFD_PRECISION*[NUM_FACES];
      _send_data_by_surface = new CMFD_PRECISION*[NUM_FACES];
      _boundary_volumes = new CMFD_PRECISION**[NUM_FACES];
      _boundary_reaction = new CMFD_PRECISION**[NUM_FACES];
      _boundary_diffusion = new CMFD_PRECISION**[NUM_FACES];
      _boundary_surface_currents = new CMFD_PRECISION**[NUM_FACES];

      _old_boundary_flux = new CMFD_PRECISION**[NUM_FACES];

      int start = 0;
      int ext = 0;
      for (int s=0; s < NUM_FACES; s++) {

        _domain_data_by_surface[s] = &_inter_domain_data[start];
        _send_data_by_surface[s] = &_send_domain_data[start];
        _boundary_volumes[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        _boundary_reaction[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        _boundary_diffusion[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        _boundary_surface_currents[s] = new CMFD_PRECISION*[num_per_side[s % 3]];

        _old_boundary_flux[s] = new CMFD_PRECISION*[num_per_side[s % 3]];

        for (int idx=0; idx < num_per_side[s % 3]; idx++) {

          _boundary_volumes[s][idx] = &_inter_domain_data[start];
          _boundary_reaction[s][idx] = &_inter_domain_data[start+1];
          _boundary_diffusion[s][idx] = &_inter_domain_data[start+ncg+1];
          _boundary_surface_currents[s][idx] =
            &_inter_domain_data[start+2*ncg+1];

          _old_boundary_flux[s][idx] = &_inter_domain_data[comm_data_size+ext];

          ext += ncg;
          start += storage_per_cell;
        }
      }

      /* Allocate memory for split current communication */
      int ns = NUM_FACES + NUM_EDGES;
      int vec_size = ns*ncg*sizeof(CMFD_PRECISION);
      int split_current_size = ncg * ns * num_boundary_cells;
      _send_split_current_data = new CMFD_PRECISION[split_current_size];
      _receive_split_current_data = new CMFD_PRECISION[split_current_size];
      //NOTE Rank 0 is at a corner
      log_printf(INFO_ONCE, "CMFD corner current comm. storage per domain = "
                 "%6.2f MB", 4 * split_current_size * sizeof(CMFD_PRECISION)
                 / float(1e6));

      _send_split_currents_array = new CMFD_PRECISION*[NUM_FACES];
      _receive_split_currents_array = new CMFD_PRECISION*[NUM_FACES];
      _off_domain_split_currents = new CMFD_PRECISION**[NUM_FACES];
      _received_split_currents = new CMFD_PRECISION**[NUM_FACES];

      start = 0;
      for (int s=0; s < NUM_FACES; s++) {

        _send_split_currents_array[s] =
          &_send_split_current_data[start];
        _receive_split_currents_array[s] =
          &_receive_split_current_data[start];
        _off_domain_split_currents[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        _received_split_currents[s] = new CMFD_PRECISION*[num_per_side[s % 3]];

        for (int idx=0; idx < num_per_side[s % 3]; idx++) {

          _off_domain_split_currents[s][idx] =
            &_send_split_current_data[start];
          memset(_off_domain_split_currents[s][idx], 0, vec_size);
          _received_split_currents[s][idx] =
            &_receive_split_current_data[start];
          memset(_received_split_currents[s][idx], 0, vec_size);

          start += ns*ncg;
        }
      }

      /* Allocate memory for communication of on-domain quantities */
      _send_volumes = new CMFD_PRECISION**[NUM_FACES];
      _send_reaction = new CMFD_PRECISION**[NUM_FACES];
      _send_diffusion = new CMFD_PRECISION**[NUM_FACES];
      _send_currents = new CMFD_PRECISION**[NUM_FACES];

      start = 0;
      for (int s=0; s < NUM_FACES; s++) {
        _send_volumes[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        _send_reaction[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        _send_diffusion[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        _send_currents[s] = new CMFD_PRECISION*[num_per_side[s % 3]];
        for (int idx=0; idx < num_per_side[s % 3]; idx++) {
          _send_volumes[s][idx] = &_send_domain_data[start];
          _send_reaction[s][idx] = &_send_domain_data[start+1];
          _send_diffusion[s][idx] = &_send_domain_data[start+ncg+1];
          _send_currents[s][idx] = &_send_domain_data[start+2*ncg+1];
          start += storage_per_cell;
        }
      }

      /* Calculate the starting and ending indexes of on-domain CMFD cells */
      int x_start = _accumulate_lmx[_domain_communicator->_domain_idx_x];
      int x_end = x_start + _local_num_x;
      int y_start = _accumulate_lmy[_domain_communicator->_domain_idx_y];
      int y_end = y_start + _local_num_y;
      int z_start = _accumulate_lmz[_domain_communicator->_domain_idx_z];
      int z_end = z_start + _local_num_z;

      _boundary_index_map.resize(NUM_FACES);

      /* Map connecting cells on x-surfaces */
      int global_ind;
      for (int y=0; y < _local_num_y; y++) {
        for (int z=0; z < _local_num_z; z++) {
          if (x_start > 0) {
            global_ind = ((z_start + z) * _num_y + y + y_start) *
                            _num_x + x_start - 1;
            _boundary_index_map.at(SURFACE_X_MIN)[global_ind] = z * _local_num_y
                                                              + y;
          }
          if (x_end < _num_x) {
            global_ind = ((z_start + z) * _num_y + y + y_start) *
                            _num_x + x_end;

            _boundary_index_map.at(SURFACE_X_MAX)[global_ind] = z * _local_num_y
                                                                + y;
          }
        }
      }

      /* Map connecting cells on y-surfaces */
      for (int x=0; x < _local_num_x; x++) {
        for (int z=0; z < _local_num_z; z++) {
          if (y_start > 0) {
            global_ind = ((z_start + z) * _num_y + y_start-1) *
                            _num_x + x + x_start;
            _boundary_index_map.at(SURFACE_Y_MIN)[global_ind] = z * _local_num_x
                                                                + x;
          }
          if (y_end < _num_y) {
            global_ind = ((z_start + z) * _num_y + y_end)
                          * _num_x + x + x_start;
            _boundary_index_map.at(SURFACE_Y_MAX)[global_ind] = z * _local_num_x
                                                                + x;
          }
        }
      }

      /* Map connecting cells on z-surfaces */
      for (int x=0; x < _local_num_x; x++) {
        for (int y=0; y < _local_num_y; y++) {
          if (z_start > 0) {
            global_ind = ((z_start-1) * _num_y + y + y_start) *
                          _num_x + x + x_start;
            _boundary_index_map.at(SURFACE_Z_MIN)[global_ind] = y * _local_num_x
                                                                + x;
          }
          if (z_end < _num_z) {
            global_ind = (z_end * _num_y + y + y_start) *
                            _num_x + x + x_start;
            _boundary_index_map.at(SURFACE_Z_MAX)[global_ind] = y * _local_num_x
                                                                + x;
          }
        }
      }
    }
#endif
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the CMFD mesh objects. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initialize the CMFD lattice and compute mesh dimensions, considering
 *        both uniform/non-uniform and 2D/3D cases.
 * @param offset the offset point of the CMFD Lattice
 * @param is_2D whether CMFD will be used in a 2D simulation (true) or 3D
 */
void Cmfd::initializeLattice(Point* offset, bool is_2D) {

  /* Deal with 2D case, set all widths Z to 1 */
  if (is_2D || _width_z > FLT_INFINITY) {
    _num_z = 1;
    _local_num_z = 1;
    _width_z = 1.0;
    _cell_width_z = 1.0;
    _cell_widths_z.resize(_num_z);
    _cell_widths_z[0] = _cell_width_z;
    setBoundary(SURFACE_Z_MIN, REFLECTIVE);
    setBoundary(SURFACE_Z_MAX, REFLECTIVE);
  }

  /* Handle use of X,Y and Z symmetries */
  if (_geometry->getSymmetry(0)) {

    // Compute current width of CMFD mesh
    double width_x = 0;
    for (int i=0; i<_cell_widths_x.size(); i++)
      width_x = width_x + _cell_widths_x[i];

    // If CMFD mesh was meant for full geometry, adapt it
    if (std::abs(width_x - _width_x * 2) < FLT_EPSILON) {
      if (_cell_widths_x.size() % 2 == 0)
        _cell_widths_x.resize(_cell_widths_x.size() / 2);
      else {
        _cell_widths_x.resize(_cell_widths_x.size() / 2 + 1);
        _cell_widths_x[_cell_widths_x.size() - 1] /= 2;
      }
    }
  }
  if (_geometry->getSymmetry(1)) {

    // Compute current width of CMFD mesh
    double width_y = 0;
    for (int i=0; i<_cell_widths_y.size(); i++)
      width_y = width_y + _cell_widths_y[i];

    // If CMFD mesh was meant for full geometry, adapt it
    if (std::abs(width_y - _width_y * 2) < FLT_EPSILON) {
      if (_cell_widths_y.size() % 2 == 0)
        _cell_widths_y.resize(_cell_widths_y.size() / 2);
      else {
        _cell_widths_y.resize(_cell_widths_y.size() / 2 + 1);
        _cell_widths_y[_cell_widths_y.size() - 1] /= 2;
      }
    }
  }
  if (_geometry->getSymmetry(2)) {

    // Compute current width of CMFD mesh
    double width_z = 0;
    for (int i=0; i<_cell_widths_z.size(); i++)
      width_z = width_z + _cell_widths_z[i];

    // If CMFD mesh was meant for full geometry, adapt it
    if (std::abs(width_z - _width_z * 2) < FLT_EPSILON) {
      if (_cell_widths_z.size() % 2 == 0)
        _cell_widths_z.resize(_cell_widths_z.size() / 2);
      else {
        _cell_widths_z.resize(_cell_widths_z.size() / 2 + 1);
        _cell_widths_z[_cell_widths_z.size() - 1] /= 2;
      }
    }
  }

  if (_non_uniform) {
    setNumX(_cell_widths_x.size());
    setNumY(_cell_widths_y.size());
    setNumZ(_cell_widths_z.size());
  }
  else {
    _cell_width_x = _width_x / _num_x;
    _cell_width_y = _width_y / _num_y;
    _cell_width_z = _width_z / _num_z;

    _cell_widths_x.resize(_num_x, _cell_width_x);
    _cell_widths_y.resize(_num_y, _cell_width_y);
    _cell_widths_z.resize(_num_z, _cell_width_z);
  }

  _accumulate_x.resize(_num_x+1, 0.0);
  _accumulate_y.resize(_num_y+1, 0.0);
  _accumulate_z.resize(_num_z+1, 0.0);

  for (int i=0; i<_num_x; i++)
    _accumulate_x[i+1] = _accumulate_x[i] + _cell_widths_x[i];

  for (int i=0; i<_num_y; i++)
    _accumulate_y[i+1] = _accumulate_y[i] + _cell_widths_y[i];

  for (int i=0; i<_num_z; i++)
    _accumulate_z[i+1] = _accumulate_z[i] + _cell_widths_z[i];

  if (fabs(_width_x - _accumulate_x[_num_x]) > FLT_EPSILON ||
      fabs(_width_y - _accumulate_y[_num_y]) > FLT_EPSILON ||
      fabs(_width_z - _accumulate_z[_num_z]) > FLT_EPSILON)
    log_printf(ERROR, "The sum of non-uniform mesh widths are not consistent "
               "with geometry dimensions. width_x = %20.17E, width_y = %20.17E"
               ", width_z = %20.17E, sum_x = %20.17E, sum_y = %20.17E, sum_z ="
               " %20.17E, diff_x = %20.17E, diff_y = %20.17E, diff_z = %20.17E"
               ", FLT_EPSILON = %20.17E", _width_x, _width_y, _width_z,
               _accumulate_x[_num_x], _accumulate_y[_num_y],
               _accumulate_z[_num_z], fabs(_width_x - _accumulate_x[_num_x]),
               fabs(_width_y - _accumulate_y[_num_y]),
               fabs(_width_z - _accumulate_z[_num_z]), FLT_EPSILON);

  /* Delete old lattice if it exists */
  if (_lattice != NULL)
    delete _lattice;

  /* Initialize the lattice */
  _lattice = new Lattice();
  _lattice->setNumX(_num_x);
  _lattice->setNumY(_num_y);
  _lattice->setNumZ(_num_z);

  if (_non_uniform)
    _lattice->setWidths(_cell_widths_x, _cell_widths_y, _cell_widths_z);
  else
    _lattice->setWidth(_cell_width_x, _cell_width_y, _cell_width_z);
  _lattice->setOffset(offset->getX(), offset->getY(), offset->getZ());
  _lattice->computeSizes();
}


/**
 * @brief Initializes a backup CMFD solver.
 * @details This backup solver is not necessary to run simulations, but may be
 *          used if the regular solver fails and the user wants to try another
 *          group structure without restarting the simulation.
 */
void Cmfd::initializeBackupCmfdSolver() {

  /* Initialize new CMFD object */
  _backup_cmfd = new Cmfd();
  _backup_cmfd->useAxialInterpolation(_use_axial_interpolation);
  _backup_cmfd->setLatticeStructure(_num_x, _num_y, _num_z);
  _backup_cmfd->setKNearest(_k_nearest);
  _backup_cmfd->setSORRelaxationFactor(_SOR_factor);
  _backup_cmfd->setCMFDRelaxationFactor(_relaxation_factor);
  _backup_cmfd->useFluxLimiting(_flux_limiting);

  /* Set one-group group structure */
  if (_backup_group_structure.size() == 0) {

    std::vector<int> all_groups;
    for (int e=0; e < _num_moc_groups; e++)
      all_groups.push_back(e+1);
    _backup_group_structure.push_back(all_groups);

    _cmfd_group_to_backup_group = new int[_num_cmfd_groups];
    for (int e=0; e < _num_cmfd_groups; e++)
      _cmfd_group_to_backup_group[e] = 0;
  }
  _backup_cmfd->setGroupStructure(_backup_group_structure);

  /* Set CMFD mesh boundary conditions */
  for (int i=0; i < 6; i++)
    _backup_cmfd->setBoundary(i, _boundaries[i]);

  /* Set CMFD mesh dimensions */
  _backup_cmfd->setWidthX(_width_x);
  _backup_cmfd->setWidthY(_width_y);
  _backup_cmfd->setWidthZ(_width_z);

  /* Initialize CMFD Maps */
  _backup_cmfd->initializeCellMap();

  /* Initialize the CMFD lattice */
  _backup_cmfd->initializeLattice(_lattice->getOffset());
  _backup_cmfd->setGeometry(_geometry);

#ifdef MPIx
  if (_domain_communicator != NULL) {

    _backup_cmfd->setNumDomains(_domain_communicator->_num_domains_x,
                                _domain_communicator->_num_domains_y,
                                _domain_communicator->_num_domains_z);
    _backup_cmfd->setDomainIndexes(_domain_communicator->_domain_idx_x,
                                   _domain_communicator->_domain_idx_y,
                                   _domain_communicator->_domain_idx_z);
  }
#endif

  /* Initialize the backup CMFD solver */
  _backup_cmfd->initialize();

  /* Initialize the CMFD energy group structure */
  _backup_cmfd->setSourceConvergenceThreshold(_source_convergence_threshold);
  _backup_cmfd->setNumMOCGroups(_num_moc_groups);
  _backup_cmfd->initializeGroupMap();

  /* Give CMFD number of FSRs and FSR property arrays */
  _backup_cmfd->setSolve3D(_SOLVE_3D);
  _backup_cmfd->setNumFSRs(_num_FSRs);
  _backup_cmfd->setFSRVolumes(_FSR_volumes);
  _backup_cmfd->setFSRMaterials(_FSR_materials);
  _backup_cmfd->setFSRFluxes(_FSR_fluxes);
  _backup_cmfd->setFSRSources(_FSR_sources);
  _backup_cmfd->setQuadrature(_quadrature);
  if (_flux_moments != NULL)
    _backup_cmfd->setFluxMoments(_flux_moments);

  /* Add FSRs to cells */
  _backup_cmfd->setCellFSRs(&_cell_fsrs);

  /* Initialize the backup CMFD solver */
  _backup_cmfd->initialize();
  _backup_cmfd->setConvergenceData(_convergence_data);
}


/**
 * @brief Copies the current from the regular to the backup CMFD solver.
 * @details The currents are condensed to the backup solver's energy structure
 *          when transfered as well.
 */
void Cmfd::copyCurrentsToBackup() {

  /* Clear currents */
  _backup_cmfd->zeroCurrents();

  /* Get the number of backup groups */
  int nbg = _backup_group_structure.size();

  /* Get the local current array */
  Vector* backup_currents = _backup_cmfd->getLocalCurrents();

  /* Copy on-node surface currents */
#pragma omp parallel for
  for (int i=0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

    for (int f=0; f < NUM_FACES; f++) {

      for (int e=0; e < _num_cmfd_groups; e++) {

        /* Sum group contributions and add to currents */
        int bg =  _cmfd_group_to_backup_group[e];
        CMFD_PRECISION val =
          _surface_currents->getValue(i, f * _num_cmfd_groups + e);
        backup_currents->incrementValue(i, f * nbg + bg, val);
      }
    }
  }

#ifdef MPIx
  /* Copy off-node surface currents */
  if (_domain_communicator != NULL) {

    CMFD_PRECISION*** off_node_currents =
      _backup_cmfd->getBoundarySurfaceCurrents();

    for (int surface=0; surface < NUM_FACES; surface++) {

      /* Extract arrays on surface */
      CMFD_PRECISION** boundary_currents = _boundary_surface_currents[surface];
      CMFD_PRECISION** backup_currents = off_node_currents[surface];

      /* Loop over all CMFD cells on the current surface */
      std::map<int, int>::iterator it;
      for (it=_boundary_index_map.at(surface).begin();
           it != _boundary_index_map.at(surface).end(); ++it) {

        int idx = it->second;

        /* Loop over cell faces */
        for (int f=0; f < NUM_FACES; f++) {

          /* Loop over CMFD coarse energy groups */
          for (int e = 0; e < _num_cmfd_groups; e++) {
            int bg =  _cmfd_group_to_backup_group[e];
            backup_currents[idx][f*nbg + bg] +=
              boundary_currents[idx][f*_num_cmfd_groups+e];
          }
        }
      }
    }
  }
#endif
}


/**
 * @brief Returns the width of a given surface
 * @param surface A surface index, from 0 to NUM_FACES - 1
 * @param global_ind global index of a CMFD cell
 * @return The surface width
 */
CMFD_PRECISION Cmfd::getSurfaceWidth(int surface, int global_ind) {

  CMFD_PRECISION width;

  int ix = global_ind % _num_x;
  int iy = (global_ind % (_num_x * _num_y)) / _num_x;
  int iz = global_ind / (_num_x * _num_y);

  if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
    return _cell_widths_y[iy] * _cell_widths_z[iz];
  else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX)
    return _cell_widths_x[ix] * _cell_widths_z[iz];
  else
    return _cell_widths_x[ix] * _cell_widths_y[iy];
}


/**
 * @brief Returns the width of the surface perpendicular to a given surface
 * @param surface A surface index, from 0 to NUM_FACES - 1
 * @param global_ind The CMFD cell global index
 * @return The perpendicular surface width
 */
CMFD_PRECISION Cmfd::getPerpendicularSurfaceWidth(int surface, int global_ind) {

  int ix = global_ind % _num_x;
  int iy = (global_ind % (_num_x * _num_y)) / _num_x;
  int iz = global_ind / (_num_x * _num_y);

  if (surface == SURFACE_X_MIN || surface == SURFACE_X_MAX)
    return _cell_widths_x[ix];
  else if (surface == SURFACE_Y_MIN || surface == SURFACE_Y_MAX)
    return _cell_widths_y[iy];
  else
    return _cell_widths_z[iz];
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
#ifndef THREED
  _SOLVE_3D = solve_3D;
#endif
}


/**
 * @brief Sets the azimuthal spacings.
 * @param azim_spacings An array of azimuthal spacings for each azimuthal angle.
 * @param num_azim the number of azimuthal angles.
 */
void Cmfd::setAzimSpacings(const std::vector<double>& azim_spacings,
                           int num_azim) {

  if (_azim_spacings != NULL)
    delete [] _azim_spacings;

  _azim_spacings = new double[num_azim/4];

  for (int a=0; a < num_azim/4; a++)
    _azim_spacings[a] = double(azim_spacings[a]);
}


/**
 * @brief Sets the polar spacings.
 * @param polar_spacings A 2D array of polar spacings for each azimuthal and
 *        polar angle combination.
 * @param num_azim the number of azimuthal angles.
 * @param num_polar the number of polar angles.
 */
void Cmfd::setPolarSpacings(const std::vector< std::vector<double> >&
                            polar_spacings, int num_azim, int num_polar) {

  if (_polar_spacings != NULL) {
    for (int a=0; a < num_azim/4; a++)
      delete [] _polar_spacings[a];
    delete [] _polar_spacings;
  }

  _polar_spacings = new double*[num_azim/4];
  for (int a=0; a < num_azim/4; a++)
    _polar_spacings[a] = new double[num_polar/2];

  for (int a=0; a < num_azim/4; a++) {
    for (int p=0; p < num_polar/2; p++)
      _polar_spacings[a][p] = double(polar_spacings[a][p]);
  }
}


/**
 * @brief Set the value of the k effective for the CMFD solver. This is meant
 *        for research / debugging purposes.
 * @param k_eff the k_eff value to set.
 */
void Cmfd::setKeff(double k_eff) {
  _k_eff = k_eff;
}


/**
 * @brief Set the backup CMFD solver's group structure. It is necessarily
 *        coarser than and must align with the regular CMFD group structure.
 * @param group_indices the indices of the CMFD groups in the MOC groups
 */
void Cmfd::setBackupGroupStructure(std::vector< std::vector<int> >
                                   group_indices) {

  /* Assign the number of backup energy groups */
  _num_backup_groups = group_indices.size();

  /* Initialize mappings */
  _backup_group_structure = group_indices;
  _cmfd_group_to_backup_group = new int[_num_cmfd_groups];
  for (int e=0; e < _num_cmfd_groups; e++)
    _cmfd_group_to_backup_group[e] = -1;

  /* Check that the mapping is valid and assign CMFD groups to backup groups */
  int cmfd_group = -1;
  int moc_group = 0;
  for (int i=0; i < group_indices.size(); i++) {
    for (int j=0; j < group_indices.at(i).size(); j++) {
      if (group_indices.at(i).at(j) != moc_group + 1) {
        log_printf(ERROR, "Invalid backup group structure: indices must be "
                   "monotonic and include all MOC groups.");
      }
      if (moc_group >= _group_indices[cmfd_group+1]) {
        cmfd_group++;
        _cmfd_group_to_backup_group[cmfd_group] = i;
      }

      if (i != _cmfd_group_to_backup_group[cmfd_group])
        log_printf(ERROR, "Invalid backup group structure: indices of backup "
                   "group structure must align with boundaries of CMFD group "
                   "structure.");

      moc_group++;
    }
  }

  /* Ensure that every CMFD group has a backup group */
  for (int e=0; e < _num_cmfd_groups; e++) {
    if (_cmfd_group_to_backup_group[e] == -1)
      log_printf(ERROR, "Invalid backup group structure: failed to find "
                 "matching index for CMFD group %d", e);
  }
}


/**
 * @brief A function that prints a summary of the CMFD input parameters.
 */
void Cmfd::printInputParamsSummary() {

  if (_flux_update_on)
    log_printf(NORMAL, "CMFD acceleration: ON");
  else
    log_printf(NORMAL, "CMFD acceleration: OFF (no MOC flux update)");

  if (_flux_update_on) {
    // Print CMFD relaxation information
    if (std::abs(_SOR_factor - 1) > FLT_EPSILON)
      log_printf(NORMAL, "CMFD inner linear solver SOR factor: %f",
                 _SOR_factor);
    log_printf(NORMAL, "CMFD corrected diffusion coef. relaxation factor: %f",
               _relaxation_factor);

    // Print CMFD interpolation techniques
    if (_centroid_update_on)
      log_printf(NORMAL, "CMFD K-nearest scheme: %d neighbors", _k_nearest);
    if (_use_axial_interpolation == 1)
      log_printf(NORMAL, "CMFD axial interpolation with axially averaged "
                 "update ratios");
    else if (_use_axial_interpolation == 2)
      log_printf(NORMAL, "CMFD axial interpolation with update ratios evaluated"
                 " at centroid Z-coordinate");

    // Print other CMFD modifications
    if (_flux_limiting)
      log_printf(INFO_ONCE, "CMFD corrected diffusion coef. bounded by "
                 "regular diffusion coef.");
    if (_balance_sigma_t)
      log_printf(INFO_ONCE, "CMFD total cross sections adjusted for matching "
                 "MOC reaction rates");
  }

  // Print CMFD space and energy mesh information
  log_printf(NORMAL, "CMFD Mesh: %d x %d x %d", _num_x, _num_y, _num_z);
  if (_flux_update_on) {
    if (_num_cmfd_groups != _num_moc_groups) {
      log_printf(NORMAL, "CMFD Group Structure:");
      log_printf(NORMAL, "\t MOC Group \t CMFD Group");
      for (int g=0; g < _num_moc_groups; g++)
        log_printf(NORMAL, "\t %d \t\t %d", g+1, getCmfdGroup(g)+1);
    }
    else
      log_printf(NORMAL, "CMFD and MOC group structures match");
  }
}


/**
 * @brief Report the physical time use by major components of the CMFD solver.
 */
void Cmfd::printTimerReport() {

  std::string msg_string;

  /* Get the total CMFD time */
  double tot_time = _timer->getSplit("Total CMFD time");
  msg_string = "  Total CMFD computation time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);

  /* Get the total XS collapse time */
  double xs_collapse_time = _timer->getSplit("Total collapse time");
  msg_string = "    XS collapse time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), xs_collapse_time);

  /* Get the total matrix construction time */
  double matrix_construction_time = _timer->getSplit("Matrix construction time");
  msg_string = "    Matrix construction time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), matrix_construction_time);

#ifdef MPIx
  /* Get the MPI communication time */
  double comm_time = _timer->getSplit("CMFD MPI communication time");
  msg_string = "    MPI communication time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), comm_time);
#endif

  /* Get the total solver time */
  double solver_time = _timer->getSplit("Total solver time");
  msg_string = "    Total CMFD solver time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), solver_time);

  /* Get the total MOC flux update time */
  double update_time = _timer->getSplit("Total MOC flux update time");
  msg_string = "    Total flux update time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), update_time);
}


/**
 * @brief Forms a full copy of the surface currents on every surface
 * @details The copy contains all surface currents including edge and corner
 *          currents explicitly. It is stored in the _full_surface_currents
 *          vector for use in debugging and diagnostics.
 */
void Cmfd::copyFullSurfaceCurrents() {

  /* Allocate full surface currents if necessary */
  if (_full_surface_currents == NULL)
    _full_surface_currents = new Vector(_cell_locks, _local_num_x,
                                        _local_num_y, _local_num_z,
                                        _num_cmfd_groups * NUM_SURFACES);

  /* Clear the currently saved surface currents */
  _full_surface_currents->clear();

  /* Copy surface currents from surface faces */
  for (int i=0; i < _local_num_x * _local_num_y * _local_num_z; i++) {
    for (int s=0; s < NUM_FACES; s++) {
      for (int g=0; g < _num_cmfd_groups; g++) {
        FP_PRECISION current =
          _surface_currents->getValue(i, s * _num_cmfd_groups + g);
        _full_surface_currents->incrementValue(i, s * _num_cmfd_groups + g,
                                               current);
      }
    }
  }

  /* Copy surface currents from edges and corners */
  std::map<int, CMFD_PRECISION>::iterator it;
  for (it = _edge_corner_currents.begin();
       it != _edge_corner_currents.end(); ++it) {
    int key = it->first;
    int cell = key / (_num_cmfd_groups * NUM_SURFACES);
    int surf_group = key % (_num_cmfd_groups * NUM_SURFACES);
    _full_surface_currents->incrementValue(cell, surf_group, it->second);
  }
}


/**
 * @brief Computes the neutron balance over each CMFD cell for both MOC and CMFD
 * @details This routine can be used once the CMFD matrices have been formed to
 *          compute the neutron balance in the CMFD cell. With regards to MOC,
 *          it loops over all fsrs in the cell to compute all reaction rates
 *          and currents.
 * //NOTE : Expect a neutron imbalance : - in CMFD if there is CMFD relaxation
 *          - in MOC (and CMFD) if using the newly computed source, which
 *          is not converged with regards to group-to-group scattering
 *          - in MOC at the boundaries, as the incoming currents are not
 *          tallied (except at ite 0, they are null)
 *          - in MOC at the reflective boundaries when tracks hit edges and
 *          corners as the contributions are double tallied (not a bug, only
 *          a problem for this routine, see NOTE for where to modify code)
 * @param pre_split whether edge currents are not split (default true)
 * @param moc_balance whether to check the MOC balance over the cell or
 *        compare the MOC and CMFD imbalance (due to group to group scattering)
 */
void Cmfd::checkNeutronBalance(bool pre_split, bool moc_balance) {

  /* Print a few warnings on routine limitations (exhaustive in docstring) */
  if (_geometry->isDomainDecomposed() && pre_split)
    log_printf(WARNING_ONCE, "MOC neutron balance is currently not checked "
               "correctly for domain decomposed geometries.");
  if (_relaxation_factor != 1.0)
    log_printf(WARNING_ONCE, "CMFD relaxation factor is not 1.0, expect CMFD "
               "vs MOC neutron imbalance after iteration 0.");
  if (moc_balance)
    log_printf(WARNING_ONCE, "MOC neutron imbalance is expected at the "
               "boundaries past the first iteration as incoming currents are "
               "not tallied.");

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
  /* Compute neutron production */
  matrixMultiplication(_M, _old_flux, &m_phi);

  /* Compute neutron transfer and loss */
  matrixMultiplication(_A, _old_flux, &a_phi);
  CMFD_PRECISION* a_phi_array = a_phi.getArray();

#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    int* coupling_sizes = NULL;
    int** coupling_indexes = NULL;
    CMFD_PRECISION** coupling_coeffs = NULL;
    CMFD_PRECISION** coupling_fluxes = NULL;
    int offset = 0;
    for (int color=0; color < 2; color++) {

      getCouplingTerms(_domain_communicator, color, coupling_sizes,
                       coupling_indexes, coupling_coeffs, coupling_fluxes,
                       _old_flux->getArray(), offset);

#pragma omp parallel for collapse(2)
      for (int iz=0; iz < _local_num_z; iz++) {
        for (int iy=0; iy < _local_num_y; iy++) {
          for (int ix=(iy+iz+color+offset)%2; ix < _local_num_x; ix+=2) {
            int cell = (iz*_local_num_y + iy)*_local_num_x + ix;
            bool on_surface = (iz==0) || (iz==_local_num_z-1) || (iy==0) ||
                 (iy==_local_num_y-1) || (ix==0) || (ix==_local_num_x-1);

            for (int g=0; g < _num_cmfd_groups; g++) {

              int row = cell * _num_cmfd_groups + g;

              if (on_surface) {
                int row_surf = _domain_communicator->mapLocalToSurface[cell];
                for (int i = 0; i < coupling_sizes[row_surf]; i++) {
                  int idx = coupling_indexes[row_surf][i] * _num_cmfd_groups
                       + g;
                  int domain =
                       _domain_communicator->domains[color][row_surf][i];
                  CMFD_PRECISION flux = coupling_fluxes[domain][idx];
                  a_phi_array[row] += coupling_coeffs[row_surf][i] * flux;
                }
              }
            }
          }
        }
      }
    }
  }
#endif

  int num_imbalanced = 0;
  double max_imbalance = 0.0;
  double max_imbalance_moc, max_imbalance_cmfd;
  int max_imbalance_cell = -1;
  int max_imbalance_grp = -1;

  /* Compute MOC balance */
  /* Loop over CMFD cells */
  for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

    bool imbalance_reported = false;
    int x = (i % (_local_num_x * _local_num_y)) % _local_num_x;
    int y = (i % (_local_num_x * _local_num_y)) / _local_num_x;
    int z = i / (_local_num_x * _local_num_y);

    Material* cell_material = _materials[i];

    /* Loop over CMFD coarse energy groups */
    for (int e = 0; e < _num_cmfd_groups; e++) {

      /* Initialize tallies */
      double total = 0.0;
      double in_scattering = 0.0;
      double fission = 0.0;

      /* Loop over FSRs in CMFD cell */
      for (int j = 0; j < _cell_fsrs.at(i).size(); j++) {

        long fsr_id = _cell_fsrs.at(i).at(j);
        Material* fsr_material = _FSR_materials[fsr_id];
        FP_PRECISION volume = _FSR_volumes[fsr_id];
        FP_PRECISION* scat = fsr_material->getSigmaS();
        FP_PRECISION* flux = &_FSR_fluxes[fsr_id*_num_moc_groups];

        /* Loop over MOC energy groups within this CMFD coarse group */
        double chi = 0.0;
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++)
          chi += fsr_material->getChiByGroup(h+1);

        /* Use new fluxes to compute the source terms */
        if (!moc_balance) {
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
        }
        /* Use the old MOC source to check neutron balance */
        else {
          in_scattering = 0.0;  // for convenience
          fission += _FSR_sources[fsr_id*_num_moc_groups + e] * volume * FOUR_PI;
        }

        /* Calculate total reaction rate in this CMFD coarse group */
        for (int h = _group_indices[e]; h < _group_indices[e+1]; h++) {
          double tot = fsr_material->getSigmaTByGroup(h+1);
          total += tot * flux[h] * volume;
        }
      }

      /* Calculate net current out of the cell */
      double net_current = 0.0;

      /* Use currents before splitting edges/corners if requested */
      if (pre_split) {

        /* Create arrays of cell indexes and bounds */
        int cell_limits[3] = {_local_num_x, _local_num_y, _local_num_z};
        int cell_ind[3];
        cell_ind[0] = i % _local_num_x;
        cell_ind[1] = (i / _local_num_x) % _local_num_y;
        cell_ind[2] = i / (_local_num_x * _local_num_y);

        /* Tally current from all surfaces including edges and corners */
        for (int s=0; s < NUM_SURFACES; s++) {

          /* Compute index and vector direction */
          int idx = s * _num_cmfd_groups + e;
          int direction[3];
          convertSurfaceToDirection(s, direction);

          /* Compute the next CMFD cell from the cell indexes and direction */
          int cmfd_cell_next = 0;
          int cell_next_ind[3];
          for (int d=0; d < 3; d++)
            cell_next_ind[d] = cell_ind[d] + direction[d];

          cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x
                         + cell_next_ind[2] * (_local_num_x * _local_num_y);

          /* Compute the opposite direction vector */
          int op_direction[3];
          for (int d=0; d < 3; d++)
            op_direction[d] = -1 * direction[d];

          /* Determine if the next CMFD cell is within the bounds */
          for (int d=0; d < 3; d++)
            if (cell_next_ind[d] < 0 || cell_next_ind[d] >= cell_limits[d])
              cmfd_cell_next = -1;

          /* If the cell is outside the bounds, handle boundaries */
          if (cmfd_cell_next == -1) {

            /* Booleans for determining surface boundary type */
            bool vacuum = false;
            bool reflective = false;
            bool transmit_avail = false;

            int transmit_direction[3] = {0,0,0};

            /* Loop over all directions to handle boundaries */
            for (int d=0; d < 3; d++) {
              if (cell_next_ind[d] < 0 || cell_next_ind[d] >= cell_limits[d]) {

                /* Form the surface for each direction */
                int partial_direction[3] = {0,0,0};
                partial_direction[d] = direction[d];
                int partial_surface =
                  convertDirectionToSurface(partial_direction);

                /* Look at the boundary type in this direction */
                if (_boundaries[partial_surface] == VACUUM) {
                  vacuum = true;
                }
                else if (_boundaries[partial_surface] == REFLECTIVE) {
                  reflective = true;
                  op_direction[d] *= -1;
                }
              }

              /* For non-boundary surfaces, save the direction */
              else if (direction[d] != 0) {
                transmit_avail = true;
                transmit_direction[d] = direction[d];
              }
            }

            /* For vacuum boundaries, tally the leakage */
            if (vacuum)
              net_current += _full_surface_currents->getValue(i, idx);

            /* For reflective boundaries, find the appropriate cell to
               deliver current if available */
            else if (reflective && transmit_avail) {
              for (int d=0; d < 3; d++)
                cell_next_ind[d] = cell_ind[d] + transmit_direction[d];
              cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] *
                              _local_num_x + cell_next_ind[2] *
                              (_local_num_x * _local_num_y);
            }
          }

          /* Transmit current to available cells */
          if (cmfd_cell_next != -1) {
            int surface_next = convertDirectionToSurface(op_direction);
            int idx_next = surface_next * _num_cmfd_groups + e;
            net_current += _full_surface_currents->getValue(i, idx) -
              _full_surface_currents->getValue(cmfd_cell_next, idx_next);
          }
        }
      }

      /* Use post-split edges/corner currents if requested */
      else {

        /* Tally current only over surface faces */
        for (int s = 0; s < NUM_FACES; s++) {
          int idx = s * _num_cmfd_groups + e;
          int cmfd_cell_next = getCellNext(i, s, false);
          int surface_next = (s + NUM_FACES / 2) % NUM_FACES;
          int idx_next = surface_next * _num_cmfd_groups + e;

          /* Out of domain currents */
          if (cmfd_cell_next == -1) {

            int i_global = getGlobalCMFDCell(i);
            int global_cmfd_cell_next = getCellNext(i_global, s);

            /* Currents at the geometry outer boundaries */
            if ((global_cmfd_cell_next == -1 && _boundaries[s] == VACUUM) ||
                (moc_balance))
              net_current += _surface_currents->getValue(i, idx);
            //FIXME For moc_balance, this is correct at iteration 0 only

            /* Currents from other domains */
            else if (global_cmfd_cell_next != -1) {
              int idx_off =
                   _boundary_index_map.at(s)[global_cmfd_cell_next];
              net_current += _surface_currents->getValue(i, idx) -
                   _boundary_surface_currents[s][idx_off][surface_next*
                                                         _num_cmfd_groups+e];
            }
          }
          /* In-domain currents */
          else
            net_current += _surface_currents->getValue(i, idx) -
                _surface_currents->getValue(cmfd_cell_next, idx_next);
        }
      }

      /* Compute balance in given cell and group */
      double moc_balance = in_scattering + fission - total - net_current;

      double cmfd_balance = m_phi.getValue(i, e) / _k_eff -
            a_phi.getValue(i, e);

      double tmp_imbalance = std::abs(moc_balance);
      if (!moc_balance)
        tmp_imbalance = std::abs(moc_balance - cmfd_balance);

      if (tmp_imbalance > max_imbalance) {
        max_imbalance = tmp_imbalance;
        max_imbalance_moc = moc_balance;
        max_imbalance_cmfd = cmfd_balance;
        max_imbalance_cell = i;
        max_imbalance_grp = e;
      }

      /* Select appropriate level of numerical tolerance */
      //NOTE Track angular fluxes being in single precision, double precision
      //     neutron balance is not expected.
#ifdef SINGLE
      FP_PRECISION tol = 5e-6;
#else
      FP_PRECISION tol = 5e-7;
#endif

      /* Report an abnormal neutron imbalance */
      if (tmp_imbalance > tol) {

        /* Output an imbalance only once for all groups on log level INFO */
        if (!imbalance_reported) {
          num_imbalanced++;
          log_printf(INFO, "Neutron imbalance in cell (%d, %d, %d) for CMFD "
                     "group %d = MOC %g, CMFD %g, diff %g", x, y, z, e,
                     moc_balance, cmfd_balance, moc_balance - cmfd_balance);
        }
        /* Output an imbalance for every group at log level DEBUG */
        else
          log_printf(DEBUG, "Neutron imbalance in cell (%d, %d, %d) for CMFD "
                     "group %d = MOC %g, CMFD %g, diff %g", x, y, z, e,
                     moc_balance, cmfd_balance, moc_balance - cmfd_balance);
        imbalance_reported = true;
      }
    }
  }

  /* Output max imbalance in CMFD mesh */
  int x = (max_imbalance_cell % (_local_num_x * _local_num_y)) % _local_num_x;
  int y = (max_imbalance_cell % (_local_num_x * _local_num_y)) / _local_num_x;
  int z = max_imbalance_cell / (_local_num_x * _local_num_y);
  if (moc_balance) {
    log_printf(NODAL, "Maximum neutron imbalance MOC %.2e (CMFD %.2e) at cell "
               "%i (%d %d %d) and group %d.", max_imbalance_moc,
               max_imbalance_cmfd, max_imbalance_cell, x, y, z,
               max_imbalance_grp);
    log_printf(NODAL, "%d CMFD cells report a MOC neutron imbalance",
               num_imbalanced);
  }
  else {
    log_printf(NODAL, "Maximum neutron imbalance between MOC and CMFD : %.2e "
               "(MOC %.2e CMFD %.2e) at cell %i (%d %d %d) and group %d.",
               max_imbalance_moc - max_imbalance_cmfd, max_imbalance_moc,
               max_imbalance_cmfd, max_imbalance_cell, x, y, z,
               max_imbalance_grp);
    log_printf(NODAL, "%d CMFD cells report a neutron imbalance between MOC "
               "and CMFD", num_imbalanced);
  }
}


/**
 * @brief Returns the color of a CMFD cell in the red/black SOR solver
 * @param cmfd_cell The cmfd cell's global ID
 */
int Cmfd::getCellColor(int cmfd_cell) {
  int ix = cmfd_cell % _num_x;
  int iy = (cmfd_cell % (_num_x * _num_y)) / _num_x;
  int iz = cmfd_cell / (_num_x * _num_y);
  int color = (ix + iy + iz) % 2;
  return color;
}


/**
 * @brief Packs reaction rates and currents into buffers for communication.
 * @details Buffer description is found in ghostCellExchange's docstring
 */
#ifdef MPIx
void Cmfd::packBuffers() {

  int current_idx[6] = {0,0,0,0,0,0};
  bool found_surfaces[NUM_FACES];

  for (int z=0; z < _local_num_z; z++) {
    for (int y=0; y < _local_num_y; y++) {
      for (int x=0; x < _local_num_x; x++) {
        for (int s=0; s < NUM_FACES; s++)
          found_surfaces[s] = false;

        /* Check that cell is at a boundary */
        if (x == 0)
          found_surfaces[SURFACE_X_MIN] = true;
        if (x == _local_num_x-1)
          found_surfaces[SURFACE_X_MAX] = true;
        if (y == 0)
          found_surfaces[SURFACE_Y_MIN] = true;
        if (y == _local_num_y-1)
          found_surfaces[SURFACE_Y_MAX] = true;
        if (z == 0)
          found_surfaces[SURFACE_Z_MIN] = true;
        if (z == _local_num_z-1)
          found_surfaces[SURFACE_Z_MAX] = true;

        /* Fill buffers with tallies */
        for (int s=0; s < NUM_FACES; s++) {
          if (found_surfaces[s]) {
            int idx = current_idx[s];
            int cell_id = ((z * _local_num_y) + y) * _local_num_x + x;
            _send_volumes[s][idx][0] = _volume_tally[cell_id][0];
            for (int e=0; e < _num_cmfd_groups; e++) {
              _send_reaction[s][idx][e] = _reaction_tally[cell_id][e];
              _send_diffusion[s][idx][e] = _diffusion_tally[cell_id][e];
              for (int f=0; f < NUM_FACES; f++) {
                _send_currents[s][idx][f*_num_cmfd_groups + e] =
                  _surface_currents->getValue(cell_id, f*_num_cmfd_groups+e);
              }
            }
            current_idx[s]++;
          }
        }
      }
    }
  }
}


/**
 * @brief Exchanges ghost cell buffers in 3D cartesian (i.e., 6 directions)
 * @details comm The cartesian MPI domain communicator object that is
 *          configured for the CMFD exchange
 *          send_buffers A 2D array of floating point data. The outer dimension
 *          corresponds to each face of the domain, while the inner dimension
 *          is the serialized buffer corresponding to the number of 2D cells
 *          to exchange times the number of energy groups.
 *          recv_buffers A 2D array of floating point data. The outer dimension
 *          corresponds to each face of the domain,  while the inner dimension
 *          is the serialized buffer corresponding to the number of 2D cells to
 *          exchange times the number of energy groups.
 */
void Cmfd::ghostCellExchange() {

  packBuffers();

  MPI_Request requests[2*NUM_FACES];

  MPI_Datatype precision;
  if (sizeof(CMFD_PRECISION) == 4)
    precision = MPI_FLOAT;
  else
    precision = MPI_DOUBLE;

  int storage_per_cell = ((2 + NUM_FACES) * _num_cmfd_groups + 1);

  int sizes[NUM_FACES];
  for (int coord=0; coord < 3; coord++) {
    for (int d=0; d<2; d++) {

      int dir = 2*d - 1;
      int surf = coord + 3*d;
      int op_surf = surf - 3*dir;
      int source, dest;

      // Figure out serialized buffer length for this face
      int size = 0;
      if (surf == SURFACE_X_MIN) {
        size = _local_num_y * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_X_MAX) {
        size = _local_num_y * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_Y_MIN) {
        size = _local_num_x * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_Y_MAX) {
        size = _local_num_x * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_Z_MIN) {
        size = _local_num_x * _local_num_y * storage_per_cell;
      }
      else if (surf == SURFACE_Z_MAX) {
        size = _local_num_x * _local_num_y * storage_per_cell;
      }

      sizes[surf] = size;

      /* Get ranks of source and destination domains, using the cartesian
         structure of the domain decomposition */
      MPI_Cart_shift(_domain_communicator->_MPI_cart, coord, dir, &source,
                     &dest);

      // Post send
      MPI_Isend(_send_data_by_surface[surf], size, precision,
          dest, 0, _domain_communicator->_MPI_cart, &requests[2*surf]);

      // Post receive
      MPI_Irecv(_domain_data_by_surface[op_surf], size, precision,
          source, 0, _domain_communicator->_MPI_cart, &requests[2*surf+1]);
    }
  }

  // Block for communication round to complete
  MPI_Waitall(12, requests, MPI_STATUSES_IGNORE);
}


/**
 * @brief Communicate split (at corners and edges) currents (respectively edge
 *        and face currents) to other domains.
 * @param faces whether the currents are for edges or faces, for unpacking
 */
void Cmfd::communicateSplits(bool faces) {

  //TODO: Form into communicateEdgeCurrents and communicateFaceCurrents
  // 1. communicate edge currents use array of length NUM_EDGES, called after
  //    vertex splits
  // 2. communicate face currents use array of length NUM_FACES, called after
  //    edge splits
  // NOTE: only communicate currents that are saved OFF DOMAIN at each step
  // NOTE: communicateFaceCurrents will use currents formed from vertex splits

  MPI_Request requests[2*NUM_FACES];

  MPI_Datatype precision;
  if (sizeof(CMFD_PRECISION) == 4)
    precision = MPI_FLOAT;
  else
    precision = MPI_DOUBLE;

  int storage_per_cell = (NUM_FACES + NUM_EDGES) * _num_cmfd_groups;

  int sizes[NUM_FACES];
  for (int coord=0; coord < 3; coord++) {
    for (int d=0; d<2; d++) {

      int dir = 2*d-1;
      int surf = coord + 3*d;
      int op_surf = surf - 3*dir;
      int source, dest;

      // Figure out serialized buffer length for this face
      int size = 0;
      if (surf == SURFACE_X_MIN) {
        size = _local_num_y * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_X_MAX) {
        size = _local_num_y * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_Y_MIN) {
        size = _local_num_x * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_Y_MAX) {
        size = _local_num_x * _local_num_z * storage_per_cell;
      }
      else if (surf == SURFACE_Z_MIN) {
        size = _local_num_x * _local_num_y * storage_per_cell;
      }
      else if (surf == SURFACE_Z_MAX) {
        size = _local_num_x * _local_num_y * storage_per_cell;
      }

      sizes[surf] = size;

      MPI_Cart_shift(_domain_communicator->_MPI_cart, coord, dir, &source, &dest);

      // Post send
      MPI_Isend(_send_split_currents_array[surf], size, precision,
          dest, 0, _domain_communicator->_MPI_cart, &requests[2*surf]);

      // Post receive
      MPI_Irecv(_receive_split_currents_array[op_surf], size, precision,
          source, 0, _domain_communicator->_MPI_cart, &requests[2*surf+1]);
    }
  }

  // Block for communication round to complete
  MPI_Waitall(12, requests, MPI_STATUSES_IGNORE);

  unpackSplitCurrents(faces);
}


/**
 * @brief Unpacks communicated split current data
 * @param faces Whether to split the currents onto surface faces
 */
void Cmfd::unpackSplitCurrents(bool faces) {

  int current_idx[6] = {0,0,0,0,0,0};
  bool found_surfaces[NUM_FACES];

  /* Loop over all CMFD cells */
  for (int z=0; z < _local_num_z; z++) {
    for (int y=0; y < _local_num_y; y++) {
      for (int x=0; x < _local_num_x; x++) {

        /* Look for boundaries touching the CMFD cell */
        for (int s=0; s < NUM_FACES; s++)
          found_surfaces[s] = false;
        if (x == 0)
          found_surfaces[SURFACE_X_MIN] = true;
        if (x == _local_num_x-1)
          found_surfaces[SURFACE_X_MAX] = true;
        if (y == 0)
          found_surfaces[SURFACE_Y_MIN] = true;
        if (y == _local_num_y-1)
          found_surfaces[SURFACE_Y_MAX] = true;
        if (z == 0)
          found_surfaces[SURFACE_Z_MIN] = true;
        if (z == _local_num_z-1)
          found_surfaces[SURFACE_Z_MAX] = true;

        /* Handle all boundaries */
        for (int s=0; s < NUM_FACES; s++) {
          if (found_surfaces[s]) {

            /* Convert the (x,y,z) indexes to a cell ID and boundary index */
            int cell_id = ((z * _local_num_y) + y) * _local_num_x + x;
            int idx = current_idx[s];

            /* Copy the appropriate face or edge information */
            if (faces) {

              /* Treat CMFD cell face currents */
              for (int f=0; f < NUM_FACES; f++) {
                for (int g=0; g < _num_cmfd_groups; g++) {

                  /* Get the face current value */
                  CMFD_PRECISION value =
                    _received_split_currents[s][idx][f * _num_cmfd_groups + g];

                  /* Treat nonzero values */
                  if (fabs(value) > FLUX_EPSILON)
                    _surface_currents->incrementValue(cell_id,
                                                      f * _num_cmfd_groups + g,
                                                      value);
                }
              }
            }
            else {

              /* Treat CMFD cell edge currents */
              for (int e=NUM_FACES; e < NUM_EDGES+NUM_FACES; e++) {

                int surf_idx = cell_id * NUM_SURFACES * _num_cmfd_groups + e *
                  _num_cmfd_groups;

                for (int g=0; g < _num_cmfd_groups; g++) {

                  /* Get the edge current value */
                  CMFD_PRECISION value =
                    _received_split_currents[s][idx][e * _num_cmfd_groups + g];

                  /* Treat nonzero values */
                  if (fabs(value) > FLUX_EPSILON) {

                    int new_ind = surf_idx + g;

                    /* Add the contribution */
                    _edge_corner_currents[new_ind] += value;
                  }
                }
              }
            }

            /* Increment the boundary index */
            current_idx[s]++;
          }
        }
      }
    }
  }
}
#endif


/**
 * @brief Converts a global CMFD cell ID into its local ID
 * @details Marked for deletion, but still used thoroughly.
 * @param cmfd_cell The global CMFD cell ID
 * @return The local CMFD cell ID, -1 if not in the domain.
 */
int Cmfd::getLocalCMFDCell(int cmfd_cell) {

  int x_start = 0;
  int y_start = 0;
  int z_start = 0;
  int x_end = _num_x;
  int y_end = _num_y;
  int z_end = _num_z;
  if (_geometry->isDomainDecomposed()) {
    if (_domain_communicator != NULL) {
      x_start = _accumulate_lmx[_domain_communicator->_domain_idx_x];
      x_end = x_start + _local_num_x;
      y_start = _accumulate_lmy[_domain_communicator->_domain_idx_y];
      y_end = y_start + _local_num_y;
      z_start = _accumulate_lmz[_domain_communicator->_domain_idx_z];
      z_end = z_start + _local_num_z;
    }
  }

  int ix = (cmfd_cell % (_num_x * _num_y)) % _num_x;
  int iy = (cmfd_cell % (_num_x * _num_y)) / _num_x;
  int iz = cmfd_cell / (_num_x * _num_y);

  int local_cmfd_cell;
  if (ix >= x_start && ix < x_end && iy >= y_start && iy < y_end &&
      iz >= z_start && iz < z_end)
    local_cmfd_cell = ((iz - z_start) * _local_num_y + iy - y_start)
                      * _local_num_x + ix - x_start;
  else
    local_cmfd_cell = -1;

  return local_cmfd_cell;
}


/**
 * @brief Converts a 3 integer vector direction into a surface
 * @details The direction is a tuplet with each value taking either
 *          +1 (positive directed), 0 (neutral), or -1 (negative directed)
 * @param direction The integer vector describing the direction
 * @return The surface associated with traveling the provided direction from
 *         the origin of the cell
 */
int Cmfd::convertDirectionToSurface(int* direction) {
  int surface = 0;
  int num_crossings = std::abs(direction[0]) + std::abs(direction[1]) +
    std::abs(direction[2]);
  if (num_crossings == 1) {
    for (int i=0; i < 3; i++) {
      int present = std::abs(direction[i]);
      int fwd = (direction[i] + 1) / 2;
      surface += present * (3 * fwd + i);
    }
  }
  else if (num_crossings == 2) {
    surface += NUM_FACES;
    int ind1 = 0;
    int ind2 = 0;
    if (direction[0] == 0) {
      ind1 = direction[1];
      ind2 = direction[2];
      surface += 8;
    }
    else if (direction[1] == 0) {
      ind1 = direction[0];
      ind2 = direction[2];
      surface += 4;
    }
    else if (direction[2] == 0) {
      ind1 = direction[0];
      ind2 = direction[1];
    }
    ind1 = (ind1 + 1) / 2;
    ind2 = (ind2 + 1) / 2;
    surface += 2 * ind2 + ind1;
  }
  else if (num_crossings == 3) {
    surface += NUM_FACES + NUM_EDGES;
    int fwd[3];
    for (int i=0; i < 3; i++)
      fwd[i] = (direction[i] + 1) / 2;
    surface += 4 * fwd[0] + 2 * fwd[1] + fwd[2];
  }
  else {
    log_printf(ERROR, "Invalid number of surface crossings");
  }
  return surface;
}


/**
 * @brief Converts a surface into a 3 integer vector direction
 * @details The direction is a tuplet with each value taking either
 *          +1 (positive directed), 0 (neutral, or -1 (negative directed)
 * @param surface The surface of interest
 * @param direction The integer vector describing the direction
 */
void Cmfd::convertSurfaceToDirection(int surface, int* direction) {
  direction[0] = 0;
  direction[1] = 0;
  direction[2] = 0;
  if (surface < NUM_FACES) {
    int ind = surface % 3;
    int dir = 2 * (surface/3) - 1;
    direction[ind] = dir;
  }
  else if (surface < NUM_FACES + NUM_EDGES) {
    surface -= NUM_FACES;
    int group = surface / 4;
    int skipped = 2 - group;
    surface = surface % 4;
    int ind[2];
    ind[0] = surface % 2;
    ind[1] = (surface - ind[0]) / 2;
    int n = 0;
    for (int i=0; i < 3; i++) {
      if (i != skipped) {
        direction[i] = 2 * ind[n] - 1;
        n++;
      }
    }
  }
  else if (surface < NUM_SURFACES) {
    surface -= NUM_FACES + NUM_EDGES;
    direction[0] = 2 * (surface / 4) - 1;
    direction[1] = 2 * ((surface / 2) % 2) - 1;
    direction[2] = 2 * (surface % 2) - 1;
  }
  else {
    log_printf(ERROR, "Invalid surface ID %d", surface);
  }
}


/**
 * @brief Returns the surface name associated with the 3 integer vector
 *        direction
 * @details The direction is a tuplet with each value taking either
 *          +1 (positive directed), 0 (neutral, or -1 (negative directed)
 * @param direction The integer vector describing the direction
 * @return A string containing the surface name
 */
std::string Cmfd::getSurfaceNameFromDirection(int* direction) {
  std::string str = "SURFACE";
  std::string variables = "XYZ";
  for (int i=0; i < 3; i++) {
    if (direction[i] != 0) {
      str += "_";
      str += variables.at(i);
      if (direction[i] < 0)
        str += "_MIN";
      else
        str += "_MAX";
    }
  }
  return str;
}


/**
 * @brief Returns the surface name associated with a surface
 * @param surface The surface of interest
 * @return A string containing the surface name
 */
std::string Cmfd::getSurfaceNameFromSurface(int surface) {
  int direction[3];
  convertSurfaceToDirection(surface, direction);
  return getSurfaceNameFromDirection(direction);
}


/**
 * @brief A debugging tool that prints all prolongation factors to file
 */
void Cmfd::printProlongationFactors() {

  /* Loop over CMFD groups */
  for (int e = 0; e < _num_cmfd_groups; e++) {

    /* Create arrays for spatial data */
    CMFD_PRECISION log_ratios[_num_x * _num_y * _num_z];
    for (int i = 0; i < _num_x * _num_y * _num_z; i++)
      log_ratios[i] = 0.0;
    for (int i = 0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

      double old_flux = _old_flux->getValue(i, e);
      double new_flux = _new_flux->getValue(i, e);
      int cell_id = getGlobalCMFDCell(i);
      log_ratios[cell_id] = log(new_flux/old_flux);
    }

#ifdef MPIx
    if (_geometry->isDomainDecomposed()) {
      CMFD_PRECISION temp_log_ratios[_num_x * _num_y * _num_z];

      /* Select appropriate floating point size for transfer */
      MPI_Datatype mpi_precision;
      if (sizeof(CMFD_PRECISION) == 4)
        mpi_precision = MPI_FLOAT;
      else
        mpi_precision = MPI_DOUBLE;

      for (int i = 0; i < _num_x * _num_y * _num_z; i++)
        temp_log_ratios[i] = log_ratios[i];
      MPI_Allreduce(temp_log_ratios, log_ratios, _num_x * _num_y * _num_z,
                    mpi_precision, MPI_SUM, _geometry->getMPICart());
    }
#endif

    /* Print prolongation factors distribution to file */
    if (_geometry->isRootDomain()) {
      long long iter = _moc_iteration;
      long long group = e;
      std::string fname = "pf_group_";
      std::string group_num = std::to_string(group);
      std::string iter_num = std::to_string(iter);
      fname += group_num;
      fname += "_iter_";
      fname += iter_num;
      fname += ".txt";
      std::ofstream out(fname);
      out<< std::setprecision(5);

      out << "[NORMAL]  Spatial distribution of prolongation factors:"
          << std::endl;
      for (int z=0; z < _num_z; z++) {
        out << " -------- z = " << z << " ----------" << std::endl;
        for (int y=0; y < _num_y; y++) {
          for (int x=0; x < _num_x; x++) {
            int ind = (z * _num_y + y) * _num_x + x;
            out << log_ratios[ind] << " ";
          }
          out << std::endl;
        }
      }
      out.close();
    }
  }
}


/**
 * @brief This function tallies the current impinging on the domain from
 *        starting fluxes
 * @details Incoming currents are tallied for use in diagnostics, debugging,
 *          and adjusting sigma-t to enforce consistency with the MOC solution,
 *          if requested
 * @param point The point where the fluxes enter the geometry
 * @param delta_x The a small x-nudge in the direction of travel
 * @param delta_y The a small y-nudge in the direction of travel
 * @param delta_z The a small z-nudge in the direction of travel
 * @param track_flux The angular fluxes impinging on the domain
 * @param weight The weight of the Track
 */
void Cmfd::tallyStartingCurrent(Point* point, double delta_x, double delta_y,
                                double delta_z, float* track_flux,
                                double weight) {

  /* Check for non-zero current */
  bool non_zero = false;
  for (int e=0; e < _num_moc_groups; e++) {
    if (fabs(track_flux[e]) > 0) {
      non_zero = true;
      break;
    }
  }
  if (!non_zero)
    return;

  /* Create local coordinate */
  LocalCoords coords;
  coords.setUniverse(_geometry->getRootUniverse());
  coords.setX(point->getX());
  coords.setY(point->getY());
  coords.setZ(point->getZ());

  /* Find the CMFD cell */
  coords.adjustCoords(delta_x, delta_y, delta_z);
  int cell = findCmfdCell(&coords);
  coords.adjustCoords(-delta_x, -delta_y, -delta_z);

  /* Check the CMFD cell */
  if (cell == -1)
    log_printf(ERROR, "Failed to find starting CMFD cell for track start "
               "point");
  int cell_x = cell % _local_num_x;
  int cell_y = (cell % (_local_num_x * _local_num_y)) / _local_num_x;
  int cell_z = cell / (_local_num_x * _local_num_y);
  int bounds[3];
  bool singular[3] = {_local_num_x == 1, _local_num_y == 1, _local_num_z == 1};
  bounds[0] = -1 * (cell_x == 0) + (cell_x == _local_num_x-1);
  bounds[1] = -1 * (cell_y == 0) + (cell_y == _local_num_y-1);
  bounds[2] = -1 * (cell_z == 0) + (cell_z == _local_num_z-1);
  if ((bounds[0] == 0 && !singular[0]) && (bounds[1] == 0 && !singular[1]) &&
      (bounds[2] == 0 && !singular[2]))
    log_printf(ERROR, "Track start point not on a boundary CMFD cell. "
               "Cell = %d (%d, %d, %d) from Track: (%3.2f, %3.2f, %3.2f) "
               "adjusted (%3.2e, %3.2e, %3.2e)", cell, cell_x, cell_y, cell_z,
               point->getX(), point->getY(), point->getZ(), delta_x, delta_y,
               delta_z);

  CMFD_PRECISION currents[_num_cmfd_groups]
       __attribute__ ((aligned(VEC_ALIGNMENT)));
  memset(currents, 0, _num_cmfd_groups * sizeof(CMFD_PRECISION));

  /* Tally currents to each CMFD group locally */
  for (int e=0; e < _num_moc_groups; e++) {

    /* Get the CMFD group */
    int cmfd_group = getCmfdGroup(e);

    /* Increment the surface group */
    currents[cmfd_group] += track_flux[e] * weight;
  }

  /* Tally starting currents to cell */
  _starting_currents->incrementValues(cell, 0, _num_cmfd_groups - 1, currents);

}


/**
 * @brief Records net currents (leakage) on every CMFD cell for every group
 */
void Cmfd::recordNetCurrents() {

#pragma omp parallel for
  for (int i=0; i < _local_num_x * _local_num_y * _local_num_z; i++) {

    for (int e=0; e < _num_cmfd_groups; e++)
      _net_currents->incrementValue(i, e,
                                    -1 * _starting_currents->getValue(i,e));

    /* Compute cell indexes */
    int cell_ind[3];
    cell_ind[0] = i % _local_num_x;
    cell_ind[1] = (i / _local_num_x) % _local_num_y;
    cell_ind[2] = i / (_local_num_x * _local_num_y);

    /* Tally current from all surfaces including edges and corners */
    for (int s=0; s < NUM_SURFACES; s++) {

      /* Check if edge/corner exists */
      if (s >= NUM_FACES) {
        int idx = i * NUM_SURFACES * _num_cmfd_groups + s * _num_cmfd_groups;
        std::map<int, CMFD_PRECISION>::iterator it =
          _edge_corner_currents.find(idx);
        if (it == _edge_corner_currents.end())
          continue;
      }

      /* Compute index and vector direction */
      int direction[3];
      convertSurfaceToDirection(s, direction);

      /* Compute the next CMFD cell from the cell indexes and direction */
      int cmfd_cell_next = 0;
      int cell_next_ind[3];
      for (int d=0; d < 3; d++)
        cell_next_ind[d] = cell_ind[d] + direction[d];

      cmfd_cell_next = cell_next_ind[0] + cell_next_ind[1] * _local_num_x
                     + cell_next_ind[2] * (_local_num_x * _local_num_y);
      if (cell_next_ind[0] < 0 || cell_next_ind[0] >= _local_num_x ||
          cell_next_ind[1] < 0 || cell_next_ind[1] >= _local_num_y ||
          cell_next_ind[2] < 0 || cell_next_ind[2] >= _local_num_z)
        cmfd_cell_next = -1;

      /* Tally net currents */
      if (s < NUM_FACES) {
        int idx = s * _num_cmfd_groups;
        for (int e=0; e < _num_cmfd_groups; e++) {
          double current = 1 * _surface_currents->getValue(i, idx+e);
          _net_currents->incrementValue(i, e, current);
        }

        if (cmfd_cell_next != -1) {
          for (int e=0; e < _num_cmfd_groups; e++) {
            double current = -1 * _surface_currents->getValue(i, idx+e);
            _net_currents->incrementValue(cmfd_cell_next, e, current);
          }
        }
      }
      else {
        int idx = i * NUM_SURFACES * _num_cmfd_groups + s * _num_cmfd_groups;
        for (int e=0; e < _num_cmfd_groups; e++) {
          double current = _edge_corner_currents.at(idx+e);
          _net_currents->incrementValue(i, e, current);
        }

        if (cmfd_cell_next != -1) {
          for (int e=0; e < _num_cmfd_groups; e++) {
            double current = -1 * _edge_corner_currents.at(idx+e);
            _net_currents->incrementValue(cmfd_cell_next, e, current);
          }
        }
      }
    }
  }
}

/**
 * @brief Set width of non-uniform meshes in x y z directions.
 * @details An example of how this may be called from Python illustrated below:
 *
 * @code
 *          cmfd.setWidths([[1,2,3], [4,5,6,7], [3.3,2.4]])
 * @endcode
 *
 * @param widths A vector of 3 vectors for the x y z sizes of non-uniform meshes
 */
void Cmfd::setWidths(std::vector< std::vector<double> > widths) {

  if (widths.size() == 3)
    _cell_widths_z = widths[2];
  else if (widths.size() == 2)
    _cell_widths_z.push_back(1.0);
  else
    log_printf(ERROR, "CMFD lattice widths must have dimension 2 or 3.");

  _non_uniform = true;
  _cell_widths_x = widths[0];
  _cell_widths_y = widths[1];
}


/**
 * @brief Print information about CMFD cell dimensions.
 * @details For debug use.
 */
void Cmfd::printCmfdCellSizes() {
  int i;
  printf("non_uniform=%d, \nNum_XYZ: %2d, %2d, %2d\n", _non_uniform,
         _num_x, _num_y, _num_z);
  printf("Num_Local_XYZ: %2d, %2d, %2d\n", _local_num_x,
         _local_num_y, _local_num_z);
  printf("width_XYZ: %f, %f, %f\n", _width_x,_width_y,_width_z);
  printf("cell_width_XYZ: %f, %f, %f\n", _cell_width_x,
         _cell_width_y,_cell_width_z);
  printf("cell_widths_XYZ:\n");
  for (i=0; i<_num_x; i++)
    printf("i=%d, %f; ",i, _cell_widths_x[i]);
  printf("\n");
  for (i=0; i<_num_y; i++)
    printf("i=%d, %f; ",i, _cell_widths_y[i]);
  printf("\n");
  for (i=0; i<_num_z; i++)
    printf("i=%d, %f; ",i, _cell_widths_z[i]);
  printf("\n");

  printf("accumulates_XYZ:\n");
  for (i=0; i<_num_x+1; i++)
    printf("i=%d, %f; ",i, _accumulate_x[i]);
  printf("\n");
  for (i=0; i<_num_y+1; i++)
    printf("i=%d, %f; ",i, _accumulate_y[i]);
  printf("\n");
  for (i=0; i<_num_z+1; i++)
    printf("i=%d, %f; ",i, _accumulate_z[i]);
  printf("\n");
}


/**
 * @brief Create a string with information about the CMFD solver.
 * @details For pretty printing in Python API
 */
std::string Cmfd::toString() {

  std::stringstream message;
  message << "CMFD acceleration at " << (void*)this << std::endl;
  message << "Mesh in XYZ: [" << _num_x << ", " << _num_y << ", " << _num_z;
  message << "]" << std::endl;
  message << "Condensing " << _num_moc_groups << " MOC groups to " <<
             _num_cmfd_groups << " CMFD groups";
  return message.str();
}
