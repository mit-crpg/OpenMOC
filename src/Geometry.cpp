#include "Geometry.h"

/**
 * @brief Resets the auto-generated unique IDs for Materials, Surfaces,
 *        Cells and Universes/Lattices to 10000.
 */
void reset_auto_ids() {
  reset_material_id();
  reset_surface_id();
  reset_cell_id();
  reset_universe_id();
}


/**
 * @brief Constructor initializes an empty Geometry.
 */
Geometry::Geometry() {

  /* Initialize CMFD object to NULL */
  _cmfd = NULL;
  _overlaid_mesh = NULL;
  _root_universe = NULL;
  _domain_decomposed = false;
  _num_domains_x = 1;
  _num_domains_y = 1;
  _num_domains_z = 1;
  _num_modules_x = 1;
  _num_modules_y = 1;
  _num_modules_z = 1;
  _domain_index_x = 1;
  _domain_index_y = 1;
  _domain_index_z = 1;
  _symmetries.resize(3, false);
  _domain_FSRs_counted = false;
  _contains_FSR_centroids = false;
  _twiddle = false;
  _loaded_from_file = false;
}


/**
 * @brief Destructor clears FSR to Cells and Materials maps.
 */
Geometry::~Geometry() {

  /* Free all materials */
  if (_loaded_from_file) {
    std::map<int, Material*> materials = _root_universe->getAllMaterials();
    std::map<int, Material*>::iterator iter;
    for (iter = materials.begin(); iter != materials.end(); ++iter)
      delete iter->second;
  }

  /* Free all surfaces */
  if (_loaded_from_file) {
    std::map<int, Surface*> surfaces = getAllSurfaces();
    std::map<int, Surface*>::iterator iter_s;
    for (iter_s = surfaces.begin(); iter_s != surfaces.end(); ++iter_s) {
      delete iter_s->second;
    }
  }

  /* Free all cells and universes */
  if (_loaded_from_file) {
    std::map<int, Cell*> cells = getAllCells();
    std::map<int, Universe*> universes = getAllUniverses();
    std::map<int, Cell*>::iterator iter_c;
    std::map<int, Universe*>::iterator iter_u;
    for (iter_c = cells.begin(); iter_c != cells.end(); ++iter_c)
      delete iter_c->second;
    for (iter_u = universes.begin(); iter_u != universes.end(); ++iter_u)
      delete iter_u->second;
  }

  /* Free FSR maps if they were initialized */
  if (_FSR_keys_map.size() != 0) {

    fsr_data **values = _FSR_keys_map.values();
    for (long i=0; i<_FSR_keys_map.size(); i++)
      delete values[i];
    delete [] values;

    ExtrudedFSR **extruded_fsrs = _extruded_FSR_keys_map.values();
    for (int i=0; i<_extruded_FSR_keys_map.size(); i++) {
      if (extruded_fsrs[i]->_mesh != NULL)
        delete [] extruded_fsrs[i]->_mesh;
      delete [] extruded_fsrs[i]->_fsr_ids;
      delete [] extruded_fsrs[i]->_materials;
      delete extruded_fsrs[i]->_coords;
      delete extruded_fsrs[i];
    }
    delete [] extruded_fsrs;

    _FSR_keys_map.clear();
    _extruded_FSR_keys_map.clear();
    _FSRs_to_keys.clear();
    _FSRs_to_material_IDs.clear();
  }
  if (_overlaid_mesh != NULL)
    delete _overlaid_mesh;
}


/**
 * @brief Returns the total width in the x-direction of the Geometry in cm.
 * @return the total width of the Geometry in the x-direction (cm)
 */
double Geometry::getWidthX() {
  return (getMaxX() - getMinX());
}


/**
 * @brief Returns the total width in the y-direction of the Geometry in cm.
 * @return the total width of the Geometry in the y-direction (cm)
 */
double Geometry::getWidthY() {
  return (getMaxY() - getMinY());
}


/**
 * @brief Returns the total width in the z-direction of the Geometry in cm.
 * @return the total width of the Geometry in the z-direction (cm)
 */
double Geometry::getWidthZ() {
  return (getMaxZ() - getMinZ());
}


/**
 * @brief Return the minimum x-coordinate contained by the Geometry.
 * @return the minimum x-coordinate (cm)
 */
double Geometry::getMinX() {
  if (_domain_decomposed) {
    double geometry_min_x = _root_universe->getMinX();
    double geometry_max_x = _root_universe->getMaxX();
    double domain_width_x = (geometry_max_x - geometry_min_x) / _num_domains_x;
    return geometry_min_x + _domain_index_x * domain_width_x;
  }
  else {
    return _root_universe->getMinX();
  }
}


/**
 * @brief Return the maximum x-coordinate contained by the Geometry.
 * @return the maximum x-coordinate (cm)
 */
double Geometry::getMaxX() {
  if (_domain_decomposed) {
    double geometry_min_x = _root_universe->getMinX();
    double geometry_max_x = _root_universe->getMaxX();
    double domain_width_x = (geometry_max_x - geometry_min_x) / _num_domains_x;
    int reverse_index_x = _num_domains_x - _domain_index_x - 1;
    return geometry_max_x - reverse_index_x * domain_width_x;
  }
  else {
    return _root_universe->getMaxX();
  }
}


/**
 * @brief Return the minimum y-coordinate contained by the Geometry.
 * @return the minimum y-coordinate (cm)
 */
double Geometry::getMinY() {
  if (_domain_decomposed) {
    double geometry_min_y = _root_universe->getMinY();
    double geometry_max_y = _root_universe->getMaxY();
    double domain_width_y = (geometry_max_y - geometry_min_y) / _num_domains_y;
    return geometry_min_y + _domain_index_y * domain_width_y;
  }
  else {
    return _root_universe->getMinY();
  }
}


/**
 * @brief Return the maximum y-coordinate contained by the Geometry.
 * @return the maximum y-coordinate (cm)
 */
double Geometry::getMaxY() {
  if (_domain_decomposed) {
    double geometry_min_y = _root_universe->getMinY();
    double geometry_max_y = _root_universe->getMaxY();
    double domain_width_y = (geometry_max_y - geometry_min_y) / _num_domains_y;
    int reverse_index_y = _num_domains_y - _domain_index_y - 1;
    return geometry_max_y - reverse_index_y * domain_width_y;
  }
  else {
    return _root_universe->getMaxY();
  }
}


/**
 * @brief Return the minimum z-coordinate contained by the Geometry.
 * @return the minimum z-coordinate (cm)
 */
double Geometry::getMinZ() {
  if (_domain_decomposed) {
    double geometry_min_z = _root_universe->getMinZ();
    double geometry_max_z = _root_universe->getMaxZ();
    double domain_width_z = (geometry_max_z - geometry_min_z) / _num_domains_z;
    return geometry_min_z + _domain_index_z * domain_width_z;
  }
  else {
    return _root_universe->getMinZ();
  }
}


/**
 * @brief Return the maximum z-coordinate contained by the Geometry.
 * @return the maximum z-coordinate (cm)
 */
double Geometry::getMaxZ() {
  if (_domain_decomposed) {
    double geometry_min_z = _root_universe->getMinZ();
    double geometry_max_z = _root_universe->getMaxZ();
    double domain_width_z = (geometry_max_z - geometry_min_z) / _num_domains_z;
    int reverse_index_z = _num_domains_z - _domain_index_z - 1;
    return geometry_max_z - reverse_index_z * domain_width_z;
  }
  else {
    return _root_universe->getMaxZ();
  }
}


/**
 * @brief Returns the boundary conditions at the minimum x-coordinate in the
 *        Geometry / domain if the geometry is domain-decomposed.
 * @return the boundary conditions for the minimum x-coordinate in the domain
 */
boundaryType Geometry::getMinXBoundaryType() {
  if (_domain_decomposed && _domain_index_x > 0)
   return INTERFACE;
  else
    return _root_universe->getMinXBoundaryType();
}


/**
 * @brief Returns the boundary conditions at the maximum x-coordinate in the
 *        Geometry / domain if the geometry is domain-decomposed.
 * @return the boundary conditions for the maximum x-coordinate in the domain
 */
boundaryType Geometry::getMaxXBoundaryType() {
  if (_domain_decomposed && _domain_index_x < _num_domains_x-1)
    return INTERFACE;
  else
    return _root_universe->getMaxXBoundaryType();
}


/**
 * @brief Returns the boundary conditions at the minimum y-coordinate in the
 *        Geometry / domain if the geometry is domain-decomposed.
 * @return the boundary conditions for the minimum y-coordinate in the domain
 */
boundaryType Geometry::getMinYBoundaryType() {
  if (_domain_decomposed && _domain_index_y > 0)
    return INTERFACE;
  else
    return _root_universe->getMinYBoundaryType();
}


/**
 * @brief Returns the boundary conditions at the maximum y-coordinate in the
 *        Geometry / domain if the geometry is domain-decomposed.
 * @return the boundary conditions for the maximum y-coordinate in the domain
 */
boundaryType Geometry::getMaxYBoundaryType() {
  if (_domain_decomposed && _domain_index_y < _num_domains_y-1)
    return INTERFACE;
  else
    return _root_universe->getMaxYBoundaryType();
}


/**
 * @brief Returns the boundary conditions at the minimum z-coordinate in the
 *        Geometry / domain if the geometry is domain-decomposed.
 * @return the boundary conditions for the minimum z-coordinate in the domain
 */
boundaryType Geometry::getMinZBoundaryType() {
  if (_domain_decomposed && _domain_index_z > 0)
    return INTERFACE;
  else
    return _root_universe->getMinZBoundaryType();
}


/**
 * @brief Returns the boundary conditions at the maximum z-coordinate in the
 *        Geometry / domain if the geometry is domain-decomposed.
 * @return the boundary conditions for the maximum z-coordinate in the domain
 */
boundaryType Geometry::getMaxZBoundaryType() {
  if (_domain_decomposed && _domain_index_z < _num_domains_z-1)
    return INTERFACE;
  else
    return _root_universe->getMaxZBoundaryType();
}


/**
 * @brief Returns the number of source regions in the Geometry domain.
 * @return number of FSRs
 */
long Geometry::getNumFSRs() {
  return _FSRs_to_keys.size();
}


/**
 * @brief Returns the number of source regions in the entire Geometry.
 * @return number of FSRs
 */
long Geometry::getNumTotalFSRs() {
  long domain_fsrs =  _FSRs_to_keys.size();
  long total_fsrs = domain_fsrs;
#ifdef MPIx
  if (_domain_decomposed)
    MPI_Allreduce(&domain_fsrs, &total_fsrs, 1, MPI_LONG, MPI_SUM, _MPI_cart);
#endif
  return total_fsrs;
}


/**
 * @brief Returns the number of energy groups for each Material's nuclear data.
 * @return the number of energy groups
 */
int Geometry::getNumEnergyGroups() {

  std::map<int, Material*> materials = getAllMaterials();

  if (materials.size() == 0)
    log_printf(ERROR, "Unable to return the number of energy groups from "
               "the Geometry since it does not contain any Materials");

  int num_groups = materials.begin()->second->getNumEnergyGroups();
  std::map<int, Material*>::iterator iter;

  for (iter = materials.begin(); iter != materials.end(); ++iter) {
    if (iter->second->getNumEnergyGroups() != num_groups)
      log_printf(ERROR, "Unable to return the number of energy groups from "
                 "the Geometry since it contains different numbers of groups: "
                 "%d and %d", num_groups, iter->second->getNumEnergyGroups());
  }

  return num_groups;
}


/**
 * @brief Returns the number of Materials in the Geometry.
 * @return the number of Materials
 */
int Geometry::getNumMaterials() {

  std::map<int, Material*> all_materials;

  if (_all_materials.size() == 0)
    all_materials = getAllMaterials();
  else
    all_materials = _all_materials;

  int num_materials = all_materials.size();
  return num_materials;
}


/**
 * @brief Returns the number of Cells in the Geometry.
 * @return the number of Cells
 */
int Geometry::getNumCells() {

  int num_cells = 0;

  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = _root_universe->getAllCells();
    num_cells = all_cells.size();
  }

  return num_cells;
}


/**
 * @brief Return a std::map container of Surface IDs (keys) with Surfaces
 *        pointers (values).
 * @return a std::map of Surfaces indexed by Surface ID in the geometry
 */
std::map<int, Surface*> Geometry::getAllSurfaces() {

  Cell* cell;
  Surface* surf;
  std::map<int, Surface*> all_surfs;
  std::map<int, Halfspace*> surfs;
  std::map<int, Cell*>::iterator c_iter;
  std::map<int, Halfspace*>::iterator s_iter;

  /* Gather all surfaces from all cells */
  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = getAllCells();

    for (c_iter = all_cells.begin(); c_iter != all_cells.end(); ++c_iter) {
      cell = (*c_iter).second;
      surfs = cell->getSurfaces();

      for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
        surf = (*s_iter).second->_surface;
        all_surfs[surf->getId()] = surf;
      }
    }
  }

  return all_surfs;
}


/**
 * @brief Return a std::map container of Material IDs (keys) with Materials
 *        pointers (values).
 * @return a std::map of Materials indexed by Material ID in the geometry
 */
std::map<int, Material*> Geometry::getAllMaterials() {

  std::map<int, Material*> all_materials;
  Cell* cell;
  Material* material;

  /* Gather all materials from all cells */
  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = _root_universe->getAllCells();
    std::map<int, Cell*>::iterator iter;

    for (iter = all_cells.begin(); iter != all_cells.end(); ++iter) {
      cell = (*iter).second;

      if (cell->getType() == MATERIAL) {
        material = cell->getFillMaterial();

        if (material != NULL)
          all_materials[material->getId()] = material;
      }
    }
  }

  return all_materials;
}


/**
 * @brief Modify scattering and total cross sections to study MOC stability.
 */
void Geometry::manipulateXS() {

  std::map<int, Material*> all_materials = getAllMaterials();
  std::map<int, Material*>::iterator mat_iter;

  for (mat_iter = all_materials.begin(); mat_iter != all_materials.end();
       ++mat_iter) {

      Material* mat = mat_iter->second;
      int ng = mat->getNumEnergyGroups();
      for (int g=0; g < ng; g++) {
        double sigma_t = mat->getSigmaTByGroup(g+1);

        /* Set all negative scatter cross sections to 0 */
        for (int gp=0; gp < ng; gp++)
          if (mat->getSigmaSByGroup(g+1, gp+1) < 0)
            mat->setSigmaSByGroup(0.0, g+1, gp+1);

        /* Gather all non-negative scattering terms */
        double sigma_s = 0.0;
        for (int gp=0; gp < ng; gp++)
          sigma_s += mat->getSigmaSByGroup(g+1, gp+1);

        /* Fudge */
        sigma_s *= 1.1;
        if (sigma_s > sigma_t)
          mat->setSigmaTByGroup(sigma_s, g+1);
      }
  }
}


/**
 * @brief Loads an array of SPH factors into the geometry's domains.
 * @details This method is called by compute_sph_factors in the 'materialize'
 *          module but may also be called by the user in Python if needed:
 *
 * @code
 *          geometry.loadSPHFactors(sph_factors.flatten(),
 *                                  double(sph_to_domains), "cell")
 * @endcode
 * @param sph_factors 1D array of SPH factors (group dependence in inner loop)
 * @param num_domains_groups number of domains times number of groups
 * @param sph_to_domain_ids map to link sph_factors array into domains, must be
          doubles //FIXME write a typemap to allow ints
 * @param num_domains total number of domains (may be larger than the number of
          domains in the OpenMOC geometry)
 * @param domain_type type of domain (material or cell) containing the cross
 *        sections that the SPH factors apply to
 */
void Geometry::loadSPHFactors(double* sph_factors, int num_domains_groups,
                              double* sph_to_domain_ids, int num_sph_domains,
                              const char* domain_type) {

  /* Check type of domain */
  if (strcmp(domain_type, "material") != 0 && strcmp(domain_type, "cell") != 0)
    log_printf(ERROR, "Domain type %s is not supported for loading SPH factor",
               &domain_type[0]);

  int num_groups = getNumEnergyGroups();

  /* If by material, loop on materials */
  if (strcmp(domain_type, "material") == 0) {
    std::map<int, Material*> all_materials = getAllMaterials();
    std::map<int, Material*>::iterator iter;

    for (int i=0; i<num_sph_domains; i++) {

      /* Find material */
      iter = all_materials.find(int(sph_to_domain_ids[i]));
      if (iter == all_materials.end()) {
        log_printf(WARNING, "SPH material %d is not in geometry",
                   int(sph_to_domain_ids[i]));
        continue;
      }
      Material* mat = iter->second;

      /* Use sph factors */
      double* sph = &sph_factors[i*num_groups];

      for (int g=1; g<=num_groups; g++) {

        mat->setSigmaTByGroup(mat->getSigmaTByGroup(g) * sph[g-1], g);
        mat->setSigmaAByGroup(mat->getSigmaAByGroup(g) * sph[g-1], g);
        mat->setSigmaFByGroup(mat->getSigmaFByGroup(g) * sph[g-1], g);
        mat->setNuSigmaFByGroup(mat->getNuSigmaFByGroup(g) * sph[g-1], g);
        for (int g2=1; g2<=num_groups; g2++) {
          mat->setSigmaSByGroup(mat->getSigmaSByGroup(g, g2) * sph[g-1], g, g2);
        }
      }
    }
  }
  /* SPH factors by cells */
  if (strcmp(domain_type, "cell") == 0) {
    std::map<int, Cell*> all_cells = getAllMaterialCells();
    std::map<int, Cell*>::iterator iter;

    for (int i=0; i<num_sph_domains; i++) {

      /* Find cell then material to use SPH factor on */
      iter = all_cells.find(int(sph_to_domain_ids[i]));
      if (iter == all_cells.end()) {
        log_printf(WARNING, "SPH cell %d is not in geometry",
                   int(sph_to_domain_ids[i]));
        continue;
      }
      Material* mat = iter->second->getFillMaterial();

      /* Use sph factors */
      double* sph = &sph_factors[i*num_groups];

      for (int g=1; g<=num_groups; g++) {

        mat->setSigmaTByGroup(mat->getSigmaTByGroup(g) * sph[g-1], g);
        mat->setSigmaAByGroup(mat->getSigmaAByGroup(g) * sph[g-1], g);
        mat->setSigmaFByGroup(mat->getSigmaFByGroup(g) * sph[g-1], g);
        mat->setNuSigmaFByGroup(mat->getNuSigmaFByGroup(g) * sph[g-1], g);
        for (int g2=1; g2<=num_groups; g2++) {
          mat->setSigmaSByGroup(mat->getSigmaSByGroup(g, g2) * sph[g-1], g, g2);
        }
      }
    }
  }
}


/**
 * @brief Return a std::map container of Cell IDs (keys) with Cells
 *        pointers (values).
 * @return a std::map of Cells indexed by Cell ID in the geometry
 */
std::map<int, Cell*> Geometry::getAllCells() {

  std::map<int, Cell*> all_cells;

  if (_root_universe != NULL)
    all_cells = _root_universe->getAllCells();

  return all_cells;
}


/**
 * @brief Return a std::map container of Cell IDs (keys) with Cells
 *        pointers (values).
 * @return a std::map of Cells indexed by Cell ID in the geometry
 */
std::map<int, Cell*> Geometry::getAllMaterialCells() {

  std::map<int, Cell*> all_material_cells;
  Cell* cell;

  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = _root_universe->getAllCells();
    std::map<int, Cell*>::iterator iter;

    for (iter = all_cells.begin(); iter != all_cells.end(); ++iter) {
      cell = (*iter).second;

      if (cell->getType() == MATERIAL)
        all_material_cells[cell->getId()] = cell;
    }
  }

  return all_material_cells;
}


/**
 * @brief Return a std::map container of Universe IDs (keys) with Universes
 *        pointers (values).
 * @return a std::map of Universes indexed by Universe ID in the geometry
 */
std::map<int, Universe*> Geometry::getAllUniverses() {

  std::map<int, Universe*> all_universes;

  if (_root_universe != NULL)
    all_universes = _root_universe->getAllUniverses();

  return all_universes;
}


/**
 * @brief Returns the Universe at the root node in the CSG tree.
 * @return the root Universe
 */
Universe* Geometry::getRootUniverse() {
  return _root_universe;
}


/**
 * @brief Returns a pointer to the CMFD object.
 * @return A pointer to the CMFD object
 */
Cmfd* Geometry::getCmfd() {
  return _cmfd;
}


/**
 * @brief Returns whether the Geometry is domain decomposed
 * @return If the domain is decomposed (true) or not (false)
 */
bool Geometry::isDomainDecomposed() {
  return _domain_decomposed;
}


/**
 * @brief Returns whether this MPI domain is the root domain
 * @return If this domain is the root domain (true) or not (false)
 */
bool Geometry::isRootDomain() {

  bool first_domain = true;
  if (_domain_decomposed)
    if (_domain_index_x != 0 || _domain_index_y != 0 || _domain_index_z != 0)
      first_domain = false;
  return first_domain;
}


/**
 * @brief Returns the domain indexes of the current domain
 * @param indexes A pointer to the array to be filled with the indexes
 */
void Geometry::getDomainIndexes(int* indexes) {
  indexes[0] = _domain_index_x;
  indexes[1] = _domain_index_y;
  indexes[2] = _domain_index_z;
}


/**
 * @brief Returns the number of domains in each direction.
 * @param structure A pointer to the array to be filled with the domain numbers
 */
void Geometry::getDomainStructure(int* structure) {
  structure[0] = _num_domains_x;
  structure[1] = _num_domains_y;
  structure[2] = _num_domains_z;
}


/**
 * @brief Sets the root Universe for the CSG tree.
 * @param root_universe the root Universe of the CSG tree.
 */
void Geometry::setRootUniverse(Universe* root_universe) {
  _root_universe = root_universe;
}


//FIXME: add setDefaultDomainDecomposition() function


/**
 * @brief Sets how many modular track laydown domains are in each MPI domain
 * @param num_x The number of modular domains in the x-direction per MPI domain
 * @param num_y The number of modular domains in the y-direction per MPI domain
 * @param num_z The number of modular domains in the z-direction per MPI domain
 */
void Geometry::setNumDomainModules(int num_x, int num_y, int num_z) {
  _num_modules_x = num_x;
  _num_modules_y = num_y;
  _num_modules_z = num_z;

  log_printf(NORMAL, "Using a [%d %d %d] inner domain structure for track "
             "laydown.", num_x, num_y, num_z);
}


/**
 * @brief Get the number of modular domains in the x-direction per MPI domain
 * @return _num_modules_x number of modular domains in the x-direction in domain
 */
int Geometry::getNumXModules() {
  return _num_modules_x;
}


/**
 * @brief Get the number of modular domains in the y-direction per MPI domain
 * @return _num_modules_y number of modular domains in the y-direction in domain
 */
int Geometry::getNumYModules() {
  return _num_modules_y;
}


/**
 * @brief Get the number of modular domains in the z-direction per MPI domain
 * @return _num_modules_z number of modular domains in the z-direction in domain
 */
int Geometry::getNumZModules() {
  return _num_modules_z;
}


/**
 * @brief Take into account domain symmetries to reduce the problem domain.
 * @param X_symmetry whether the domain is symmetric in X
 * @param Y_symmetry whether the domain is symmetric in Y
 * @param Z_symmetry whether the domain is symmetric in Z
 */
void Geometry::useSymmetry(bool X_symmetry, bool Y_symmetry, bool Z_symmetry) {

  if (_root_universe->getCells().size() > 1)
    log_printf(ERROR, "To take advantage of the problem symmetries, use a root"
               " universe AND a root cell to contain the CSG.");

#ifdef ONLYVACUUMBC
  if (X_symmetry || Y_symmetry || Z_symmetry)
    log_printf(ERROR, "Using symmetries requires reflective boundary conditions"
               ", re-compile without the ONLYVACUUMBC flag.");
#endif

  // Keep track of symmetries used
  _symmetries[0] = X_symmetry;
  _symmetries[1] = Y_symmetry;
  _symmetries[2] = Z_symmetry;

  if (X_symmetry) {
    // Get center plane
    double mid_x = (_root_universe->getMaxX() + _root_universe->getMinX()) / 2;
    XPlane* symX = new XPlane(mid_x);
    symX->setBoundaryType(REFLECTIVE);

    // Add plane to root cell
    Cell* root_cell = _root_universe->getCells().begin()->second;
    root_cell->addSurface(+1, symX);
    log_printf(NORMAL, "Using X symmetry to restrict domain to [%.3f %.3f] cm",
                mid_x, _root_universe->getMaxX());
  }

  if (Y_symmetry) {
    // Get center plane
    double mid_y = (_root_universe->getMaxY() + _root_universe->getMinY()) / 2;
    YPlane* symY = new YPlane(mid_y);
    symY->setBoundaryType(REFLECTIVE);

    // Add plane to root cell
    Cell* root_cell = _root_universe->getCells().begin()->second;
    root_cell->addSurface(+1, symY);
    log_printf(NORMAL, "Using Y symmetry to restrict domain to [%.3f %.3f] cm",
               mid_y, _root_universe->getMaxY());
  }

  if (Z_symmetry) {
    // Get center planes
    double mid_z = (_root_universe->getMaxZ() + _root_universe->getMinZ()) / 2;
    ZPlane* symZ = new ZPlane(mid_z);
    symZ->setBoundaryType(REFLECTIVE);

    // Add plane to root cell
    Cell* root_cell = _root_universe->getCells().begin()->second;
    root_cell->addSurface(+1, symZ);
    log_printf(NORMAL, "Using Z symmetry to restrict domain to [%.3f %.3f] cm",
               mid_z, _root_universe->getMaxZ());
  }

  // Reset boundaries to trigger boundary calculation again
  _root_universe->resetBoundaries();
}


/**
 * @brief Get the symmetries used to restrict the domain
 * @return a boolean indicating if the symmetry along this axis is used
 */
bool Geometry::getSymmetry(int axis) {
  return _symmetries[axis];
}


#ifdef MPIx
/**
 * @brief Domain decomposes the Geometry with MPI and modular ray tracing
 * @param nx The number of MPI domains in the x-direction
 * @param ny The number of MPI domains in the y-direction
 * @param nz The number of MPI domains in the z-direction
 * @param comm The MPI communicator to be used to communicate between domains
 */
void Geometry::setDomainDecomposition(int nx, int ny, int nz, MPI_Comm comm) {

  /* Calculate number of domains and get the number of MPI ranks */
  int num_domains = nx*ny*nz;
  int num_ranks;
  MPI_Comm_size(comm, &num_ranks);

  /* Check that the number of domains equals the number of ranks */
  log_set_ranks(comm);
  if (num_ranks != num_domains)
    log_printf(ERROR, "Number of ranks is %d and number of domains is %d",
               num_ranks, num_domains);


  /* Check that the root universe has been set */
  if (_root_universe == NULL)
    log_printf(ERROR, "The root universe must be set before domain "
                      "decomposition.");

  /* Check that the Geometry needs to be decomposed */
  if (num_domains > 1) {

    /* Make note of the domain decomposition */
    _domain_decomposed = true;
    _num_domains_x = nx;
    _num_domains_y = ny;
    _num_domains_z = nz;

    /* Create the MPI Communicator */
    int dims[3] = {nx, ny, nz};
    int wrap[3] = {false, false, false};
    int ret = MPI_Cart_create(comm, 3, dims, wrap, true, &_MPI_cart);
    log_set_ranks(_MPI_cart);

    /* Determine the domain indexes */
    int rank;
    int cart_coords[3];
    MPI_Comm_rank(_MPI_cart, &rank);
    MPI_Cart_coords(_MPI_cart, rank, 3, cart_coords);
    _domain_index_x = cart_coords[0];
    _domain_index_y = cart_coords[1];
    _domain_index_z = cart_coords[2];

    /* Set the bounds in a 1x1x1 Lattice object */
    _domain_bounds = new Lattice();
    _domain_bounds->setNumX(1);
    _domain_bounds->setNumY(1);
    _domain_bounds->setNumZ(1);
    double width_x = getWidthX();
    double width_y = getWidthY();
    double width_z = getWidthZ();
    _domain_bounds->setWidth(width_x, width_y, width_z);
    double offset_x = width_x / 2.0 + getMinX();
    double offset_y = width_y / 2.0 + getMinY();
    double offset_z = width_z / 2.0 + getMinZ();
    _domain_bounds->setOffset(offset_x, offset_y, offset_z);
    _domain_bounds->computeSizes();
    log_printf(NORMAL, "Successfully set %d x %d x %d domain decomposition",
                        nx, ny, nz);
  }
}


/**
 * @brief Returns the MPI communicator to communicate between MPI domains
 * @return The MPI communicator
 */
MPI_Comm Geometry::getMPICart() {
  if (!_domain_decomposed)
    log_printf(ERROR, "Tried to get MPI Cart but domain decomposition has not "
               "yet been set");
  return _MPI_cart;
}


/**
 * @brief Returns the rank of the specified neighboring domain.
 * @param offset_x The shift in the x-direction (number of domains)
 * @param offset_y The shift in the y-direction (number of domains)
 * @param offset_z The shift in the z-direction (number of domains)
 * @return The rank of the neighboring domain
 */
int Geometry::getNeighborDomain(int offset_x, int offset_y, int offset_z) {

  int neighbor_rank = -1;
  int neighbor_coords[3];
  neighbor_coords[0] = offset_x + _domain_index_x;
  neighbor_coords[1] = offset_y + _domain_index_y;
  neighbor_coords[2] = offset_z + _domain_index_z;
  if (_domain_decomposed)
    if (neighbor_coords[0] >= 0 && neighbor_coords[0] < _num_domains_x &&
        neighbor_coords[1] >= 0 && neighbor_coords[1] < _num_domains_y &&
        neighbor_coords[2] >= 0 && neighbor_coords[2] < _num_domains_z)
      MPI_Cart_rank(_MPI_cart, neighbor_coords, &neighbor_rank);

  return neighbor_rank;
}
#endif


/**
 * @brief Sets the pointer to a CMFD object used for acceleration.
 * @param cmfd a pointer to the CMFD object
 */
void Geometry::setCmfd(Cmfd* cmfd) {
  _cmfd = cmfd;
}


/**
 * @brief Sets a global overlaid mesh with the given mesh height
 * @details The global overlaid mesh is overlaid across the entire Geometry
 * @param axial_mesh_height The desired height of axial mesh cells
 * @param num_x number of divisions in the X direction
 * @param num_y number of divisions in the Y direction
 * @param num_radial_domains number of radial domains
 * @param radial_domains array with the indexes of each domain in X and Y
 */
void Geometry::setOverlaidMesh(double axial_mesh_height, int num_x, int num_y,
                               int num_radial_domains, int* radial_domains) {

  /* Get the global Geometry boundaries */
  double min_x = _root_universe->getMinX();
  double max_x = _root_universe->getMaxX();
  double min_y = _root_universe->getMinY();
  double max_y = _root_universe->getMaxY();
  double min_z = _root_universe->getMinZ();
  double max_z = _root_universe->getMaxZ();

  /* Create the lattice */
  _overlaid_mesh = new Lattice();
  int real_num_x = 1;
  int real_num_y = 1;
  if (num_x > 0 && num_y > 0) {
    if (radial_domains == NULL) {
      real_num_x = num_x;
      real_num_y = num_y;
    }
    else {
      for (int i=0; i < num_radial_domains; i++) {
        if (radial_domains[2*i] == _domain_index_x &&
            radial_domains[2*i+1] == _domain_index_y) {
          real_num_x = num_x;
          real_num_y = num_y;
        }
      }
    }
  }
  num_x = real_num_x;
  num_y = real_num_y;
  _overlaid_mesh->setNumX(num_x);
  _overlaid_mesh->setNumY(num_y);

  /* Determine actual axial mesh spacing from desired spacing */
  double total_width_z = max_z - min_z;
  int num_cells_z = total_width_z / axial_mesh_height;
  axial_mesh_height = total_width_z / num_cells_z;
  double mesh_width_x = (max_x - min_x) / num_x;
  double mesh_width_y = (max_y - min_y) / num_y;
  _overlaid_mesh->setNumZ(num_cells_z);
  _overlaid_mesh->setWidth(mesh_width_x, mesh_width_y, axial_mesh_height);

  /* Create the center point */
  Point offset;
  offset.setX(min_x + (max_x - min_x)/2.0);
  offset.setY(min_y + (max_y - min_y)/2.0);
  offset.setZ(min_z + (max_z - min_z)/2.0);
  _overlaid_mesh->setOffset(offset.getX(), offset.getY(), offset.getZ());

  _overlaid_mesh->computeSizes();

  log_printf(NORMAL, "Set global axial mesh of width %6.4f cm",
             axial_mesh_height);
}


/**
 * @brief Clears all boundary conditions from the Geometry
 */
void Geometry::clearBoundaries() {

  /* Extract all surfaces from the Geometry */
  std::map<int, Surface*> all_surfaces = getAllSurfaces();
  std::map<int, Surface*>::iterator it;

  /* Iterate over all surfaces */
  for (it = all_surfaces.begin(); it != all_surfaces.end(); ++it) {

    /* Remove boundary conditions */
    Surface* surface = it->second;
    if (surface->getBoundaryType() != BOUNDARY_NONE)
      surface->setBoundaryType(BOUNDARY_NONE);
  }
}


/**
 * @brief Find the Cell that this LocalCoords object is in at the lowest level
 *        of the nested Universe hierarchy.
 * @details This method assumes that the LocalCoords has been initialized
 *          with coordinates and a Universe ID. The method will recursively
 *          find the Cell on the lowest level of the nested Universe hierarchy
 *          by building a linked list of LocalCoords from the LocalCoord
 *          passed in as an argument down to the lowest level Cell found. In
 *          the process it will set the coordinates at each level of the
 *          hierarchy for each LocalCoord in the linked list for the Lattice
 *          or Universe that it is in. If the LocalCoords is outside the bounds
 *          of the Geometry or on the boundaries this method will return NULL;
 *          otherwise it will return a pointer to the Cell that is found by the
 *          recursive Geometry::findCell(...) method.
 * @param coords pointer to a LocalCoords object
 * @return returns a pointer to a Cell if found, NULL if no Cell found
 */
Cell* Geometry::findCellContainingCoords(LocalCoords* coords) {

  Universe* univ = coords->getUniverse();
  Cell* cell;

  /* Check if the coords are inside the geometry bounds */
  if (univ->getId() == _root_universe->getId()) {
    if (!withinBounds(coords))
      return NULL;
  }

  if (univ->getType() == SIMPLE)
    cell = univ->findCell(coords);
  else
    cell = static_cast<Lattice*>(univ)->findCell(coords);

  return cell;
}


/**
 * @brief Find the first Cell of a Track segment with a starting Point that is
 *        represented by the LocalCoords method parameter.
 * @details This method assumes that the LocalCoords has been initialized
 *          with coordinates and a Universe ID. This method will move the
 *          initial starting point by a small amount along the direction of
 *          the Track in order to ensure that the track starts inside of a
 *          distinct FSR rather than on the boundary between two of them.
 *          The method will recursively find the LocalCoords by building a
 *          linked list of LocalCoords from the LocalCoords passed in as an
 *          argument down to the Cell found in the lowest level of the nested
 *          Universe hierarchy. In the process, the method will set the
 *          coordinates at each level in the nested Universe hierarchy for
 *          each LocalCoord in the linked list for the Lattice or Universe
 *          that it is in.
 * @param coords pointer to a LocalCoords object
 * @param azim the azimuthal angle for a trajectory projected from the LocalCoords
 * @param polar the polar angle for a trajectory projected from the LocalCoords
 * @return returns a pointer to a cell if found, NULL if no cell found
*/
Cell* Geometry::findFirstCell(LocalCoords* coords, double azim, double polar) {

  double delta_x = cos(azim) * sin(polar) * TINY_MOVE;
  double delta_y = sin(azim) * sin(polar) * TINY_MOVE;
  double delta_z = cos(polar) * TINY_MOVE;
  coords->adjustCoords(delta_x, delta_y, delta_z);
  return findCellContainingCoords(coords);
}


/**
 * @brief Find the Material for a flat source region ID.
 * @details  This method finds the fsr_id within the
 *           _FSR_to_material_IDs map and returns the corresponding
 *           pointer to the Material object.
 * @param fsr_id a FSR id
 * @return a pointer to the Material that this FSR is in
 */
Material* Geometry::findFSRMaterial(long fsr_id) {

  std::map<int, Material*> all_materials;

  if (_all_materials.size() == 0)
    all_materials = getAllMaterials();
  else
    all_materials = _all_materials;

  int mat_id = _FSRs_to_material_IDs.at(fsr_id);
  if (all_materials.find(mat_id) == all_materials.end())
      log_printf(ERROR, "Failed to find FSR Material for FSR %ld with "
                 "Material ID %d", fsr_id, mat_id);

  return all_materials[mat_id];
}


/**
 * @brief Finds the next Cell for a LocalCoords object along a trajectory
 *        defined by two angles : azimuthal from 0 to 2Pi, polar from 0 to Pi.
 * @details The method will update the LocalCoords passed in as an argument
 *          to be the one at the boundary of the next Cell crossed along the
 *          given trajectory. It will do this by finding the minimum distance
 *          to the surfaces at all levels of the coords hierarchy.
 *          If the LocalCoords is outside the bounds of the Geometry or on
 *          the boundaries this method will return NULL; otherwise it will
 *          return a pointer to the Cell that the LocalCoords will reach
 *          next along its trajectory.
 * @param coords pointer to a LocalCoords object
 * @param azim the azimuthal angle of the trajectory
 * @param polar the polar angle of the trajectory
 * @return a pointer to a Cell if found, NULL if no Cell found
 */
Cell* Geometry::findNextCell(LocalCoords* coords, double azim, double polar) {

  Cell* cell = NULL;
  double dist;
  double min_dist = std::numeric_limits<double>::infinity();

  /* Save angles in case of a rotated cell */
  double old_azim = azim;
  double old_polar = polar;

  /* Get highest level coords */
  coords = coords->getHighestLevel();

  /* If the current coords is outside the domain, return NULL */
  if (_domain_decomposed) {
    Point* point = coords->getPoint();
    if (!_domain_bounds->containsPoint(point))
      return NULL;
  }

  /* If the current coords is not in any Cell, return NULL */
  if (coords->getLowestLevel()->getCell() == NULL)
    return NULL;

  /* If the current coords is inside a Cell, look for next Cell */
  else {

    /* Descend universes/coord until at the lowest level.
     * At each universe/lattice level get distance to next
     * universe or lattice cell. Compare to get new min_dist. */
    while (coords != NULL) {

      /* If we reach a LocalCoord in a Lattice, find the distance to the
       * nearest lattice cell boundary */
      if (coords->getType() == LAT) {
        Lattice* lattice = coords->getLattice();
        dist = lattice->minSurfaceDist(coords->getPoint(), azim, polar);
      }
      /* If we reach a LocalCoord in a Universe, find the distance to the
       * nearest cell surface */
      else {
        Cell* cell = coords->getCell();
        dist = cell->minSurfaceDist(coords->getPoint(), azim, polar);

        /* Apply rotation to direction. Position has already been modified by
           universe->findCell() in findCellContainingCoords() */
        if (cell->isRotated()) {

          double* matrix = cell->getRotationMatrix();
          double uvw[3] = {sin(polar)*cos(azim), sin(polar)*sin(azim),
                           cos(polar)};
          double rot_uvw[3] = {matrix[0]*uvw[0] + matrix[1]*uvw[1] +
                               matrix[2]*uvw[2], matrix[3]*uvw[0] +
                               matrix[4]*uvw[1] + matrix[5]*uvw[2],
                               matrix[6]*uvw[0] + matrix[7]*uvw[1] +
                               matrix[8]*uvw[2]};
          polar = acos(rot_uvw[2]);
          azim = atan2(rot_uvw[1], rot_uvw[0]);
          if (azim < 0)
            azim += 2 * M_PI;
        }
      }

      /* Recheck min distance */
      min_dist = std::min(dist, min_dist);

      /* Descend one level, exit if at the lowest level of coordinates */
      if (coords->getNext() == NULL)
        break;
      else
        coords = coords->getNext();
    }

    /* Reset coords direction in case there was a rotated cell */
    coords = coords->getHighestLevel();
    coords->prune();
    azim = old_azim;
    polar = old_polar;

    /* Check for distance to an overlaid mesh */
    if (_overlaid_mesh != NULL) {
      dist = _overlaid_mesh->minSurfaceDist(coords->getPoint(), azim, polar);
      min_dist = std::min(dist, min_dist);
    }

    /* Check for distance to nearest CMFD mesh cell boundary */
    if (_cmfd != NULL) {
      Lattice* lattice = _cmfd->getLattice();
      dist = lattice->minSurfaceDist(coords->getPoint(), azim, polar);
      min_dist = std::min(dist, min_dist);
    }

    /* Check for distance to nearest domain boundary */
    bool domain_boundary = false;
    if (_domain_decomposed) {
      dist = _domain_bounds->minSurfaceDist(coords->getPoint(), azim, polar);
      if (dist - min_dist < ON_SURFACE_THRESH) {
        min_dist = dist;
        domain_boundary = true;
      }
    }

    /* Move point and get next cell */
    double delta_x = cos(azim) * sin(polar) * (min_dist + TINY_MOVE);
    double delta_y = sin(azim) * sin(polar) * (min_dist + TINY_MOVE);
    double delta_z = cos(polar) * (min_dist + TINY_MOVE);
    coords->adjustCoords(delta_x, delta_y, delta_z);

    if (domain_boundary)
      return NULL;
    else
      return findCellContainingCoords(coords);
  }
}


/**
 * @brief Find and return the ID of the flat source region that a given
 *        LocalCoords object resides within.
 * @param coords a LocalCoords object pointer
 * @return the FSR ID for a given LocalCoords object
 */
long Geometry::findFSRId(LocalCoords* coords) {

  long fsr_id;
  LocalCoords* curr = coords;
  curr = coords->getLowestLevel();

  /* Generate unique FSR key */
  int thread_id = omp_get_thread_num();
  std::string& fsr_key = _fsr_keys[thread_id];
  getFSRKeyFast(coords, fsr_key);

  /* If FSR has not been encountered, update FSR maps and vectors */
  if (!_FSR_keys_map.contains(fsr_key)) {

    /* Try to get a clean copy of the fsr_id, adding the FSR data
       if necessary where -1 indicates the key was already added */
    fsr_id = _FSR_keys_map.insert_and_get_count(fsr_key, NULL);
    if (fsr_id == -1) {
      fsr_data volatile* fsr;
      do {
        fsr = _FSR_keys_map.at(fsr_key);
      } while (fsr == NULL);
      fsr_id = fsr->_fsr_id;
    }
    else {

      /* Add FSR information to FSR key map and FSR_to vectors */
      fsr_data* fsr = new fsr_data;
      fsr->_fsr_id = fsr_id;
      _FSR_keys_map.update(fsr_key, fsr);
      Point* point = new Point();
      point->setCoords(coords->getHighestLevel()->getX(),
                       coords->getHighestLevel()->getY(),
                       coords->getHighestLevel()->getZ());

      /* Get the cell that contains coords */
      Cell* cell = findCellContainingCoords(curr);
      fsr->_point = point;
      fsr->_mat_id = cell->getFillMaterial()->getId();

      /* If CMFD acceleration is on, add FSR CMFD cell to FSR data */
      if (_cmfd != NULL)
        fsr->_cmfd_cell = _cmfd->findCmfdCell(coords->getHighestLevel());
    }
  }

  /* If FSR has already been encountered, get the fsr id from map */
  else {
    fsr_data volatile* fsr;
    do {
      fsr = _FSR_keys_map.at(fsr_key);
    } while (fsr == NULL);

    fsr_id = fsr->_fsr_id;
  }

  return fsr_id;
}


/**
 * @brief Finds and returns a pointer to the axially extruded flat source
 *        region that a given LocalCoords object resides within.
 * @param coords a LocalCoords object pointer
 * @return the ID of an extruded FSR for a given LocalCoords object
 */
int Geometry::findExtrudedFSR(LocalCoords* coords) {

  int fsr_id;
  LocalCoords* curr = coords;
  curr = coords->getLowestLevel();

  /* Generate unique FSR key */
  std::string fsr_key = getFSRKey(coords);

  /* If FSR has not been encountered, update FSR maps and vectors */
  if (!_extruded_FSR_keys_map.contains(fsr_key)) {

    /* Try to get a clean copy of the fsr_id, adding the FSR data
       if necessary where -1 indicates the key was already added */
    fsr_id = _extruded_FSR_keys_map.insert_and_get_count(fsr_key, NULL);
    if (fsr_id == -1) {
      ExtrudedFSR volatile* fsr;
      int count = 0;
      do {
        fsr = _extruded_FSR_keys_map.at(fsr_key);
        count++;
        if (count > 1e8)
          log_printf(ERROR, "Application stuck in finding extruded FSR");
      } while (fsr == NULL);
      fsr_id = fsr->_fsr_id;
    }
    else {

      /* Add FSR information to FSR key map and FSR_to vectors */
      ExtrudedFSR* fsr = new ExtrudedFSR;
      fsr->_fsr_id = fsr_id;
      fsr->_num_fsrs = 0;
      fsr->_coords = new LocalCoords(0, 0, 0, true);
      coords->copyCoords(fsr->_coords);
      _extruded_FSR_keys_map.update(fsr_key, fsr);
    }
  }

  /* If FSR has already been encountered, get the fsr id from map */
  else {
    ExtrudedFSR volatile* fsr;
    int count = 0;
    do {
      fsr = _extruded_FSR_keys_map.at(fsr_key);
      count++;
      if (count > 1e8)
        log_printf(ERROR, "Application stuck in finding extruded FSR");
    } while (fsr == NULL);

    fsr_id = fsr->_fsr_id;
  }

  return fsr_id;
}


/**
 * @brief Return the ID of the flat source region that a given
 *        LocalCoords object resides within.
 * @param coords a LocalCoords object pointer
 * @param err_check whether to fail instead of returning -1 if not found
 * @return the FSR ID for a given LocalCoords object
 */
long Geometry::getFSRId(LocalCoords* coords, bool err_check) {

  long fsr_id = 0;
  std::string fsr_key;

  try {
    /* Generate unique FSR key */
    int thread_id = omp_get_thread_num();
    std::string& fsr_key = _fsr_keys[thread_id];
    getFSRKeyFast(coords, fsr_key);
    if (!_FSR_keys_map.contains(fsr_key) && !err_check)
      return -1;
    fsr_id = _FSR_keys_map.at(fsr_key)->_fsr_id;
  }
  catch(std::exception &e) {
    if (err_check) {
        log_printf(ERROR, "Could not find FSR ID with key: %s. Try creating "
                  "geometry with finer track spacing", fsr_key.c_str());
    }
    else {
      return -1;
    }
  }

  return fsr_id;
}


/**
 * @brief Returns the rank of the domain containing the given coordinates.
 * @param coords The coordinates used to search the domain rank
 * @return The rank of the domain containing the coordinates
 */
int Geometry::getDomainByCoords(LocalCoords* coords) {

  int domain = 0;
#ifdef MPIx
  if (_domain_decomposed) {
    int domain_idx[3];
    double width_x = _root_universe->getMaxX() - _root_universe->getMinX();
    domain_idx[0] = (coords->getPoint()->getX() - _root_universe->getMinX())
                    * _num_domains_x / width_x;
    double width_y = _root_universe->getMaxY() - _root_universe->getMinY();
    domain_idx[1] = (coords->getPoint()->getY() - _root_universe->getMinY())
                    * _num_domains_y / width_y;
    double width_z = _root_universe->getMaxZ() - _root_universe->getMinZ();
    domain_idx[2] = (coords->getPoint()->getZ() - _root_universe->getMinZ())
                    * _num_domains_z / width_z;

    MPI_Cart_rank(_MPI_cart, domain_idx, &domain);
  }
#endif
  return domain;
}


/**
 * @brief Returns a map from cells to FSRs.
 * @return a map from cells to FSRs contained in those cells
 */
std::map<Cell*, std::vector<long> > Geometry::getCellsToFSRs() {

  std::map<int, Cell*> all_cells = _root_universe->getAllCells();
  std::map<int, Cell*>::iterator cell_iter;
  Cell* fsr_cell;
  std::map<Cell*, std::vector<long> > cells_to_fsrs;
  size_t num_FSRs = _FSR_keys_map.size();

  for (cell_iter = all_cells.begin(); cell_iter != all_cells.end();
       ++cell_iter) {

    /* Search for this Cell in all FSRs */
    for (long r=0; r < num_FSRs; r++) {
      fsr_cell = findCellContainingFSR(r);
      if (cell_iter->first == fsr_cell->getId())
        cells_to_fsrs[cell_iter->second].push_back(r);
    }
  }

  return cells_to_fsrs;
}


/**
 * @brief Return the global ID of the flat source region that a given
 *        LocalCoords object resides within.
 * @param coords a LocalCoords object pointer
 * @param err_check whether to fail instead of returning -1 if not found
 * @return the FSR ID for a given LocalCoords object
 */
long Geometry::getGlobalFSRId(LocalCoords* coords, bool err_check) {

  /* Check if the Geometry is domain decomposed */
  if (!_domain_decomposed) {
    return getFSRId(coords, err_check);
  }
  else {
    long temp_fsr_id = 0;
    long global_fsr_id = 0;
#ifdef MPIx
    int domain = getDomainByCoords(coords);
    int rank;
    MPI_Comm_rank(_MPI_cart, &rank);
    if (domain == rank)
      temp_fsr_id = getFSRId(coords);
    MPI_Allreduce(&temp_fsr_id, &global_fsr_id, 1, MPI_LONG, MPI_SUM,
                  _MPI_cart);

    /* Count FSRs on each domain if not already counted */
    if (!_domain_FSRs_counted)
      countDomainFSRs();

    /* Add FSR count from prior domains */
    for (long i=0; i < domain; i++)
      global_fsr_id += _num_domain_FSRs.at(i);

#endif
    return global_fsr_id;
  }
}


/**
 * @brief Return the characteristic point for a given FSR ID.
 * @param fsr_id the FSR ID
 * @return the FSR's characteristic point
 */
Point* Geometry::getFSRPoint(long fsr_id) {

  Point* point;

  try {
    std::string& key = _FSRs_to_keys[fsr_id];
    point = _FSR_keys_map.at(key)->_point;
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not find characteristic point in FSR: %d", fsr_id);
  }

  return point;
}


/**
 * @brief Return the centroid for a given FSR ID.
 * @param fsr_id the FSR ID
 * @return the FSR's centroid
 */
Point* Geometry::getFSRCentroid(long fsr_id) {

  if (fsr_id < _FSR_keys_map.size())
    return _FSRs_to_centroids[fsr_id];
  else
    log_printf(ERROR, "Could not find centroid in FSR: %d.", fsr_id);
  return NULL;
}


/**
 * @brief Returns whether any FSR centroids have been set
 */
bool Geometry::containsFSRCentroids() {
  return _contains_FSR_centroids;
}


#ifdef MPIx
/**
 * @brief Counts the number of FSRs in each MPI domain
 */
void Geometry::countDomainFSRs() {

  /* Gather the number of FSRs into an array */
  int num_domains = _num_domains_x * _num_domains_y * _num_domains_z;
  long num_domains_array[num_domains];
  long my_fsrs = getNumFSRs();
  MPI_Allgather(&my_fsrs, 1, MPI_LONG, num_domains_array, 1, MPI_LONG,
                _MPI_cart);

  /* Convert to a vector */
  _num_domain_FSRs.resize(num_domains);
  for (int i=0; i < num_domains; i++)
    _num_domain_FSRs.at(i) = num_domains_array[i];
  _domain_FSRs_counted = true;
}


/**
 * @brief Finds the local FSR ID and domain rank from a global FSR ID.
 * @param global_fsr_id The global unique identifier of an FSR
 * @param local_fsr_id The unique identifier of an FSR on its current domain
 * @param domain The rank of the domain containing the FSR
 */
void Geometry::getLocalFSRId(long global_fsr_id, long &local_fsr_id,
                             int &domain) {

  /* Count FSRs on each domain if not already counted */
  if (!_domain_FSRs_counted)
    countDomainFSRs();

  /* Determine the local domain where the global FSR resides */
  long cum_fsrs = 0;
  domain = -1;
  for (int i=0; i < _num_domains_x * _num_domains_y * _num_domains_z; i++) {
    if (cum_fsrs + _num_domain_FSRs.at(i) > global_fsr_id) {
      domain = i;
      break;
    }
    else {
      cum_fsrs += _num_domain_FSRs.at(i);
    }
  }

  /* Ensure a domain was found with the FSR ID */
  if (domain == -1)
    log_printf(ERROR, "No domain was found with the global FSR ID %ld. The "
               "total number of FSRs in the Geometry is %ld.", global_fsr_id,
               getNumTotalFSRs());

  local_fsr_id = global_fsr_id - cum_fsrs;
}
#endif


/**
 * @brief Returns the FSR centroid of a global FSR.
 * @param global_fsr_id The global unique identifier of the FSR
 * @return A vector containing the coordinates of the FSR centroid
 */
std::vector<double> Geometry::getGlobalFSRCentroidData(long global_fsr_id) {
  double xyz[3];
  if (!_domain_decomposed) {
    Point* centroid = getFSRCentroid(global_fsr_id);
    xyz[0] = centroid->getX();
    xyz[1] = centroid->getY();
    xyz[2] = centroid->getZ();
  }
#ifdef MPIx
  else {

    /* Determine the domain and local FSR ID */
    long fsr_id;
    int domain;
    getLocalFSRId(global_fsr_id, fsr_id, domain);

    /* Get the FSR centroid in the correct domain */
    int rank;
    MPI_Comm_rank(_MPI_cart, &rank);
    double temp_xyz[3];
    if (rank == domain) {
      Point* centroid = getFSRCentroid(fsr_id);
      temp_xyz[0] = centroid->getX();
      temp_xyz[1] = centroid->getY();
      temp_xyz[2] = centroid->getZ();
    }
    else {
      temp_xyz[0] = 0;
      temp_xyz[1] = 0;
      temp_xyz[2] = 0;
    }

    /* Broadcast the centroid */
    MPI_Allreduce(temp_xyz, xyz, 3, MPI_DOUBLE, MPI_SUM, _MPI_cart);
  }
#endif

  /* Convert centroid data into a vector */
  std::vector<double> data(3);
  for (int i=0; i<3; i++)
    data.at(i) = xyz[i];
  return data;
}


/**
 * @brief Return the CMFD cell for a given FSR ID
 * @param fsr_id the FSR ID
 * @return the CMFD cell
 */
int Geometry::getCmfdCell(long fsr_id) {
  return _FSRs_to_CMFD_cells[fsr_id];
}


/**
 * @brief Reserves memory for threaded constructs used by the Geometry class.
 * @param num_threads the number of threads
 */
void Geometry::setNumThreads(int num_threads) {

  /* Allocate fsr_keys for enough threads */
  reserveKeyStrings(num_threads);

  /* Allocate threaded constructs of the ExtrudedFSR maps */
  _FSR_keys_map.setNumThreads(num_threads);
  _extruded_FSR_keys_map.setNumThreads(num_threads);
}


/**
 * @brief Reserves memory for the FSR key strings for each thread
 * @param num_threads the number of threads
 */
void Geometry::reserveKeyStrings(int num_threads) {
  int string_size = 255;
  _fsr_keys.resize(num_threads);
  for (int i=0; i<num_threads; ++i) {
    _fsr_keys[i].reserve(string_size);
  }
}


/**
 * @brief Generate a string FSR "key" that identifies an FSR by its
 *        unique hierarchical lattice/universe/cell structure.
 * @details Since not all FSRs will reside on the absolute lowest universe
 *          level and Cells might overlap other cells, it is important to
 *          have a method for uniquely identifying FSRs. This method
 *          creates a unique FSR key by constructing a structured string
 *          that describes the hierarchy of lattices/universes/cells.
 * @param coords a LocalCoords object pointer
 * @param key a reference to the FSR key
 */
void Geometry::getFSRKeyFast(LocalCoords* coords, std::string& key) {

  LocalCoords* curr = coords->getHighestLevel();

  /* Assess the string size of the key */
  Point* point = curr->getPoint();
  int total_size = 0;
  if (_cmfd != NULL) {
    total_size += getNumDigits(_cmfd->getLattice()->getLatX(point));
    total_size += getNumDigits(_cmfd->getLattice()->getLatY(point));
    total_size += getNumDigits(_cmfd->getLattice()->getLatZ(point));
    total_size += 9;
  }
  if (_overlaid_mesh != NULL) {
    total_size += getNumDigits(_overlaid_mesh->getLatX(point));
    total_size += getNumDigits(_overlaid_mesh->getLatY(point));
    total_size += getNumDigits(_overlaid_mesh->getLatZ(point));
    total_size += 4;
  }
  while (curr != NULL) {
    if (curr->getType() == LAT) {
      total_size += getNumDigits(curr->getLattice()->getId());
      total_size += getNumDigits(curr->getLatticeX());
      total_size += getNumDigits(curr->getLatticeY());
      total_size += getNumDigits(curr->getLatticeZ());
      total_size += 6;
    }
    else {
      total_size += getNumDigits(curr->getUniverse()->getId());
      total_size += 2;
    }
    if (curr->getNext() == NULL)
      break;
    else
      curr = curr->getNext();
  }
  total_size += getNumDigits(curr->getCell()->getId());
  total_size += 1;
  int version_num = coords->getVersionNum();
  if (version_num != 0) {
    total_size += getNumDigits(version_num);
    total_size += 2;
  }

  /* Resize key */
  if (total_size > 255)
    log_printf(ERROR, "Found key exceeding the 255 character threshold");
  key.resize(total_size);
  curr = curr->getHighestLevel();

  /* If CMFD is on, get CMFD latice cell and write to key */
  int ind = 0;
  if (_cmfd != NULL) {
    key.replace(ind, 5, "CMFD(");
    ind += 5;
    printToString(key, ind, _cmfd->getLattice()->getLatX(point));
    key.at(ind) = ',';
    ind++;
    printToString(key, ind, _cmfd->getLattice()->getLatY(point));
    key.at(ind) = ',';
    ind++;
    printToString(key, ind, _cmfd->getLattice()->getLatZ(point));
    key.replace(ind, 2, "):");
    ind += 2;
  }

  /* If a global overlaid mesh is present, get the axial mesh cell */
  if (_overlaid_mesh != NULL) {
    key.at(ind) = 'A';
    ind++;
    printToString(key, ind, _overlaid_mesh->getLatZ(point));
    key.at(ind) = 'R';
    ind++;
    printToString(key, ind, _overlaid_mesh->getLatX(point));
    key.at(ind) = ',';
    ind++;
    printToString(key, ind, _overlaid_mesh->getLatY(point));
    key.at(ind) = ':';
    ind++;
  }

  /* Descend the linked list hierarchy until the lowest level has
   * been reached */
  while (curr != NULL) {

    /* write lattice to key */
    if (curr->getType() == LAT) {
      key.at(ind) = 'L';
      ind++;
      printToString(key, ind, curr->getLattice()->getId());
      key.at(ind) = '(';
      ind++;
      printToString(key, ind, curr->getLatticeX());
      key.at(ind) = ',';
      ind++;
      printToString(key, ind, curr->getLatticeY());
      key.at(ind) = ',';
      ind++;
      printToString(key, ind, curr->getLatticeZ());
      key.replace(ind, 2, "):");
      ind += 2;
    }
    else {

      /* write universe ID to key */
      key.at(ind) = 'U';
      ind++;
      printToString(key, ind, curr->getUniverse()->getId());
      key.at(ind) = ':';
      ind++;
    }

    /* If lowest coords reached break; otherwise get next coords */
    if (curr->getNext() == NULL)
      break;
    else
      curr = curr->getNext();
  }

  /* write cell id to key */
  key.at(ind) = 'C';
  ind++;
  printToString(key, ind, curr->getCell()->getId());

  /* write version number to key */
  if (version_num != 0) {
    key.replace(ind, 2, ":V");
    ind += 2;
    printToString(key, ind, version_num);
  }
}


//FIXME Find a better way to do this, without a function call
/**Using std::stringstream would be more clear.
 * @brief Get the number of digits in base 10 of a number
 * @param number the number of interest
 * @return the number of digits in base 10 of a number
 */
int Geometry::getNumDigits(int number) {
  if (number < 0)
    log_printf(ERROR, "Trying to get the digits of negative number %d", number);
  int ref = 10;
  int num_digits = 1;
  while (number >= ref) {
    num_digits++;
    ref *= 10;
  }
  return num_digits;
}


//FIXME Find a better way to do this, without a function call
/**Using std::stringstream would be more clear.
 * @brief Print a number to a given String.
 * @param str the string to print to
 * @param index the last index in that string
 * @param value the number to print
 * @return the number of digits in base 10 of a number
 */
void Geometry::printToString(std::string& str, int& index, int value) {

  char digits[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

  int num_digits = getNumDigits(value);
  for (int i=0; i < num_digits; i++) {
    int base = 1;
    for (int j=0; j < num_digits - i - 1; j++)
      base *= 10;
    char digit = digits[value / base];
    str.at(index+i) = digit;
    value -= (value / base) * base;
  }
  index += num_digits;
}


/**
 * @brief Generate a string FSR "key" for the FSR where the point reside in. A
          string FSR "key" identifies an FSR by its unique hierarchical
 *        lattice/universe/cell structure.
 * @details Since not all FSRs will reside on the absolute lowest universe
 *          level and Cells might overlap other cells, it is important to
 *          have a method for uniquely identifying FSRs. This method
 *          creates a unique FSR key by constructing a structured string
 *          that describes the hierarchy of lattices/universes/cells.
 * @param coords a LocalCoords object pointer
 * @return the FSR key
 */
std::string Geometry::getFSRKey(LocalCoords* coords) {

  std::stringstream key;
  LocalCoords* curr = coords->getHighestLevel();
  std::ostringstream curr_level_key;

  /* If CMFD is on, get CMFD latice cell and write to key */
  if (_cmfd != NULL) {
      curr_level_key << _cmfd->getLattice()->getLatX(curr->getPoint());
      key << "CMFD(" << curr_level_key.str() << ",";
      curr_level_key.str(std::string());
      curr_level_key << _cmfd->getLattice()->getLatY(curr->getPoint());
      key << curr_level_key.str() << ",";
      curr_level_key.str(std::string());
      curr_level_key << _cmfd->getLattice()->getLatZ(curr->getPoint());
      key << curr_level_key.str() << "):";
  }

  /* If a global overlaid mesh is present, record mesh cells */
  if (_overlaid_mesh != NULL) {
      curr_level_key.str(std::string());
      curr_level_key << _overlaid_mesh->getLatZ(curr->getPoint());
      key << "A" << curr_level_key.str() << ":";
      curr_level_key.str(std::string());
      curr_level_key << _overlaid_mesh->getLatX(curr->getPoint());
      key << "R" << curr_level_key.str() << ",";
      curr_level_key.str(std::string());
      curr_level_key << _overlaid_mesh->getLatY(curr->getPoint());
      key << curr_level_key.str() << ":";
  }

  /* Descend the linked list hierarchy until the lowest level has
   * been reached */
  while (curr != NULL) {

    /* Clear string stream */
    curr_level_key.str(std::string());

    if (curr->getType() == LAT) {

      /* Write lattice ID and lattice cell to key */
      curr_level_key << curr->getLattice()->getId();
      key << "L" << curr_level_key.str() << "(";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeX();
      key << curr_level_key.str() << ",";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeY();
      key << curr_level_key.str() << ",";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeZ();
      key << curr_level_key.str() << "):";
    }
    else {
      /* write universe ID to key */
      curr_level_key << curr->getUniverse()->getId();
      key << "U" << curr_level_key.str() << " : ";
    }

    /* If lowest coords reached break; otherwise get next coords */
    if (curr->getNext() == NULL)
      break;
    else
      curr = curr->getNext();
  }

  /* clear string stream */
  curr_level_key.str(std::string());

  /* write cell id to key */
  curr_level_key << curr->getCell()->getId();
  key << "C" << curr_level_key.str();

  /* write version number to key */
  int version_num = coords->getVersionNum();
  if (version_num != 0) {
    curr_level_key.str(std::string());
    curr_level_key << version_num;
    key << ":V" << curr_level_key.str();
  }

  return key.str();
}


/**
 * @brief Return a pointer to an ExtrudedFSR by its extruded FSR ID
 * @param extruded_fsr_id the extruded FSR ID
 * @return a pointer to the ExtrudedFSR
 */
ExtrudedFSR* Geometry::getExtrudedFSR(int extruded_fsr_id) {
  return _extruded_FSR_lookup[extruded_fsr_id];
}


/**
 * @brief Subdivides all Cells in the Geometry into rings and angular sectors
 *        aligned with the z-axis.
 * @details This method is called by the Geometry::initializeFlatSourceRegions()
 *          method but may also be called by the user in Python if needed:
 *
 * @code
 *          geometry.subdivideCells()
 * @endcode
 */
void Geometry::subdivideCells() {

  /* Compute equivalent radius with the same area as the Geometry */
  /* This is used as the maximum radius for all ringified Cells */
  double width_x = _root_universe->getMaxX() - _root_universe->getMinX();
  double width_y = _root_universe->getMaxY() - _root_universe->getMinY();
  double max_radius = sqrt(width_x * width_x + width_y * width_y) / 2;

  /* Recursively subdivide Cells into rings and sectors */
  _root_universe->subdivideCells(max_radius);
}


/**
 * @brief Compute the number of flat source regions in the Geometry and
 *        initialize CMFD.
 * @details This method is intended to be called by the user before initiating
 *          source iteration. This method first subdivides all Cells by calling
 *          the Geometry::subdivideCells() method. Then it initializes the CMFD
 *          object.
 */
void Geometry::initializeFlatSourceRegions() {

  log_printf(NORMAL, "Initializing flat source regions...");

  /* Subdivide Cells into sectors and rings */
  subdivideCells();

  /* Build collections of neighbor Cells for optimized ray tracing */
  //FIXME
  //_root_universe->buildNeighbors();

  /* Create map of Material IDs to Material pointers */
  _all_materials = getAllMaterials();

  /* Initialize absorption XS if CMFD absent */
  if (_cmfd == NULL) {
    std::map<int, Material*>::iterator iter;
    for (iter = _all_materials.begin(); iter != _all_materials.end(); ++iter)
      iter->second->getSigmaA();
  }

  /* Initialize CMFD */
  if (_cmfd != NULL)
    initializeCmfd();
}


/**
 * @brief This method performs ray tracing to create Track segments within each
 *        flat source region in the Geometry.
 * @details This method starts at the beginning of a Track and finds successive
 *          intersection points with FSRs as the Track crosses through the
 *          Geometry and creates segment structs and adds them to the Track.
 * @param track a pointer to a track to segmentize
 * @param z_coord the axial height at which the 2D plane of the geometry is
 *        formed
 */
void Geometry::segmentize2D(Track* track, double z_coord) {

  /* Track starting Point coordinates and azimuthal angle */
  double x0 = track->getStart()->getX();
  double y0 = track->getStart()->getY();
  double z0 = z_coord;
  double phi = track->getPhi();
  double delta_x, delta_y, delta_z;

  /* Length of each segment */
  double length;
  Material* material;
  int fsr_id;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords start(x0, y0, z0, true);
  LocalCoords end(x0, y0, z0, true);
  start.setUniverse(_root_universe);
  end.setUniverse(_root_universe);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&end, phi);
  Cell* prev;

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL)
    log_printf(ERROR, "Could not find a Cell containing the start Point "
               "of this Track: %s", track->toString().c_str());

  /* While the end of the segment's LocalCoords is still within the Geometry,
   * move it to the next Cell, create a new segment, and add it to the
   * Geometry */
  while (curr != NULL) {

    end.copyCoords(&start);

    /* Find the next Cell along the Track's trajectory */
    prev = curr;
    curr = findNextCell(&end, phi);

    /* Checks that segment does not have the same start and end Points */
    if (fabs(start.getX() - end.getX()) < FLT_EPSILON
        && fabs(start.getY() - end.getY()) < FLT_EPSILON)
      log_printf(ERROR, "Created 2D segment with same start and end point: "
                 "x = %f, y = %f, z = %f", start.getX(), start.getY(),
                 start.getZ());

    /* Find the segment length, Material and FSR ID */
    length = double(end.getPoint()->distanceToPoint(start.getPoint()));
    material = prev->getFillMaterial();
    fsr_id = findFSRId(&start);

    /* Create a new Track segment */
    segment new_segment;
    new_segment._material = material;
    new_segment._length = length;
    new_segment._region_id = fsr_id;

    log_printf(DEBUG, "segment start x = %f, y = %f; end x = %f, y = %f",
               start.getX(), start.getY(), end.getX(), end.getY());

    /* Save indices of CMFD Mesh surfaces that the Track segment crosses */
    if (_cmfd != NULL) {

      /* Find cmfd cell that segment lies in */
      int cmfd_cell = _cmfd->findCmfdCell(&start);

      /* Reverse nudge from surface to determine whether segment start or end
       * points lie on a CMFD surface. */
      delta_x = cos(phi) * TINY_MOVE;
      delta_y = sin(phi) * TINY_MOVE;
      start.adjustCoords(-delta_x, -delta_y);
      end.adjustCoords(-delta_x, -delta_y);

      /* Calculate CMFD surfaces */
      int cmfd_surfaces[2];
      cmfd_surfaces[0] = _cmfd->findCmfdSurface(cmfd_cell, &end);
      cmfd_surfaces[1] = _cmfd->findCmfdSurface(cmfd_cell, &start);

      /* Ensure surfaces are x-y surfaces (no z-crossings) */
      /* Note: this code takes advantage of the numeric representation of
         surfaces to find a mapping that removes z-surfaces */
      for (int d=0; d<2; d++) {
        int local_surface = cmfd_surfaces[d] % NUM_SURFACES;
        if (local_surface == 2 || local_surface == 5) {
            cmfd_surfaces[d] = -1;
        }
        else if (local_surface > 9) {
          int cell = cmfd_surfaces[d] / NUM_SURFACES;
          int half_surf = local_surface / 2;
          if (local_surface > 17) {
            int quart_surf = half_surf / 2;
            local_surface = 2 + quart_surf + (half_surf == 2*quart_surf);
            cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface;
          }
          else {
            local_surface = (half_surf > 6) + 3 *
                (local_surface != 2*half_surf);
            cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface;
          }
        }
      }

      /* Save CMFD surfaces */
      new_segment._cmfd_surface_fwd = cmfd_surfaces[0];
      new_segment._cmfd_surface_bwd = cmfd_surfaces[1];

      /* Re-nudge segments from surface. */
      start.adjustCoords(delta_x, delta_y);
      end.adjustCoords(delta_x, delta_y);
    }

    /* Calculate the local centroid of the segment if available */
    //FIXME Consider reversing nudge
    Point* starting_point = start.getHighestLevel()->getPoint();
    new_segment._starting_position[0] = starting_point->getX();
    new_segment._starting_position[1] = starting_point->getY();
    if (_contains_FSR_centroids) {
      Point* centroid = getFSRCentroid(fsr_id);
      double x_start = starting_point->getX() - centroid->getX();
      double y_start = starting_point->getY() - centroid->getY();
      new_segment._starting_position[0] = x_start;
      new_segment._starting_position[1] = y_start;
    }

    /* Add the segment to the Track */
    track->addSegment(&new_segment);
  }

  log_printf(DEBUG, "Created %d segments for Track: %s",
             track->getNumSegments(), track->toString().c_str());

  /* Truncate the linked list for the LocalCoords */
  start.prune();
  end.prune();
}


/**
 * @brief This method performs ray tracing to create Track segments within each
 *        flat source region in the Geometry.
 * @details This method starts at the beginning of a Track and finds successive
 *          intersection points with FSRs as the Track crosses through the
 *          Geometry and creates segment structs and adds them to the Track.
 * @param track a pointer to a track to segmentize
 * @param OTF_setup whether this routine is called during OTF ray tracing setup
 */
void Geometry::segmentize3D(Track3D* track, bool OTF_setup) {

  /* Track starting Point coordinates and azimuthal angle */
  double x0 = track->getStart()->getX();
  double y0 = track->getStart()->getY();
  double z0 = track->getStart()->getZ();
  double phi = track->getPhi();
  double theta = track->getTheta();

  /* Length of each segment */
  double length;
  long fsr_id;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords start(x0, y0, z0, true);
  LocalCoords end(x0, y0, z0, true);
  start.setUniverse(_root_universe);
  end.setUniverse(_root_universe);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&end, phi, theta);
  Cell* prev;

  /* Vector to fill coordinates if necessary */
  std::vector<LocalCoords*> fsr_coords;

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL)
    log_printf(ERROR, "Could not find a Cell containing the start Point "
               "of this Track3D: %s", track->toString().c_str());

  /* While the end of the segment's LocalCoords is still within the Geometry,
   * move it to the next Cell, create a new segment, and add it to the
   * Geometry */
  while (curr != NULL) {

    end.copyCoords(&start);

    /* Find the next Cell along the Track's trajectory */
    prev = curr;
    curr = findNextCell(&end, phi, theta);

    /* Checks to make sure that new Segment does not have the same start
     * and end Points */
    if (fabs(start.getX() - end.getX()) < FLT_EPSILON &&
        fabs(start.getY() - end.getY()) < FLT_EPSILON &&
        fabs(start.getZ() - end.getZ()) < FLT_EPSILON) {
      log_printf(ERROR, "Created a Track3D segment with the same start and end "
                 "point: x = %f, y = %f, z = %f", start.getX(),
                 start.getY(), start.getZ());
    }

    /* Find the segment length and its region's id */
    length = end.getPoint()->distanceToPoint(start.getPoint());
    long fsr_id = findFSRId(&start);

    /* Create a new Track segment */
    segment new_segment;
    new_segment._material = prev->getFillMaterial();
    new_segment._length = length;
    new_segment._region_id = fsr_id;

    log_printf(DEBUG, "segment start x = %f, y = %f, z = %f; "
               "end x = %f, y = %f, z = %f",
               start.getX(), start.getY(), start.getZ(),
               end.getX(), end.getY(), end.getZ());

    /* Save indices of CMFD Mesh surfaces that the Track segment crosses */
    if (_cmfd != NULL && !OTF_setup) {

      /* Find cmfd cell that segment lies in */
      int cmfd_cell = _cmfd->findCmfdCell(&start);

      /* Reverse nudge from surface to determine whether segment start or end
       * points lie on a cmfd surface. */
      double delta_x = cos(phi) * sin(theta) * TINY_MOVE;
      double delta_y = sin(phi) * sin(theta) * TINY_MOVE;
      double delta_z = cos(theta) * TINY_MOVE;
      start.adjustCoords(-delta_x, -delta_y, -delta_z);
      end.adjustCoords(-delta_x, -delta_y, -delta_z);

      new_segment._cmfd_surface_fwd =
        _cmfd->findCmfdSurface(cmfd_cell, &end);
      new_segment._cmfd_surface_bwd =
        _cmfd->findCmfdSurface(cmfd_cell, &start);

      /* Re-nudge segments from surface. */
      start.adjustCoords(delta_x, delta_y, delta_z);
      end.adjustCoords(delta_x, delta_y, delta_z);
    }

    /* For regular 3D tracks, get starting position relative to FSR centroid */
    if (!OTF_setup) {
      Point* starting_point = start.getHighestLevel()->getPoint();
      new_segment._starting_position[0] = starting_point->getX();
      new_segment._starting_position[1] = starting_point->getY();
      new_segment._starting_position[2] = starting_point->getZ();
      if (_contains_FSR_centroids) {
        Point* centroid = getFSRCentroid(fsr_id);
        double x_start = starting_point->getX() - centroid->getX();
        double y_start = starting_point->getY() - centroid->getY();
        double z_start = starting_point->getZ() - centroid->getZ();
        new_segment._starting_position[0] = x_start;
        new_segment._starting_position[1] = y_start;
        new_segment._starting_position[2] = z_start;
      }
    }

    /* Add the segment to the Track */
    track->addSegment(&new_segment);
  }

  log_printf(DEBUG, "Created %d segments for Track3D: %s",
             track->getNumSegments(), track->toString().c_str());

  /* Truncate the linked list for the LocalCoords */
  start.prune();
  end.prune();
}


/**
 * @brief This method performs ray tracing to create extruded track segments
 *        within each flat source region in the implicit 2D superposition plane
 *        of the Geometry by 2D ray tracing across input heights that encompass
 *        all radial geometric details.
 * @details This method starts at the beginning of an extruded track and finds
 *          successive intersection points with FSRs as the extruded track
 *          crosses radially through the Geometry at defined z-coords. The
 *          minimum distance to intersection of all z-coords is chosen leading
 *          to implicitly capturing all geometric radial detail at the defined
 *          z-heights, saving the lengths and region IDs to the extruded track
 *          and initializing ExtrudedFSR structs in the traversed FSRs.
 * @param flattened_track a pointer to a 2D track to segmentize into regions of
 *        extruded FSRs
 * @param z_coords a vector of axial heights in the root geometry at which
 *        the Geometry is segmentized radially
 */
void Geometry::segmentizeExtruded(Track* flattened_track,
    std::vector<double> z_coords) {

  /* Track starting Point coordinates and azimuthal angle */
  double x0 = flattened_track->getStart()->getX();
  double y0 = flattened_track->getStart()->getY();
  double z0 = z_coords[0];
  double phi = flattened_track->getPhi();
  double delta_x, delta_y, delta_z;

  /* Length of each segment */
  double length;
  int min_z_ind;
  int region_id;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords start(x0, y0, z0, true);
  LocalCoords end(x0, y0, z0, true);
  start.setUniverse(_root_universe);
  end.setUniverse(_root_universe);

  /* Create two localCoords to check results */
  LocalCoords test_ext_coords(0,0,0,true);
  LocalCoords test_start_coords(0,0,0,true);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&end, phi);

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL) {
    int dom = _domain_index_x + _domain_index_y * _num_domains_x +
      _domain_index_z * _num_domains_x * _num_domains_y;
    log_printf(ERROR, "Could not find a Cell containing the start Point "
               "of this Track: %s on domain %d with bounds [%f, %f] x [%f, %f]"
               " x [%f, %f]", flattened_track->toString().c_str(), dom,
               getMinX(), getMaxX(), getMinY(), getMaxY(), getMinZ(),
               getMaxZ());
  }

  /* While the end of the segment's LocalCoords is still within the Geometry,
   * move it to the next Cell, create a new segment, and add it to the
   * Geometry */
  int find_cell_count = 0;
  while (curr != NULL) {

    /* Check if stuck in loop */
    find_cell_count++;
    if (find_cell_count > 1e6)
      log_printf(ERROR, "Caught in infinite loop finding next cell");

    /* Records the minimum length to a 2D intersection */
    double min_length = std::numeric_limits<double>::infinity();
    region_id = -1;
    min_z_ind = -1;

    /* Copy end coordinates to start */
    end.copyCoords(&start);

    /* Loop over all z-heights to find shortest 2D intersection */
    for (int i=0; i < z_coords.size(); i++) {

      /* Change z-height and copy starting coordinates to end */
      start.setZ(z_coords[i]);
      start.prune();
      findCellContainingCoords(&start);
      start.copyCoords(&end);

      /* Find the next Cell along the Track's trajectory */
      curr = findNextCell(&end, phi);

      /* Checks that segment does not have the same start and end Points */
      if (fabs(start.getX() - end.getX()) < FLT_EPSILON &&
          fabs(start.getY() - end.getY()) < FLT_EPSILON)
        log_printf(ERROR, "Created segment with same start and end "
                   "point: x = %f, y = %f", start.getX(), start.getY());

      /* Find the segment length and extruded FSR */
      length = double(end.getPoint()->distanceToPoint(start.getPoint()));

      /* Check if the segment length is the smallest found */
      if (length < min_length) {
        min_length = length;
        min_z_ind = i;
      }
    }

    /* Traverse across shortest segment */
    start.prune();
    start.setZ(z_coords[min_z_ind]);
    findCellContainingCoords(&start);

    bool found_coordinate = false;

    int next_version = 0;
    for (int v=0; v < MAX_VERSION_NUM; v++) {

      /* Find FSR using starting coordinate */
      start.setVersionNum(v);
      region_id = findExtrudedFSR(&start); //FIXME
      std::string fsr_key = getFSRKey(&start);

      /* Get the coordinate of the extruded FSR */
      LocalCoords* volatile retrieved_coords = NULL;
      do {
        retrieved_coords = _extruded_FSR_keys_map.at(fsr_key)->_coords;
      } while (retrieved_coords == NULL);
      LocalCoords* ext_coords = retrieved_coords;

      /* Create coordinate copies */
      ext_coords->copyCoords(&test_ext_coords);
      start.copyCoords(&test_start_coords);

      /* Check to see that this point contains the cell of every axial level */
      bool coords_contained = true;
      for (int i=0; i < z_coords.size(); i++) {

        /* Check the FSR key at this level */
        test_start_coords.setZ(z_coords[i]);
        test_start_coords.prune();
        test_start_coords.setVersionNum(0);
        findCellContainingCoords(&test_start_coords);
        fsr_key = getFSRKey(&test_start_coords);

        test_ext_coords.setZ(z_coords[i]);
        test_ext_coords.prune();
        test_ext_coords.setVersionNum(0);
        findCellContainingCoords(&test_ext_coords);
        std::string ext_fsr_key = getFSRKey(&test_ext_coords);

        /* Check that FSR keys match */
        if (fsr_key != ext_fsr_key) {
          coords_contained = false;
          break;
        }
      }

      /* Check if we found a valid coordinate */
      if (coords_contained) {
        found_coordinate = true;
        break;
      }

      /* Reset the starting coordinate */
      next_version++;
    }

    if (next_version >= MAX_VERSION_NUM)
      log_printf(ERROR, "Exceeded the maximum version number in 2D extruded "
                 "FSRs");

    /* Move the coordinates to the next intersection */
    start.copyCoords(&end);
    curr = findNextCell(&end, phi, M_PI_2);

    log_printf(DEBUG, "segment start x = %f, y = %f; end x = %f, y = %f",
               start.getX(), start.getY(), end.getX(), end.getY());

    /* Check that the region ID is valid */
    if (region_id == -1)
      log_printf(ERROR, "Failed to find a valid FSR during axial extruded "
                 "segmentation");

    /* Create a new 2D Track segment with extruded region ID */
    segment new_segment;
    new_segment._length = min_length;
    new_segment._region_id = region_id;

    /* Save indices of CMFD Mesh surfaces that the Track segment crosses */
    if (_cmfd != NULL) {

      /* Find cmfd cell that segment lies in */
      int cmfd_cell = _cmfd->findCmfdCell(&start);

      /* Reverse nudge from surface to determine whether segment start or end
       * points lie on a CMFD surface. */
      delta_x = cos(phi) * TINY_MOVE;
      delta_y = sin(phi) * TINY_MOVE;
      start.adjustCoords(-delta_x, -delta_y, 0);
      end.adjustCoords(-delta_x, -delta_y, 0);

      /* Calculate CMFD surfaces */
      int cmfd_surfaces[2];
      cmfd_surfaces[0] = _cmfd->findCmfdSurface(cmfd_cell, &end);
      cmfd_surfaces[1] = _cmfd->findCmfdSurface(cmfd_cell, &start);

      /* Ensure surfaces are x-y surfaces (no z-crossings) */
      /* Note: this code takes advantage of the numeric representation of
         surfaces to find a mapping that removes z-surfaces */
      for (int d=0; d<2; d++) {
        int local_surface = cmfd_surfaces[d] % NUM_SURFACES;
        if (local_surface == 2 || local_surface == 5) {
            cmfd_surfaces[d] = -1;
        }
        else if (local_surface > 9) {
          int cell = cmfd_surfaces[d] / NUM_SURFACES;
          int half_surf = local_surface / 2;
          if (local_surface > 17) {
            int quart_surf = half_surf / 2;
            local_surface = 2 + quart_surf + (half_surf == 2*quart_surf);
            cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface;
          }
          else {
            local_surface = (half_surf > 6) + 3 *
                (local_surface != 2*half_surf);
            cmfd_surfaces[d] = cell * NUM_SURFACES + local_surface;
          }
        }
      }

      /* Save CMFD surfaces */
      new_segment._cmfd_surface_fwd = cmfd_surfaces[0];
      new_segment._cmfd_surface_bwd = cmfd_surfaces[1];

      /* Re-nudge segments from surface. */
      start.adjustCoords(delta_x, delta_y, 0);
      end.adjustCoords(delta_x, delta_y, 0);
    }

    /* Add the segment to the 2D track */
    flattened_track->addSegment(&new_segment);
  }

  /* Truncate the linked list for the LocalCoords */
  start.prune();
  end.prune();
}


/**
 * @brief Fixes the FSR map size so that the map is static and fast.
 */
void Geometry::fixFSRMaps() {
  _FSR_keys_map.setFixedSize();
  _extruded_FSR_keys_map.setFixedSize();
}


/**
 * @brief Rays are shot vertically through each ExtrudedFSR struct to calculate
 *        the axial mesh and initialize 3D FSRs.
 * @details From a 2D point within each FSR, a temporary 3D track is created
 *          starting at the bottom of the geometry and extending vertically to
 *          the top of the geometry. These tracks are segmented using the
 *          segmentize3D routine to calculate the distances between axial
 *          intersections forming the axial meshes if necessary and
 *          initializing the 3D FSRs as new regions are traversed.
 * @param global_z_mesh A global z mesh used for ray tracing. If the vector's
 *        length is zero, z meshes are local and need to be created for every
 *        ExtrudedFSR. The global_z_mesh is local to the domain when using
 *        domain decomposition.
 */
void Geometry::initializeAxialFSRs(std::vector<double> global_z_mesh) {

  log_printf(NORMAL, "Initializing 3D FSRs in axially extruded regions...");

  /* Determine the extent of the axial geometry */
  double min_z = getMinZ();
  double max_z = getMaxZ();

  /* Extract list of extruded FSRs */
  ExtrudedFSR** extruded_FSRs = _extruded_FSR_keys_map.values();

  std::string msg = "initializing 3D FSRs";
  Progress progress(_extruded_FSR_keys_map.size(), msg, 0.1, this, true);

  /* Re-allocate the FSR keys map with the new anticipated size */
  int anticipated_size = 2 * _extruded_FSR_keys_map.size();
  if (_overlaid_mesh != NULL)
    anticipated_size *= _overlaid_mesh->getNumZ();
  else
    if (_cmfd != NULL)
      anticipated_size *= _cmfd->getLocalNumZ();
  if (anticipated_size > _FSR_keys_map.bucket_count())
    _FSR_keys_map.realloc(anticipated_size);
  long total_number_fsrs_in_stack = 0;

  /* Loop over extruded FSRs */
#pragma omp parallel for
  for (int i=0; i < _extruded_FSR_keys_map.size(); i++) {

    progress.incrementCounter();

    /* Extract coordinates of extruded FSR */
    ExtrudedFSR* extruded_FSR = extruded_FSRs[i];
    double x0 = extruded_FSR->_coords->getX();
    double y0 = extruded_FSR->_coords->getY();

    /* Determine if there is a global mesh or local meshes should be created */
    if (global_z_mesh.size() > 0) {

      /* Allocate materials in extruded FSR */
      size_t num_regions = global_z_mesh.size() - 1;
      extruded_FSR->_num_fsrs = num_regions;
      extruded_FSR->_materials = new Material*[num_regions];
      extruded_FSR->_fsr_ids = new long[num_regions];

      /* Loop over all regions in the global mesh */
#pragma omp critical
      {
        for (int n=0; n < num_regions; n++) {

          /* Set the axial coordinate at the midpoint of mesh boundaries */
          double midpt = (global_z_mesh[n] + global_z_mesh[n+1]) / 2;
          LocalCoords coord(x0, y0, midpt);
          coord.setUniverse(_root_universe);

          /* Get the FSR ID and material */
          Cell* cell = findCellContainingCoords(&coord);
          long fsr_id = findFSRId(&coord);
          Material* material = cell->getFillMaterial();

          /* Set the FSR ID and material */
          extruded_FSR->_fsr_ids[n] = fsr_id;
          extruded_FSR->_materials[n] = material;
        }
      }
    }
    else {

      /* Create vertical track in the extruded FSR */
      Track3D track;
      track.setValues(x0, y0, min_z, x0, y0, max_z, 0, 0);

      /* Shoot vertical track through the geometry to initialize 3D FSRs */
      segmentize3D(&track, true);

      /* Extract segments from track */
      int num_segments = track.getNumSegments();
      segment* segments = track.getSegments();

      /* Allocate materials and mesh in extruded FSR */
      extruded_FSR->_num_fsrs = (size_t) num_segments;
      extruded_FSR->_materials = new Material*[num_segments];
      extruded_FSR->_fsr_ids = new long[num_segments];
      extruded_FSR->_mesh = new double[num_segments+1];

      /* Initialize values in extruded FSR */
      for (int s=0; s < num_segments; s++) {
        extruded_FSR->_materials[s] = segments[s]._material;
        extruded_FSR->_fsr_ids[s] = segments[s]._region_id;
      }

      /* Initialize z mesh */
      double level = min_z;
      extruded_FSR->_mesh[0] = level;
      for (int s=0; s < num_segments; s++) {
        level += segments[s]._length;
        if (std::abs(level - max_z) < 1e-12)
          level = max_z;
        extruded_FSR->_mesh[s+1] = level;
      }
    }
    /* Keep track of the number of FSRs in extruded FSRs */
#pragma omp atomic update
    total_number_fsrs_in_stack += extruded_FSR->_num_fsrs;
  }

  delete [] extruded_FSRs;
#ifdef MPIx
  if (_domain_decomposed)
    MPI_Barrier(_MPI_cart);
#endif

  // Output the extruded FSR storage requirement
  float size = total_number_fsrs_in_stack * (sizeof(double) + sizeof(long) +
             sizeof(Material*)) + _extruded_FSR_keys_map.size() * (sizeof(
             _extruded_FSR_keys_map.keys()[0]) + sizeof(ExtrudedFSR) +
             (LOCAL_COORDS_LEN + 1) * sizeof(LocalCoords));
  float max_size = size;
#ifdef MPIX
    if (isDomainDecomposed())
      MPI_Allreduce(&size, &max_size, 1, MPI_FLOAT, MPI_MAX, _MPI_cart);
#endif

  log_printf(INFO_ONCE, "Max storage for extruded FSRs per domain = %.2f MB",
             max_size / float(1e6));

  /* Re-order FSR IDs so they are sequential in the axial direction */
  reorderFSRIDs();
  log_printf(NORMAL, "Finished initializing 3D FSRs");
}


/**
 * @brief Reorders FSRs so that they are contiguous in the axial direction.
 */
void Geometry::reorderFSRIDs() {

  /* Extract list of extruded FSRs */
  ExtrudedFSR** extruded_FSRs = _extruded_FSR_keys_map.values();

  std::string msg = "reordering FSR IDs";
  Progress progress(_extruded_FSR_keys_map.size(), msg, 0.1, this, true);

  /* Get the FSR data objects */
  long curr_id = 0;
  fsr_data** value_list = _FSR_keys_map.values();
  fsr_data** fsr_data_objects = new fsr_data*[_FSR_keys_map.size()];

  /* Create a mapping of old to new IDs */
  long* id_mapping = new long[_FSR_keys_map.size()];
  bool* id_remapped = new bool[_FSR_keys_map.size()];
#pragma omp parallel for
  for (long i=0; i < _FSR_keys_map.size(); i++) {
    long id = value_list[i]->_fsr_id;
    fsr_data_objects[id] = value_list[i];
    id_mapping[i] = i;
    id_remapped[i] = false;
  }

  /* Loop over extruded FSRs */
  long count = 0;
  for (int i=0; i < _extruded_FSR_keys_map.size(); i++) {

    progress.incrementCounter();

    /* Extract coordinates of extruded FSR */
    ExtrudedFSR* extruded_FSR = extruded_FSRs[i];

    /* Get the number of FSRs in this axial region */
    size_t num_local_fsrs = extruded_FSR->_num_fsrs;

    /* Re-assign the IDs of all axial FSRs */
    for (int j=0; j < num_local_fsrs; j++) {

      long previous_id = extruded_FSR->_fsr_ids[j];
      if (!id_remapped[previous_id]) {
        id_mapping[previous_id] = count;
        id_remapped[previous_id] = true;
        count++;
      }
      long new_id = id_mapping[previous_id];

      fsr_data_objects[previous_id]->_fsr_id = new_id;
      extruded_FSR->_fsr_ids[j] = new_id;
    }
  }

  delete [] extruded_FSRs;
  delete [] value_list;
  delete [] fsr_data_objects;
  delete [] id_mapping;
  delete [] id_remapped;
}


/**
 * @brief Initialize key and material ID vectors for lookup by FSR ID.
 * @details This function initializes and sets reverse lookup vectors by FSR ID.
 *      This is called after the FSRs have all been identified and allocated
 *      during segmentation. This function must be called after
 *      Geometry::segmentize() has completed. It should not be called if tracks
 *      are loaded from a file.
 */
void Geometry::initializeFSRVectors() {

  /* Get keys and values from map */
  log_printf(NORMAL, "Initializing FSR lookup vectors");
  std::string *key_list = _FSR_keys_map.keys();
  fsr_data **value_list = _FSR_keys_map.values();

  /* Allocate vectors */
  size_t num_FSRs = _FSR_keys_map.size();
  _FSRs_to_keys = std::vector<std::string>(num_FSRs);
  _FSRs_to_centroids = std::vector<Point*>(num_FSRs, NULL);
  _FSRs_to_material_IDs = std::vector<int>(num_FSRs);
  _FSRs_to_CMFD_cells = std::vector<int>(num_FSRs);
  _contains_FSR_centroids = false;

  /* Fill vectors key and material ID information */
#pragma omp parallel for
  for (long i=0; i < num_FSRs; i++) {
    std::string key = key_list[i];
    fsr_data* fsr = value_list[i];
    long fsr_id = fsr->_fsr_id;
    _FSRs_to_keys.at(fsr_id) = key;
    _FSRs_to_material_IDs.at(fsr_id) = fsr->_mat_id;
  }

  /* Add cmfd information serially */
  if (_cmfd != NULL) {
    for (long i=0; i < num_FSRs; i++) {
      fsr_data* fsr = value_list[i];
      long fsr_id = fsr->_fsr_id;
      _cmfd->addFSRToCell(fsr->_cmfd_cell, fsr_id);
      _FSRs_to_CMFD_cells.at(fsr_id) = fsr->_cmfd_cell;
    }
  }

  /* Output approximate storage for various FSR maps, locks, volumes... */
  long size = num_FSRs * (sizeof(fsr_data) + sizeof(omp_lock_t) +
       sizeof(FP_PRECISION) + sizeof(fsr_data*) + 2 * sizeof(key_list[0]) +
       sizeof(Point*) + sizeof(Point) + 2*sizeof(int));
  long max_size = size;
#ifdef MPIX
  if (isDomainDecomposed())
    MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX, _MPI_cart);
#endif
  log_printf(INFO_ONCE, "Max FSR, maps and data, storage per domain = %.2f MB",
             max_size / float(1e6));

  /* Check if extruded FSRs are present */
  size_t num_extruded_FSRs = _extruded_FSR_keys_map.size();
  if (num_extruded_FSRs > 0) {

    /* Allocate extruded FSR lookup vector and fill with extruded FSRs by ID */
    _extruded_FSR_lookup = std::vector<ExtrudedFSR*>(num_extruded_FSRs);
    ExtrudedFSR **extruded_value_list = _extruded_FSR_keys_map.values();
#pragma omp parallel for
    for (int i=0; i < num_extruded_FSRs; i++) {
      long fsr_id = extruded_value_list[i]->_fsr_id;
      _extruded_FSR_lookup[fsr_id] = extruded_value_list[i];
    }

    delete [] extruded_value_list;
  }

  /* Delete key and value lists */
  delete[] key_list;
  delete[] value_list;
}


/**
 * @brief Determines the fissionability of each Universe within this Geometry.
 * @details A Universe is determined fissionable if it contains a Cell
 *          filled by a Material with a non-zero fission cross-section. Note
 *          that this method recurses through all Universes at each level in
 *          the nested Universe hierarchy. Users should only call this method
 *          without a parameter (the default) from Python as follows to ensure
 *          that the recursion starts from the uppermost Universe level:
 *
 * @code
 *          geometry.computeFissionability()
 * @endcode
 *
 * @param univ the Universe of interest (default is NULL)
 */
void Geometry::computeFissionability(Universe* univ) {

  bool fissionable = false;

  Material* material;
  std::map<int, Material*> materials;
  std::map<int, Material*>::iterator mat_iter;

  Universe* universe;
  std::map<int, Universe*> universes;
  std::map<int, Universe*>::iterator univ_iter;

  /* If no Universe was passed in as an argument, then this is the first
   * recursive call from a user via Python, so get the base Universe */
  if (univ == NULL)
    univ = _root_universe;

  /* If a Universe was passed in as an argument, then this is a recursive
   * call with a Universe at a lower level in the nested Universe hierarchy */
  if (univ->getType() == SIMPLE) {
    materials = univ->getAllMaterials();
    universes = univ->getAllUniverses();
  }

  else
    universes = static_cast<Lattice*>(univ)->getAllUniverses();

  /* Loop over the nested Universes first to ensure that fissionability
   * is set at each nested Universe level */
  for (univ_iter=universes.begin(); univ_iter != universes.end(); ++univ_iter) {
    universe = univ_iter->second;

    /* Recursively check whether this nested Universe is fissionable */
    computeFissionability(universe);

    if (universe->isFissionable())
      fissionable = true;
  }

  /* Loop over the Materials in this Universe at this level */
  for (mat_iter=materials.begin(); mat_iter != materials.end(); ++mat_iter) {
    material = mat_iter->second;

    /* Check whether this Material is fissionable or not */
    if (material->isFissionable())
      fissionable = true;
  }

  /* Set this Universe's fissionability based on the nested Universes
   * and Materials within it */
  univ->setFissionability(fissionable);
}


/**
 * @brief Get the material, cell or FSR IDs on a 2D spatial grid.
 * @details This is a helper method for the openmoc.plotter module.
 *          This method may also be called by the user in Python if needed.
 *          A user must initialize NumPy arrays with the x and y grid
 *          coordinates input to this function. This function then fills
 *          a NumPy array with the domain IDs for each coordinate. An example
 *          of how this function might be called in Python is as follows:
 *
 * @code
 *          grid_x = numpy.arange(-2., +2., 100)
 *          grid_y = numpy.arange(-2., +2., 100)
 *          domain_ids = geometry.getSpatialDataOnGrid(
 *              grid_x, grid_y, 20., 'xy', 'material')
 * @endcode
 *
 * @param dim1 a numpy array of the first dimension's coordinates
 * @param dim2 a numpy array of the second dimension's coordinates
 * @param offset The coordinate at which the plane is located
 * @param plane The plane for which data is gathered ('xy', 'xz', 'yz')
 * @param domain_type the type of domain ('fsr', 'material', 'cell')
 * @return a NumPy array or list of the domain IDs
 */
std::vector<long> Geometry::getSpatialDataOnGrid(std::vector<double> dim1,
                                                 std::vector<double> dim2,
                                                 double offset,
                                                 const char* plane,
                                                 const char* domain_type) {

  /* Instantiate a vector to hold the domain IDs */
  std::vector<long> domains(dim1.size() * dim2.size());

  /* Allocate fsr keys in case the track generator has not */
  setNumThreads(omp_get_max_threads());

  /* Extract the source region IDs */
#pragma omp parallel for collapse(2)
  for (int i=0; i < dim1.size(); i++) {
    for (int j=0; j < dim2.size(); j++) {

      Cell* cell;
      LocalCoords point(0, 0, 0, true);

      /* Find the Cell containing this point */
      if (strcmp(plane, "xy") == 0)
        point.getPoint()->setCoords(dim1[i], dim2[j], offset);
      else if (strcmp(plane, "xz") == 0)
        point.getPoint()->setCoords(dim1[i], offset, dim2[j]);
      else if (strcmp(plane, "yz") == 0)
        point.getPoint()->setCoords(offset, dim1[i], dim2[j]);
      else
        log_printf(ERROR, "Unable to extract spatial data for "
                          "unsupported plane %s", plane);

      point.setUniverse(_root_universe);
      cell = _root_universe->findCell(&point);
      domains[i+j*dim1.size()] = -1;

      /* Extract the ID of the domain of interest */
      if (withinGlobalBounds(&point) && cell != NULL) {
        if (strcmp(domain_type, "fsr") == 0)
          domains[i+j*dim1.size()] = getGlobalFSRId(&point, false);
        else if (strcmp(domain_type, "material") == 0)
          domains[i+j*dim1.size()] = cell->getFillMaterial()->getId();
        else if (strcmp(domain_type, "cell") == 0)
          domains[i+j*dim1.size()] = cell->getId();
        else
          log_printf(ERROR, "Unable to extract spatial data for "
                            "unsupported domain type %s", domain_type);
      }
    }
  }

  /* Return the domain IDs */
  return domains;
}


/**
 * @brief Converts this Geometry's attributes to a character array.
 * @details This method calls the toString() method for all Materials,
 *          Surfaces, Cell, Universes and Lattices contained by the Geometry.
 * @return a character array of this Geometry's class attributes
 */
std::string Geometry::toString() {

  std::stringstream string;

  std::map<int, Cell*> all_cells = _root_universe->getAllCells();
  std::map<int, Universe*> all_universes = _root_universe->getAllUniverses();

  std::map<int, Cell*>::iterator cell_iter;
  std::map<int, Universe*>::iterator univ_iter;

  string << "\n\tCells:\n\t\t";
  for (cell_iter = all_cells.begin(); cell_iter != all_cells.end(); ++cell_iter)
    string << cell_iter->second->toString() << "\n\t\t";

  string << "\n\tUniverses:\n\t\t";
  for (univ_iter = all_universes.begin();
       univ_iter != all_universes.end(); ++univ_iter)
    string << univ_iter->second->toString() << "\n\t\t";

  std::string formatted_string = string.str();
  formatted_string.erase(formatted_string.end()-3);

  return formatted_string;
}


/**
 * @brief Prints a string representation of all of the Geometry's attributes to
 *        the console.
 * @details This method calls the printString() method for all Materials,
 *          Surfaces, Cell, Universes and Lattices contained by the Geometry.
 */
void Geometry::printString() {
  log_printf(RESULT, "%s", toString().c_str());
}


/**
 * @brief Prints FSR layout to file
 * @details This provides a way to get the functionality of the
 *              plot_flat_source_regions Python function without Python
 * @param plane The "xy", "xz", or "yz" plane in which to extract flat source
 *        regions
 * @param gridsize The number of points to plot in each direction
 * @param offset The offset of the plane in the third dimension
 * @param bounds_x a two valued array for the plotted x-limits
 * @param bounds_y a two valued array for the plotted y-limits
 * @param bounds_z a two valued array for the plotted z-limits
 */
void Geometry::printFSRsToFile(const char* plane, int gridsize, double offset,
                               double* bounds_x, double* bounds_y,
                               double* bounds_z) {

  /* Get geometry min and max */
  double min_x = _root_universe->getMinX();
  double max_x = _root_universe->getMaxX();
  double min_y = _root_universe->getMinY();
  double max_y = _root_universe->getMaxY();
  double min_z = _root_universe->getMinZ();
  double max_z = _root_universe->getMaxZ();

  if (bounds_x != NULL) {
    min_x = bounds_x[0];
    max_x = bounds_x[1];
  }
  if (bounds_y != NULL) {
    min_y = bounds_y[0];
    max_y = bounds_y[1];
  }
  if (bounds_z != NULL) {
    min_z = bounds_z[0];
    max_z = bounds_z[1];
  }

  /* Create coordinate vectors */
  std::vector<double> dim1(gridsize);
  std::vector<double> dim2(gridsize);

  /* Determine minimum and maximum values */
  double dim1_min = -1;
  double dim1_max = -1;
  double dim2_min = -1;
  double dim2_max = -1;

  /* x-y plane */
  if (strcmp(plane, "xy") == 0) {
    dim1_min = min_x;
    dim1_max = max_x;
    dim2_min = min_y;
    dim2_max = max_y;
  }

  /* x-z plane */
  else if (strcmp(plane, "xz") == 0) {
    dim1_min = min_x;
    dim1_max = max_x;
    dim2_min = min_z;
    dim2_max = max_z;
  }

  /* y-z plane */
  else if (strcmp(plane, "yz") == 0) {
    dim1_min = min_y;
    dim1_max = max_y;
    dim2_min = min_z;
    dim2_max = max_z;
  }

  else {
    log_printf(ERROR, "Plane type %s unrecognized", plane);
  }

  /* Create grid */
  double width1 = (dim1_max - dim1_min) / (gridsize + 1);
  double width2 = (dim2_max - dim2_min) / (gridsize + 1);
  for (int i=0; i < gridsize; i++) {
    dim1.at(i) = dim1_min + (i+1) * width1;
    dim2.at(i) = dim2_min + (i+1) * width2;
  }

  /* Retrieve data */
  log_printf(NORMAL, "Getting FSR layout on domains");
  std::vector<long> domain_data = getSpatialDataOnGrid(dim1, dim2, offset,
                                                       plane, "fsr");

  long* fsr_array = new long[domain_data.size()];
#pragma omp parallel for
  for (int i=0; i < domain_data.size(); i++) {
    fsr_array[i] = domain_data.at(i) + 1;
  }

#ifdef MPIx
  log_printf(NORMAL, "Communicating FSR layout accross domains");
  long* reduced_fsr_array = new long[domain_data.size()];
  MPI_Allreduce(fsr_array, reduced_fsr_array, domain_data.size(),
                MPI_LONG, MPI_SUM, _MPI_cart);
#pragma omp parallel for
  for (int i=0; i < domain_data.size(); i++)
    fsr_array[i] = reduced_fsr_array[i];
  delete [] reduced_fsr_array;
#endif

  /* Print to file */
  log_printf(NORMAL, "Printing FSRs to file");
  if (isRootDomain()) {
    std::ofstream out("fsr-printout.txt");
    out << "[HEADER] FSR printout" << std::endl;
    out << "[HEADER] Plane = " << plane << std::endl;
    out << "[HEADER] Offset = " << offset << std::endl;
    out << "[HEADER] Bounds = (" << dim1_min << ", " << dim1_max << ") x ("
        << dim2_min << ", " << dim2_max << ")" << std::endl;
    out << "[HEADER] Gridsize = " << gridsize << std::endl;
    for (int i=0; i < domain_data.size(); i++) {
      out << fsr_array[i] << " ";
    }
  }

  delete [] fsr_array;
}


/**
 * @brief This is a method that initializes the CMFD Lattice and sets
 *          CMFD parameters.
 */
void Geometry::initializeCmfd() {

  /* Get the global Geometry boundary conditions */
  boundaryType min_x_bound = _root_universe->getMinXBoundaryType();
  boundaryType max_x_bound = _root_universe->getMaxXBoundaryType();
  boundaryType min_y_bound = _root_universe->getMinYBoundaryType();
  boundaryType max_y_bound = _root_universe->getMaxYBoundaryType();
  boundaryType min_z_bound = _root_universe->getMinZBoundaryType();
  boundaryType max_z_bound = _root_universe->getMaxZBoundaryType();

  /* Get the global Geometry boundaries */
  double min_x = _root_universe->getMinX();
  double max_x = _root_universe->getMaxX();
  double min_y = _root_universe->getMinY();
  double max_y = _root_universe->getMaxY();
  double min_z = _root_universe->getMinZ();
  double max_z = _root_universe->getMaxZ();

  /* Set CMFD mesh boundary conditions */
  _cmfd->setBoundary(SURFACE_X_MIN, min_x_bound);
  _cmfd->setBoundary(SURFACE_Y_MIN, min_y_bound);
  _cmfd->setBoundary(SURFACE_Z_MIN, min_z_bound);
  _cmfd->setBoundary(SURFACE_X_MAX, max_x_bound);
  _cmfd->setBoundary(SURFACE_Y_MAX, max_y_bound);
  _cmfd->setBoundary(SURFACE_Z_MAX, max_z_bound);

  /* Set CMFD mesh dimensions */
  _cmfd->setWidthX(max_x - min_x);
  _cmfd->setWidthY(max_y - min_y);
  _cmfd->setWidthZ(max_z - min_z);

  /* Initialize the CMFD lattice */
  Point offset;
  offset.setX(min_x + (max_x - min_x)/2.0);
  offset.setY(min_y + (max_y - min_y)/2.0);
  if (std::abs(min_z + (max_z - min_z)/2.0) < FLT_INFINITY)
    offset.setZ(min_z + (max_z - min_z)/2.0);
  else
    offset.setZ(0.);

  _cmfd->setGeometry(this);
  _cmfd->initializeLattice(&offset);

#ifdef MPIx
  if (_domain_decomposed) {

    /* Check that CMFD mesh is compatible with domain decomposition */
    _cmfd->setNumDomains(_num_domains_x, _num_domains_y, _num_domains_z);
    _cmfd->setDomainIndexes(_domain_index_x, _domain_index_y, _domain_index_z);
  }
#endif
  /* Initialize CMFD Maps */
  _cmfd->initializeCellMap();
}


/**
 * @brief This is a method that initializes the initial spectrum calculator.
 */
void Geometry::initializeSpectrumCalculator(Cmfd* spectrum_calculator) {

  /* Setup the CMFD lattice with the domain dimensions */
  spectrum_calculator->setLatticeStructure(_num_domains_x, _num_domains_y,
                                           _num_domains_z);

  /* Get the global Geometry boundary conditions */
  boundaryType min_x_bound = _root_universe->getMinXBoundaryType();
  boundaryType max_x_bound = _root_universe->getMaxXBoundaryType();
  boundaryType min_y_bound = _root_universe->getMinYBoundaryType();
  boundaryType max_y_bound = _root_universe->getMaxYBoundaryType();
  boundaryType min_z_bound = _root_universe->getMinZBoundaryType();
  boundaryType max_z_bound = _root_universe->getMaxZBoundaryType();

  /* Get the global Geometry boundaries */
  double min_x = _root_universe->getMinX();
  double max_x = _root_universe->getMaxX();
  double min_y = _root_universe->getMinY();
  double max_y = _root_universe->getMaxY();
  double min_z = _root_universe->getMinZ();
  double max_z = _root_universe->getMaxZ();

  /* Set spectrum caclulator boundary conditions */
  spectrum_calculator->setBoundary(SURFACE_X_MIN, min_x_bound);
  spectrum_calculator->setBoundary(SURFACE_Y_MIN, min_y_bound);
  spectrum_calculator->setBoundary(SURFACE_Z_MIN, min_z_bound);
  spectrum_calculator->setBoundary(SURFACE_X_MAX, max_x_bound);
  spectrum_calculator->setBoundary(SURFACE_Y_MAX, max_y_bound);
  spectrum_calculator->setBoundary(SURFACE_Z_MAX, max_z_bound);

  /* Set spectrum calculator dimensions */
  spectrum_calculator->setWidthX(max_x - min_x);
  spectrum_calculator->setWidthY(max_y - min_y);
  spectrum_calculator->setWidthZ(max_z - min_z);

  /* Initialize CMFD Maps */
  spectrum_calculator->initializeCellMap();

  /* Initialize the CMFD lattice */
  Point offset;
  offset.setX(min_x + (max_x - min_x)/2.0);
  offset.setY(min_y + (max_y - min_y)/2.0);
  offset.setZ(min_z + (max_z - min_z)/2.0);

  spectrum_calculator->setGeometry(this);
  spectrum_calculator->initializeLattice(&offset);

#ifdef MPIx
  if (_domain_decomposed) {
    spectrum_calculator->setNumDomains(_num_domains_x, _num_domains_y,
                                       _num_domains_z);
    spectrum_calculator->setDomainIndexes(_domain_index_x, _domain_index_y,
                                          _domain_index_z);
  }
#endif

  /* Add FSRs to domain cell */
  for (long r=0; r < getNumFSRs(); r++)
    spectrum_calculator->addFSRToCell(0, r);
}


/**
 * @brief Returns a pointer to the map that maps FSR keys to FSR IDs.
 * @return pointer to _FSR_keys_map map of FSR keys to FSR IDs
 */
ParallelHashMap<std::string, fsr_data*>& Geometry::getFSRKeysMap() {
  return _FSR_keys_map;
}


/**
 * @brief Returns a pointer to the map that maps FSR keys to extruded FSRs.
 * @return pointer to _FSR_keys_map map of FSR keys to extruded FSRs
 */
ParallelHashMap<std::string, ExtrudedFSR*>& Geometry::getExtrudedFSRKeysMap() {
  return _extruded_FSR_keys_map;
}


/**
 * @brief Returns the vector that maps FSR IDs to FSR key hashes.
 * @return _FSR_keys_map map of FSR keys to FSR IDs
 */
std::vector<std::string>& Geometry::getFSRsToKeys() {
  return _FSRs_to_keys;
}


/**
 * @brief Returns the vector that maps FSR IDs to extruded FSRs.
 * @return _extruded_FSR_lookup map of FSR keys to extruded FSRs
 */
std::vector<ExtrudedFSR*>& Geometry::getExtrudedFSRLookup() {
  return _extruded_FSR_lookup;
}


/**
 * @brief Returns a vector indexed by flat source region IDs which contains
 *        the corresponding Material IDs.
 * @return an integer vector of FSR-to-Material IDs indexed by FSR ID
 */
std::vector<int>& Geometry::getFSRsToMaterialIDs() {
  return _FSRs_to_material_IDs;
}


/**
 * @brief Return a vector indexed by flat source region IDs which contains
 *        pointers to the corresponding Centroid information.
 * @return an array of centroid pointers indexed by FSR ID
 */
std::vector<Point*>& Geometry::getFSRsToCentroids() {
  return _FSRs_to_centroids;
}


/**
 * @brief Return a vector indexed by flat source region IDs which contains
 *        the corresponding CMFD cell.
 * @return an integer vector of FSR to CMFD cell IDs indexed by FSR ID
 */
std::vector<int>& Geometry::getFSRsToCMFDCells() {
  return _FSRs_to_CMFD_cells;
}


/**
 * @brief Determines whether a point is within the bounding box of the domain.
 * @param coords a populated LocalCoords linked list
 * @return boolean indicating whether the coords is within the domain
 */
bool Geometry::withinBounds(LocalCoords* coords) {

  double x = coords->getX();
  double y = coords->getY();
  double z = coords->getZ();

  if (x <= getMinX() || x >= getMaxX() || y <= getMinY() || y >= getMaxY()
      || z <= getMinZ() || z >= getMaxZ())
    return false;
  else
    return true;
}


/**
 * @brief Determines whether a point is within the bounding box of the Geometry.
 * @param coords a populated LocalCoords linked list
 * @return boolean indicating whether the coords is within the geometry
 */
bool Geometry::withinGlobalBounds(LocalCoords* coords) {

  double x = coords->getX();
  double y = coords->getY();
  double z = coords->getZ();

  if (x <= _root_universe->getMinX() || x >= _root_universe->getMaxX() ||
      y <= _root_universe->getMinY() || y >= _root_universe->getMaxY() ||
      z <= _root_universe->getMinZ() || z >= _root_universe->getMaxZ())
    return false;
  else
    return true;
}


/**
 * @brief Finds the Cell containing a given fsr ID.
 * @param fsr_id an FSR ID.
 */
Cell* Geometry::findCellContainingFSR(long fsr_id) {

  std::string& key = _FSRs_to_keys[fsr_id];
  Point* point = _FSR_keys_map.at(key)->_point;
  LocalCoords coords(point->getX(), point->getY(), point->getZ(), true);
  coords.setUniverse(_root_universe);
  Cell* cell = findCellContainingCoords(&coords);

  return cell;
}


/**
 * @brief Sets the centroid for an FSR.
 * @details The _FSR_keys_map stores a hash of a std::string representing
 *          the Lattice/Cell/Universe hierarchy for a unique region
 *          and the associated FSR data. _centroid is a point that represents
 *          the numerical centroid of an FSR computed using all segments
 *          contained in the FSR. This method is used by the TrackGenerator
 *          to set the centroid after segments have been created. It is
 *          important to note that this method is a helper function for the
 *          TrackGenerator and should not be explicitly called by the user.
 * @param fsr a FSR id
 * @param centroid a Point representing the FSR centroid
 */
void Geometry::setFSRCentroid(long fsr, Point* centroid) {
  _contains_FSR_centroids = true;
  std::string& key = _FSRs_to_keys[fsr];
  _FSR_keys_map.at(key)->_centroid = centroid;
  _FSRs_to_centroids[fsr] = centroid;
}


/**
 * @brief Sets the boolean keeping track of FSR centroids generation to false.
 * @details The FSR numbering may change if trying to use the same geometry
 *          in a restart run. The centroids generated with the previous FSR
 *          numbering should not be used when ray tracing again.
 */
void Geometry::resetContainsFSRCentroids() {
  _contains_FSR_centroids = false;
}


/**
 * @brief Returns a vector of z-coords defining a superposition of all axial
 *        boundaries in the Geometry.
 * @details The Geometry is traversed to retrieve all Z-planes and implicit
 *          z-boundaries, such as lattice boundaries. The levels of all these
 *          z-boundaries are rounded and added to a set containing no
 *          duplicates, creating a mesh.
 * @param include_overlaid_mesh whether to include an overlaid mesh in the
 *        set of unique z-coords
 * @return a vector of z-coords
 */
std::vector<double> Geometry::getUniqueZHeights(bool include_overlaid_mesh) {

  /* Get the bounds of the geometry */
  double min_z = getMinZ();
  double max_z = getMaxZ();

  /* Initialize set for axial mesh */
  std::set<double> unique_mesh;
  unique_mesh.insert(min_z);
  unique_mesh.insert(max_z);

  /* Initialize vector of unvisited universes and add the root universe */
  std::vector<Universe*> universes;
  universes.push_back(_root_universe);

  /* Initialize vector of offsets */
  std::vector<double> offsets;
  offsets.push_back(0.0);

  /* Cycle through known universes */
  while (!universes.empty()) {

    /* Get the last universe and explore it */
    Universe* curr_universe = universes.back();
    universes.pop_back();

    /* Get the z-offset of the universe */
    double z_offset = offsets.back();
    offsets.pop_back();

    /* Store a vector of the z_heights before rounding */
    std::vector<double> z_heights;

    /* Check if universe is actually a lattice */
    universeType type = curr_universe->getType();
    if (type == LATTICE) {

      /* Get lattice dimensions */
      Lattice* lattice = static_cast<Lattice*>(curr_universe);
      int nx = lattice->getNumX();
      int ny = lattice->getNumY();
      int nz = lattice->getNumZ();

      /* Calculate z-intersections */
      const std::vector<double>& accumulatez = lattice->getAccumulateZ();
      /* Get offset of the lattice */
      double offset = z_offset + lattice->getMinZ();
      for (int k=0; k<nz+1; k++) {
        double z_height = accumulatez[k] + offset;
        z_heights.push_back(z_height);
      }

      /* Add universes to unvisited universes vector */
      for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
          for (int k=0; k<nz; k++) {
            Universe* new_universe = lattice->getUniverse(i, j, k);
            universes.push_back(new_universe);
            offsets.push_back(z_offset);
          }
        }
      }
    }

    /* Otherwise check if universe is simple, contains cells */
    else if (type == SIMPLE) {

      /* Get all cells in the universe */
      std::map<int, Cell*> cells = curr_universe->getCells();

      /* Cycle through all cells */
      std::map<int, Cell*>::iterator cell_iter;
      for (cell_iter = cells.begin(); cell_iter != cells.end(); ++cell_iter) {

        /* Get surfaces bounding the cell */
        std::map<int, Halfspace*> surfaces =
          cell_iter->second->getSurfaces();

        /* Cycle through all surfaces and add them to the set */
        std::map<int, Halfspace*>::iterator surf_iter;
        for (surf_iter = surfaces.begin(); surf_iter != surfaces.end();
            ++surf_iter) {

          /* Extract surface type */
          Surface* surface = surf_iter->second->_surface;
          surfaceType surf_type = surface->getSurfaceType();

          /* Treat surface types */
          if (surf_type == PLANE) {

            /* Extract plane paramters */
            Plane* plane = static_cast<Plane*>(surface);
            double A = plane->getA();
            double B = plane->getB();
            double C = plane->getC();
            double D = plane->getD();

            /* Check if there is a z-component */
            if (fabs(C) > FLT_EPSILON) {

              /* Check if plane has a continuous varying slope */
              if (fabs(A) > FLT_EPSILON || fabs(B) > FLT_EPSILON)
                log_printf(ERROR, "Continuous axial variation found in the "
                          "Geometry during axial on-the-fly ray tracing. "
                          "Axial on-the-fly ray tracing only supports "
                          "geometries that are capable of having an axially "
                          "extruded representation");

              /* Otherwise, surface is a z-plane */
              else
                z_heights.push_back(-D/C + z_offset);
            }
          }

          /* Treat explicit z-planes */
          else if (surf_type == ZPLANE) {
            ZPlane* zplane = static_cast<ZPlane*>(surface);
            z_heights.push_back(zplane->getZ() + z_offset);
          }
        }

        /* Add min and max z-height to the cell */
        double z_limits[2];
        z_limits[0] = cell_iter->second->getMinZ();
        z_limits[1] = cell_iter->second->getMaxZ();
        for (int i=0; i < 2; i++) {
          if (std::abs(z_limits[i]) != std::numeric_limits<double>::infinity())
            z_heights.push_back(z_limits[i] + z_offset);
        }

        /* See if cell is filled with universes or lattices */
        cellType cell_type = cell_iter->second->getType();
        if (cell_type == FILL) {
          Universe* new_universe = cell_iter->second->getFillUniverse();
          universes.push_back(new_universe);
          offsets.push_back(z_offset);
        }
      }
    }

    /* Add rounded z-heights to the set of unique z-heights */
    for (int i=0; i < z_heights.size(); i++) {

      /* Round z-height */
      z_heights[i] = floor(z_heights[i] / ON_SURFACE_THRESH)
        * ON_SURFACE_THRESH;

      /* Add the rounded z-height to the set */
      if (z_heights[i] > min_z && z_heights[i] < max_z)
        unique_mesh.insert(z_heights[i]);
    }
  }

  /* Include overlaid mesh heights if requested */
  if (include_overlaid_mesh && _overlaid_mesh != NULL) {
    int num_z = _overlaid_mesh->getNumZ();
    double dz = (max_z - min_z) / num_z;
    for (int i=1; i < num_z; i++)
      unique_mesh.insert(min_z + i*dz);
  }

  /* Get a vector of the unique z-heights in the Geometry */
  std::vector<double> unique_heights;
  std::set<double>::iterator iter;
  for (iter = unique_mesh.begin(); iter != unique_mesh.end(); ++iter)
    unique_heights.push_back(static_cast<double>(*iter));

  std::sort(unique_heights.begin(), unique_heights.end());
  return unique_heights;
}


/**
 * @brief Returns a vector of z-coords defining potential unique radial planes
 *        in the Geometry.
 * @details The Geometry is traversed to retrieve all Z-planes and implicit
 *          z-boundaries, such as lattice boundaries. The mid points of this
 *          mesh are then used to construct a vector of all potential unique
 *          radial planes and returned to the user.
 * @return a vector of z-coords
 */
std::vector<double> Geometry::getUniqueZPlanes() {

  /* Get a vector of all unique z-heights in the Geometry */
  std::vector<double> unique_heights = getUniqueZHeights();

  /* Use the midpoints to construct all possible unique radial planes */
  std::vector<double> unique_z_planes;
  for (int i=1; i < unique_heights.size(); i++) {
    double mid = (unique_heights[i-1] + unique_heights[i]) / 2;
    unique_z_planes.push_back(mid);
  }

  /* Output the unique Z heights for debugging */
  std::stringstream string;
  string << unique_heights.size() - 1 << " unique Z domains with bounds: ";
  for (int i=0; i < unique_heights.size(); i++)
    string << unique_heights[i] << " ";
  string << "(cm)";
  log_printf(INFO, "%s", string.str().c_str());

  return unique_z_planes;
}


/**
 * @brief Prints all Geometry and Material details to a Geometry restart file.
 * @param filename The name of the file where the data is printed
 */
void Geometry::dumpToFile(std::string filename) {

  FILE* out;
  out = fopen(filename.c_str(), "w");
  if (out == NULL)
    log_printf(ERROR, "Geometry file %s cannot be written. Wrong folder?",
               &filename[0]);

  /* Print number of energy groups */
  int num_groups = getNumEnergyGroups();
  fwrite(&num_groups, sizeof(int), 1, out);

  /* Print all material information */
  std::map<int, Material*> all_materials = getAllMaterials();
  int num_materials = all_materials.size();
  fwrite(&num_materials, sizeof(int), 1, out);
  std::map<int, Material*>::iterator material_iter;
  for (material_iter = all_materials.begin();
      material_iter != all_materials.end(); ++material_iter) {
    int key = material_iter->first;
    Material* mat = material_iter->second;
    int id = mat->getId();
    char* name = mat->getName();

    /* Print key and general material information */
    fwrite(&key, sizeof(int), 1, out);
    fwrite(&id, sizeof(int), 1, out);
    if (strcmp(name, "") == 0) {
      int length = 0;
      fwrite(&length, sizeof(int), 1, out);
    }
    else {
      int length = std::char_traits<char>::length(name);
      fwrite(&length, sizeof(int), 1, out);
      fwrite(name, sizeof(char), length, out);
    }

    /* Print total cross-section */
    for (int g=0; g < num_groups; g++) {
      double value = mat->getSigmaTByGroup(g+1);
      fwrite(&value, sizeof(double), 1, out);
    }

    /* Print fission cross-section */
    for (int g=0; g < num_groups; g++) {
      double value = mat->getSigmaFByGroup(g+1);
      fwrite(&value, sizeof(double), 1, out);
    }

    /* Print nu * fisison cross-section */
    for (int g=0; g < num_groups; g++) {
      double value = mat->getNuSigmaFByGroup(g+1);
      fwrite(&value, sizeof(double), 1, out);
    }

    /* Print neutron emission spectrum (chi) */
    for (int g=0; g < num_groups; g++) {
      double value = mat->getChiByGroup(g+1);
      fwrite(&value, sizeof(double), 1, out);
    }

    /* Print absorption cross section */
    //NOTE This can be used to transfer another XS, like U238 absorption
    FP_PRECISION* xs = mat->getSigmaA();
    for (int g=0; g < num_groups; g++) {
      double value = xs[g];
      fwrite(&value, sizeof(double), 1, out);
    }

    /* Print scattering cross-section */
    for (int g=0; g < num_groups; g++) {
      for (int gp=0; gp < num_groups; gp++) {
        double value = mat->getSigmaSByGroup(g+1, gp+1);
        fwrite(&value, sizeof(double), 1, out);
      }
    }
  }

  /* Print root universe ID */
  int root_id = _root_universe->getId();
  fwrite(&root_id, sizeof(int), 1, out);

  /* Retrieve all surfaces, cells, and universes */
  std::map<int, Surface*> all_surfaces = getAllSurfaces();
  std::map<int, Cell*> all_cells = _root_universe->getAllCells();
  std::map<int, Universe*> all_universes = _root_universe->getAllUniverses();

  std::map<int, Surface*>::iterator surface_iter;
  std::map<int, Cell*>::iterator cell_iter;
  std::map<int, Universe*>::iterator univ_iter;

  /* Print all surface information */
  int num_surfaces = all_surfaces.size();
  fwrite(&num_surfaces, sizeof(int), 1, out);
  for (surface_iter = all_surfaces.begin(); surface_iter != all_surfaces.end();
      ++surface_iter) {

    /* Get key, value pair and general surface information */
    int key = surface_iter->first;
    Surface* value = surface_iter->second;
    int id = value->getId();
    char* name = value->getName();
    surfaceType st = value->getSurfaceType();
    boundaryType bt = value->getBoundaryType();

    /* Print key and general surface information */
    fwrite(&key, sizeof(int), 1, out);
    fwrite(&id, sizeof(int), 1, out);
    if (strcmp(name, "") == 0) {
      int length = 0;
      fwrite(&length, sizeof(int), 1, out);
    }
    else {
      int length = std::char_traits<char>::length(name);
      fwrite(&length, sizeof(int), 1, out);
      fwrite(name, sizeof(char), length, out);
    }
    fwrite(&st, sizeof(surfaceType), 1, out);
    fwrite(&bt, sizeof(boundaryType), 1, out);

    /* Treat specific surface types */
    if (st == PLANE) {
      Plane* pl = static_cast<Plane*>(value);
      double a = pl->getA();
      double b = pl->getB();
      double c = pl->getC();
      double d = pl->getD();
      fwrite(&a, sizeof(double), 1, out);
      fwrite(&b, sizeof(double), 1, out);
      fwrite(&c, sizeof(double), 1, out);
      fwrite(&d, sizeof(double), 1, out);
    }
    else if (st == ZCYLINDER) {
      ZCylinder* zcyl = static_cast<ZCylinder*>(value);
      double x = zcyl->getX0();
      double y = zcyl->getY0();
      double radius = zcyl->getRadius();
      fwrite(&x, sizeof(double), 1, out);
      fwrite(&y, sizeof(double), 1, out);
      fwrite(&radius, sizeof(double), 1, out);
    }
    else if (st == XPLANE) {
      XPlane* xpl = static_cast<XPlane*>(value);
      double x = xpl->getX();
      fwrite(&x, sizeof(double), 1, out);
    }
    else if (st == YPLANE) {
      YPlane* ypl = static_cast<YPlane*>(value);
      double y = ypl->getY();
      fwrite(&y, sizeof(double), 1, out);
    }
    else if (st == ZPLANE) {
      ZPlane* zpl = static_cast<ZPlane*>(value);
      double z = zpl->getZ();
      fwrite(&z, sizeof(double), 1, out);
    }
    else {
      log_printf(ERROR, "Unsupported surface type for surface ID: %d, name:",
                 " %s", id, name);
    }
  }

  /* Print all cell information */
  int num_cells = all_cells.size();
  fwrite(&num_cells, sizeof(int), 1, out);
  for (cell_iter = all_cells.begin(); cell_iter != all_cells.end();
      ++cell_iter) {

    /* Get key, value pair and general cell information */
    int key = cell_iter->first;
    Cell* cell = cell_iter->second;
    int id = cell->getId();
    char* name = cell->getName();
    cellType ct = cell->getType();

    /* Print key and general cell information */
    fwrite(&key, sizeof(int), 1, out);
    fwrite(&id, sizeof(int), 1, out);
    if (strcmp(name, "") == 0) {
      int length = 0;
      fwrite(&length, sizeof(int), 1, out);
    }
    else {
      int length = std::char_traits<char>::length(name);
      fwrite(&length, sizeof(int), 1, out);
      fwrite(name, sizeof(char), length, out);
    }
    fwrite(&ct, sizeof(cellType), 1, out);

    /* Print cell fill information */
    if (ct == MATERIAL) {
      Material* mat = cell->getFillMaterial();
      int mat_id = mat->getId();
      fwrite(&mat_id, sizeof(int), 1, out);
    }
    else if (ct == FILL) {
      Universe* univ = cell->getFillUniverse();
      int univ_id = univ->getId();
      fwrite(&univ_id, sizeof(int), 1, out);
    }

    /* Print cell rotations */
    bool rot = cell->isRotated();
    fwrite(&rot, sizeof(bool), 1, out);
    if (rot) {
      double rotation[3];
      rotation[0] = cell->getPhi("radians");
      rotation[1] = cell->getTheta("radians");
      rotation[2] = cell->getPsi("radians");
      fwrite(rotation, sizeof(double), 3, out);
    }

    /* Print cell translations */
    bool trans = cell->isTranslated();
    fwrite(&trans, sizeof(bool), 1, out);
    if (trans) {
      double* translation = cell->getTranslation();
      fwrite(translation, sizeof(double), 3, out);
    }

    /* Print ring / sector information */
    int num_rings = cell->getNumRings();
    int num_sectors = cell->getNumSectors();
    fwrite(&num_rings, sizeof(int), 1, out);
    fwrite(&num_sectors, sizeof(int), 1, out);

    /* Print parent cell */
    bool has_parent = cell->hasParent();
    fwrite(&has_parent, sizeof(bool), 1, out);
    if (has_parent) {
      int parent_id = cell->getParent()->getId();
      fwrite(&parent_id, sizeof(int), 1, out);
    }

    /* Print region and halfspaces */
    Region* region = cell->getRegion();
    std::vector<Region*> all_nodes;
    std::vector<Region*>::iterator node_iter;

    if (region != NULL) {
      /* Add local Region, head of the CSG Region tree, to the printed nodes */
      all_nodes = region->getAllNodes();
      all_nodes.insert(all_nodes.begin(), region);
    }

    int num_nodes = all_nodes.size();
    fwrite(&num_nodes, sizeof(int), 1, out);

    for (node_iter = all_nodes.begin(); node_iter != all_nodes.end();
         ++node_iter) {

      int region_type = (*node_iter)->getRegionType();
      int num_subnodes = (*node_iter)->getAllNodes().size();
      fwrite(&region_type, sizeof(int), 1, out);
      fwrite(&num_subnodes, sizeof(int), 1, out);

      if (region_type == HALFSPACE) {
        int surface_id =
             static_cast<Halfspace*>(*node_iter)->getSurface()->getId();
        int halfspace = static_cast<Halfspace*>(*node_iter)->getHalfspace();
        fwrite(&surface_id, sizeof(int), 1, out);
        fwrite(&halfspace, sizeof(int), 1, out);
      }
    }

    //FIXME Print neighbors or decide to re-compute them
  }

  /* Print all universe information */
  int num_universes = all_universes.size();
  fwrite(&num_universes, sizeof(int), 1, out);
  for (univ_iter = all_universes.begin();
       univ_iter != all_universes.end(); ++univ_iter) {

    /* Get key, value pair and general universe information */
    int key = univ_iter->first;
    Universe* universe = univ_iter->second;
    int id = universe->getId();
    char* name = universe->getName();
    universeType ut = universe->getType();

    /* Print key and general universe information */
    fwrite(&key, sizeof(int), 1, out);
    fwrite(&id, sizeof(int), 1, out);
    if (strcmp(name, "") == 0) {
      int length = 0;
      fwrite(&length, sizeof(int), 1, out);
    }
    else {
      int length = std::char_traits<char>::length(name);
      fwrite(&length, sizeof(int), 1, out);
      fwrite(name, sizeof(char), length, out);
    }
    fwrite(&ut, sizeof(universeType), 1, out);

    if (ut == SIMPLE) {
      /* Print all cells in the universe */
      std::map<int, Cell*> cells = universe->getCells();
      int num_universe_cells = cells.size();
      fwrite(&num_universe_cells, sizeof(int), 1, out);
      for (cell_iter = cells.begin(); cell_iter != cells.end(); ++cell_iter) {
        int cell_id = cell_iter->first;
        fwrite(&cell_id, sizeof(int), 1, out);
      }
    }
    else if (ut == LATTICE) {
      /* Print lattice information */
      Lattice* lattice = static_cast<Lattice*>(universe);
      int num_x = lattice->getNumX();
      int num_y = lattice->getNumY();
      int num_z = lattice->getNumZ();
      double width_x = lattice->getWidthX();
      double width_y = lattice->getWidthY();
      double width_z = lattice->getWidthZ();
      double* offset = lattice->getOffset()->getXYZ();
      fwrite(&num_x, sizeof(int), 1, out);
      fwrite(&num_y, sizeof(int), 1, out);
      fwrite(&num_z, sizeof(int), 1, out);
      fwrite(&width_x, sizeof(double), 1, out);
      fwrite(&width_y, sizeof(double), 1, out);
      fwrite(&width_z, sizeof(double), 1, out);
      fwrite(offset, sizeof(double), 3, out);
      bool non_uniform = lattice->getNonUniform();
      fwrite(&non_uniform, sizeof(bool), 1, out);
      if (non_uniform) {
        const std::vector<double> widths_x = lattice->getWidthsX();
        const std::vector<double> widths_y = lattice->getWidthsY();
        const std::vector<double> widths_z = lattice->getWidthsZ();
        fwrite(&widths_x[0], sizeof(double), num_x, out);
        fwrite(&widths_y[0], sizeof(double), num_y, out);
        fwrite(&widths_z[0], sizeof(double), num_z, out);
      }

      /* Get universes */
      Universe* universes[num_x * num_y * num_z];
      for (int i=0; i < num_x; i++) {
        for (int j=0; j < num_y; j++) {
          for (int k =0; k < num_z; k++) {
            int idx = (num_z-1-k) * num_x * num_y + (num_y-1-j) * num_x + i;
            universes[idx] = lattice->getUniverse(i, j, k);
          }
        }
      }
      for (int i=0; i < num_x * num_y * num_z; i++) {
        int universe_id = universes[i]->getId();
        fwrite(&universe_id, sizeof(int), 1, out);
      }
    }
  }

  /* Close the output file */
  fclose(out);
}


/**
 * @brief Loads all Geometry and Material details from a Geometry restart file
 * @param filename The name of the file from which the data is loaded
 * @param twiddle Whether the bytes are inverted (BGQ) or not
 */
void Geometry::loadFromFile(std::string filename, bool twiddle) {

  _twiddle = twiddle;
  _loaded_from_file = true;

  FILE* in;
  in = fopen(filename.c_str(), "r");

  if (_root_universe != NULL)
    delete _root_universe;

  log_printf(NORMAL, "Reading Geometry from %s", &filename[0]);
  if (in == NULL)
    log_printf(ERROR, "Geometry file %s was not found.", &filename[0]);

  std::map<int, Surface*> all_surfaces;
  std::map<int, Cell*> all_cells;
  std::map<int, Universe*> all_universes;
  std::map<int, Material*> all_materials;

  std::map<int, int> fill_cell_universes;
  std::map<int, int> cell_parent;
  std::map<int, int*> lattice_universes;

  /* Read number of energy groups */
  int num_groups;
  int ret = twiddleRead(&num_groups, sizeof(int), 1, in);

  /* Read all material infromation */
  int num_materials;
  ret = twiddleRead(&num_materials, sizeof(int), 1, in);
  for (int i=0; i < num_materials; i++) {

    /* Get key, value pair and cross section information */
    int key, id;
    int length;
    char* str;
    const char* name;
    ret = twiddleRead(&key, sizeof(int), 1, in);
    ret = twiddleRead(&id, sizeof(int), 1, in);
    ret = twiddleRead(&length, sizeof(int), 1, in);
    if (length > 0) {
      str = new char[length+1];
      ret = twiddleRead(str, sizeof(char), length, in);
      str[length] = '\0';
      name = str;
    }
    else {
      name = "";
    }

    /* Create Material */
    all_materials[key] = new Material(id, name);
    if (strcmp(name, "") != 0)
      delete [] name;
    Material* mat = all_materials[key];
    mat->setNumEnergyGroups(num_groups);

    /* Set total cross-section */
    double value;
    for (int g=0; g < num_groups; g++) {
      ret = twiddleRead(&value, sizeof(double), 1, in);
      mat->setSigmaTByGroup(value, g+1);
    }

    /* Set fission cross-section */
    for (int g=0; g < num_groups; g++) {
      ret = twiddleRead(&value, sizeof(double), 1, in);
      mat->setSigmaFByGroup(value, g+1);
    }

    /* Set nu * fisison cross-section */
    for (int g=0; g < num_groups; g++) {
      ret = twiddleRead(&value, sizeof(double), 1, in);
      mat->setNuSigmaFByGroup(value, g+1);
    }

    /* Set neutron emission spectrum (chi) */
    for (int g=0; g < num_groups; g++) {
      ret = twiddleRead(&value, sizeof(double), 1, in);
      mat->setChiByGroup(value, g+1);
    }

    /* Set absorption cross section */
    for (int g=0; g < num_groups; g++) {
      ret = twiddleRead(&value, sizeof(double), 1, in);
      mat->setSigmaAByGroup(value, g+1);
    }

    /* Set scattering cross-section */
    for (int g=0; g < num_groups; g++) {
      for (int gp=0; gp < num_groups; gp++) {
        ret = twiddleRead(&value, sizeof(double), 1, in);
        mat->setSigmaSByGroup(value, g+1, gp+1);
      }
    }
  }

  /* Read root universe ID */
  int root_id;
  ret = twiddleRead(&root_id, sizeof(int), 1, in);

  /* Read all surface information */
  int num_surfaces;
  ret = twiddleRead(&num_surfaces, sizeof(int), 1, in);
  for (int i=0; i < num_surfaces; i++) {

    /* Get key, value pair and general surface information */
    int key, id;
    int length;
    char* str;
    const char* name;
    ret = twiddleRead(&key, sizeof(int), 1, in);
    ret = twiddleRead(&id, sizeof(int), 1, in);
    ret = twiddleRead(&length, sizeof(int), 1, in);
    if (length > 0) {
      str = new char[length+1];
      ret = twiddleRead(str, sizeof(char), length, in);
      str[length] = '\0';
      name = str;
    }
    else {
      name = "";
    }
    surfaceType st;
    boundaryType bt;
    ret = twiddleRead(&st, sizeof(surfaceType), 1, in);
    ret = twiddleRead(&bt, sizeof(boundaryType), 1, in);


    /* Treat specific surface types */
    if (st == PLANE) {
      double a, b, c, d;
      ret = twiddleRead(&a, sizeof(double), 1, in);
      ret = twiddleRead(&b, sizeof(double), 1, in);
      ret = twiddleRead(&c, sizeof(double), 1, in);
      ret = twiddleRead(&d, sizeof(double), 1, in);
      all_surfaces[key] = new Plane(a, b, c, d, id, name);
    }
    else if (st == ZCYLINDER) {
      double x, y, radius;
      ret = twiddleRead(&x, sizeof(double), 1, in);
      ret = twiddleRead(&y, sizeof(double), 1, in);
      ret = twiddleRead(&radius, sizeof(double), 1, in);
      all_surfaces[key] = new ZCylinder(x, y, radius, id, name);
    }
    else if (st == XPLANE) {
      double x;
      ret = twiddleRead(&x, sizeof(double), 1, in);
      all_surfaces[key] = new XPlane(x, id, name);
    }
    else if (st == YPLANE) {
      double y;
      ret = twiddleRead(&y, sizeof(double), 1, in);
      all_surfaces[key] = new YPlane(y, id, name);
    }
    else if (st == ZPLANE) {
      double z;
      ret = twiddleRead(&z, sizeof(double), 1, in);
      all_surfaces[key] = new ZPlane(z, id, name);
    }
    else {
      log_printf(ERROR, "Unsupported surface type %s", name);
    }
    if (strcmp(name, "") != 0)
      delete [] name;

    /* Check that the key and ID match */
    if (key != id) {
      std::string str = all_surfaces[key]->toString();
      log_printf(ERROR, "Surface key %d does not match it's corresponding ID "
                        "%d for surface:\n%s", key, id, str.c_str());
    }

    /* Set boundary */
    all_surfaces[key]->setBoundaryType(bt);
  }

  /* Read all cell information */
  int num_cells;
  ret = twiddleRead(&num_cells, sizeof(int), 1, in);
  for (int i=0; i < num_cells; i++) {

    /* Get key, value pair and general cell information */
    int key, id;
    int length;
    char* str;
    const char* name;
    ret = twiddleRead(&key, sizeof(int), 1, in);
    ret = twiddleRead(&id, sizeof(int), 1, in);
    ret = twiddleRead(&length, sizeof(int), 1, in);
    if (length > 0) {
      str = new char[length+1];
      ret = twiddleRead(str, sizeof(char), length, in);
      str[length] = '\0';
      name = str;
    }
    else {
      name = "";
    }
    cellType ct;
    ret = twiddleRead(&ct, sizeof(cellType), 1, in);

    /* Create the cell */
    all_cells[key] = new Cell(id, name);
    if (strcmp(name, "") != 0)
      delete [] name;

    /* Fill the cell */
    if (ct == MATERIAL) {
      int mat_id;
      ret = twiddleRead(&mat_id, sizeof(int), 1, in);
      all_cells[key]->setFill(all_materials[mat_id]);
    }
    else if (ct == FILL) {
      int univ_id;
      ret = twiddleRead(&univ_id, sizeof(int), 1, in);
      fill_cell_universes[key] = univ_id;
    }

    /* Read cell rotations */
    bool rot;
    ret = twiddleRead(&rot, sizeof(bool), 1, in);
    if (rot) {
      double rotation[3];
      ret = twiddleRead(rotation, sizeof(double), 3, in);
      all_cells[key]->setRotation(rotation, 3, "radians");
    }

    /* Read cell translations */
    bool trans;
    ret = twiddleRead(&trans, sizeof(bool), 1, in);
    if (trans) {
      double translation[3];
      ret = twiddleRead(translation, sizeof(double), 3, in);
      all_cells[key]->setTranslation(translation, 3);
    }

    /* Read ring / sector information */
    int num_rings, num_sectors;
    ret = twiddleRead(&num_rings, sizeof(int), 1, in);
    ret = twiddleRead(&num_sectors, sizeof(int), 1, in);
    all_cells[key]->setNumRings(num_rings);
    all_cells[key]->setNumSectors(num_sectors);

    /* Read parent cell */
    bool has_parent;
    ret = twiddleRead(&has_parent, sizeof(bool), 1, in);
    if (has_parent) {
      int parent_id;
      ret = twiddleRead(&parent_id, sizeof(int), 1, in);
      cell_parent[key] = parent_id;
    }

    /* Read region */
    int num_nodes;
    ret = twiddleRead(&num_nodes, sizeof(int), 1, in);

    /* Vector to store the number of subnodes for each node */
    std::vector<int> i_subnodes;
    std::vector<int>::iterator iter;

    /* This loop on the number of nodes reproduces the CSG tree of the region
       , going down to the leaves to add Halfspaces, and adding logical nodes
       (Intersection, Union, Complement) at the other levels */
    for (int n=0; n < num_nodes; n++) {
      int region_type, last_region_type;
      int num_subnodes;
      twiddleRead(&region_type, sizeof(int), 1, in);
      twiddleRead(&num_subnodes, sizeof(int), 1, in);

      /* Keep number of subnodes at all levels */
      if (num_subnodes > 0)
        i_subnodes.push_back(num_subnodes+1);

      /* Remove zero values from subnode vector, and go up in region tree */
      for (iter=i_subnodes.begin(); iter<i_subnodes.end(); ) {
        if ((*iter) == 0) {
          i_subnodes.erase(iter);
          all_cells[key]->goUpOneRegionLogical();
        }
        else
          iter++;
      }

      /* Add surface */
      if (region_type == HALFSPACE) {
        int surface_id;
        int halfspace;
        ret = twiddleRead(&surface_id, sizeof(int), 1, in);
        ret = twiddleRead(&halfspace, sizeof(int), 1, in);
        all_cells[key]->addSurfaceInRegion(halfspace, all_surfaces[surface_id]);
      }

      /* Add Region logical node which will contain surfaces */
      else
        all_cells[key]->addLogicalNode(region_type);

      /* Remove one from all subnode levels */
      std::for_each(i_subnodes.begin(), i_subnodes.end(), [](int& d) {d-=1;});
    }

    /* Check that the key and ID match */
    if (key != id) {
      std::string str = all_cells[key]->toString();
      log_printf(ERROR, "Cell key %d does not match its corresponding ID "
                        "%d for cell:\n%s", key, id, str.c_str());
    }
  }

  /* Read all universe information */
  int num_universes;
  ret = twiddleRead(&num_universes, sizeof(int), 1, in);
  for (int i=0; i < num_universes; i++) {

    /* Get key, value pair and general universe information */
    int key, id;
    int length;
    char* str;
    const char* name;
    ret = twiddleRead(&key, sizeof(int), 1, in);
    ret = twiddleRead(&id, sizeof(int), 1, in);
    ret = twiddleRead(&length, sizeof(int), 1, in);
    if (length > 0) {
      str = new char[length+1];
      ret = twiddleRead(str, sizeof(char), length, in);
      str[length] = '\0';
      name = str;
    }
    else {
      name = "";
    }
    universeType ut;
    ret = twiddleRead(&ut, sizeof(universeType), 1, in);

    if (ut == SIMPLE) {

      /* Read all cells in the universe */
      all_universes[key] = new Universe(id, name);
      int num_universe_cells;
      ret = twiddleRead(&num_universe_cells, sizeof(int), 1, in);
      for (int c=0; c < num_universe_cells; c++) {
        int cell_id;
        ret = twiddleRead(&cell_id, sizeof(int), 1, in);
        all_universes[key]->addCell(all_cells[cell_id]);
      }
    }
    else if (ut == LATTICE) {

      /* Read lattice information */
      int num_x, num_y, num_z;
      double width_x, width_y, width_z;
      double offset[3];
      ret = twiddleRead(&num_x, sizeof(int), 1, in);
      ret = twiddleRead(&num_y, sizeof(int), 1, in);
      ret = twiddleRead(&num_z, sizeof(int), 1, in);
      ret = twiddleRead(&width_x, sizeof(double), 1, in);
      ret = twiddleRead(&width_y, sizeof(double), 1, in);
      ret = twiddleRead(&width_z, sizeof(double), 1, in);
      ret = twiddleRead(offset, sizeof(double), 3, in);

      std::vector<double> widths_x(num_x), widths_y(num_y), widths_z(num_z);
      bool non_uniform = false;
      ret = twiddleRead(&non_uniform, sizeof(bool), 1, in);

      /* Read widths vectors if the lattice is non-uniform */
      if (non_uniform) {
        ret = twiddleRead(&widths_x[0], sizeof(double), num_x, in);
        ret = twiddleRead(&widths_y[0], sizeof(double), num_y, in);
        ret = twiddleRead(&widths_z[0], sizeof(double), num_z, in);
      }

      /* Create lattice */
      Lattice* new_lattice = new Lattice(id, name);
      all_universes[key] = new_lattice;
      new_lattice->setNumX(num_x);
      new_lattice->setNumY(num_y);
      new_lattice->setNumZ(num_z);
      if (non_uniform) {
        new_lattice->setWidths(widths_x, widths_y, widths_z);
        new_lattice->setWidth(1, 1, 1);
      }
      else
        new_lattice->setWidth(width_x, width_y, width_z);

      new_lattice->setNonUniform(non_uniform);
      new_lattice->setOffset(offset[0], offset[1], offset[2]);

      /* Get universes */
      lattice_universes[key] = new int[num_x*num_y*num_z];
      for (int j=0; j < num_x * num_y * num_z; j++) {
        int universe_id;
        ret = twiddleRead(&universe_id, sizeof(int), 1, in);
        lattice_universes[key][j] = universe_id;
      }
    }
    if (strcmp(name, "") != 0)
      delete [] name;

    /* Check that the key and ID match */
    if (key != id) {
      std::string str = all_universes[key]->toString();
      log_printf(ERROR, "Universe key %d does not match it's corresponding ID "
                        "%d for surface:\n%s", key, id, str.c_str());
    }
  }

  /* Set universe fills in cells */
  std::map<int, int>::iterator id_iter;
  for (id_iter = fill_cell_universes.begin();
       id_iter != fill_cell_universes.end(); ++id_iter)
    all_cells[id_iter->first]->setFill(all_universes[id_iter->second]);

  /* Set parent cells */
  for (id_iter = cell_parent.begin(); id_iter != cell_parent.end(); ++id_iter)
    all_cells[id_iter->first]->setParent(all_cells[id_iter->second]);

  /* Set lattice universes */
  std::map<int, int*>::iterator lattice_iter;
  for (lattice_iter = lattice_universes.begin();
       lattice_iter != lattice_universes.end(); ++lattice_iter) {
    int id = lattice_iter->first;
    int* array = lattice_iter->second;
    Lattice* lattice = static_cast<Lattice*>(all_universes[id]);
    int num_x = lattice->getNumX();
    int num_y = lattice->getNumY();
    int num_z = lattice->getNumZ();
    Universe* universes[num_x * num_y * num_z];
    for (int i=0; i < num_x * num_y * num_z; i++) {
      universes[i] = all_universes[array[i]];
    }
    lattice->setUniverses(num_z, num_y, num_x, universes);
    delete [] lattice_iter->second;
  }

  /* Set root universe */
  _root_universe = all_universes[root_id];

  /* Close the input file */
  fclose(in);

  log_printf(NORMAL, "Read complete");
}


/**
 * @brief Read an integer array from file.
 * @param ptr the integer array to fill with the data read
 * @param size the size of each element to read (here size(int))
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(int* ptr, size_t size, size_t nmemb,
                             FILE* stream) {
  size_t ret = fread(ptr, size, nmemb, stream);
  if (_twiddle)
    for (int i=0; i < nmemb; i++)
      ptr[i] = __builtin_bswap32(ptr[i]);
  return ret;
}


/**
 * @brief Read a boolean array from file.
 * @param ptr the boolean array to fill with the data read
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(bool* ptr, size_t size, size_t nmemb, FILE* stream) {
  size_t ret = fread(ptr, size, nmemb, stream);
  return ret;
}


/**
 * @brief Read an array of universeType from file.
 * @param ptr the array to fill with the data read
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(universeType* ptr, size_t size, size_t nmemb,
                             FILE* stream) {
  size_t ret = fread(ptr, size, nmemb, stream);
  int* arr = reinterpret_cast<int*>(ptr);
  if (_twiddle)
    for (int i=0; i < nmemb; i++)
      arr[i] = __builtin_bswap32(arr[i]);
  return ret;
}


/**
 * @brief Read an array of cellType from file.
 * @param ptr the array to fill with the read data
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(cellType* ptr, size_t size, size_t nmemb,
                             FILE* stream) {
  size_t ret = fread(ptr, size, nmemb, stream);
  int* arr = reinterpret_cast<int*>(ptr);
  if (_twiddle)
    for (int i=0; i < nmemb; i++)
      arr[i] = __builtin_bswap32(arr[i]);
  return ret;
}


/**
 * @brief Read an array of surfaceType from file.
 * @param ptr the array to fill with the data read
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(surfaceType* ptr, size_t size, size_t nmemb,
                             FILE* stream) {
  size_t ret = fread(ptr, size, nmemb, stream);
  int* arr = reinterpret_cast<int*>(ptr);
  if (_twiddle)
    for (int i=0; i < nmemb; i++)
      arr[i] = __builtin_bswap32(arr[i]);
  return ret;
}


/**
 * @brief Read an array of boundaryType from file.
 * @param ptr the array to fill with the data read
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(boundaryType* ptr, size_t size, size_t nmemb,
                             FILE* stream) {
  size_t ret = fread(ptr, size, nmemb, stream);
  int* arr = reinterpret_cast<int*>(ptr);
  if (_twiddle)
    for (int i=0; i < nmemb; i++)
      arr[i] = __builtin_bswap32(arr[i]);
  return ret;
}


/**
 * @brief Read an array of char from file.
 * @param ptr the array to fill with the data read
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(char* ptr, size_t size, size_t nmemb, FILE* stream) {
  size_t ret = fread(ptr, size, nmemb, stream);
  return ret;
}


/**
 * @brief Read an array of double from file.
 * @param ptr the array to fill with the data read
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(double* ptr, size_t size, size_t nmemb, FILE* stream) {
  long* arr = reinterpret_cast<long*>(ptr);
  size_t ret = fread(arr, size, nmemb, stream);
  if (_twiddle)
    for (int i=0; i < nmemb; i++)
      arr[i] = __builtin_bswap64(arr[i]);
  return ret;
}


/**
 * @brief Read an array of long int from file.
 * @param ptr the array to fill with the data read
 * @param size the size of each element to read
 * @param nmemb the number of elements to read
 * @param stream the file to read from
 * @return the return status of the read operation
 */
size_t Geometry::twiddleRead(long* ptr, size_t size, size_t nmemb, FILE* stream) {
  long* arr = ptr;
  size_t ret = fread(arr, size, nmemb, stream);
  if (_twiddle)
    for (int i=0; i < nmemb; i++)
      arr[i] = __builtin_bswap64(arr[i]);
  return ret;
}
