/*
 @file    MCSolver.cpp
 @brief   utility functions for monte carlo neutron simulation
 @author  Luke Eure
 @date    January 12 2016
*/

#include "MCSolver.h"

/*
 @brief   constructor for MCSolver
*/
MCSolver::MCSolver() : Solver(NULL) {}


/*
 @brief   deconstructor for MCSolver
*/
MCSolver::~MCSolver() {}


/*
 @brief   set the geometry and root universe for the solver
 @param   geometry the geometry of the solver
*/
void MCSolver::setGeometry(Geometry* geometry) {
  _geometry = geometry;
  _root_universe = geometry->getRootUniverse();
}

/**
 * @brief Initializes the MCSolver fluxes and tally regions
 * @details The Geometry hierarchy is traversed to find all unique regions
 *      within the Geometry
 */
void MCSolver::initialize() {

  // Retrieve simulation parameters from the Geometry
  _num_FSRs = _geometry->getNumCells();
  _num_groups = _geometry->getNumEnergyGroups();
  _num_materials = _geometry->getNumMaterials();

  // Initialize a vector of universes and cells encountered and their
  // associated levels in the Geometry hierarch
  std::vector<Universe*> universes;
  std::vector<int> levels;
  std::vector<Cell*> cell_list;

  // Track the number of passes through a Lattice position with a hash map
  ParallelHashMap<LocalCoords*, int> lattice_pos;

  // Initialize the encountered universes with the root universe at level 0
  universes.push_back(_geometry->getRootUniverse());
  levels.push_back(0);

  // Track the current level within the Geometry hierarchy and the LocalCoords
  int level = -1;
  LocalCoords* coords = NULL;

  // Iterate through all universes in the Geometry hierarchy
  while (universes.size() > 0) {

    // Pick the next universe from the list
    Universe* universe = universes.back();
    universes.pop_back();

    // Prune LocalCoords as needed when the tree is ascended
    for (int i=levels.back(); i <= level; i++) {
      coords = coords->getPrev();
      if (i == level)
        coords->prune();
    }
    level = levels.back();
    levels.pop_back();

    // Get lattice position if necessary
    if (coords != NULL) {
      if (coords->getType() == LAT) {

        // Get the dimensions of the Lattice
        Lattice* lattice = coords->getLattice();
        int nx = lattice->getNumX();
        int ny = lattice->getNumY();
        int nz = lattice->getNumZ();

        // Use the number of passes thorugh the Lattice to calculate the
        // indexes
        int pos = lattice_pos.at(coords);
        int id = nx*ny*nz - pos - 1;
        int x = id / (ny * nz);
        int y = (id - x*ny*nz) / nz;
        int z = id - x*ny*nz - y*nz;

        // Set the indexes in the lattice and increment the number of passes
        coords->setLatticeX(x);
        coords->setLatticeY(y);
        coords->setLatticeZ(z);
        lattice_pos.update(coords, pos+1);
      }
    }

    // Check if the universe is a lattice or simple
    if (universe->getType() == LATTICE) {

      Lattice* lattice = static_cast<Lattice*>(universe);

      // Create a LocalCoords object with the current Lattice
      LocalCoords* lat_coords = new LocalCoords(0,0,0);
      lat_coords->setLattice(lattice);
      lat_coords->setUniverse(universe);
      lat_coords->setType(LAT);

      // Add the LocalCoords to the linked list
      coords->setNext(lat_coords);
      lat_coords->setPrev(coords);
      coords = lat_coords;
      lattice_pos.insert(coords, 0);

      // Add all universes in the Lattice to the list of encountered universes
      for (int i=0; i < lattice->getNumX(); i++) {
        for (int j=0; j < lattice->getNumY(); j++) {
          for (int k=0; k < lattice->getNumZ(); k++) {
            universes.push_back(lattice->getUniverse(i,j,k));
            levels.push_back(level+1);
            cell_list.push_back(NULL);
          }
        }
      }
    }
    else {

      // Create a new LocalCoords object with the current universe
      LocalCoords* univ_coords = new LocalCoords(0,0,0);
      univ_coords->setUniverse(universe);
      univ_coords->setType(UNIV);

      // Set the cell associated with the universe at the previous level and
      // add the LocalCoords to the linked list
      if (coords != NULL) {
        coords->setCell(cell_list.back());
        cell_list.pop_back();
        coords->setNext(univ_coords);
      }
      univ_coords->setPrev(coords);
      coords = univ_coords;

      // loop over all cells in the universe
      std::map<int, Cell*> cells = universe->getCells();
      std::map<int, Cell*>::iterator iter;
      for (iter = cells.begin(); iter != cells.end(); iter++) {

        // Extract the cell and its type
        Cell* cell = iter->second;
        cellType type = cell->getType();
        coords->setCell(cell);

        // Check cell types
        if (type == MATERIAL) {
          // Add the unique region to the FSR map
          _geometry->findFSRId(coords->getHighestLevel());
        }
        else if (type == FILL) {
          // Add nested universes to the list of encountered universes
          universes.push_back(cell->getFillUniverse());
          cell_list.push_back(cell);
          levels.push_back(level+1);
        }
      }
    }
  }

  // Prune the LocalCoords to
  LocalCoords* prev_coords = coords;
  while (prev_coords != NULL) {
    coords = prev_coords;
    prev_coords = coords->getPrev();
  }
  coords->prune();

  // Initialize the FSR and flux arrays
  initializeFluxArrays();
  _geometry->initializeFSRVectors();
  _scalar_flux = new FP_PRECISION[_geometry->getNumFSRs() *
                  _geometry->getNumEnergyGroups()];
}


/*
 @brief   generates and transports neutron histories, calculates the mean
      crow distance
 @param   n_histories number of neutron histories to run
 @param   mat a Material object containing information
      about the material
 @param   flux a Flux object containing information about the flux
 @param   num_batches the number of batches to be tested
 @param   num_groups the number of neutron energy groups
*/
void MCSolver::computeEigenvalue(int n_histories, int num_batches,
                                 int num_groups) {

  // initialize fsid
  _geometry->initializeFSRs();

  // create arrays for tallies and fissions
  std::vector <Tally> tallies(5);
  Fission fission_banks;

  bool first_round = true;

  for (int batch=1; batch <= num_batches; ++batch) {

    // clear flux data
    for (int i=0; i<_num_groups; ++i) {
      for (int j=0; j<_geometry->getNumFSRs(); ++j) {
        _scalar_flux(j,i) = 0.0;
      }
    }

    // assign new fission locations to old fission locations
    fission_banks.newBatch();

    // clear tallies for leaks absorptions and fissions
    tallies[LEAKS].clear();
    tallies[ABSORPTIONS].clear();
    tallies[FISSIONS].clear();

    // simulate neutron behavior
    for (int i=0; i<n_histories; ++i) {
      
      /*
      trackSingleNeutron();
      std::cout << "finishes neutron\n";
      exit(0);
      */ 
      i = 6039;
      std::cout << "Trans " << batch << " " << i << std::endl;
      transportNeutron(tallies, first_round, &fission_banks,
          num_groups, i);// batch, 5383, 2, false);
    }

    // give results
    _k_eff = tallies[FISSIONS].getCount() /
      (tallies[LEAKS].getCount() + tallies[ABSORPTIONS].getCount());

    double sumStandardDev = tallies[LEAKS].getStandardDeviation()
      + tallies[FISSIONS].getStandardDeviation()
      + tallies[ABSORPTIONS].getStandardDeviation();
    std::cout << "\nFor batch " << batch << ", k = " << _k_eff
      << " with standard deviation of " << sumStandardDev << std::endl;

    double fissions;
    fissions = tallies[FISSIONS].getCount();
    std::cout << "fissions: " << fissions/fissions << std::endl;
    std::cout << "leaks: " << tallies[LEAKS].getCount()/fissions
      << std::endl;
    std::cout << "absorptions: " << tallies[ABSORPTIONS].getCount()/fissions
      << std::endl;

    first_round = false;
  }
  double mean_crow_distance = tallies[CROWS].getCount()
    / tallies[NUM_CROWS].getCount();
  std::cout << "Mean crow fly distance = " << mean_crow_distance << std::endl;

}


/*
 @brief   function that generates a neutron and measures how
          far it travels before being absorbed.
 @details   A neutron is created in the bounding box using
            sample_location() for the first batch and sample_fission_site()
            for the rest of the batches. It moves a distance determined by
            sample_distance(). It is then either absorbed or
            scattered as determined by sample_interaction(). When
            it is absorbed, its distance from its starting point
            is appended to crow_distances. If the absorption creates
            a fission event, the number of neutrons emited is sampled.
            The location of the fission event is added to a list of fission
            events.
            of the bounding box
 @param   tallies a dictionary containing tallies of crow distances,
          leakages, absorptions, and fissions
 @param   flux a Flux object containing information about the flux
 @param   old_fission_banks containing the old fission bank
 @param   new_fission_banks containing the new fission bank
 @param   num_groups the number of neutron energy groups
 @param   neutron_num an int used for indexing neutrons
*/
void MCSolver::transportNeutron(std::vector <Tally> &tallies,
    bool first_round, Fission* fission_banks, int num_groups, int neutron_num) {

  const double BOUNDARY_ERROR = 1e-8;

  Point* neutron_position = new Point();

  // new way to sample neutron and set its direction
  Neutron neutron(neutron_num);
  neutron.sampleDirection();

  // get and set neutron starting poinit
  if (first_round)
    sampleLocation(&neutron);
  else {
    fission_banks->sampleSite(&neutron);
  }

  // set neutron_position pointer
  neutron.getPositionVector(neutron_position);

  // set neutron group
  Material* cell_mat;
  Cell* cell_obj;
  int group;

  // get cell material
  LocalCoords* neutron_coord_position = new LocalCoords(
      neutron_position->getX(), neutron_position->getY(),
      neutron_position->getZ());
  neutron_coord_position->setUniverse(_root_universe);
  cell_obj = _geometry->findCellContainingCoords(neutron_coord_position);
  cell_mat = cell_obj->getFillMaterial();

  //std::cout << "neutron start position " << neutron_coord_position->toString()
   // << std::endl;

  std::vector <double> chi(num_groups);
  for (int g=0; g<num_groups; ++g) {
    chi[g] = cell_mat->getChiByGroup(g+1);
  }
  group = neutron.sampleNeutronEnergyGroup(chi);
  neutron.setGroup(group);

  // follow neutron while it's alive
  while (neutron.alive()) {

    // use a LocalCoords for the start and end of each segment
    double neutron_distance;

    //sample a distance to travel in the material
    neutron_distance =
      -log(neutron.arand()) / cell_mat->getSigmaTByGroup(group+1);

    // track neutron until collision or escape
    while (neutron_distance > 0.00000) {

      double x0 = neutron.getPosition(0);
      double y0 = neutron.getPosition(1);
      double z0 = neutron.getPosition(2);

      // calculate phi and ensure it's in the correct quadrant
      double phi = atan(neutron.getDirection(1)/neutron.getDirection(0));
      if (neutron.getDirection(0)<0)
        phi += M_PI;
      if (phi < 0)
        phi += 2*M_PI;

      LocalCoords start(x0, y0, z0);
      LocalCoords end(x0, y0, z0);
      start.setUniverse(_root_universe);
      end.setUniverse(_root_universe);
      start.setPhi(phi);
      end.setPhi(phi);

      // find the Cell containing the Track starting Point
      Cell* curr = _geometry->findFirstCell(&end);
      Cell* prev;

      neutron_coord_position->setX(neutron.getPosition(0));
      neutron_coord_position->setY(neutron.getPosition(1));
      neutron_coord_position->setZ(neutron.getPosition(2));
      neutron_coord_position->setUniverse(_root_universe);
      _geometry->findFirstCell(neutron_coord_position);

        //std::cout << "\n\nstart position: " << neutron_position->getX()
          //<< " " << neutron_position->getY() << std::endl;

      // if starting Point was outside the bounds of the Geometry
      if (curr == NULL)
        log_printf(ERROR,
            "Could not find a material-filled Cell containing"
            " the start Point of this Track: %s");

      // While the end of the segment's LocalCoords is still
      // within the Geometry move it to the next Cell, create
      // a new segment, and add it to the Geometry
      while (curr != NULL) {

        end.copyCoords(&start);
        end.setPhi(phi);

        // Find the next Cell along the Track's trajectory
        prev = curr;

        if (neutron_num == 193) {
          std::cout << "pre-position: " << start.getX() << " " << start.getY()
            << "\n\n";
          std::cout << "phi: " << phi << std::endl << std::endl;
        }


        curr = _geometry->findNextCell(&end);
          std::cout << "phi: " << phi << std::endl;
          std::cout << "prev: " << prev->getName() << std::endl;
          if (curr != NULL)
          std::cout << "curr: " << curr->getName() << std::endl;
          std::cout << "position: " << start.getX() << " " << start.getY()
            << "\n\n";

        // Checks that segment does not have the same
        // start and end Points
        if (start.getX() == end.getX()
            && start.getY() == end.getY())
          log_printf(ERROR,
              "Created segment with same start and end point:"
              " x = %f, y = %f",start.getX(), start.getY());

        // Find the segment length, Material and
        // FSR ID of the last cell
        FP_PRECISION length;
        length = 
          FP_PRECISION(end.getPoint()->distanceToPoint(start.getPoint()));
        int fsr_id = _geometry->findFSRId(&start);
        if (fsr_id > 80) {
          std::cout << "ERROR " << std::endl;
          exit(1);
        }

        // if neutron's path doesn't end in this cell
        if (length < neutron_distance) {

          // add distance travelled to flux, shorten distance and
          // move neutron
          _scalar_flux(fsr_id, neutron.getGroup()) += length;
          neutron_distance -= length;
          neutron.setPosition(0, end.getX());
          neutron.setPosition(1, end.getY());
          neutron.setPosition(2, end.getZ());
          std::cout << "continue neutron position after move "
            << neutron.getPosition(0)
            << " " << neutron.getPosition(1) << std::endl;
        }

        // if the neutron's path ends in this cell
        else {
          _scalar_flux(fsr_id, neutron.getGroup())
            += neutron_distance;
          neutron.move(neutron_distance);
          std::cout << " end in cell neutron position after move "
            << neutron.getPosition(0)
            << " " << neutron.getPosition(1) << std::endl;
          neutron_distance = 0;

          // break out of while (curr!= null)
          curr = NULL;
        }
      } // while curr != null

      neutron.move(-TINY_MOVE);
          std::cout << "continue neutron position after tiny move "
            << neutron.getPosition(0)
            << " " << neutron.getPosition(1) << std::endl;

      // find out which boundary the neutron is on
      std::vector <int> box_lim_bound;
      box_lim_bound.clear();
      if (std::abs(neutron.getPosition(0)
            - _geometry->getMinX()) < BOUNDARY_ERROR)
        box_lim_bound.push_back(0);
      if (std::abs(neutron.getPosition(0)
            - _geometry->getMaxX()) < BOUNDARY_ERROR)
        box_lim_bound.push_back(1);
      if (std::abs(neutron.getPosition(1)
            - _geometry->getMinY()) < BOUNDARY_ERROR)
        box_lim_bound.push_back(2);
      if (std::abs(neutron.getPosition(1)
            - _geometry->getMaxY()) < BOUNDARY_ERROR)
        box_lim_bound.push_back(3);


      // put boundaryTypes into vector for easy iteration
      std::vector <boundaryType> bound_types (4);
      bound_types[0] = _geometry->getMinXBoundaryType();
      bound_types[1] = _geometry->getMaxXBoundaryType();
      bound_types[2] = _geometry->getMinYBoundaryType();
      bound_types[3] = _geometry->getMaxYBoundaryType();

      // put boundary locations into vector for easy iteration
      std::vector <double> bound_locations (4);
      bound_locations[0] = _geometry->getMinX();
      bound_locations[1] = _geometry->getMaxX();
      bound_locations[2] = _geometry->getMinY();
      bound_locations[3] = _geometry->getMaxY();

      // check boundary conditions on all hit surfaces
      for (int sur_side=0; sur_side <4; ++sur_side) {
        int axis = sur_side/2;
        int side = sur_side%2;

        // if sur_side is in box_lim_bound
        if (std::find(box_lim_bound.begin(),
              box_lim_bound.end(),sur_side)
            != box_lim_bound.end()) {

          // if the neutron is reflected
          if (bound_types[sur_side] == 1) {
            neutron.reflect(axis);
          }

          // if the neutron escapes
          if (bound_types[sur_side] == 0) {
            neutron.kill();
            neutron_distance = 0;
            tallies[LEAKS] += 1;
            break;
          }
        }
      } // check boundary conditions
    } // while distance > 0

    // check interaction
    if (neutron.alive()) {
      neutron_coord_position->setX(neutron_position->getX());
      neutron_coord_position->setY(neutron_position->getY());
      neutron_coord_position->setZ(neutron_position->getZ());
      neutron_coord_position->setUniverse(_root_universe);
      cell_obj =
        _geometry->findCellContainingCoords(neutron_coord_position);
      cell_mat = cell_obj->getFillMaterial();

      // calculate sigma_s for a group in order to sample an interaction
      std::vector <double> sigma_s_group;
      double sum_sigma_s_group=0;
      for (int g=1; g<=cell_mat->getNumEnergyGroups(); ++g) {
        sigma_s_group.push_back(cell_mat->getSigmaSByGroup(group+1, g));
        sum_sigma_s_group += cell_mat->getSigmaSByGroup(group+1, g);
      }

      // calculate sigma_a in order to sample an interaction
      double sigma_a;
      sigma_a =
        cell_mat->getSigmaTByGroup(group+1) - sum_sigma_s_group;

      // sample an interaction
      int neutron_interaction =
        (int) (neutron.arand() < (sigma_a
              / cell_mat->getSigmaTByGroup(group+1)));

      // scattering event
      if (neutron_interaction == 0) {

        // sample scattered direction
        neutron.sampleDirection();

        // sample new energy group
        int new_group = neutron.sampleScatteredGroup(sigma_s_group,
          group);

        // set new group
        neutron.setGroup(new_group);

        // findfirstcell
        neutron_coord_position->setX(neutron.getPosition(0));
        neutron_coord_position->setY(neutron.getPosition(1));
        neutron_coord_position->setZ(neutron.getPosition(2));
        neutron_coord_position->setUniverse(_root_universe);
        _geometry->findFirstCell(neutron_coord_position);
        
      }

      // absorption event
      else {

        // tally absorption
        tallies[ABSORPTIONS] += 1;

        // sample for fission event
        group = neutron.getGroup();
        neutron_coord_position->setX(neutron_position->getX());
        neutron_coord_position->setY(neutron_position->getY());
        neutron_coord_position->setZ(neutron_position->getZ());
        neutron_coord_position->setUniverse(_root_universe);
        cell_obj =
          _geometry->findCellContainingCoords(neutron_coord_position);
        cell_mat = cell_obj->getFillMaterial();
        neutron.getPositionVector(neutron_position);

        // sample whether fission event occurs
        int fission_occurs =
          neutron.arand() < cell_mat->getSigmaFByGroup(group+1)
          / sigma_a;

        // fission event
        if (fission_occurs == 1) {

          // sample number of neutrons released during fission
          double nu = cell_mat->getNuSigmaFByGroup(1)
            / cell_mat->getSigmaFByGroup(1);
          int lower = (int) nu;
          int add = (int) (neutron.arand() < nu -lower);
          int num_neutrons = lower + add;
          for (int i=0; i<num_neutrons; ++i) {
            fission_banks->add(neutron_position);
            tallies[FISSIONS] += 1;
          }
        }

        // end neutron history
        neutron.kill();

      } // absorption event
    } // if neutron alive
  } // while neutron alive

  delete neutron_coord_position;

  // tally crow distance
  double crow_distance;
  crow_distance = neutron.getDistance(neutron_position);
  tallies[CROWS] += crow_distance;
  tallies[NUM_CROWS] += 1;
}


/*
 @brief   returns k, the ration of fission events to absoprtions plus leakages
 @return  the effective k-value
*/
FP_PRECISION MCSolver::getKeff() {
  return _k_eff;
}


/*
 @brief   returns the geometry used by the solver
 @return  the geometry used by the solver
*/
Geometry* MCSolver::getGeometry() {
  return _geometry;
}


/*
 @brief   returns the flux value for neutrons of a given group in a fsr
 @param   fsr_id the id number of the fsr for which to find the flux
 @param   group the energy group for which to find the flux
 @return  the neutron flux
*/
FP_PRECISION MCSolver::getFlux(int fsr_id, int group) {
  return _scalar_flux(_geometry->getNumFSRs(), _num_groups);
}


// functions that allow MCSolver to be compativle with Solver
void MCSolver::computeFSRFissionRates(
    double* fission_rates, int num_FSRs) {}
void MCSolver::getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) {}
void MCSolver::setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) {}
void MCSolver::flattenFSRFluxes(FP_PRECISION value) {}
void MCSolver::computeFSRFissionSources() {}
void MCSolver::computeFSRScatterSources() {}
void MCSolver::initializeSourceArrays() {}
void MCSolver::addSourceToScalarFlux() {}
void MCSolver::initializeFluxArrays() {}
void MCSolver::computeFSRSources() {}
void MCSolver::normalizeFluxes() {}
void MCSolver::zeroTrackFluxes() {}
void MCSolver::transportSweep() {}
void MCSolver::storeFSRFluxes() {}
void MCSolver::computeKeff() {}
double MCSolver::computeResidual(residualType res_type){
  return 0;
}


/*
 @brief initialize the fsrs
*/
void MCSolver::initializeFSRs(Lattice* lattice) {
}


/*
 @brief   function that samples a random location within the root cell.
 @details   a position along each axis is randomly and uniformally sampled in
      the bounding box. Each positin is set as the neutron's position
      along that axis.
 @param   neutron, the Neutron whose position will be set
*/
void MCSolver::sampleLocation(Neutron* neutron) {
  Point sampled_location;
  double width = _geometry->getMaxX() - _geometry->getMinX();
  double coord = _geometry->getMinX() + width * neutron->arand();
  neutron->setPosition(0, coord);

  width = _geometry->getMaxY() - _geometry->getMinY();
  coord = _geometry->getMinY() + width * neutron->arand();
  neutron->setPosition(1, coord);
}


/*
 @brief   function that transports a neutron using geometry.segmentize()
 @param   tallies a dictionary containing tallies of crow distances,
          leakages, absorptions, and fissions
 @param   first_round a bool that is true if this is the first round
 @param   fission_banks containing the fission bank
 @param   num_groups the number of neutron energy groups
 @param   neutron_num an int used for indexing neutrons
 @param   batch the batch number
 @param   input_neutron to be used in debuggin bad neutrons
*/
void MCSolver::transportNeutronWithTrack(std::vector <Tally> &tallies,
    bool first_round, Fission* fission_banks, int num_groups, int neutron_num,
    int batch, int write_neutron, int write_batch,
     bool read, Neutron* input_neutron) {

  std::cout.precision(32);
  Point* neutron_position = new Point();

  // new way to sample neutron and set its direction
  Neutron neutron(neutron_num);

  bool debug = (input_neutron != NULL);

  // if not debugging
  if (!debug) {

    neutron.sampleDirection();

    // get and set neutron starting point
    if (first_round)
      sampleLocation(&neutron);
    else {
      fission_banks->sampleSite(&neutron);
    }
  }

  // if debugging
  else {
    std::cout << "gets to else\n";
    neutron.setDirection(0, input_neutron->getDirection(0));
    neutron.setDirection(1, input_neutron->getDirection(1));
    neutron.setDirection(2, input_neutron->getDirection(2));
    neutron.setPosition(0, input_neutron->getPosition(0));
    neutron.setPosition(1, input_neutron->getPosition(1));
    neutron.setPosition(2, input_neutron->getPosition(2));
    neutron.setGroup(input_neutron->getGroup());
    std::cout << "sets attributes\n";
  }

  if (neutron_num == write_neutron and batch == write_batch){
    std::cout << "writing neutron " << neutron_num << " in batch " << batch
      << std::endl;
      saveBadNeutron(&neutron, neutron_num, batch);
      std::cout << "saves neutron\n";
    }


  // set neutron_position pointer
  neutron.getPositionVector(neutron_position);

  // set neutron group
  Material* cell_mat;
  Cell* cell_obj;
  int group;

  // get cell material
  LocalCoords* neutron_coord_position = new LocalCoords(
      neutron_position->getX(), neutron_position->getY(),
      neutron_position->getZ());
  neutron_coord_position->setUniverse(_root_universe);
  cell_obj = _geometry->findCellContainingCoords(neutron_coord_position);
  cell_mat = cell_obj->getFillMaterial();

  std::vector <double> chi(num_groups);
  for (int g=0; g<num_groups; ++g) {
    chi[g] = cell_mat->getChiByGroup(g+1);
  }
  group = neutron.sampleNeutronEnergyGroup(chi);
  neutron.setGroup(group);

  // follow neutron while it's alive
  while (neutron.alive()) {

    // use a LocalCoords for the start and end of each segment
    double neutron_distance;

    //sample a distance to travel in the material
    neutron_distance =
      -log(neutron.arand()) / cell_mat->getSigmaTByGroup(group+1);

    bool first_boundary_check = true;
    bool completed_track;

    // track neutron until collision or escape
    while (neutron_distance > 0.00000) {

      // check if the neutron is on a boundary
      std::vector <int> box_lim_bound;
      box_lim_bound.clear();
      if (std::abs(neutron.getPosition(0) - _geometry->getMinX())
          < ON_SURFACE_THRESH)
        box_lim_bound.push_back(0);
      if (std::abs(neutron.getPosition(0) - _geometry->getMaxX())
          < ON_SURFACE_THRESH)
        box_lim_bound.push_back(1);
      if (std::abs(neutron.getPosition(1) - _geometry->getMinY())
          < ON_SURFACE_THRESH)
        box_lim_bound.push_back(2);
      if (std::abs(neutron.getPosition(1) - _geometry->getMaxY())
          < ON_SURFACE_THRESH)
        box_lim_bound.push_back(3);
  
      if (!first_boundary_check) {
        // throw error if neutron did not complete segment
        if (box_lim_bound.empty() and completed_track)
          log_printf(ERROR,
              "track complete but no surface found %s");
        if (!box_lim_bound.empty() and !completed_track)
          log_printf(ERROR,
              "surface found but track not completed %s");
      }
      first_boundary_check = false;

      // put boundaryTypes into vector for easy iteration
      std::vector <boundaryType> bound_types (4);
      bound_types[0] = _geometry->getMinXBoundaryType();
      bound_types[1] = _geometry->getMaxXBoundaryType();
      bound_types[2] = _geometry->getMinYBoundaryType();
      bound_types[3] = _geometry->getMaxYBoundaryType();

      // put boundary locations into vector for easy iteration
      std::vector <double> bound_locations (4);
      bound_locations[0] = _geometry->getMinX();
      bound_locations[1] = _geometry->getMaxX();
      bound_locations[2] = _geometry->getMinY();
      bound_locations[3] = _geometry->getMaxY();

      // check boundary conditions on all hit surfaces
      for (int sur_side=0; sur_side <4; ++sur_side) {
        int axis = sur_side/2;
        int side = sur_side%2;

        // if sur_side is in box_lim_bound
        if (std::find(box_lim_bound.begin(), box_lim_bound.end(),sur_side)
            != box_lim_bound.end()) {

          // if the neutron is reflected
          if (bound_types[sur_side] == 1) {

            neutron.reflect(axis);

            //neutron.move(TINY_MOVE);

            neutron_coord_position->setX(neutron.getPosition(0));
            neutron_coord_position->setY(neutron.getPosition(1));
            neutron_coord_position->setZ(neutron.getPosition(2));
            neutron_coord_position->setUniverse(_root_universe);
            Cell* test = _geometry->findFirstCell(neutron_coord_position);
  
          }

          // if the neutron escapes
          if (bound_types[sur_side] == 0) {
            neutron.kill();
            neutron_distance = 0;
            tallies[LEAKS] += 1;
            break;
          }
        }
      } // check boundary conditions
      

      double x0 = neutron.getPosition(0);
      double y0 = neutron.getPosition(1);
      double z0 = neutron.getPosition(2);

      // calculate phi and ensure it's in the correct quadrant
      double phi = atan(neutron.getDirection(1)/neutron.getDirection(0));
      if (neutron.getDirection(0)<0)
        phi += M_PI;
      if (phi < 0)
        phi += 2*M_PI;

      Track* track = new Track();
      track->setValues(x0, y0, 0, 0, 0, 0, phi);
      _geometry->segmentize(track);
      segment* segments = track->getSegments();

if (debug) { 
      std::cout << "--------------------------------\n";
      std::cout << "x0 yo zo: " << x0 << " " << y0 << " " << z0 << std::endl;
      std::cout << "phi " << phi << std::endl;
      std::cout << "neutron direction: " << neutron.getDirection(0) << " "
        << neutron.getDirection(1) << " "
        << neutron.getDirection(2) << std::endl;
      std::cout << "segment lengths: " << segments[0]._length
        << ",  " << segments[1]._length << " " << segments[2]._length
        << std::endl;
      std::cout << "num segs " << track->getNumSegments() << std::endl;
      std::cout << "track: " << track->toString() << std::endl;
    }

      double length;
      double length_3d; // length in three dimensions
      completed_track = false; // did neutron complete segment?
      int fsr_id;

      for (int i=0; i < track->getNumSegments(); i++) {

        length = segments[i]._length;
        double denom = sqrt(neutron.getDirection(0)*neutron.getDirection(0)
                            + neutron.getDirection(1)*neutron.getDirection(1));

        length_3d = length/denom;

        fsr_id = segments[i]._region_id;

        if (not neutron_distance > 0.00000)
          break;

        // if the neutron passes through this region
        if (length_3d - neutron_distance < TINY_MOVE) {
          neutron.move(length_3d);
          _scalar_flux(fsr_id, neutron.getGroup()) += length_3d;
          neutron_distance -= length_3d;
/*          std::cout << "\n1 continues through cell " << neutron.getPosition(0) 
            << " " << neutron.getPosition(1) << std::endl;
         std::cout << "neutron radius: "
            << sqrt(neutron.getPosition(0)*neutron.getPosition(0)
                    + neutron.getPosition(1)*neutron.getPosition(1))
            << std::endl;
   */       
          // if the neutron has completed this segment
          if (i == track->getNumSegments() - 1) {
            completed_track = true;
          }
        }

        // if the neutron's path ends in this region
        else {
          neutron.move(neutron_distance);
          _scalar_flux(fsr_id, neutron.getGroup()) += neutron_distance;
          neutron_distance = 0;
          break;
        }
      } // for the Track

      if (completed_track) {
//        std::cout << "neutron has completed the segment" << std::endl;
        //neutron.move(-TINY_MOVE);
      }

    } // while distance > 0

    // check interaction
    if (neutron.alive()) {
      neutron_coord_position->setX(neutron_position->getX());
      neutron_coord_position->setY(neutron_position->getY());
      neutron_coord_position->setZ(neutron_position->getZ());
      neutron_coord_position->setUniverse(_root_universe);
      cell_obj = _geometry->findCellContainingCoords(neutron_coord_position);
      cell_mat = cell_obj->getFillMaterial();

      // calculate sigma_s for a group in order to sample an interaction
      std::vector <double> sigma_s_group;
      double sum_sigma_s_group=0;
      for (int g=1; g<=cell_mat->getNumEnergyGroups(); ++g) {
        sigma_s_group.push_back(cell_mat->getSigmaSByGroup(group+1, g));
        sum_sigma_s_group += cell_mat->getSigmaSByGroup(group+1, g);
      }

      // calculate sigma_a in order to sample an interaction
      double sigma_a = cell_mat->getSigmaTByGroup(group+1) - sum_sigma_s_group;

      // sample an interaction
      int neutron_interaction = (int) (neutron.arand() < (sigma_a
              / cell_mat->getSigmaTByGroup(group+1)));

      // scattering event
      if (neutron_interaction == 0) {

        // sample scattered direction
        neutron.sampleDirection();

        // sample new energy group
        int new_group = neutron.sampleScatteredGroup(sigma_s_group, group);

        // set new group
        neutron.setGroup(new_group);

        // findfirstcell
        neutron_coord_position->setX(neutron.getPosition(0));
        neutron_coord_position->setY(neutron.getPosition(1));
        neutron_coord_position->setZ(neutron.getPosition(2));
        neutron_coord_position->setUniverse(_root_universe);
        Cell* test = _geometry->findFirstCell(neutron_coord_position);

        // if findFirstCell nudged the neutron out of the boundary nudge it
        // backwards then use findFirstCell again
        if (test == NULL) {
          neutron.move(-TINY_MOVE);
          neutron_coord_position->setX(neutron.getPosition(0));
          neutron_coord_position->setY(neutron.getPosition(1));
          neutron_coord_position->setZ(neutron.getPosition(2));
          neutron_coord_position->setUniverse(_root_universe);
          _geometry->findFirstCell(neutron_coord_position);
          neutron.move(TINY_MOVE);
        }
      }

      // absorption event
      else {

        // tally absorption
        tallies[ABSORPTIONS] += 1;

        // sample for fission event
        group = neutron.getGroup();
        neutron_coord_position->setX(neutron.getPosition(0));
        neutron_coord_position->setY(neutron.getPosition(1));
        neutron_coord_position->setZ(neutron.getPosition(2));
        neutron_coord_position->setUniverse(_root_universe);
        cell_obj = _geometry->findCellContainingCoords(neutron_coord_position);
        cell_mat = cell_obj->getFillMaterial();
        neutron.getPositionVector(neutron_position);

        // sample whether fission event occurs
        int fission_occurs =
          neutron.arand() < cell_mat->getSigmaFByGroup(group+1) / sigma_a;

        // fission event
        if (fission_occurs == 1) {

          // sample number of neutrons released during fission
          double nu = cell_mat->getNuSigmaFByGroup(1)
            / cell_mat->getSigmaFByGroup(1);
          int lower = (int) nu;
          int add = (int) (neutron.arand() < nu -lower);
          int num_neutrons = lower + add;
          for (int i=0; i<num_neutrons; ++i) {
            fission_banks->add(neutron_position);
            tallies[FISSIONS] += 1;
          }
        }

        // end neutron history
        neutron.kill();

      } // absorption event
    } // if neutron alive
  } // while neutron alive

  delete neutron_coord_position;

  // tally crow distance
  double crow_distance;
  crow_distance = neutron.getDistance(neutron_position);
  tallies[CROWS] += crow_distance;
  tallies[NUM_CROWS] += 1;
}


/**
  * brief prints neutron info to binary file so it can be run again easily
  */
void MCSolver::saveBadNeutron(Neutron* neutron, int neutron_num, int batch) {
  FILE* out;
  out = fopen("neutron", "wb");
  std::cout << "writes file\n";

  double x = neutron->getPosition(0);
  double y = neutron->getPosition(1);
  double z = neutron->getPosition(2);
  fwrite(&x, sizeof(double), 1, out);
  fwrite(&y, sizeof(double), 1, out);
  fwrite(&z, sizeof(double), 1, out);
  std::cout << "x pos written " << x << std::endl;
  std::cout << "y pos written " << y << std::endl;
  std::cout << "z pos written " << z << std::endl;

  double x_hat = neutron->getDirection(0);
  double y_hat = neutron->getDirection(1);
  double z_hat = neutron->getDirection(2);
  fwrite(&x_hat, sizeof(double), 1, out);
  fwrite(&y_hat, sizeof(double), 1, out);
  fwrite(&z_hat, sizeof(double), 1, out);

  int group = neutron->getGroup();
  fwrite(&group, sizeof(int), 1, out);
  fwrite(&neutron_num, sizeof(int), 1, out);
  fwrite(&batch, sizeof(batch), 1, out);

  fclose(out);
}

/**
  * brief run simulation of a single neutron for debugging. Neutron is read
  * from binary file
  */
void MCSolver::trackSingleNeutron() {


  FILE* in;
  in = fopen("neutron", "r");

  int ret;

  // read neutron position
  double x, y, z;
  ret = fread(&x, sizeof(double), 1, in);
  ret = fread(&y, sizeof(double), 1, in);
  ret = fread(&z, sizeof(double), 1, in);

  std::cout << "x pos read " << x << std::endl;
  std::cout << "y pos read " << y << std::endl;
  std::cout << "z pos read " << z << std::endl;

  // read neutron direction
  double x_hat, y_hat, z_hat;
  ret = fread(&x_hat, sizeof(double), 1, in);
  ret = fread(&y_hat, sizeof(double), 1, in);
  ret = fread(&z_hat, sizeof(double), 1, in);
  
  // read neutron group, number, and batch number
  int group, neutron_num, batch;
  ret = fread(&group, sizeof(int), 1, in);
  ret = fread(&neutron_num, sizeof(int), 1, in);
  ret = fread(&batch, sizeof(int), 1, in);

  // create neutron and set its attributes
  Neutron badNeutron(neutron_num);
  badNeutron.setPosition(0, x);
  badNeutron.setPosition(1, y);
  badNeutron.setPosition(2, z);
  badNeutron.setDirection(0, x);
  badNeutron.setDirection(1, y);
  badNeutron.setDirection(2, z);
  badNeutron.setGroup(group);

  fclose(in);

  // create structures needed to transport
  std::vector <Tally> tallies(5);
  Fission fission_banks;
  bool first_round = false;
  
  std::cout << "ready to transport\n";

  // transport
  transportNeutronWithTrack(tallies, first_round, &fission_banks,
      _geometry->getNumEnergyGroups(), neutron_num, batch, neutron_num, batch,
      false, &badNeutron);

}
