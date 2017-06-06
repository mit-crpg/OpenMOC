#include "Mesh.h"

#ifndef SWIG

/**
 * @brief The Mesh constructor
 * @details If no lattice is given, a default lattice can be constructed with
 *          the Mesh::createLattice function.
 * @param solver The solver from which scalar fluxes and cross-sections are
 *        extracted
 * @param lattice An optional parameter for the lattice across which reaction
 *        rates are tallied
 */
Mesh::Mesh(Solver* solver, Lattice* lattice) {
  _solver = solver;
  _lattice = lattice;
  _lattice_allocated = false;
}


/**
 * @brief The Mesh destrcutor deletes its lattice if the lattice was allocated
 *        internally
 */
Mesh::~Mesh() {
  if (_lattice_allocated)
    delete _lattice;
}


/**
 * @brief Creates an internal lattice over which to tally reaction rates with
 *        the user-input dimensions
 * @param num_x the number of mesh cells in the x-direction
 * @param num_y the number of mesh cells in the y-direction
 * @param num_z the number of mesh cells in the z-direction
 */
void Mesh::createLattice(int num_x, int num_y, int num_z) {

  /* Delete the current lattice if currently allocated */
  if (_lattice_allocated)
    delete _lattice;

  /* Allocate a new lattice */
  _lattice = new Lattice();
  _lattice->setNumX(num_x);
  _lattice->setNumY(num_y);
  _lattice->setNumZ(num_z);
  _lattice_allocated = true;

  /* Get the root universe */
  Geometry* geometry = _solver->getGeometry();
  Universe* root_universe = geometry->getRootUniverse();

  /* Determine the geometry widths in each direction */
  double width_x = (root_universe->getMaxX() - root_universe->getMinX())/num_x;
  double width_y = (root_universe->getMaxY() - root_universe->getMinY())/num_y;
  double width_z = (root_universe->getMaxZ() - root_universe->getMinZ())/num_z;

  /* Determine the center-point of the geometry */
  double offset_x = (root_universe->getMinX() + root_universe->getMaxX()) / 2;
  double offset_y = (root_universe->getMinY() + root_universe->getMaxY()) / 2;
  double offset_z = (root_universe->getMinZ() + root_universe->getMaxZ()) / 2;

  /* Create the Mesh lattice */
  _lattice->setWidth(width_x, width_y, width_z);
  _lattice->setOffset(offset_x, offset_y, offset_z);
}


/**
 * @brief Tallies reaction rates of the given type over the Mesh lattice
 * @param rx The type of reaction to tally
 * @return The reaction rates in a 1D vector indexed by the lattice cell IDs
 */
std::vector<double> Mesh::getReactionRates(RxType rx) {

  /* Check that the Mesh contains a lattice */
  if (_lattice == NULL)
    log_printf(ERROR, "A Lattice must be set or created to get reaction rates "
                      "form a Mesh object");

  /* Extract fluxes and geometry information */
  Geometry* geometry = _solver->getGeometry();
  NEW_PRECISION* volumes = _solver->getTrackGenerator()->getFSRVolumesBuffer();
  NEW_PRECISION* fluxes = _solver->getFluxesArray();
  long num_fsrs = geometry->getNumFSRs();

  /* Create a 1D array of reaction rates with the appropriate size */
  std::vector<double> rx_rates;
  int size = _lattice->getNumX() * _lattice->getNumY() * _lattice->getNumZ();
  rx_rates.resize(size);

  /* Extract the number of groups */
  int num_groups = geometry->getNumEnergyGroups();

  /* Create temporary array for cross-sections */
  NEW_PRECISION temp_array[num_groups];

  /* Loop over all flat source regions */
  for (long r=0; r < num_fsrs; r++) {

    /* Determine the FSR material and which Mesh cell contains it */
    Material* mat = geometry->findFSRMaterial(r);
    Point* pt = geometry->getFSRPoint(r);
    int lat_cell = _lattice->getLatticeCell(pt);

    /* Determine the volume and cross-sections of the FSR */
    NEW_PRECISION volume = volumes[r];
    NEW_PRECISION* xs_array;
    switch (rx) {
      case FISSION_RX:
        xs_array = mat->getSigmaF();
        break;
      case TOTAL_RX:
        xs_array = mat->getSigmaT();
        break;
      case ABSORPTION_RX:
        {
          xs_array = temp_array;
          NEW_PRECISION* scattering = mat->getSigmaS();
          for (int g=0; g < num_groups; g++) {
            xs_array[g] = 0.0;
            for (int gp=0; gp < num_groups; gp++) {
              xs_array[g] += scattering[gp*num_groups + g];
            }
          }
        }
        break;
      case FLUX_RX:
        xs_array = temp_array;
        for (int g=0; g < num_groups; g++)
          xs_array[g] = 1.0;
        break;
      default:
        log_printf(ERROR, "Unrecgonized reaction type in Mesh object");
    }

    /* Tally the reaction rates summed across all groups */
    for (int g=0; g < num_groups; g++) {
      double xs = xs_array[g];
      rx_rates.at(lat_cell) += fluxes[r*num_groups + g] * volume * xs;
    }
  }

  /* If domain decomposed, do a reduction */
#ifdef MPIx
  if (geometry->isDomainDecomposed()) {
    MPI_Comm comm = geometry->getMPICart();
    double* rx_rates_array = &rx_rates[0];
    double* rx_rates_send = new double[size];
    for (int i=0; i < size; i++)
      rx_rates_send[i] = rx_rates_array[i];
    MPI_Allreduce(rx_rates_send, rx_rates_array, size, MPI_DOUBLE, MPI_SUM,
                  comm);
    delete [] rx_rates_send;
  }
#endif

  return rx_rates;
}


/**
 * @brief Tallies reaction rates of the given type over the Mesh lattice
 * @param rx The type of reaction to tally
 * @return The reaction rates in a 3D vector indexed by the lattice cell
 *         x, y, and z indexes
 */
Vector3D Mesh::getFormattedReactionRates(RxType rx) {

  /* Extract reaction rates */
  Vector3D rx_rates;
  std::vector<double> rx_rates_array = getReactionRates(rx);

  /* Format reaction rates into a 3D array */
  int num_x = _lattice->getNumX();
  for (int i=0; i < num_x; i++) {
    int num_y = _lattice->getNumY();
    std::vector<std::vector<double> > vector_2D;
    for (int j=0; j < num_y; j++) {
      int num_z = _lattice->getNumZ();
      std::vector<double> vector_1D;
      for (int k=0; k < num_z; k++) {
        int idx = k * num_x * num_y + j * num_x + i;
        vector_1D.push_back(rx_rates_array[idx]);
      }
      vector_2D.push_back(vector_1D);
    }
    rx_rates.push_back(vector_2D);
  }
  return rx_rates;
}

#endif
