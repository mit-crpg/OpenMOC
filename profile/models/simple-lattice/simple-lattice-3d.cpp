#include "../../../src/CPULSSolver.h"
#include "../../../src/Mesh.h"
#include "../../../src/log.h"
#include <array>
#include <iostream>

int main(int argc, char* argv[]) {

#ifdef MPIx
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  /* Define simulation parameters */
#ifdef OPENMP
  int num_threads = omp_get_num_threads();
#else
  int num_threads = 1;
#endif
  double azim_spacing = 0.1;
  int num_azim = 4;
  double polar_spacing = 0.5;
  int num_polar = 2;
  double tolerance = 1e-5;
  int max_iters = 1000;

  int fuel_rings = 1;
  int moderator_rings = 1;
  int fuel_sectors = 1;
  int moderator_sectors = 1;

  /* Define material properties */
  log_printf(NORMAL, "Defining material properties...");

  const size_t num_groups = 7;
  std::map<std::string, std::array<double, num_groups> > sigma_a;
  std::map<std::string, std::array<double, num_groups> > nu_sigma_f;
  std::map<std::string, std::array<double, num_groups> > sigma_f;
  std::map<std::string, std::array<double, num_groups*num_groups> > sigma_s;
  std::map<std::string, std::array<double, num_groups> > chi;
  std::map<std::string, std::array<double, num_groups> > sigma_t;

  /* Define water cross-sections */
  sigma_a["Water"] = std::array<double, num_groups> {6.0105E-4, 1.5793E-5,
      3.3716E-4, 0.0019406, 0.0057416, 0.015001, 0.037239};
  nu_sigma_f["Water"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_f["Water"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_s["Water"] = std::array<double, num_groups*num_groups>
      {0.0444777, 0.1134, 7.2347E-4, 3.7499E-6, 5.3184E-8, 0.0, 0.0,
      0.0, 0.282334, 0.12994, 6.234E-4, 4.8002E-5, 7.4486E-6, 1.0455E-6,
      0.0, 0.0, 0.345256, 0.22457, 0.016999, 0.0026443, 5.0344E-4,
      0.0, 0.0, 0.0, 0.0910284, 0.41551, 0.063732, 0.012139,
      0.0, 0.0, 0.0, 7.1437E-5, 0.139138, 0.51182, 0.061229,
      0.0, 0.0, 0.0, 0.0, 0.0022157, 0.699913, 0.53732,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.13244, 2.4807};
  chi["Water"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_t["Water"] = std::array<double, num_groups> {0.159206, 0.41297,
    0.59031, 0.58435, 0.718, 1.25445, 2.65038};

  /* Define UO2 cross-sections */
  sigma_a["UO2"] = std::array<double, num_groups> {0.0080248, 0.0037174,
    0.026769, 0.096236, 0.03002, 0.11126, 0.28278};
  nu_sigma_f["UO2"] = std::array<double, num_groups> {0.02005998, 0.002027303,
    0.01570599, 0.04518301, 0.04334208, 0.2020901, 0.5257105};
  sigma_f["UO2"] = std::array<double, num_groups> {0.00721206, 8.19301E-4,
    0.0064532, 0.0185648, 0.0178084, 0.0830348, 0.216004};
  sigma_s["UO2"] = std::array<double, num_groups*num_groups>
      {0.127537, 0.042378, 9.4374E-6, 5.5163E-9, 0.0, 0.0, 0.0,
      0.0, 0.324456, 0.0016314, 3.1427E-9, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.45094, 0.0026792, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.452565, 0.0055664, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.2525E-4, 0.271401, 0.010255, 1.0021E-8,
      0.0, 0.0, 0.0, 0.0, 0.0012968, 0.265802, 0.016809,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0085458, 0.27308};
  chi["UO2"] = std::array<double, num_groups> {0.58791, 0.41176, 3.3906E-4,
    1.1761E-7, 0.0, 0.0, 0.0};
  sigma_t["UO2"] = std::array<double, num_groups> {0.177949, 0.329805,
    0.480388, 0.554367, 0.311801, 0.395168, 0.564406};

  /* Create materials */
  log_printf(NORMAL, "Creating materials...");
  std::map<std::string, Material*> materials;
  std::map<std::string, std::array<double, num_groups> >::iterator it;
  int id_num = 0;
  for (it = sigma_t.begin(); it != sigma_t.end(); it++) {

    std::string name = it->first;
    materials[name] = new Material(id_num, name.c_str());
    materials[name]->setNumEnergyGroups(num_groups);
    id_num++;

    materials[name]->setSigmaT(sigma_t[name].data(), num_groups);
    materials[name]->setSigmaS(sigma_s[name].data(), num_groups*num_groups);
    materials[name]->setSigmaF(sigma_f[name].data(), num_groups);
    materials[name]->setNuSigmaF(nu_sigma_f[name].data(), num_groups);
    materials[name]->setChi(chi[name].data(), num_groups);
  }

  /* Create surfaces */
  XPlane xmin(-2.0);
  XPlane xmax( 2.0);
  YPlane ymin(-2.0);
  YPlane ymax( 2.0);
  ZPlane zmin(-10.0);
  ZPlane zmax( 10.0);

  xmin.setBoundaryType(REFLECTIVE);
  ymin.setBoundaryType(REFLECTIVE);
  zmin.setBoundaryType(VACUUM);
  xmax.setBoundaryType(REFLECTIVE);
  ymax.setBoundaryType(REFLECTIVE);
  zmax.setBoundaryType(VACUUM);

  ZCylinder large_pin(0.0, 0.0, 0.4);
  ZCylinder medium_pin(0.0, 0.0, 0.3);
  ZCylinder small_pin(0.0, 0.0, 0.2);

  /* Create cells */
  log_printf(NORMAL, "Creating cells...");

  Cell large_fuel;
  large_fuel.setFill(materials["UO2"]);
  large_fuel.addSurface(-1, &large_pin);
  large_fuel.setNumRings(fuel_rings);
  large_fuel.setNumSectors(fuel_sectors);

  Cell large_moderator;
  large_moderator.setFill(materials["Water"]);
  large_moderator.addSurface(+1, &large_pin);
  large_moderator.setNumRings(moderator_rings);
  large_moderator.setNumSectors(moderator_sectors);

  Cell medium_fuel;
  medium_fuel.setFill(materials["UO2"]);
  medium_fuel.addSurface(-1, &medium_pin);
  medium_fuel.setNumRings(fuel_rings);
  medium_fuel.setNumSectors(fuel_sectors);

  Cell medium_moderator;
  medium_moderator.setFill(materials["Water"]);
  medium_moderator.addSurface(+1, &medium_pin);
  medium_moderator.setNumRings(moderator_rings);
  medium_moderator.setNumSectors(moderator_sectors);

  Cell small_fuel;
  small_fuel.setFill(materials["UO2"]);
  small_fuel.addSurface(-1, &small_pin);
  small_fuel.setNumRings(fuel_rings);
  small_fuel.setNumSectors(fuel_sectors);

  Cell small_moderator;
  small_moderator.setFill(materials["Water"]);
  small_moderator.addSurface(+1, &small_pin);
  small_moderator.setNumRings(moderator_rings);
  small_moderator.setNumSectors(moderator_sectors);

  Cell root_cell;
  root_cell.addSurface(+1, &xmin);
  root_cell.addSurface(-1, &xmax);
  root_cell.addSurface(+1, &ymin);
  root_cell.addSurface(-1, &ymax);
  root_cell.addSurface(+1, &zmin);
  root_cell.addSurface(-1, &zmax);

  /* Create universes */
  log_printf(NORMAL, "Creating universes...");

  Universe pin1;
  Universe pin2;
  Universe pin3;
  Universe root_universe;

  pin1.addCell(&large_fuel);
  pin1.addCell(&large_moderator);
  pin2.addCell(&medium_fuel);
  pin2.addCell(&medium_moderator);
  pin3.addCell(&small_fuel);
  pin3.addCell(&small_moderator);
  root_universe.addCell(&root_cell);

  /* Create lattice */
  log_printf(NORMAL, "Creating simple 4 x 4 lattice...");

  Lattice lattice;
  lattice.setWidth(1.0, 1.0, 20.0);

  Universe* lattice_matrix[4*4*1];
  {
    int mold[4*4] = {1, 2, 1, 2,
                     2, 3, 2, 3,
                     1, 2, 1, 2,
                     2, 3, 2, 3};
    std::map<int, Universe*> names = {{1, &pin1}, {2, &pin2}, {3, &pin3}};
    for (int z=0; z < 1; z++)
      for (int n=0; n < 4*4; n++)
        lattice_matrix[z*4*4 + n] = names[mold[n]];

    lattice.setUniverses(1, 4, 4, lattice_matrix);
  }
  root_cell.setFill(&lattice);

  /* Create CMFD mesh */
  log_printf(NORMAL, "Creating Cmfd mesh...");

  Cmfd cmfd;
  cmfd.useAxialInterpolation(true);
  cmfd.setLatticeStructure(4,4,4);
  std::vector<std::vector<int> > cmfd_group_structure;
  cmfd_group_structure.resize(2);
  for (int g=0; g<3; g++)
    cmfd_group_structure.at(0).push_back(g+1);
  for (int g=3; g<7; g++)
    cmfd_group_structure.at(1).push_back(g+1);
  cmfd.setGroupStructure(cmfd_group_structure);
  cmfd.setKNearest(9);
  cmfd.useFluxLimiting(false);
  //setCentroidUpdateOn is the prerequisites for k-nearest update. 
  cmfd.setCentroidUpdateOn(false);
  //useAxialInterpolation sets the axial flux update to be 2nd order interpolation.
  cmfd.useAxialInterpolation(true);
  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");

  Geometry geometry;
  geometry.setRootUniverse(&root_universe);
#ifdef MPIx
  geometry.setDomainDecomposition(1,1,1, MPI_COMM_WORLD);
#endif
  //geometry.setNumDomainModules(2,2,1);
  geometry.setCmfd(&cmfd);
  geometry.initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");

  TrackGenerator3D track_generator(&geometry, num_azim, num_polar,
                                   azim_spacing, polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setSegmentFormation(OTF_STACKS);
  track_generator.generateTracks();

  /* Run simulation */
  CPUSolver solver(&track_generator);
  solver.setVerboseIterationReport();
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  Lattice mesh_lattice;
  Mesh mesh(&solver);
  mesh.createLattice(4, 4, 10);
  Vector3D rx_rates = mesh.getFormattedReactionRates(FISSION_RX);

  int my_rank = 0;
#ifdef MPIx
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  if (my_rank == 0) {
    for (int k=0; k < rx_rates.at(0).at(0).size(); k++) {
      std::cout << " -------- z = " << k << " ----------" << std::endl;
      for (int j=0; j < rx_rates.at(0).size(); j++) {
        for (int i=0; i < rx_rates.size(); i++) {
          std::cout << rx_rates.at(i).at(j).at(k) << " ";
        }
        std::cout << std::endl;
      }
    }
  }

  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif

  return 0;
}
