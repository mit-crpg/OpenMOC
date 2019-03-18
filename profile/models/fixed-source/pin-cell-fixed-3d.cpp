#include "../../../src/CPUSolver.h"
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
  double polar_spacing = 0.1;
  int num_polar = 6;
  double tolerance = 1e-5;
  int max_iters = 1000;

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
  log_printf(NORMAL, "Creating surfaces...");

  ZCylinder* pin = new ZCylinder(0.0, 0.0, 1.0);
  XPlane xmin(-2.0);
  XPlane xmax( 2.0);
  YPlane ymin(-2.0);
  YPlane ymax( 2.0);
  ZPlane zmin(-2.0);
  ZPlane zmax( 2.0);

  xmin.setBoundaryType(REFLECTIVE);
  ymin.setBoundaryType(REFLECTIVE);
  zmin.setBoundaryType(REFLECTIVE);
  xmax.setBoundaryType(REFLECTIVE);
  ymax.setBoundaryType(REFLECTIVE);
  zmax.setBoundaryType(REFLECTIVE);

  /* Create cells */
  log_printf(NORMAL, "Creating cells...");

  Cell* fuel = new Cell();
  fuel->setFill(materials["Water"]);
  fuel->addSurface(-1, pin);
  fuel->addSurface(+1, &zmin);
  fuel->addSurface(-1, &zmax);
  fuel->setNumSectors(8);

  Cell* moderator = new Cell();
  moderator->setFill(materials["Water"]);
  moderator->addSurface(+1, pin);
  moderator->addSurface(+1, &xmin);
  moderator->addSurface(-1, &xmax);
  moderator->addSurface(+1, &ymin);
  moderator->addSurface(-1, &ymax);
  moderator->addSurface(+1, &zmin);
  moderator->addSurface(-1, &zmax);
  moderator->setNumSectors(8);

  /* Add universes */
  log_printf(NORMAL, "Creating universes...");

  Universe* root_universe = new Universe();
  root_universe->addCell(fuel);
  root_universe->addCell(moderator);

  /* Creat the geometry */
  log_printf(NORMAL, "Creating geometry...");

  Geometry* geometry = new Geometry();
  geometry->setRootUniverse(root_universe);
  geometry->initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");

  TrackGenerator3D track_generator(geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setSegmentFormation(OTF_TRACKS);
  track_generator.generateTracks();

  /* Run simulation */
  CPUSolver solver(&track_generator);
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.setFixedSourceByCell(fuel, 1, 1.0);
  solver.computeFlux(max_iters, false);
  solver.printTimerReport();

  solver.printFluxesTemp();
  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}
