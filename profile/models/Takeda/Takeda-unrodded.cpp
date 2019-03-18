#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
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
  int num_azim = 16;
  double polar_spacing = 0.75;
  int num_polar = 6;
  double tolerance = 1e-6;
  int max_iters = 1000;
  int refines = 1;

  /* Set logging information */
  set_log_level("NORMAL");
  log_printf(TITLE, "Simulating the Takeda Unrodded Benchmark Problem...");

  /* Define material properties */
  log_printf(NORMAL, "Defining material properties...");

  const size_t num_groups = 2;
  std::map<std::string, std::array<double, num_groups> > sigma_a;
  std::map<std::string, std::array<double, num_groups> > sigma_t;
  std::map<std::string, std::array<double, num_groups*num_groups> > sigma_s;
  std::map<std::string, std::array<double, num_groups> > sigma_f;
  std::map<std::string, std::array<double, num_groups> > nu_sigma_f;
  std::map<std::string, std::array<double, num_groups> > chi;

  /* Define core cross-sections */
  sigma_a["Core"] = std::array<double, num_groups> {0.00852709, 0.158196};
  sigma_t["Core"] = std::array<double, num_groups> {0.223775, 1.03864};
  sigma_s["Core"] = std::array<double, num_groups*num_groups>
      {0.192423, 0.0228253,
        0.00, 0.880439};
  sigma_f["Core"] = std::array<double, num_groups> {0.0004, 0.1};
  nu_sigma_f["Core"] = std::array<double, num_groups> {0.00909319, 0.290183};
  chi["Core"] = std::array<double, num_groups> {1.0, 0.0};

  /* Define reflector cross-sections */
  sigma_a["Reflector"] = std::array<double, num_groups> {0.000416392,
                                                         0.0202999};
  sigma_t["Reflector"] = std::array<double, num_groups> {0.250367, 1.64482};
  sigma_s["Reflector"] = std::array<double, num_groups*num_groups>
      {0.193446, 0.0565042,
       0.00, 1.62452};
  sigma_f["Reflector"] = std::array<double, num_groups> {0.0, 0.0};
  nu_sigma_f["Reflector"] = std::array<double, num_groups> {0.0, 0.0};
  chi["Reflector"] = std::array<double, num_groups> {1.0, 0.0};

  /* Define void cross-sections */
  sigma_a["Void"] = std::array<double, num_groups> {0.0000465132, 0.00132890};
  sigma_t["Void"] = std::array<double, num_groups> {0.0128407, 0.0120676};
  sigma_s["Void"] = std::array<double, num_groups*num_groups>
      {0.01277, 0.0000240997,
       0.00, 0.0107387};
  sigma_f["Void"] = std::array<double, num_groups> {0.0, 0.0};
  nu_sigma_f["Void"] = std::array<double, num_groups> {0.0, 0.0};
  chi["Void"] = std::array<double, num_groups> {1.0, 0.0};

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

  /* Create boundaries of the geometry */
  log_printf(NORMAL, "Creating surfaces...");

  XPlane xmin(-12.5);
  XPlane xmax( 12.5);
  YPlane ymin(-12.5);
  YPlane ymax( 12.5);
  ZPlane zmin(-12.5);
  ZPlane zmax( 12.5);

  xmin.setBoundaryType(REFLECTIVE);
  ymin.setBoundaryType(REFLECTIVE);
  zmin.setBoundaryType(REFLECTIVE);
  xmax.setBoundaryType(VACUUM);
  ymax.setBoundaryType(VACUUM);
  zmax.setBoundaryType(VACUUM);

  /* Create cells in the geometry */
  log_printf(NORMAL, "Creating cells...");

  Cell* core_cell = new Cell();
  core_cell->setFill(materials["Core"]);

  Cell* reflector_cell = new Cell();
  reflector_cell->setFill(materials["Reflector"]);

  Cell* void_cell = new Cell();
  void_cell->setFill(materials["Void"]);

  /* Root Cell* */
  Cell* root_cell = new Cell(1, "root");
  root_cell->addSurface(+1, &xmin);
  root_cell->addSurface(-1, &xmax);
  root_cell->addSurface(+1, &ymin);
  root_cell->addSurface(-1, &ymax);
  root_cell->addSurface(+1, &zmin);
  root_cell->addSurface(-1, &zmax);

  /* Create Universes */
  log_printf(NORMAL, "Creating universes...");

  Universe* core = new Universe();
  Universe* reflector = new Universe();
  Universe* void_u = new Universe();
  Universe* root_universe = new Universe();

  core->addCell(core_cell);
  reflector->addCell(reflector_cell);
  void_u->addCell(void_cell);
  root_universe->addCell(root_cell);

  /* Setup Takeda core */
  log_printf(NORMAL, "Creating Takeda core...");

  Lattice* lattice = new Lattice();
  lattice->setWidth(5.0/refines, 5.0/refines, 5.0/refines);

  Universe* mold[] = {reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      void_u,    reflector};

  /* Refine lattice */
  Universe* refined_mold[5*5*5*refines*refines*refines];
  for (int x=0; x < refines; x++) {
    for (int i=0; i < 5; i++) {
      for (int y=0; y < refines; y++) {
        for (int j=0; j < 5; j++) {
          for (int z=0; z < refines; z++) {
            for (int k=0; k < 5; k++) {
              int index = 5*5*refines*refines*refines*k + 5*5*refines*refines*z
                + 5*refines*refines*j + 5*refines*y + i*refines + x;
              refined_mold[index] = mold[5*5*k + 5*j + i];
            }
          }
        }
      }
    }
  }
  lattice->setUniverses(5*refines, 5*refines, 5*refines, refined_mold);
  root_cell->setFill(lattice);

  /* Create CMFD mesh */
  log_printf(NORMAL, "Creating Cmfd mesh...");
  Cmfd* cmfd = new Cmfd();
  cmfd->setCMFDRelaxationFactor(0.7);
  cmfd->setLatticeStructure(5, 5, 5);
  cmfd->setKNearest(1);

  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");

  Geometry geometry;
  geometry.setRootUniverse(root_universe);
  geometry.setCmfd(cmfd);
  geometry.initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");

  Quadrature* quad = new EqualAnglePolarQuad();
  quad->setNumPolarAngles(num_polar);
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setSegmentFormation(OTF_TRACKS);
  std::vector<double> seg_zones {-12.5, 12.5};
  track_generator.setSegmentationZones(seg_zones);
  track_generator.setNumThreads(num_threads);
  track_generator.setQuadrature(quad);
  track_generator.generateTracks();

  /* Run simulation */
  CPULSSolver solver(&track_generator);
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}
