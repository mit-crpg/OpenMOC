#include "../../../src/CPUSolver.h"
#include "../../../src/log.h"
#include <array>
#include <iostream>

int main() {

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
  double azim_spacing = 0.1;
  int num_azim = 4;
  double tolerance = 1e-5;
  int max_iters = 1000;

  /* Set logging information */
  set_log_level("NORMAL");
  log_printf(TITLE, "Simulating a one group homogeneous infinite medium...");
  log_printf(HEADER, "The reference keff = 1.43...");

  /* Create materials */
  log_printf(NORMAL, "Creating materials...");
  Material* infinite_medium = new Material();
  infinite_medium->setNumEnergyGroups(1);
  double sigmaF[1] = {0.0414198575};
  double nuSigmaF[1] = {0.0994076580};
  double sigmaS[1] = {0.383259177};
  double chi[1] = {1.0};
  double sigmaT[1] = {0.452648699};
  infinite_medium->setSigmaF(sigmaF, 1);
  infinite_medium->setNuSigmaF(nuSigmaF, 1);
  infinite_medium->setSigmaS(sigmaS, 1);
  infinite_medium->setChi(chi, 1);
  infinite_medium->setSigmaT(sigmaT, 1);

  /* Create surfaces */
  log_printf(NORMAL, "Creating surfaces...");
  double L = 200.0;
  XPlane left(-L/2);
  XPlane right(L/2);
  YPlane top(L/2);
  YPlane bottom(-L/2);

  left.setBoundaryType(REFLECTIVE);
  right.setBoundaryType(REFLECTIVE);
  top.setBoundaryType(REFLECTIVE);
  bottom.setBoundaryType(REFLECTIVE);

  /* Create cells */
  log_printf(NORMAL, "Creating cells...");

  Cell* cell = new Cell();
  cell->setFill(infinite_medium);
  cell->addSurface(+1, &left);
  cell->addSurface(-1, &right);
  cell->addSurface(+1, &bottom);
  cell->addSurface(-1, &top);

  /* Create universes */
  log_printf(NORMAL, "Creating universes...");

  Universe* root_universe = new Universe();
  root_universe->addCell(cell);

  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.setRootUniverse(root_universe);

  /* Generate tracks */
  log_printf(NORMAL, "Initializing the track generator...");
  TrackGenerator track_generator(&geometry, num_azim, azim_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.generateTracks();

  /* Run simulation */
  CPUSolver solver(&track_generator);
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  return 0;
}
