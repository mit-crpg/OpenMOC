#include "../src/CPUSolver.h"
#include "../src/Solver.h"
#include "../src/log.h"
#include <array>
#include <iostream>

int main(){

  /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
  double track_spacing = 0.1;
  int num_azim = 4;
  double tolerance = 1e-5;
  int max_iters = 1000;

  /* Set logging information */
  set_log_level("NORMAL");
  log_printf(TITLE, "Simulating a one group homogeneous two directional"
     " gradient...");
  
  /* Create materials */
  log_printf(NORMAL, "Creating materials...");
  Material basic_material;
  basic_material.setNumEnergyGroups(1);
  double sigmaA[1] = {0.069389522};
  double sigmaF[1] = {0.0414198575};
  double nuSigmaF[1] = {0.0994076580};
  double sigmaS[1] = {0.383259177};
  double chi[1] = {1.0};
  double sigmaT[1] = {0.452648699};
  basic_material.setSigmaA(sigmaA, 1);
  basic_material.setSigmaF(sigmaF, 1);
  basic_material.setNuSigmaF(nuSigmaF, 1);
  basic_material.setSigmaS(sigmaS, 1);
  basic_material.setChi(chi, 1);
  basic_material.setSigmaT(sigmaT, 1);

  /* Create surfaces */
  log_printf(NORMAL, "Creating surfaces...");
  double L = 100.0;
  XPlane left(-L/2);
  XPlane right(L/2);
  YPlane top(L/2);
  YPlane bottom(-L/2);

  left.setBoundaryType(VACUUM);
  right.setBoundaryType(VACUUM);
  top.setBoundaryType(REFLECTIVE);
  bottom.setBoundaryType(REFLECTIVE);

  /* Create cells */
  log_printf(NORMAL, "Creating cells...");
  Cell fill;
  fill.setFill(&basic_material);

  Cell root_cell;
  root_cell.addSurface(+1, &left);
  root_cell.addSurface(-1, &right);
  root_cell.addSurface(+1, &bottom);
  root_cell.addSurface(-1, &top);

  /* Create universes */
  log_printf(NORMAL, "Creating universes...");
  Universe fill_universe;
  fill_universe.addCell(&fill);

  Universe root_universe;
  root_universe.addCell(&root_cell);

  /* Create lattice */
  log_printf(NORMAL, "Creating 100 x 1 lattice...");
  int num_cells_x = 100;
  int num_cells_y = 1;
  Lattice lattice;
  lattice.setWidth(L/num_cells_x, L/num_cells_y);
  
  Universe** matrix = new Universe*[num_cells_x * num_cells_y];
  for(int j=0; j<num_cells_y; j++)
    for(int i=0; i<num_cells_x; i++)
      matrix[(num_cells_y-1-j)*num_cells_x + i] = &fill_universe;
  lattice.setUniverses(num_cells_y, num_cells_x, matrix);
  root_cell.setFill(&lattice);

  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.setRootUniverse(&root_universe);
  geometry.initializeFlatSourceRegions();

  /* Generate tracks */
  log_printf(NORMAL, "Initializing the track generator...");
  TrackGenerator track_generator(&geometry, num_azim, track_spacing);
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
