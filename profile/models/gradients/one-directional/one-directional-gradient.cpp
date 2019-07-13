#include "../../../../src/CPUSolver.h"
#include "../../../../src/log.h"
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
  double tolerance = 1e-5;
  int max_iters = 2000;

  /* Set logging information */
  set_log_level("NORMAL");
  log_printf(TITLE, "Simulating a one group homogeneous two directional"
     " gradient...");

  /* Create materials */
  log_printf(NORMAL, "Creating materials...");
  Material* basic_material = new Material();
  basic_material->setNumEnergyGroups(1);
  double sigmaF[1] = {0.0414198575};
  double nuSigmaF[1] = {0.0994076580};
  double sigmaS[1] = {0.383259177};
  double chi[1] = {1.0};
  double sigmaT[1] = {0.452648699};
  basic_material->setSigmaF(sigmaF, 1);
  basic_material->setNuSigmaF(nuSigmaF, 1);
  basic_material->setSigmaS(sigmaS, 1);
  basic_material->setChi(chi, 1);
  basic_material->setSigmaT(sigmaT, 1);

  /* Create surfaces */
  log_printf(NORMAL, "Creating surfaces...");
  double L = 100.0;
  XPlane left(-L/2);
  XPlane right(L/2);
  YPlane top(L/2);
  YPlane bottom(-L/2);

  left.setBoundaryType(VACUUM);
  right.setBoundaryType(REFLECTIVE);
  top.setBoundaryType(REFLECTIVE);
  bottom.setBoundaryType(REFLECTIVE);

  /* Create cells */
  log_printf(NORMAL, "Creating cells...");
  Cell* fill = new Cell();
  fill->setFill(basic_material);

  Cell* root_cell = new Cell();
  root_cell->addSurface(+1, &left);
  root_cell->addSurface(-1, &right);
  root_cell->addSurface(+1, &bottom);
  root_cell->addSurface(-1, &top);

  /* Create universes */
  log_printf(NORMAL, "Creating universes...");
  Universe* fill_universe = new Universe();
  fill_universe->addCell(fill);

  Universe* root_universe = new Universe();
  root_universe->addCell(root_cell);

  /* Create lattice */
  log_printf(NORMAL, "Creating 100 x 1 lattice...");
  int num_cells_x = 100;
  int num_cells_y = 1;
  Lattice* lattice = new Lattice();
  lattice->setWidth(L/num_cells_x, L/num_cells_y);

  Universe** matrix = new Universe*[num_cells_x * num_cells_y];
  for (int j=0; j<num_cells_y; j++)
    for (int i=0; i<num_cells_x; i++)
      matrix[(num_cells_y-1-j)*num_cells_x + i] = fill_universe;
  lattice->setUniverses(1, num_cells_y, num_cells_x, matrix);
  root_cell->setFill(lattice);

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

#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}
