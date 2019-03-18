#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
#include "../../../src/log.h"
#include "../../../src/Mesh.h"
#include <array>
#include <iostream>
#include "helper-code/group-structures.h"
#include <fenv.h>

int main(int argc,  char* argv[]) {

#ifdef MPIx
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW & ~FE_OVERFLOW);
  /* Define geometry to load */
  std::string file = "tw-final-full-core-2D.geo";

  /* Define simulation parameters */
#ifdef OPENMP
  int num_threads = omp_get_num_threads();
#else
  int num_threads = 1;
#endif
 
  double azim_spacing = 0.05;
  int num_azim = 8;
  double polar_spacing = 1.5;
  int num_polar = 4;

  double tolerance = 1e-4;
  int max_iters = 550;
  
  /* Create CMFD lattice */
  Cmfd cmfd;
  cmfd.useAxialInterpolation(true);
  cmfd.setLatticeStructure(17*17, 17*17, 5);
  cmfd.setKNearest(1);
  std::vector<std::vector<int> > cmfd_group_structure =
      get_group_structure(70, 70);
  cmfd.setGroupStructure(cmfd_group_structure);
  cmfd.setCMFDRelaxationFactor(0.5);
  cmfd.setSORRelaxationFactor(1.6);
  cmfd.useFluxLimiting(true);

  /* Load the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.loadFromFile(file, true);

  //FIXME
  int num_rad_discr = 96;
  int rad_discr_domains[96*2];
  int ind = 0;
  for (int i=0; i < 17; i++) {
    int offset = std::abs(8-i);
    if (offset == 8) {
      for (int j=0; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset == 7) {
      for (int j=0; j < 5; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
      for (int j=12; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset == 6) {
      for (int j=0; j < 3; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
      for (int j=14; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset == 5 || offset == 4) {
      for (int j=0; j < 2; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
      for (int j=15; j < 17; j++) {
        rad_discr_domains[2*ind] = i;
        rad_discr_domains[2*ind+1] = j;
        ind++;
      }
    }
    else if (offset < 4) {
      rad_discr_domains[2*ind] = i;
      rad_discr_domains[2*ind+1] = 0;
      ind++;
      rad_discr_domains[2*ind] = i;
      rad_discr_domains[2*ind+1] = 16;
      ind++;
    }
  }

  geometry.setCmfd(&cmfd); //FIXME OFF /ON
  log_printf(NORMAL, "Pitch = %8.6e", geometry.getMaxX() - geometry.getMinX());
  log_printf(NORMAL, "Pitch = %8.6e", geometry.getMaxY() - geometry.getMinY());
  log_printf(NORMAL, "Height = %8.6e", geometry.getMaxZ() - geometry.getMinZ());
#ifdef MPIx
  geometry.setDomainDecomposition(17, 17, 1, MPI_COMM_WORLD);
#endif
  geometry.setOverlaidMesh(10.0, 17*17*1, 17*17*1, num_rad_discr, 
                           rad_discr_domains);
  geometry.initializeFlatSourceRegions();

  /* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");
  
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setSegmentFormation(OTF_STACKS);
  double z_arr[] = {120., 130.};
  std::vector<double> segmentation_zones(z_arr, z_arr + sizeof(z_arr) 
                                         / sizeof(double));
  track_generator.setSegmentationZones(segmentation_zones);
  track_generator.generateTracks();

  /* Find fissionable material */
  std::map<int, Material*> materials = geometry.getAllMaterials();
  for (std::map<int, Material*>::iterator it = materials.begin();
       it != materials.end(); ++it) {
    Material* mat = it->second;
    for (int g=0; g < 70; g++) {
      if (mat->getSigmaTByGroup(g+1) < 0.0)
          log_printf(ERROR, "Negative sigma-t of %6.4f", mat->getSigmaTByGroup(g+1));
    }
  }

  /* Run simulation */
  CPULSSolver solver(&track_generator); //FIXME LS / FS
  solver.stabilizeTransport(0.25);
  //solver.setResidualByReference("ref-fluxes-fc");
  solver.setNumThreads(num_threads);
  solver.setVerboseIterationReport();
  solver.setConvergenceThreshold(tolerance);
  solver.setCheckXSLogLevel(INFO);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  Lattice mesh_lattice;
  Mesh mesh(&solver);
  mesh.createLattice(17*17, 17*17, 1);
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
