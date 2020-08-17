#include "CPUSolver.h"
#include "CPULSSolver.h"
#include "log.h"
#include "Mesh.h"
#include "RunTime.h"
#include <array>
#include <iostream>

int main(int argc, char* argv[]) {

#ifdef MPIx
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  log_set_ranks(MPI_COMM_WORLD);
  if (provided < MPI_THREAD_SERIALIZED) {
    log_printf(ERROR, "Not enough thread support level in the MPI library");
  }
#endif

  int arg_index = 0;
  std::string msg_string;
  while (arg_index < argc) {
    msg_string += argv[arg_index];
    msg_string += " ";
    arg_index++;
  }

  RuntimeParameters runtime;
  runtime.setRuntimeParameters(argc, argv);

  /* Exit early if only asked to display help message */
  if (runtime._print_usage)
    return 0;

  /* stuck here for debug tools to attach */
  while (runtime._debug_flag) ;

  /* Define simulation parameters */
#ifdef OPENMP
  int num_threads = runtime._num_threads;
#else
  int num_threads = 1;
#endif

  /* Set logging information */
  if(runtime._log_filename)
    set_log_filename(runtime._log_filename);
  set_log_level(runtime._log_level);
  set_line_length(120);

  log_printf(NORMAL, "Run-time options: %s", msg_string.c_str());
  log_printf(NORMAL, "Azimuthal spacing = %f", runtime._azim_spacing);
  log_printf(NORMAL, "Azimuthal angles = %d", runtime._num_azim);
  log_printf(NORMAL, "Polar spacing = %f", runtime._polar_spacing);
  log_printf(NORMAL, "Polar angles = %d", runtime._num_polar);

  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry *geometry = new Geometry();
  if(runtime._geo_filename.empty())
    log_printf(ERROR, "No geometry file is provided");
  geometry->loadFromFile(runtime._geo_filename);
#ifdef MPIx
  geometry->setDomainDecomposition(runtime._NDx, runtime._NDy, runtime._NDz,
                                   MPI_COMM_WORLD);
#endif
  geometry->setNumDomainModules(runtime._NMx, runtime._NMy, runtime._NMz);

  if ((runtime._NCx >0 && runtime._NCy > 0 && runtime._NCz > 0) ||
      (!runtime._cell_widths_x.empty() && !runtime._cell_widths_y.empty() &&
       !runtime._cell_widths_z.empty())) {

    /* Create CMFD mesh */
    log_printf(NORMAL, "Creating CMFD mesh...");
    Cmfd* cmfd = new Cmfd();
    cmfd->setSORRelaxationFactor(runtime._SOR_factor);
    cmfd->setCMFDRelaxationFactor(runtime._CMFD_relaxation_factor);
    if(runtime._cell_widths_x.empty() || runtime._cell_widths_y.empty() ||
        runtime._cell_widths_z.empty())
      cmfd->setLatticeStructure(runtime._NCx, runtime._NCy, runtime._NCz);
    else {
      std::vector< std::vector<double> > cmfd_widths{runtime._cell_widths_x,
          runtime._cell_widths_y, runtime._cell_widths_z};
      cmfd->setWidths(cmfd_widths);
    }
    if (!runtime._CMFD_group_structure.empty())
      cmfd->setGroupStructure(runtime._CMFD_group_structure);
    cmfd->setKNearest(runtime._knearest);
    cmfd->setCentroidUpdateOn(runtime._CMFD_centroid_update_on);
    cmfd->useAxialInterpolation(runtime._use_axial_interpolation);

    geometry->setCmfd(cmfd);
  }

  geometry->initializeFlatSourceRegions();

  /* Initialize track generator and generate tracks */
  log_printf(NORMAL, "Initializing the track generator...");
  Quadrature* quad = NULL;
  switch(runtime._quadraturetype) {
    case(0): quad = new TYPolarQuad(); break;
    case(1): quad = new LeonardPolarQuad(); break;
    case(2): quad = new GLPolarQuad(); break;
    case(3): quad = new EqualWeightPolarQuad(); break;
    case(4): quad = new EqualAnglePolarQuad(); break;
  }

  quad->setNumAzimAngles(runtime._num_azim);
  quad->setNumPolarAngles(runtime._num_polar);
  TrackGenerator3D track_generator(geometry, runtime._num_azim,
                                   runtime._num_polar, runtime._azim_spacing,
                                   runtime._polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setQuadrature(quad);
  track_generator.setSegmentFormation((segmentationType)
                                      runtime._segmentation_type);
  if(!runtime._seg_zones.empty())
    track_generator.setSegmentationZones(runtime._seg_zones);
  track_generator.generateTracks();

  /* Initialize solver and run simulation */
  CPUSolver *solver = NULL;
  if(runtime._linear_solver)
    solver= new CPULSSolver(&track_generator);
  else
    solver= new CPUSolver(&track_generator);
  if(runtime._verbose_report)
    solver->setVerboseIterationReport();
  solver->setNumThreads(num_threads);
  solver->setConvergenceThreshold(runtime._tolerance);
  solver->computeEigenvalue(runtime._max_iters,
                           (residualType)runtime._MOC_src_residual_type);
  if(runtime._time_report)
    solver->printTimerReport();

  /* Extract reaction rates */
  int my_rank = 0;
#ifdef MPIx
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
  std::string rxtype[4] = {"FISSION_RX", "TOTAL_RX", "ABSORPTION_RX", "FLUX_RX"};


  for(int m=0; m< runtime._output_mesh_lattices.size(); m++) {
    Mesh mesh(solver);
    mesh.createLattice(runtime._output_mesh_lattices[m][0],
                       runtime._output_mesh_lattices[m][1],
                       runtime._output_mesh_lattices[m][2]);
    Vector3D rx_rates = mesh.getFormattedReactionRates
                        ((RxType)runtime._output_types[m]);

    if (my_rank == 0) {
      std::cout << "Output " << m << ", reaction type: "
                << rxtype[runtime._output_types[m]]
                << ", lattice: " << runtime._output_mesh_lattices[m][0] << ","
                << runtime._output_mesh_lattices[m][1] << ","
                << runtime._output_mesh_lattices[m][2] << std::endl;
      for (int k=0; k < rx_rates.at(0).at(0).size(); k++) {
        for (int j=rx_rates.at(0).size()-1; j >= 0 ; j--) {
          for (int i=0; i < rx_rates.size(); i++) {
            std::cout << rx_rates.at(i).at(j).at(k) << " ";
          }
          std::cout << std::endl;
        }
      }
    }
  }

  for(int m=0; m<runtime._non_uniform_mesh_lattices.size(); m++) {
    Mesh mesh(solver);
    Vector3D rx_rates = mesh.getNonUniformFormattedReactionRates
        (runtime._non_uniform_mesh_lattices[m],
        (RxType)runtime._output_types[m+runtime._output_mesh_lattices.size()]);
    if (my_rank == 0) {
       std::cout <<"Output " << m+runtime._output_mesh_lattices.size()
        << ", reaction type: "
        << rxtype[runtime._output_types[m+runtime._output_mesh_lattices.size()]]
        << ", lattice: "
        << runtime._non_uniform_mesh_lattices[m][0].size() << ","
        << runtime._non_uniform_mesh_lattices[m][1].size() << ","
        << runtime._non_uniform_mesh_lattices[m][2].size() << std::endl;
      for (int k=0; k < rx_rates.at(0).at(0).size(); k++) {
        for (int j=rx_rates.at(0).size()-1; j >= 0 ; j--) {
          for (int i=0; i < rx_rates.size(); i++) {
            std::cout << rx_rates.at(i).at(j).at(k) << " ";
          }
          std::cout << std::endl;
        }
      }
    }
  }

  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}
