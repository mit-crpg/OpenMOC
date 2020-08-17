#ifdef __cplusplus
#include "RunTime.h"
#endif


/**
 * @brief Process the run time options.
 * @param argc number of run time command words
 * @param argv content of run time command words
 */
int RuntimeParameters::setRuntimeParameters(int argc, char *argv[]) {

  int arg_index = 0;

  /* Parse the run time commands*/
  while (arg_index < argc) {
    if(strcmp(argv[arg_index], "-debug") == 0) {
      arg_index++;
      _debug_flag = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-domain_decompose") == 0) {
      int *pointer[] = {&_NDx, &_NDy, &_NDz};
      int i = 0;
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        *pointer[i] = atoi(p);
        i++;
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-num_domain_modules") == 0) {
      int *pointer[] = {&_NMx, &_NMy, &_NMz};
      int i = 0;
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        *pointer[i] = atoi(p);
        i++;
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-CMFD_lattice") == 0) {
      int *pointer[] = {&_NCx, &_NCy, &_NCz};
      int i = 0;
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        *pointer[i] = atoi(p);
        i++;
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-output_mesh_lattice") == 0) {
      std::vector<int> tmp;
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        tmp.push_back(atoi(p));
        buf = NULL;
      }
      _output_mesh_lattices.push_back(tmp);
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-output_type") == 0 ) {
      arg_index++;
      _output_types.push_back(atoi(argv[arg_index++]));
    }
    else if(strcmp(argv[arg_index], "-num_threads") == 0 ) {
      arg_index++;
      _num_threads = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-azim_spacing") == 0) {
      arg_index++;
      _azim_spacing = atof(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-num_azim") == 0) {
      arg_index++;
      _num_azim = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-polar_spacing") == 0) {
      arg_index++;
      _polar_spacing = atof(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-num_polar") == 0) {
      arg_index++;
      _num_polar = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-MOC_src_tolerance") == 0) {
      arg_index++;
      _tolerance = atof(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-max_iters") == 0) {
      arg_index++;
      _max_iters= atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-log_level") == 0) {
      arg_index++;
      _log_level = argv[arg_index++];
    }
    else if(strcmp(argv[arg_index], "-knearest") == 0) {
      arg_index++;
      _knearest = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-CMFD_flux_update_on") == 0) {
      arg_index++;
      _CMFD_flux_update_on = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-CMFD_centroid_update_on") == 0) {
      arg_index++;
      _CMFD_centroid_update_on = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-use_axial_interpolation") == 0) {
      arg_index++;
      _use_axial_interpolation = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-help") == 0 ||
            strcmp(argv[arg_index], "--help") == 0 ||
            strcmp(argv[arg_index], "-h") == 0 ||
            strcmp(argv[arg_index], "--h") == 0) {
      _print_usage = 1;
      break;
    }
    else if(strcmp(argv[arg_index], "-log_filename") == 0) {
      arg_index++;
      _log_filename = argv[arg_index++];
    }
    else if(strcmp(argv[arg_index], "-geo_filename") == 0) {
      arg_index++;
      _geo_filename = std::string(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-widths_x") == 0) {
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      char *ql, *qr;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        ql = strtok_r(p, "*", &qr);
        if(strcmp(qr, "")) {
          for(int i=0; i< atoi(qr); i++)
            _cell_widths_x.push_back(atof(ql));
        }
        else
          _cell_widths_x.push_back(atof(p));
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-widths_y") == 0) {
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      char *ql, *qr;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        ql = strtok_r(p, "*", &qr);
        if(strcmp(qr, "")) {
          for(int i=0; i< atoi(qr); i++)
            _cell_widths_y.push_back(atof(ql));
        }
        else
          _cell_widths_y.push_back(atof(p));
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-widths_z") == 0) {
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      char *ql, *qr;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        ql = strtok_r(p, "*", &qr);
        if(strcmp(qr, "")) {
          for(int i=0; i< atoi(qr); i++)
            _cell_widths_z.push_back(atof(ql));
        }
        else
          _cell_widths_z.push_back(atof(p));
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-ls_solver") == 0) {
      arg_index++;
      _linear_solver=atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-seg_zones") == 0) {
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *p;
      char *ql, *qr;
      while((p = strtok_r(buf, ",", &outer_ptr)) != NULL) {
        ql = strtok_r(p, "*", &qr);
        if(strcmp(qr, "")) {
          for(int i=0; i< atoi(qr); i++)
            _seg_zones.push_back(atof(ql));
        }
        else
          _seg_zones.push_back(atof(p));
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-MOC_src_residual_type") == 0) {
      arg_index++;
      _MOC_src_residual_type = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-SOR_factor") == 0) {
      arg_index++;
      _SOR_factor = atof(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-CMFD_relaxation_factor") == 0) {
      arg_index++;
      _CMFD_relaxation_factor = atof(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-segmentation_type") == 0) {
      arg_index++;
      _segmentation_type = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-CMFD_group_structure") == 0) {
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *inner_ptr = NULL;
      char *p;
      while((p = strtok_r(buf, "/", &outer_ptr)) != NULL) {
        buf = p;
        char *ql, *qr;
        std::vector<int> tmp;
        while((p = strtok_r(buf, ",", &inner_ptr)) != NULL) {
          ql = strtok_r(p, "-", &qr);
          if(strcmp(qr, "")) {
            for(int i=atoi(ql); i<= atoi(qr); i++)
              tmp.push_back(i);
          }
          else
            tmp.push_back(atoi(p));
          buf = NULL;
        }

        _CMFD_group_structure.push_back(tmp);
        buf = NULL;
      }
      arg_index++;
    }
    else if(strcmp(argv[arg_index], "-verbose_report") == 0) {
      arg_index++;
      _verbose_report = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-time_report") == 0) {
      arg_index++;
      _time_report = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-test_run") == 0) {
      arg_index++;
      _test_run = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-quadraturetype") == 0) {
      arg_index++;
      _quadraturetype = atoi(argv[arg_index++]);
    }
    else if(strcmp(argv[arg_index], "-non_uniform_output") == 0) {
      arg_index++;
      char *buf = argv[arg_index];
      char *outer_ptr = NULL;
      char *inner_ptr = NULL;
      char *p;
      std::vector<std::vector<double> > widths_offset;
      while((p = strtok_r(buf, "/", &outer_ptr)) != NULL) {
        buf = p;
        char *ql, *qr;
        std::vector<double> tmp;
        while((p = strtok_r(buf, ",", &inner_ptr)) != NULL) {
          ql = strtok_r(p, "*", &qr);
          if(strcmp(qr, "")) {
            for(int i=0; i< atoi(qr); i++)
              tmp.push_back(atof(ql));
          }
          else
            tmp.push_back(atof(p));
          buf = NULL;
        }
        widths_offset.push_back(tmp);
        buf = NULL;
      }
      _non_uniform_mesh_lattices.push_back(widths_offset);
      arg_index++;
    }
    else {
      arg_index++;
    }
  }
  int myid;
#ifdef MPIx
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
  if ((_print_usage) && (myid == 0)) {
    printf("Usage: %s [<options>], default value in (). Example commands as "
           "below\n", argv[0]);
    printf("\n");
    printf(
      "mpiexec -n 8 ./run_time_standard                                    \\\n"
      "-debug                   1                                          \\\n"
      "-log_level               NORMAL                                     \\\n"
      "-domain_decompose        2,2,2                                      \\\n"
      "-num_domain_modules      1,1,1                                      \\\n"
      "-num_threads             1                                          \\\n"
      "-log_filename            test_problem.log                           \\\n"
      "-geo_filename            test_problem.geo                           \\\n"
      "-azim_spacing            0.10                                       \\\n"
      "-num_azim                32                                         \\\n"
      "-polar_spacing           0.5                                        \\\n"
      "-num_polar               6                                          \\\n"
      "-seg_zones               -1.0,2.0,3.0                               \\\n"
      "-segmentation_type       3                                          \\\n"
      "-quadraturetype          2                                          \\\n"
      "-CMFD_group_structure    1-3/4,5/6-8,9                              \\\n"
      "-CMFD_lattice            2,3,3                                      \\\n"
      "-widths_x                1.36,1.26*15,1.36*2                        \\\n"
      "-widths_y                1.36,1.26*15,1.36*2                        \\\n"
      "-widths_z                1.19*54                                    \\\n"
      "-CMFD_flux_update_on     1                                          \\\n"
      "-knearest                2                                          \\\n"
      "-CMFD_centroid_update_on 1                                          \\\n"
      "-use_axial_interpolation 2                                          \\\n"
      "-SOR_factor              1.5                                        \\\n"
      "-CMFD_relaxation_factor  0.7                                        \\\n"
      "-ls_solver               1                                          \\\n"
      "-max_iters               100                                        \\\n"
      "-MOC_src_residual_type   1                                          \\\n"
      "-MOC_src_tolerance       1.0E-5                                     \\\n"
      "-output_mesh_lattice     -output_mesh_lattice 5,5,9 -output_type 0  \\\n"
      "-output_mesh_lattice     -output_mesh_lattice 5,5,9 -output_type 1  \\\n"
      "-non_uniform_output      1.26*3/1*3/4.*3/-1.,1.,-1. -output_type 1  \\\n"
      "-verbose_report          1                                          \\\n"
      "-time_report             1                                          \\\n"
    );

    printf("\n");
    printf("General parameters\n");
    printf("-debug                  : (0) or 1, waits in while loop for GDB to"
           " attach\n");
    printf("-log_level              : (NORMAL)\n");
    printf("-domain_decompose       : (1,1,1) domain decomposition structure\n");
    printf("-num_domain_modules     : (1,1,1) modular structure in a domain\n");
    printf("-num_threads            : (1) Number of OpenMP threads to use\n");
    printf("-log_filename           : (NULL) the file name of the log file\n");
    printf("-geo_filename           : (NULL) the file name of the geometry "
           "file\n");
    printf("\n");

    printf("Track generating parameters\n");
    printf("-azim_spacing           : (0.05)\n");
    printf("-num_azim               : (64)\n");
    printf("-polar_spacing          : (0.75)\n");
    printf("-num_polar              : (10)\n");
    printf("-seg_zones              : (null) set the segmentation zones\n");
    printf("-segmentation_type      : (3-OTF_STACKS) 0-EXPLICIT_2D, "
           "1-EXPLICIT_3D, 2-OTF_TRACKS, 3-OTF_STACKS \n");
    printf("-quadraturetype         : (2 - GAUSS_LEGENDRE) is default value\n"
           "                           0 - TABUCHI_YAMAMOTO\n"
           "                           1 - LEONARD\n"
           "                           2 - GAUSS_LEGENDRE\n"
           "                           3 - EQUAL_WEIGHT\n"
           "                           4 - EQUAL_ANGLE\n"
          );
    printf("\n");

    printf("CMFD parameters\n");
    printf("-CMFD_group_structure   : (No group condensation) set CMFD group "
           "structure with ',', '-', and '/' \n");
    printf("-CMFD_lattice           : (0,0,0) Uniform CMFD lattice structure."
           "If both CMFD_lattice and widths are set,\n"
           "                          CMFD_lattice will be overridded \n");
    printf("-widths_x               : (NULL) the widths of non-uniform CMFD "
           "meshes in x direction, use '*' for repeat\n");
    printf("-widths_y               : (NULL) the widths of non-uniform CMFD "
           "meshes in y direction, use '*' for repeat\n");
    printf("-widths_z               : (NULL) the widths of non-uniform CMFD "
           "meshes in z direction, use '*' for repeat\n");
    printf("-CMFD_flux_update_on    : (1) switch of the CMFD update\n");
    printf("-knearest               : (1) number of CMFD knearest neighbors\n");
    printf("-CMFD_centroid_update_on: (1) switch of the CMFD knearest "
           "centroid update\n");
    printf("-use_axial_interpolation: (0) option of the CMFD quadratic axial "
           "interpolation update\n"
           "                           0 - No axial interpolation\n"
           "                           1 - FSR axially averaged value\n"
           "                           2 - centroid z-coordinate evaluated "
           "value\n");
    printf("-SOR_factor             : (1.0) set CMFD SOR relaxation factor\n");
    printf("-CMFD_relaxation_factor : (1.0) set CMFD relaxation factor\n");
    printf("\n");

    printf("MOC solver parameters\n");
    printf("-ls_solver              : (1) set the linear source solver\n");
    printf("-max_iters              : (1000) Maximum number of outter "
           "iterations\n");
    printf("-MOC_src_residual_type  : (1-FISSION_SOURCE) 0-SCALAR_FLUX, "
           "1-FISSION_SOURCE, 2-TOTAL_SOURCE\n");
    printf("-MOC_src_tolerance      : (1.0E-4) MOC source convergence "
           "tolerance\n");
    printf("\n");

    printf("Output parameters\n");
    printf("-output_mesh_lattice    : (0,0,0) Uniform reaction output mesh "
           "lattice\n");
    printf("-non_uniform_output     : set the XYZ widths and offset of "
           "non_uniform lattice for reaction output, use '*' for repeat\n");
    printf("-output_type            : (0 - FISSION_RX) set the output reaction "
           "types, an output_type should always follow output lattice \n"
           "                           0 - FISSION_RX\n"
           "                           1 - NU-FISSION_RX\n"
           "                           2 - TOTAL_RX\n"
           "                           3 - ABSORPTION_RX\n"
           "                           4 - FLUX_RX\n"
          );
    printf("-verbose_report         : (1) switch of the verbose iteration "
           "report\n");
    printf("-time_report            : (1) switch of the time report\n");
    printf("-test_run               : (0) switch of the test running mode\n");

    printf("\n");
  }

  if (_print_usage) {
#ifdef MPIx
    MPI_Finalize();
#endif
  }
  return 0;
}
