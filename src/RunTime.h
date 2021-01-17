/**
 * @file RunTime.h
 * @brief Utility functions for processing run time options.
 * @details OpenMOC provide another option of input by .cpp file. With these run
 *          time options, it's possible to run various problems with different
 *          parameters via a standard built executive file, avoiding the
 *          necessity to build the code each time.
 * @author Wenbin Wu (wenbin@mit.edu)
 * @date October 12, 2018
 *
 */

#ifndef RUNTIME_H_
#define RUNTIME_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include <iostream>
#include <string.h>
#include <vector>
#endif

#ifdef MPIx
#include <mpi.h>
#endif

#ifdef SWIG
#define printf PySys_WriteStdout
#endif

/* A Vector3D is simply a 3-dimensional std::vector of doubles */
typedef std::vector<std::vector<std::vector<double> > > DoubleVector3D;

/**
 * @brief Structure for run time options.
 */
struct RuntimeParameters {
  RuntimeParameters() : _print_usage(false), _debug_flag(false), _NDx(1),
    _NDy(1), _NDz(1), _NMx(1), _NMy(1), _NMz(1), _NCx(0), _NCy(0), _NCz(0),
    _num_threads(1), _azim_spacing(0.05), _num_azim(64), _polar_spacing(0.75),
    _num_polar(10), _tolerance(1.0E-4), _max_iters(1000), _knearest(1),
    _CMFD_flux_update_on(true), _CMFD_centroid_update_on(false),
    _use_axial_interpolation(0), _log_filename(NULL), _linear_solver(true),
    _MOC_src_residual_type(1), _SOR_factor(1.0), _CMFD_relaxation_factor(1.0),
    _segmentation_type(3), _verbose_report(true), _time_report(true),
    _log_level((char*)"NORMAL"), _quadraturetype(2), _test_run(false) {}

  /* Only display help message */
  bool _print_usage;

  /* To debug or not when running, infinite while loop */
  bool _debug_flag;

  /* Level of verbosity for the execution */
  char* _log_level;

  /* Domain decomposition structure */
  int _NDx, _NDy, _NDz;

  /* Modules structure, used to define sub-domains */
  int _NMx, _NMy, _NMz;

  /* Number of OpenMP threads */
  int _num_threads;

  /* Log file name */
  char* _log_filename;

  /* Geometry file name */
  std::string _geo_filename;

  /* Space and angle quadrature parameters */
  double _azim_spacing;
  int _num_azim;
  double _polar_spacing;
  int _num_polar;

  /* Segmentation zones for 2D extruded segmentation */
  std::vector<double> _seg_zones;

  /* Segmentation type of track generation */
  int _segmentation_type;

  /* Polar quadrature type */
  int _quadraturetype;

  /* CMFD group structure */
  std::vector<std::vector<int> > _CMFD_group_structure;

  /* CMFD lattice structure, used for uniform CMFD */
  int _NCx, _NCy, _NCz;

  /** Physical dimensions of non-uniform CMFD meshes (for whole geometry) */
  std::vector<double> _cell_widths_x;
  std::vector<double> _cell_widths_y;
  std::vector<double> _cell_widths_z;

  /* Whether to update MOC fluxes with the CMFD transport acceleration */
  bool _CMFD_flux_update_on;

  /* The order of the k-nearest update */
  int _knearest;

  /* K-nearest update or conventional update */
  bool _CMFD_centroid_update_on;

  /* Whether to use axial interpolation for the CMFD update */
  int _use_axial_interpolation;

  /* CMFD linear solver SOR factor */
  double _SOR_factor;

  /* CMFD relaxation factor */
  double _CMFD_relaxation_factor;

  /* Linear source solver if true */
  bool _linear_solver;

  /* The maximum number of MOC source iterations */
  int _max_iters;

  /* Type of MOC source residual for assessing convergence */
  int _MOC_src_residual_type;

  /* MOC source residual convergence criterion */
  double _tolerance;

  /* Mesh dimension for reaction rate output on a uniform lattice */
  std::vector<std::vector<int> > _output_mesh_lattices;

  /* widths and offsets of multiple output meshes with non-uniform lattice */
  DoubleVector3D _non_uniform_mesh_lattices;

  /* Reaction types for mesh outputs, both uniform and non-uniform */
  std::vector<int> _output_types;

  /* Print a verbose report at every transport iteration */
  bool _verbose_report;

  /* Print a solver execution time report */
  bool _time_report;

  /* Whether to run the code for a test */
  bool _test_run;

  /* Setter, parses command line input */
  int setRuntimeParameters(int argc, char *argv[]);
};

#endif /* RUNTIME_H_ */
