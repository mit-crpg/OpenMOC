/*
 * Options.cpp
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */


#include "Options.h"


/**
 * Options constructor
 * @param argc the number of command line arguments from console
 * @param argv a char array of command line arguments from console
 */
Options::Options(int argc, const char **argv) {

	_input_path = "";

	/* Default track spacing */
	_track_spacing = 0.1;

	/* Default number of azimuthal angles */
	_num_azim = 16;

	/* Default number of OpenMP threads */
	_num_omp_threads = 1;

	/* Default logging level */
	_verbosity = "NORMAL";

	/* Default will not dump geometry to the terminal */
	_dump_geometry = false;

	/* Default will create plots with .png format */
	_extension = "png";

	/* Pixel dimensions for plots */
	_bit_dimension = 1000;

	/* Default will not plot materials, cells, FSRs, tracks, or segments */
	_plot_specs = false;

	/* Default will not plot fluxes */
	_plot_fluxes = false;

	/* Default will not compute pin powers */
	_compute_pin_powers = false;

	/* Default will not compute on the host */
	_compute_on_cpu = false;

	/* Default will not compute on the device */
	_compute_on_gpu = false;

	/* Blocks per kernel call on the GPU */
	_num_gpu_blocks =64;

	/* Threads per block on the GPU */
	_num_gpu_threads = 64;


	/* Iterate over the command line arguments passed in */
	for (int i = 0; i < argc; i++) {

		if (i > 0) {

			if (LAST("--inputpath") || LAST("-ip"))
				_input_path = argv[i];

			else if (LAST("--trackspacing") || LAST("-ts"))
				_track_spacing = atof(argv[i]);

			else if (LAST("--numazimuthal") || LAST("-na"))
				_num_azim = atoi(argv[i]);

			else if (LAST("--numompthreads") || LAST("-nomp"))
				_num_omp_threads = atoi(argv[i]);

			else if (LAST("--bitdimension") || LAST("-bd"))
							_bit_dimension = atoi(argv[i]);

			else if (LAST("--verbosity") || LAST("-v"))
				_verbosity = strdup(argv[i]);

			else if (THIS("--dumpgeometry") || THIS("-dg"))
				_dump_geometry = true;

			else if (LAST("--extension") || LAST("-ex"))
							_extension = argv[i];

			else if (THIS("--plotspecs") || THIS("-ps"))
				_plot_specs = true;

			else if (THIS("--plotfluxes") || THIS("-pf"))
				_plot_fluxes = true;

			else if (THIS("--computepowers") || THIS("-cp"))
				_compute_pin_powers = true;

			else if (THIS("-cpu"))
				_compute_on_cpu = true;

			else if (THIS("-gpu"))
				_compute_on_gpu = true;

			else if (LAST("--numblocks") || LAST("-B"))
				_num_gpu_blocks = atoi(argv[i]);

			else if (LAST("--numthreads") || LAST("-T"))
				_num_gpu_threads = atoi(argv[i]);

		}
	}

	_num_omp_threads = std::min(_num_omp_threads, _num_azim/4);
	omp_set_num_threads(_num_omp_threads);

	if (_input_path == "") {

		/* Checks the working directory to set the relative path for input files
		 * This is important so that default input files work when program is run
		 * from both eclipse and the console */
		if (std::string(getenv("PWD")).find("Release") != std::string::npos)
			_relative_path = "../";
		else
			_relative_path = "";

		/* Default path to the directory with input files */
		_input_path = "xml-sample/SimpleLattice/";

		/* Default geometry input file */
		_geometry_file = _relative_path + _input_path + "geometry.xml";

		/* Default material input file */
		_material_file = _relative_path + _input_path + "material.xml";

	}
	else {
		/* User-defined geometry input file */
		_geometry_file = _input_path + "geometry.xml";

		/* User-defined material input file */
		_material_file = _input_path + "material.xml";
	}
}


/**
 * Default destructor does not have anything to free
 */
Options::~Options(void) { }


/**
 * Returns a character array with the path to the geometry input file.
 * By default this will return /xml-sample/SimpleLattice/geometry.xml
 * if path is not set at runtime from the console
 * @return path to the geometry input file
 */
const char *Options::getGeometryFile() const {
    return _geometry_file.c_str();
}


/**
 * Returns a character array with the path to the material input file.
 * By default this will return /xml-sample/SimpleLattice/material.xml
 * if path is not set at runtime from the console
 * @return path to the geometry input file
 */
const char *Options::getMaterialFile() const {
    return _material_file.c_str();
}


/**
 * Returns a boolean representing whether or not to dump the geometry to the
 * console. If true, the geometry will be printed out after parsing is complete
 * @return whether or not to dump the geometry to the console
 */
bool Options::dumpGeometry() const {
	return _dump_geometry;
}


/**
 * Returns the number of azimuthal angles. By default this will return
 * 128 angles if not set at runtime from the console
 * @return the number of azimuthal angles
 */
int Options::getNumAzim() const {
    return _num_azim;
}


/**
 * Returns the pixel dimensions for plots. By default this will return
 * 1000 bits (or pixels) if not set at runtime from the console
 * @return the y dimension of plots.
 */
int Options::getBitDimension() const{
	return _bit_dimension;
}


/**
 * Returns the track spacing. By default this will return 0.05 if not
 * set at runtime from the console
 * @return the track spacing
 */
double Options::getTrackSpacing() const {
    return _track_spacing;
}


/**
 * Returns the verbosity logging level. By default this will return NORMAL
 * if not set at runtime from the console
 * @return the verbosity
 */
const char* Options::getVerbosity() const {
    return _verbosity.c_str();
}


/**
 * Returns the image files extension. By default this will return png
 * if not set at runtime from the console
 * @return the image files extension
 */
std::string Options::getExtension() const {
    return _extension.c_str();
}


std::string Options::getTrackInputFilename() {

	std::stringstream test_filename;

	int path_end_index = _material_file.find("material.xml");
	test_filename << _material_file.substr(0, path_end_index)
						<< _num_azim << "_angles_"
						<< _track_spacing << "_spacing.tracks";

    _track_input_filename = test_filename.str();

	return _track_input_filename;
}


bool Options::useTrackInputFile() {

	getTrackInputFilename();

	struct stat buf;
	if (stat(_track_input_filename.c_str(), &buf) != -1)
		return true;
	else
		return false;
}



/**
 * Returns a boolean representing whether or not to plot the specs.
 * If true, the specs will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotSpecs() const {
	return _plot_specs;
}

/**
 * Returns a boolean representing whether or not to plot the cells.
 *  If true, the cells will be plotted in a file of _extension type
 * @return whether or not to plot materials
 */
bool Options::plotFluxes() const {
	return _plot_fluxes;
}


/**
 * Returns a boolean representing whether or not to compute the powers
 * in each pin. If true, txt files with the pin powers will be created
 * in a new directory called "PinPowers"
 * @return whether or not to compute the pin powers
 */
bool Options::computePinPowers() const {
	return _compute_pin_powers;
}

/**
 * Returns a boolean representing whether or not to perform fixed source
 * iteration calculation on the CPU. By default this will return false
 * if not set at runtime from the console
 * @return whether or not to compute flux and keff on the CPU
 */
bool Options::computeOnCPU() const {
	return _compute_on_cpu;
}


/**
 * Returns a boolean representing whether or not to perform fixed source
 * iteration calculation on the GPU. By default this will return false
 * if not set at runtime from the console
 * @return whether or not to compute flux and keff on the GPU
 */
bool Options::computeOnGPU() const {
	return _compute_on_gpu;
}


/**
 * Returns the number of blocks of threads to execute on the GPU for
 * each DeviceSolver kernel call
 * @return number of GPU thread blocks
 */
int Options::getNumThreadBlocks() const {
	return _num_gpu_blocks;
}


/**
 * Returns the number of threads per block to execute on the GPU for
 * each DeviceSolver kernel call
 * @return number of GPU threads per block
 */
int Options::getNumThreadsPerBlock() const {
	return _num_gpu_threads;
}
