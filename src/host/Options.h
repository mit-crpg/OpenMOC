/*
 * Options.h
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */


#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <sstream>
#include <string.h>
#include <sys/stat.h>
#include <omp.h>
#include "log.h"
#include "configurations.h"

#define LAST(str) (strcmp(argv[i-1], (str)) == 0)
#define THIS(str) (strcmp(argv[i], (str)) == 0)


class Options {

private:

	std::string _input_path;
	std::string _relative_path;
	std::string _geometry_file;
	std::string _material_file;
	std::string _extension;

	double _track_spacing;
	int _num_azim;
	int _num_omp_threads;
	int _bit_dimension;
	std::string _verbosity;
	bool _dump_geometry;
	bool _plot_specs;
	bool _plot_fluxes;
	bool _compute_pin_powers;
	bool _compute_on_cpu;
	bool _compute_on_gpu;
	bool _use_track_inputfile;
	std::string _track_input_filename;

	int _num_gpu_blocks;
	int _num_gpu_threads;


public:

	Options(int argc, const char **argv);
    ~Options(void);
    const char *getGeometryFile() const;
    const char *getMaterialFile() const;
    bool dumpGeometry() const;
    int getNumAzim() const;
    int getBitDimension() const;
    double getTrackSpacing() const;
    const char* getVerbosity() const;
    std::string getExtension() const;
    std::string getTrackInputFilename();
    bool useTrackInputFile();
    bool plotSpecs() const;
    bool plotFluxes() const;
    bool computePinPowers() const;
    bool computeOnCPU() const;
    bool computeOnGPU() const;
    int getNumThreadBlocks() const;
    int getNumThreadsPerBlock() const;

};

#endif
