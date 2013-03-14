/*
 * gpua-moc.cu
 *
 *  Created on: Apr 24, 2012
 *      Author: wboyd
 */

/* Host code */
#include "../host/Geometry.h"
#include "../host/TrackGenerator.h"
#include "../host/Parser.h"
#include "../host/Options.h"
#include "../host/configurations.h"
#include "../host/Plotter.h"
#include "../host/log.h"
#include "../host/Solver.h"

/* Device code */
//#include <cutil.h>
#include "DeviceQuery.h"
#include "Timer.h"
#include "DeviceSolver.h"

int main(int argc, const char **argv) {

	/* Create an options class to parse command line options */
	Options opts(argc, argv);

	/* Set the verbosity */
	log_setlevel(opts.getVerbosity());

/******************************************************************************
 ***************************  Start OpenMOC  **********************************
 *****************************************************************************/

	log_printf(NORMAL, "Starting OpenMOC...");

	FP_PRECISION host_k_eff, dev_k_eff;
	Timer timer;

	/* Initialize the parser and time the parser */
	timer.startHostTimer();
	Parser parser(&opts);
	timer.stopHostTimer();
	timer.recordHostSplit("Parsing input files");

	/* Initialize the geometry with surfaces, cells & materials */
	timer.resetHostTimer();
	timer.startHostTimer();
	Geometry geometry(&parser);
	timer.stopHostTimer();
	timer.recordHostSplit("Geomery initialization");

	/* Print out geometry to console if requested at runtime*/
	if (opts.dumpGeometry())
		geometry.printString();

	Plotter plotter(&geometry, opts.getBitDimension(), opts.getExtension(),
			opts.plotSpecs(), opts.plotFluxes());

	/* Initialize the trackgenerator */
	TrackGenerator track_generator(&geometry, &plotter, opts.getNumAzim(),
				       opts.getTrackSpacing());

	/* Generate tracks */
	timer.resetHostTimer();
	timer.startHostTimer();
	track_generator.generateTracks(opts.useTrackInputFile(),
									opts.getTrackInputFilename());
	timer.stopHostTimer();
	timer.recordHostSplit("Generating tracks");

	/* Segment tracks */
	timer.resetHostTimer();
	timer.startHostTimer();
	track_generator.segmentize();
	timer.stopHostTimer();
	timer.recordHostSplit("Segmenting tracks");


	if (opts.computeOnCPU()) {
		/* Host fixed source iteration to solve for k_eff */
		Solver solver(&geometry, &track_generator, &plotter);

		timer.resetHostTimer();
		timer.startHostTimer();
		host_k_eff = solver.computeKeff(MAX_ITERATIONS);
		timer.stopHostTimer();
		timer.recordHostSplit("Fixed source iteration on host");

		/* Compute pin powers if requested at run time */
		if (opts.computePinPowers())
			solver.computePinPowers();
	}

	/* Device fixed source iteration to solve for k_eff */
	if (opts.computeOnGPU()) {

		machineContainsGPU();
		attachGPU();
		printBasicDeviceInfo();
		printDetailedDeviceInfo();

		DeviceSolver devicesolver(&geometry, &track_generator,
												&plotter, &opts);

		devicesolver.validateDeviceDataIntegrity();
		timer.resetHostTimer();
		timer.startHostTimer();
		dev_k_eff = devicesolver.computeKeff(MAX_ITERATIONS);
		timer.stopHostTimer();
		timer.recordHostSplit("Fixed source iteration on device");

		if (opts.computePinPowers())
			devicesolver.computePinPowers();
	}


	/* Print the results from fixed source iteration */
	if (opts.computeOnCPU())
		log_printf(RESULT, "Host k_eff = %f", host_k_eff);
	if (opts.computeOnGPU())
		log_printf(RESULT, "Device k_eff = %f", dev_k_eff);


	/* Print timer splits to console */
	log_printf(NORMAL, "Program complete");
	timer.printSplits();

	return 1;
}
