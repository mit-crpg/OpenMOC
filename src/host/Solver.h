/*
 * Solver.h
 *
 *  Created on: Feb 7, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "../host/Quadrature.h"
#include "../host/TrackGenerator.h"
#include "../host/FlatSourceRegion.h"



class Solver {
private:

	/* host class attributes */
	Geometry* _geom;
	Quadrature* _quad;
	FlatSourceRegion* _FSRs;
	Track** _tracks;
	int* _num_tracks;
	int _num_azim;
    int _tot_num_tracks;
	int _num_FSRs;
	FP_PRECISION* _FSRs_to_fluxes[NUM_ENERGY_GROUPS + 1];
	FP_PRECISION* _FSRs_to_powers;
	FP_PRECISION* _FSRs_to_pin_powers;
	FP_PRECISION _k_eff;
	Plotter* _plotter;

#if !STORE_PREFACTORS
	FP_PRECISION* _prefactor_array;
	int _prefactor_array_size;
	int _prefactor_max_index;
	FP_PRECISION _prefactor_spacing;
#endif

	void precomputePrefactors();
    FP_PRECISION computePrefactor(FP_PRECISION sigma_t, FP_PRECISION length,
    													FP_PRECISION sintheta);
	void initializeFSRs();

public:

	Solver(Geometry* geom, TrackGenerator* track_generator, Plotter* plotter);
	virtual ~Solver();
	void zeroTrackFluxes();
	void initFSRs();
	void zeroFSRFluxes();
	void renormalizeFluxes();
	FP_PRECISION computeFSRSources();
	void computeRatios();
	void updateFSRFluxes();
	void updateKeff();
	FP_PRECISION** getFSRtoFluxMap();
	void fixedSourceIteration(int max_iterations);
	FP_PRECISION computeKeff(int max_iterations);
	void plotFluxes();
	void checkTrackSpacing();
	void computePinPowers();

};

#endif /* SOLVER_H_ */
