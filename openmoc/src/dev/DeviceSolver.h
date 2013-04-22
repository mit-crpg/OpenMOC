/*
 * DeviceSolver.h
 *
 *  Created on: Aug 5, 2012
 *      Author: will
 */

#ifndef DEVICESOLVER_H_
#define DEVICESOLVER_H_

#include "../host/Solver.h"
#include "../host/Plotter.h"
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include "DeviceFlatSourceRegion.h"
#include "DeviceTrack.h"
#include "DeviceMaterial.h"
#include <cutil_math.h>


class DeviceSolver {

private:

	/* host pointers */
	Geometry* _geom;
	Quadrature* _quad;
	TrackGenerator* _track_generator;
	Plotter* _plotter;
	int _B;
	int _T;
	FP_PRECISION* _FSRs_to_fluxes[NUM_ENERGY_GROUPS + 1];
	FP_PRECISION* _FSRs_to_powers;
	FP_PRECISION* _FSRs_to_pin_powers;


	/* device pointers */
	dev_flatsourceregion* _FSRs;
    dev_material* _materials;
    Track** _host_tracks;
    dev_track* _dev_tracks;
	int* _num_tracks;
	int* _track_index_offsets;
    FP_PRECISION* _fission_source;
	thrust::device_vector<FP_PRECISION> _fission_source_vec;
	FP_PRECISION* _renorm_factor;
	FP_PRECISION* _k_eff;
	FP_PRECISION* _tot_abs;
	FP_PRECISION* _tot_fission;
	FP_PRECISION* _source_residual_norm;
    thrust::device_vector<FP_PRECISION> _tot_abs_vec;
	thrust::device_vector<FP_PRECISION> _tot_fission_vec;
	thrust::device_vector<FP_PRECISION> _source_residual_norm_vec;
#if !STORE_PREFACTORS
	FP_PRECISION* _prefactor_array;
#endif

public:

	DeviceSolver(Geometry* geom, TrackGenerator* track_generator,
								Plotter* plotter, Options* options);
	virtual ~DeviceSolver();

	void allocateDeviceMemory();
	void validateDeviceDataIntegrity();
	void initializeFSRs(dev_flatsourceregion* dev_FSRs);
	int computeScalarTrackIndex(int i, int j);
	void computePrefactors();

	void fixedSourceIteration(int max_iterations);
	FP_PRECISION computeKeff(int max_iterations);
	void plotFluxes();
	void computePinPowers();
};


#if STORE_PREFACTORS
__global__ void computePrefactorsOnDevice(dev_track* tracks,
										dev_material* materials);
#endif

__device__ FP_PRECISION computePrefactorOnDevice(FP_PRECISION sigma_t,
						FP_PRECISION length, FP_PRECISION sintheta);

__global__ void zeroTrackFluxesOnDevice(dev_track* tracks);

__global__ void normalizeTrackFluxesOnDevice(dev_track* tracks,
												FP_PRECISION* renorm_factor);

__global__ void initDeviceFSRsOnDevice(dev_flatsourceregion* FSRs);

__global__ void zeroFSRFluxesOnDevice(dev_flatsourceregion* FSRs);

__global__ void computeTotalFissionSourceOnDevice(dev_flatsourceregion* FSRs,
												dev_material* material,
												FP_PRECISION* fission_source);

__global__ void computeFSRSourcesOnDevice(dev_flatsourceregion* dFSRs,
										dev_material* material,
										FP_PRECISION* curr_keff,
										FP_PRECISION* source_residual_norm);

__global__ void normalizeFSRFluxesOnDevice(dev_flatsourceregion* dev_FSRs,
													FP_PRECISION* renorm_factor);

__global__ void computeAbsAndFissionOnDevice(dev_flatsourceregion* dev_FSRs,
							dev_material* material,
							FP_PRECISION* tot_abs, FP_PRECISION* tot_fission);

__global__ void normalizeFluxToVolumeOnDevice(dev_flatsourceregion* FSRs,
												dev_material* materials);

#if STORE_PREFACTORS

__global__ void propagateTrackFluxForwardOnDevice(dev_track* tracks,
										int* track_index_offsets,
										dev_flatsourceregion* FSRs);

__global__ void propagateTrackFluxReverseOnDevice(dev_track* tracks,
													int* track_index_offsets,
													dev_flatsourceregion* FSRs);

#else

__global__ void propagateTrackFluxForwardOnDevice(dev_track* dev_tracks,
												int* track_index_offsets,
												dev_material* materials,
												dev_flatsourceregion* dFSRs,
												FP_PRECISION* prefactor_array);

__global__ void propagateTrackFluxReverseOnDevice(dev_track* dev_tracks,
												int* track_index_offsets,
												dev_material* materials,
												dev_flatsourceregion* FSRs,
												FP_PRECISION* prefactor_array);
#endif

__device__ double atomicAdd(double* address, double val);


#endif /* DEVICESOLVER_H_ */
