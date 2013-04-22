#include "../host/log.h"
#include <cutil.h>

#ifndef DEVICEQUERY_H_
#define DEVICEQUERY_H_

#include "DeviceQuery.h"

bool machineContainsGPU() {
	int count;

	CUDA_SAFE_CALL(cudaGetDeviceCount(&count));

	log_printf(NORMAL, "%d CUDA-enabled devices are present", count);

	if (count > 0)
		return true;
	else
		return false;
}


void attachGPU() {
	cudaDeviceReset();
	int dev = 0;
	CUDA_SAFE_CALL(cudaSetDevice(dev));
	return;
}



void printBasicDeviceInfo() {
	int dev = 0;
	cudaDeviceProp prop;

	CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, dev));

	log_printf(NORMAL, "Device name: %s", prop.name);
	log_printf(NORMAL, "Device compute capability: %d.%d", prop.major, prop.minor);
	log_printf(NORMAL, "Device # multiprocessors: %d", prop.multiProcessorCount);
	log_printf(NORMAL, "Device clock rate: %d", prop.clockRate);

}
void printDetailedDeviceInfo() {

	int dev = 0;
	cudaDeviceProp prop;

	CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, dev));

	log_printf(NORMAL, "Device total global memory: %ld", prop.totalGlobalMem);
	log_printf(NORMAL, "Device total constant memory: %ld", prop.totalConstMem);

	log_printf(NORMAL, "Shared memory per multiprocessor: %ld", prop.sharedMemPerBlock);
	log_printf(NORMAL, "Registers per multiprocessor: %d", prop.regsPerBlock);
	log_printf(NORMAL, "Threads in warp: %d", prop.warpSize);
	log_printf(NORMAL, "Max threads per multiprocessor: %d", prop.maxThreadsPerMultiProcessor);
	log_printf(NORMAL, "Max threads per block: %d", prop.maxThreadsPerBlock);
	log_printf(NORMAL, "Max thread dimensions: [%d, %d, %d]", prop.maxThreadsDim[0],
									prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
	log_printf(NORMAL, "Max grid dimensions: [%d, %d, %d]", prop.maxGridSize[0],
									prop.maxGridSize[1], prop.maxGridSize[2]);

}


#endif /* DEVICEQUERY_H_ */
