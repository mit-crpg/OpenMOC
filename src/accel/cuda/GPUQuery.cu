#include "GPUQuery.h"


/**
 * @brief Queries a node to determine whether it contains one or more GPUs.
 * @return True if the node contains a GPU, false otherwise.
 */
bool machine_contains_gpu() {

  int count;
  cudaGetDeviceCount(&count);

  if (count > 0)
    return true;
  else
    return false;
}


/**
 * @brief Sets the primary CUDA-enabled device to be the GPU with a given ID.
 * @param id the ID for the GPU to attache
 */
void attach_gpu(int id) {

  if (!machine_contains_gpu()) {
    log_printf(WARNING, "Unable to attach GPU %d since no GPU is "
               "attached to the machine", id);
    return;
  }

  cudaDeviceReset();
  cudaSetDevice(id);

  return;
}


/**
 * @brief Prints the basic device info for the CUDA-enabled device in use.
 * @details Prints the name, compute capability, # multiprocessors and
 *          the clock rate of the device.
 */
void print_basic_gpu_info() {

  if (!machine_contains_gpu()) {
    log_printf(WARNING, "Unable to print basic GPU info since no GPU"
               " is attached to the machine");
    return;
  }

  int dev;
  cudaGetDevice(&dev);

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, dev);

  log_printf(NORMAL, "Device name: %s", prop.name);
  log_printf(NORMAL, "Device compute capability: %d.%d",prop.major, prop.minor);
  log_printf(NORMAL, "Device # multiprocessors: %d", prop.multiProcessorCount);
  log_printf(NORMAL, "Device clock rate: %d", prop.clockRate);

}


/**
 * @brief Prints the detailed device info for the CUDA-enabled device in use.
 * @details Prints the total global and constant memory, shared memory and
 *          registers per multiprocessor, # threads per warp, maximum #
 *          threads per multiprocessor, maximum # threads per block, maximum
 *          threadblock dimensions, and maximum grid dimensions.
 */
void print_detailed_gpu_info() {

  if (!machine_contains_gpu()) {
    log_printf(WARNING, "Unable to print detailed device info since no GPU"
               " is attached to the machine");
    return;
  }

  int dev;
  cudaGetDevice(&dev);

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, dev);

  log_printf(NORMAL, "Device total global memory: %ld", prop.totalGlobalMem);
  log_printf(NORMAL, "Device total constant memory: %ld", prop.totalConstMem);

  log_printf(NORMAL, "Shared memory per multiprocessor: %ld",
             prop.sharedMemPerBlock);
  log_printf(NORMAL, "Registers per multiprocessor: %d", prop.regsPerBlock);
  log_printf(NORMAL, "Threads in warp: %d", prop.warpSize);
  log_printf(NORMAL, "Max threads per multiprocessor: %d",
             prop.maxThreadsPerMultiProcessor);
  log_printf(NORMAL, "Max threads per block: %d", prop.maxThreadsPerBlock);
  log_printf(NORMAL, "Max thread dimensions: [%d, %d, %d]",
             prop.maxThreadsDim[0], prop.maxThreadsDim[1],
             prop.maxThreadsDim[2]);
  log_printf(NORMAL, "Max grid dimensions: [%d, %d, %d]", prop.maxGridSize[0],
             prop.maxGridSize[1], prop.maxGridSize[2]);
}


/**
 * @brief Returns the number of threads in a CUDA warp for the attached GPU.
 */
int get_num_threads_per_warp() {

  if (!machine_contains_gpu()) {
    log_printf(WARNING, "Unable to return the number of threads per warp "
               "since no GPU is attached to the machine");
    return 0;
  }

  int dev;
  cudaGetDevice(&dev);

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, dev);

  return prop.warpSize;
}
