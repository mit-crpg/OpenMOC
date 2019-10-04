#include<iostream>
#include "dev_exponential.h"
using namespace std;

__global__ void compute_exp_on_gpu(float* x_dev)
{
  // *x_dev = 1.-__expf(*x_dev);
  *x_dev = dev_exponential(*x_dev);
}

int main()
{
  float x;
  float* dev_x;
  while (1)
  {
    cout << "enter number:\n" << endl;
    cin >> x;
    cudaMalloc((void**)&dev_x, sizeof(float));
    cudaMemcpy((void**)dev_x, (void**)&x, sizeof(float),
        cudaMemcpyHostToDevice);
    compute_exp_on_gpu<<<1,1>>>(dev_x);
    cudaMemcpy((void**)&x, (void**)dev_x, sizeof(float),
        cudaMemcpyDeviceToHost);
    cout << "Device calculated: " << x << endl;
  }
}
