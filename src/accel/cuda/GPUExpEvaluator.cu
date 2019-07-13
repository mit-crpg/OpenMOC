#include "GPUExpEvaluator.h"


/** A boolean whether to use linear interpolation to compute exponentials */
__constant__ bool interpolate[1];

/** The maximum allowable optical length represented in the table */
__constant__ FP_PRECISION max_optical_length[1];

/** The inverse spacing for the exponential linear interpolation table */
__constant__ FP_PRECISION inverse_exp_table_spacing[1];

/** An array for the sines of the polar angle in the polar Quadrature set */
extern __constant__ FP_PRECISION sin_thetas[MAX_POLAR_ANGLES_GPU];

/** Twice the number of polar angles */
extern __constant__ int num_polar[1];


/**
 * @brief Given a pointer to an ExpEvaluator on the host and a
 *        GPUExpEvaluator on the GPU, copy all of the properties from
*         the ExpEvaluator object on the host to the GPU.
 * @details This routine is called by the GPUSolver::initializeExpEvaluator()
 *          private class method and is not intended to be called directly.
 * @param eavluator_h pointer to a ExpEvaluator on the host
 * @param evaluator_d pointer to a GPUExpEvaluator on the GPU
 */
void clone_exp_evaluator(ExpEvaluator* evaluator_h,
                         GPUExpEvaluator* evaluator_d) {

  /* Copy a boolean indicating whether or not to use the linear interpolation
   * table or the exp intrinsic function to constant memory on the device */
  bool interpolate_exp = evaluator_h->isUsingInterpolation();
  cudaMemcpyToSymbol(interpolate, (void*)&interpolate_exp,
                     sizeof(bool), 0, cudaMemcpyHostToDevice);

  if (evaluator_h->isUsingInterpolation()) {

    /* Copy inverse table spacing to constant memory on the device */
    FP_PRECISION inverse_spacing_h = 1.0 / evaluator_h->getTableSpacing();
    cudaMemcpyToSymbol(inverse_exp_table_spacing, (void*)&inverse_spacing_h,
                       sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

    /* Copy the number of table entries to constant memory on the device */
    FP_PRECISION max_optical_length_h = evaluator_h->getMaxOpticalLength();
    cudaMemcpyToSymbol(max_optical_length, (void*)&max_optical_length_h,
               sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

    /* Allocate memory for the interpolation table on the device */
    int exp_table_size_h = evaluator_h->getTableSize();
    FP_PRECISION* exp_table_h = evaluator_h->getExpTable();

    FP_PRECISION* exp_table_d;
    cudaMalloc((void**)&exp_table_d, exp_table_size_h * sizeof(FP_PRECISION));
    cudaMemcpy((void*)exp_table_d, (void*)exp_table_h,
               exp_table_size_h * sizeof(FP_PRECISION),
               cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&evaluator_d->_exp_table, (void*)&exp_table_d,
               sizeof(FP_PRECISION*), cudaMemcpyHostToDevice);
  }

  return;
}


/**
 * @brief Computes the exponential term for a optical length and polar angle.
 * @details This method computes \f$ 1 - exp(-\tau/sin(\theta_p)) \f$
 *          for some optical path length and polar angle. This method
 *          uses either a linear interpolation table (default) or the
 *          exponential intrinsic exp(...) function.
 * @param tau the optical path length (e.g., sigma_t times length)
 * @param polar the polar angle index
 * @return the evaluated exponential
 */
__device__ FP_PRECISION GPUExpEvaluator::computeExponential(FP_PRECISION tau,
                                                            int polar) {

  FP_PRECISION exponential;

  /* Evaluate the exponential using the linear interpolation table */
  if (*interpolate) {
    tau = min(tau, (*max_optical_length));
    int index = floor(tau * (*inverse_exp_table_spacing));
    index *= (*num_polar);
    exponential = (1. - (_exp_table[index + 2 * polar] * tau +
                         _exp_table[index + 2 * polar +1]));
  }

  /* Evalute the exponential using the intrinsic exp(...) function */
  else {
    #ifdef SINGLE
    exponential = 1.0 - __expf(- tau / sin_thetas[polar]);
    #else
    exponential = 1.0 - exp(- tau / sin_thetas[polar]);
    #endif
  }

  return exponential;
}
