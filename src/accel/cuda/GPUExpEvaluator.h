/**
 * @file GPUExpEvaluator.h
 * @brief The GPUExpEvaluator class.
 * @details This class is based on the ExpEvaluator class for execution
 *          on NVIDIA GPUs.
 * @date April 11, 2015.
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GPUEXPEVALUATOR_H_
#define GPUEXPEVALUATOR_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include "../../ExpEvaluator.h"
#endif


/**
 * @class GPUExpEvaluator GPUExpEvaluator.h "src/accel/cuda/ExpEvaluator.h"
 * @brief This is a class for evaluating exponentials on GPUs.
 * @details The ExpEvaluator includes different algorithms to evaluate
 *          exponentials with varying degrees of accuracy and speed. This
 *          is a helper class for the Solver and its subclasses and it not
 *          intended to be initialized as a standalone object.
 */
class GPUExpEvaluator {

private:

public:

  /** The exponential linear interpolation table */
  FP_PRECISION* _exp_table;

  __device__ FP_PRECISION computeExponential(FP_PRECISION tau, int polar);
};


void clone_exp_evaluator(ExpEvaluator* evaluator_h,
                         GPUExpEvaluator* evaluator_d);


#endif /* EXPEVALUATOR_H_ */
