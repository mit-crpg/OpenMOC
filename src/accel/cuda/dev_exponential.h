/**
 * Compute exponentials on the GPU using some of Colin Josey's finest work.
 */
#pragma once

 /* Coefficients for numerator */
__device__ const FP_PRECISION p0 = 1.0;
__device__ const FP_PRECISION p1 = 2.4172687328033081 * 1E-1;
__device__ const FP_PRECISION p2 = 6.2804790965268531 * 1E-2;
__device__ const FP_PRECISION p3 = 1.0567595009016521 * 1E-2;
__device__ const FP_PRECISION p4 = 1.0059468082903561 * 1E-3;
__device__ const FP_PRECISION p5 = 1.9309063097411041 * 1E-4;

/* Coefficients for denominator */
__device__ const FP_PRECISION d0 = 1.0;
__device__ const FP_PRECISION d1 = 7.4169266112320541 * 1E-1;
__device__ const FP_PRECISION d2 = 2.6722515319494311 * 1E-1;
__device__ const FP_PRECISION d3 = 6.1643725066901411 * 1E-2;
__device__ const FP_PRECISION d4 = 1.0590759992367811 * 1E-2;
__device__ const FP_PRECISION d5 = 1.0057980007137651 * 1E-3;
__device__ const FP_PRECISION d6 = 1.9309063097411041 * 1E-4;

/**
 * @brief Computes the exponential term for a optical length and polar angle.
 * @details This method computes \f$ 1 - exp(-x) \f$
 *          for some optical path length and polar angle. Uses the rational
            approximation of Josey from exponentials.h, but using device constants.
 * @param tau the optical path length (e.g., sigma_t times length)
 * @param polar the polar angle index
 * @return the evaluated exponential
 */
__forceinline__ __device__ FP_PRECISION dev_exponential(FP_PRECISION x)
{
  FP_PRECISION exponential;
  
  FP_PRECISION num, den;
 
  den = d6*x + d5;
  den = den*x + d4;
  den = den*x + d3;
  den = den*x + d2;
  den = den*x + d1;
  den = den*x + d0;
  den = d0 / den; 
 
  num = p5*x + p4;
  num = num*x + p3;
  num = num*x + p2;
  num = num*x + p1;
  num = num*x + p0;
 
  // TODO Redo fractional approximation for 1-exp(-x) rather than (1-exp(-x))/x
  exponential = num * den * x;
 
  return exponential;
}
