#include "dev_exponential.h"

/**
 * @brief Computes the exponential term for a optical length and polar angle.
 * @details This method computes \f$ 1 - exp(-x) \f$
 *          for some optical path length and polar angle. Uses the rational
            approximation of Josey from exponentials.h, but using device constants.
 * @param tau the optical path length (e.g., sigma_t times length)
 * @param polar the polar angle index
 * @return the evaluated exponential
 */
__device__ FP_PRECISION dev_exponential(FP_PRECISION x)
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

  // exponential = num * den * x;
  // TODO:
  exponential = 1.0 - __expf(-x);

  return exponential;
}
