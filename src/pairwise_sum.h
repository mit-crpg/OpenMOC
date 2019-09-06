/**
 * @file pairwise_sum.h
 * @brief Utility function for the accurate pairwise sum of a list of floating
 *        point numbers.
 * @author William Boyd (wboyd@mit.edu)
 * @date June 13, 2013
 */

/**
 * @brief Performs a pairwise sum of an array of numbers.
 * @details This type of summation uses a divide-and-conquer algorithm which
 *          is necessary to bound the error for summations of large sequences
 *          of numbers.
 * @param vector an array of numbers
 * @param length the length of the array
 * @return the sum of all numbers in the array
 */
template <typename T, typename L>
inline double pairwise_sum(T* vector, L length) {

  double sum = 0;

  /* Base case: if length is less than SIMD vector length, perform summation */
  if (length < VEC_LENGTH) {

#pragma omp simd reduction(+:sum)
    for (L i=0; i < length; i++)
      sum += vector[i];
  }

  else {
    L offset = length % 2;
    length = floor(length / 2);
    sum = pairwise_sum<T>(&vector[0], length) +
          pairwise_sum<T>(&vector[length], length+offset);
  }

  return sum;
}
