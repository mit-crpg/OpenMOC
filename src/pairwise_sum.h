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
template <typename T>
inline T pairwise_sum(T* vector, int length) {

    T sum = 0;

    /* Base case: if length is less than 16, perform summation */
    if (length < 16) {
        #pragma simd reduction(+:sum)
	for (int i=0; i < length; i++)
	    sum += vector[i];
    }

    else {
	int offset = length % 2;
        length = floor(length / 2);
	sum = pairwise_sum<T>(&vector[0], length) + 
	      pairwise_sum<T>(&vector[length], length+offset);	
    }

    return sum;
}
#ifdef MIC
#pragma offload_attribute(pop)
#endif
