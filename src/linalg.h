/**
 * @file linalg.h
 * @brief Utility function for performing linear algebra in CMFD solver.
 * @author Samuel Shaner (shaner@mit.edu)
 * @date August 26, 2014
 */

/**
 * @brief Copy a vector to another vector.
 * @param vector_from vector to be copied
 * @param vector_to vector to receive copied data
 * @param length length of vector
 */
template <typename T>
inline void vector_copy(T* vector_from, T* vector_to, int length){

  for (int i = 0; i < length; i++)
    vector_to[i] = vector_from[i];
}


/**
 * @brief Assign all elements in a matrix to zero.
 * @param matrix matrix to be zeroed
 * @param width width of matrix row
 * @param length length of matrix column
 */
template <typename T>
inline void matrix_zero(T** matrix, int width, int length){

  for (int i = 0; i < length; i++){
    for (int g = 0; g < width; g++)
      matrix[i][g] = 0.0;
  }
}


/**
 * @brief Assign all elements in a matrix to zero.
 * @param vector vector to be zeroed
 * @param length length of vector
 */
template <typename T>
inline void vector_zero(T* vector, int length){

  for (int i = 0; i < length; i++)
    vector[i] = 0.0;
}


/**
 * @brief Multiply matrix by vector (i.e., y = M *x).
 * @param matrix source matrix
 * @param vector_x x vector
 * @param vector_y y vector
 * @param num_blocks number of cell blocks in M matrix.
 * @param block_width number of elements in cell blocks in M matrix.
 */
template <typename T>
inline void matrix_multiplication(T** matrix, T* vector_x, 
                                  T* vector_y, int num_blocks, 
                                  int block_width){

  vector_zero(vector_y, num_blocks*block_width); 

  for (int i = 0; i < num_blocks; i++){
    for (int g = 0; g < block_width; g++){
      for (int e = 0; e < block_width; e++){
        vector_y[i*block_width+g] += matrix[i][g*block_width+e] 
            * vector_x[i*block_width+e];
      }
    }
  }
}


/**
 * @brief Scale vectgor by some scalar value.
 * @param vector vector to be scaled
 * @param scale_value value to scale vector
 * @param length vector length
 */
template <typename T>
inline void vector_scale(T* vector, T scale_value, int length){

  for (int i = 0; i < length; i++)
    vector[i] *= scale_value;
}
