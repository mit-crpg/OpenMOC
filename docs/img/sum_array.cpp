/* File sum_array.cpp */

/* Define function implementation */
double sum_array(double* input_array, int length) {

  /* Initialize sum */
  double sum = 0.;

  /* Compute sum of array elements */
  for (int i=0; i < length; i++)
    sum += input_array[i];

  return sum;
}
