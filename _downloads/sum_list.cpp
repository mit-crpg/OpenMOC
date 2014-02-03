/* File sum_list.cpp */
#include "sum_list.h"

double sum_list(double* input_list, int length) {

  /* Initialize sum */
  double sum = 0.;

  /* Compute sum of array elements */
  for (int i=0; i < length; i++)
    sum += input_list[i];

  return sum;
}
