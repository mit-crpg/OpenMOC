/* File get_rand_list.cpp */
#include "get_rand_list.h"

/* Define function implementation */
std::vector<double> get_rand_list(int length) {

  /* Allocate memory for the C++ STL vector */
  std::vector<double> output_list(length);

  /* Populate vector with random numbers */
  for (int i=0; i < length; i++)
    output_list[i] = ((double) rand()) / RAND_MAX;

  return output_list;
}
