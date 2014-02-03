/* File print_list.cpp */
#include "print_list.h"

void print_list(int length, double* list) {
  printf("Printing a Python list from C/C++\n");

  /* Loop over each list value and print to the screen */
  for (int i=0; i < length; i++)
    printf("list[%d] = %f\n", i, list[i]);
}
