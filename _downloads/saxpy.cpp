/* File saxpy.cpp */
#include "saxpy.h"

/* Define global variables */
int length;
double a;
double* x;
double* y;

void set_array_length(int n) {
  length = n;
}

void initialize_data() {

  /* Allocate memory for arrays */
  x = (double*)malloc(length*sizeof(double));
  y = (double*)malloc(length*sizeof(double));

  /* Initialize data with random numbers in [0,1] */
  a = float(rand()) / RAND_MAX; 

  for (int i=0; i < length; i++) {
    x[i] = float(rand()) / RAND_MAX;
    y[i] = float(rand()) / RAND_MAX;
  }
}

void free_data() {
  free(x);
  free(y);
}

void print_data() {
  printf("a = %f\n", a);

  for (int i=0; i < length; i++)
    printf("x[%d] = %f\ty[%d] = %f\n", i, x[i], i, y[i]);
}

void saxpy() {
  for (int i=0; i < length; i++)
    y[i] = a * x[i] + y[i];
}
