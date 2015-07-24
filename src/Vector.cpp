#include "Vector.h"

Vector::Vector(int num_x, int num_y, int num_groups) {

  setNumX(num_x);
  setNumY(num_y);
  setNumGroups(num_groups);  
  _num_rows = _num_x*_num_y*_num_groups;

  /* Initialize array and set all to 0.0 */
  _array = new FP_PRECISION[_num_rows];
  setAll(0.0);

  /* Allocate memory for OpenMP locks for each Vector cell */ 
  _cell_locks = new omp_lock_t[_num_x*_num_y];

  /* Loop over all Vector cells to initialize OpenMP locks */
  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_x*_num_y; r++)
    omp_init_lock(&_cell_locks[r]);
}


Vector::~Vector() {

  if (_array != NULL)
    delete [] _array;

  if (_cell_locks != NULL)
    delete [] _cell_locks;  
}


void Vector::incrementValue(int cell, int group, FP_PRECISION val) {

  if (cell >= _num_x*_num_y || cell < 0)
    log_printf(ERROR, "Unable to increment Vector value for cell %i"
               " which is not between 0 and %i", cell, _num_x*_num_y-1);
  else if (group >= _num_groups || group < 0)
    log_printf(ERROR, "Unable to increment Vector value for group %i"
               " which is not between 0 and %i", group, _num_groups-1);

  /* Atomically increment the Vector value from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell]);

  _array[cell*_num_groups + group] += val;

  /* Release Vector cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell]);
}


void Vector::setAll(FP_PRECISION val) {
  std::fill_n(_array, _num_rows, val);
}


void Vector::setValue(int cell, int group, FP_PRECISION val) {

  if (cell >= _num_x*_num_y || cell < 0)
    log_printf(ERROR, "Unable to set Vector value for cell %i"
               " which is not between 0 and %i", cell, _num_x*_num_y-1);
  else if (group >= _num_groups || group < 0)
    log_printf(ERROR, "Unable to set Vector value for group %i"
               " which is not between 0 and %i", group, _num_groups-1);

  /* Atomically set the Vector value from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell]);

  _array[cell*_num_groups + group] = val;

  /* Release Vector cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell]);
}


void Vector::clear() {
  setAll(0.0);
}


void Vector::scaleByValue(FP_PRECISION val) {

  #pragma omp parallel for schedule(guided)
  for (int i=0; i < _num_rows; i++)
    _array[i] *= val;
}


std::string Vector::toString() {

  std::stringstream string;
  string << std::setprecision(6);

  string << std::endl;
  string << "Vector" << std::endl;
  string << " Num rows: " << _num_rows << std::endl;

  for (int row=0; row < _num_rows; row++)
    string << " ( " << row << "): " << _array[row] << std::endl;

  string << "End Vector" << std::endl;

  return string.str();
}


void Vector::printString() {
  log_printf(NORMAL, toString().c_str());
}


void Vector::copyTo(Vector* vector) {
  std::copy(_array, _array + _num_rows, vector->getArray());
}


FP_PRECISION Vector::getValue(int cell, int group) {
  return _array[cell*_num_groups + group];
}


FP_PRECISION* Vector::getArray() {
  return _array;
}


int Vector::getNumX() {
  return _num_x;
}


int Vector::getNumY() {
  return _num_y;
}


int Vector::getNumGroups() {
  return _num_groups;
}


int Vector::getNumRows() {
  return _num_rows;
}


FP_PRECISION Vector::getSum() {
  return pairwise_sum(_array, _num_rows);
}


void Vector::setNumX(int num_x) {

  if (num_x < 1)
    log_printf(ERROR, "Unable to set Vector num x to non-positive value %i",
               num_x);

  _num_x = num_x;
}


void Vector::setNumY(int num_y) {

  if (num_y < 1)
    log_printf(ERROR, "Unable to set Vector num y to non-positive value %i",
               num_y);

  _num_y = num_y;
}


void Vector::setNumGroups(int num_groups) {

  if (num_groups < 1)
    log_printf(ERROR, "Unable to set Vector num groups to non-positive value"
               " %i", num_groups);

  _num_groups = num_groups;
}
