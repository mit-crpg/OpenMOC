#include "Vector.h"

/**
 * @brief Constructor initializes Vector object as a floating point array
 *        and sets the vector dimensions.
 * @details The vector is ordered by cell (as opposed to by group) on the
 *          outside to be consistent with the Matrix object. Locks are used to
 *          make the vector object thread-safe against concurrent writes the
 *          same value. One lock locks out multiple rows of the vector at a
 *          time representing multiple groups in the same cell.
 * @param cell_locks OpenMP locks for atomic cell operations.
 * @param num_x The number of cells in the x direction.
 * @param num_y The number of cells in the y direction.
 * @param num_z The number of cells in the z direction.
 * @param num_groups The number of energy groups in each cell.
 */
Vector::Vector(omp_lock_t* cell_locks, int num_x, int num_y, int num_z,
               int num_groups) {

  //TODO Change initialization of Vector to clean Valgrind log
  setNumX(num_x);
  setNumY(num_y);
  setNumZ(num_z);
  setNumGroups(num_groups);
  _num_rows = _num_x*_num_y*_num_z*_num_groups;

  /* Initialize array and set all to 0.0 */
  _array = new CMFD_PRECISION[_num_rows];
  setAll(0.0);

  /* Set OpenMP locks for each Vector cell */
  if (cell_locks == NULL)
    log_printf(ERROR, "Unable to create a Vector without an array of cell "
               "locks");

  _cell_locks = cell_locks;
}


/**
 * @brief Destructor deletes the arrays used to represent the vector.
 */
Vector::~Vector() {

  if (_array != NULL)
    delete [] _array;
}


/**
 * @brief Increment a value in the vector.
 * @details This method takes a cell and group and floating
 *          point value. The cell and group are used to compute the
 *          row in the vector. If a value exists for the row,
 *          the value is incremented by val; otherwise, it is set to val.
 * @param cell The cell location.
 * @param group The group location.
 * @param val The value used to increment the row location.
 */
void Vector::incrementValue(int cell, int group, CMFD_PRECISION val) {

  if (cell >= _num_x*_num_y*_num_z || cell < 0)
    log_printf(ERROR, "Unable to increment Vector value for cell %d"
               " which is not between 0 and %d", cell, _num_x*_num_y*_num_z-1);
  else if (group >= _num_groups || group < 0)
    log_printf(ERROR, "Unable to increment Vector value for group %d"
               " which is not between 0 and %d", group, _num_groups-1);

  /* Atomically increment the Vector value from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell]);

  _array[cell*_num_groups + group] += val;

  /* Release Vector cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell]);
#ifdef INTEL
#pragma omp flush
#endif
}


/**
 * @brief Increment values in the vector.
 * @details This method takes a cell, first group, last group, and floating
 *          point value. The cell and groups are used to compute the
 *          rows in the vector. If values exist for the rows,
 *          the values are incremented by vals; otherwise, they are set.
 * @param cell The cell location.
 * @param group_first The first group location to increment.
 * @param group_last The last group location to increment.
 * @param vals The values used to increment the row locations.
 *        NOTE: vals at group_first must be aligned with cache boundaries.
 */
void Vector::incrementValues(int cell, int group_first, int group_last,
                             CMFD_PRECISION* vals) {

  if (cell >= _num_x*_num_y*_num_z || cell < 0)
    log_printf(ERROR, "Unable to increment Vector values for cell %d"
               " which is not between 0 and %d", cell, _num_x*_num_y*_num_z-1);
  else if (group_first >= _num_groups || group_first < 0)
    log_printf(ERROR, "Unable to increment Vector values for first group %d"
               " which is not between 0 and %d", group_first, _num_groups-1);
  else if (group_last >= _num_groups || group_last < 0)
    log_printf(ERROR, "Unable to increment Vector values for last group %d"
               " which is not between 0 and %d", group_last, _num_groups-1);
  else if (group_first > group_last)
    log_printf(ERROR, "Unable to increment Vector values with first group %d"
               " greater than last group %d", group_first, group_last);

  /* Atomically increment the Vector values from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell]);

#pragma omp simd aligned(vals)
  for (int g=group_first; g <= group_last; g++)
    _array[cell*_num_groups + g] += vals[g-group_first];

  /* Release Vector cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell]);
#ifdef INTEL
#pragma omp flush
#endif
}


/**
 * @brief Fill vector with a value.
 * @param val value to use to fill
 */
void Vector::setAll(CMFD_PRECISION val) {
  std::fill_n(_array, _num_rows, val);
}


/**
 * @brief Set a value in the vector.
 * @details This method takes a cell and group and floating
 *          point value. The cell and group are used to compute the
 *          row and column in the vector. The location of the corresponding
 *          row is set to val.
 * @param cell The cell location.
 * @param group The group location.
 * @param val The value used to set the row location.
 */
void Vector::setValue(int cell, int group, CMFD_PRECISION val) {

  if (cell >= _num_x*_num_y*_num_z || cell < 0)
    log_printf(ERROR, "Unable to set Vector value for cell %d"
               " which is not between 0 and %d", cell, _num_x*_num_y*_num_z-1);
  else if (group >= _num_groups || group < 0)
    log_printf(ERROR, "Unable to set Vector value for group %d"
               " which is not between 0 and %d", group, _num_groups-1);

  /* Atomically set the Vector value from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell]);

  _array[cell*_num_groups + group] = val;

  /* Release Vector cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell]);
#ifdef INTEL
#pragma omp flush
#endif
}


/**
 * @brief Set values in the vector.
 * @details This method takes a cell, first group, last group, and floating
 *          point value. The cell and groups are used to compute the
 *          rows in the vector. If a value exist for the rows,
 *          the values are overwritten.
 * @param cell The cell location.
 * @param group_first The first group location to set.
 * @param group_last The last group location to set.
 * @param vals The values used to set the row locations.
 */
void Vector::setValues(int cell, int group_first, int group_last,
                       CMFD_PRECISION* vals) {

  if (cell >= _num_x*_num_y*_num_z || cell < 0)
    log_printf(ERROR, "Unable to set Vector values for cell %d"
               " which is not between 0 and %d", cell, _num_x*_num_y*_num_z-1);
  else if (group_first >= _num_groups || group_first < 0)
    log_printf(ERROR, "Unable to set Vector values for first group %d"
               " which is not between 0 and %d", group_first, _num_groups-1);
  else if (group_last >= _num_groups || group_last < 0)
    log_printf(ERROR, "Unable to set Vector values for last group %d"
               " which is not between 0 and %d", group_last, _num_groups-1);
  else if (group_first > group_last)
    log_printf(ERROR, "Unable to set Vector values with first group %d"
               " greater than last group %d", group_first, group_last);

  /* Atomically set the Vector values from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell]);

#pragma omp simd
  for (int g=group_first; g <= group_last; g++)
    _array[cell*_num_groups + g] = vals[g-group_first];

  /* Release Vector cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell]);
#ifdef INTEL
#pragma omp flush
#endif
}


/**
 * @brief Clear all values in the vector.
 */
void Vector::clear() {
  setAll(0.0);
}


/**
 * @brief Scales the vector by a given value.
 * @param val The value to scale the vector by.
 */
void Vector::scaleByValue(CMFD_PRECISION val) {

#pragma omp parallel for schedule(guided)
  for (int i=0; i < _num_rows; i++)
    _array[i] *= val;
}


/**
 * @brief Print the vector object to the log file.
 */
void Vector::printString() {

  std::stringstream string;
  string << std::setprecision(6);

  string << std::endl;
  string << "Vector" << std::endl;
  string << " Num rows: " << _num_rows << std::endl;

  for (int row=0; row < _num_rows; row++)
    string << " ( " << row << "): " << _array[row] << std::endl;

  string << "End Vector" << std::endl;

  std::cout << string.str() << std::endl;
  //log_printf(NORMAL, string.str().c_str());
}


/**
 * @brief Copy the values from the current vector to an input vector.
 * @param vector The vector to copy values to.
 */
void Vector::copyTo(Vector* vector) {
  std::copy(_array, _array + _num_rows, vector->getArray());
}


/**
 * @brief Get the array describing the vector.
 * @return The array describing the vector.
 */
CMFD_PRECISION* Vector::getArray() {
  return _array;
}


/**
 * @brief Get the number of cells in the x dimension.
 * @return The number of cells in the x dimension.
 */
int Vector::getNumX() {
  return _num_x;
}


/**
 * @brief Get the number of cells in the y dimension.
 * @return The number of cells in the y dimension.
 */
int Vector::getNumY() {
  return _num_y;
}


/**
 * @brief Get the number of cells in the z dimension.
 * @return The number of cells in the z dimension.
 */
int Vector::getNumZ() {
  return _num_z;
}


/**
 * @brief Get the number of groups in each cell.
 * @return The number of groups in each cell.
 */
int Vector::getNumGroups() {
  return _num_groups;
}


/**
 * @brief Get the number of rows in the vector.
 * @return The number of rows in the vector.
 */
int Vector::getNumRows() {
  return _num_rows;
}


/**
 * @brief Get the sum of all the values in the vector.
 * @return The sum of all the values in the vector.
 */
double Vector::getSum() {
  return pairwise_sum(_array, _num_rows);
}


/**
 * @brief Get the number of negative values in the vector.
 * @return The number of negative values in the vector.
 */
long  Vector::getNumNegativeValues() {

  long num_negative_values = 0;
#pragma omp simd reduction(+:num_negative_values)
  for (long i=0; i<_num_rows; i++) 
    num_negative_values += (_array[i] < 0);

  return num_negative_values;
}


/**
 * @brief Set the number of cells in the x dimension.
 * @param num_x The number of cells in the x dimension.
 */
void Vector::setNumX(int num_x) {

  if (num_x < 1)
    log_printf(ERROR, "Unable to set Vector num x to non-positive value %d",
               num_x);

  _num_x = num_x;
}


/**
 * @brief Set the number of cells in the y dimension.
 * @param num_y The number of cells in the y dimension.
 */
void Vector::setNumY(int num_y) {

  if (num_y < 1)
    log_printf(ERROR, "Unable to set Vector num y to non-positive value %d",
               num_y);

  _num_y = num_y;
}


/**
 * @brief Set the number of cells in the z dimension.
 * @param num_z The number of cells in the z dimension.
 */
void Vector::setNumZ(int num_z) {

  if (num_z < 1)
    log_printf(ERROR, "Unable to set Vector num z to non-positive value %d",
               num_z);

  _num_z = num_z;
}


/**
 * @brief Set the number of groups in each cell.
 * @param num_groups The number of groups in each cell.
 */
void Vector::setNumGroups(int num_groups) {

  if (num_groups < 1)
    log_printf(ERROR, "Unable to set Vector num groups to non-positive value"
               " %d", num_groups);

  _num_groups = num_groups;
}


/**
 * @brief Return the array of cell locks for atomic cell operations.
 * @return an array of cell locks
 */
omp_lock_t* Vector::getCellLocks() {
  return _cell_locks;
}
