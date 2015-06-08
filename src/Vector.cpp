#include "Vector.h"

Vector::Vector(int num_x, int num_y, int num_z, int num_groups){

  _num_x = num_x;
  _num_y = num_y;
  _num_z = num_z;
  _num_groups = num_groups;
  
  _num_rows = _num_x*_num_y*_num_z*_num_groups;

  _array = new FP_PRECISION[_num_rows];
  std::fill_n(_array, _num_rows, 0.0);

}


Vector::~Vector(){

  if (_array != NULL)
    delete [] _array;
}


void Vector::incrementValue(int row, FP_PRECISION val){
  log_printf(DEBUG, "incrementing row: %i val: %f", row, val);
  _array[row] += val;
}


void Vector::incrementValueByCoords(int x, int y, int z, int g, FP_PRECISION val){
  int row = ((z*_num_y + y)*_num_x + x)*_num_groups + g;
  incrementValue(row, val);
}


void Vector::incrementValueByCell(int cell, int g, FP_PRECISION val){
  int row = cell*_num_groups + g;
  incrementValue(row, val);
}


void Vector::setAll(FP_PRECISION val){
  std::fill_n(_array, _num_rows, val);
}


void Vector::setValue(int row, FP_PRECISION val){
  _array[row] = val;
}


void Vector::setValueByCoords(int x, int y, int z, int g, FP_PRECISION val){
  int row = ((z*_num_y + y)*_num_x + x)*_num_groups + g;
  setValue(row, val);
}


void Vector::setValueByCell(int cell, int g, FP_PRECISION val){
  int row = cell*_num_groups + g;
  setValue(row, val);
}


void Vector::clear(){
  setAll(0.0);
}


void Vector::scaleByValue(FP_PRECISION val){

  for (int i=0; i < _num_rows; i++)
    _array[i] *= val;
}


std::string Vector::toString(){

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


void Vector::printString(){
  log_printf(NORMAL, toString().c_str());
}


void Vector::copyTo(Vector* vector){
  std::copy(_array, _array + _num_rows, vector->getArray());
}


FP_PRECISION Vector::getValue(int row){
  return _array[row];
}


FP_PRECISION Vector::getValueByCoords(int x, int y, int z, int g){
  int row = ((z*_num_y + y)*_num_x + x)*_num_groups + g;
  return getValue(row);
}


FP_PRECISION Vector::getValueByCell(int cell, int g){
  int row = cell*_num_groups + g;
  return getValue(row);
}


FP_PRECISION* Vector::getArray(){
  return _array;
}


int Vector::getNumX(){
  return _num_x;
}


int Vector::getNumY(){
  return _num_y;
}


int Vector::getNumZ(){
  return _num_z;
}


int Vector::getNumGroups(){
  return _num_groups;
}


int Vector::getNumRows(){
  return _num_rows;
}


void Vector::random(){
  for (int i=0; i < _num_rows; i++)
    _array[i] = static_cast<FP_PRECISION>(rand()) / static_cast<FP_PRECISION>(RAND_MAX);
}


FP_PRECISION Vector::sum(){
  return pairwise_sum(_array, _num_rows);
}
