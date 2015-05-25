#include "Vector.h"

Vector::Vector(int num_x, int num_y, int num_z, int num_groups){

  _num_x = num_x;
  _num_y = num_y;
  _num_z = num_z;
  _num_groups = num_groups;
  
  _num_rows = _num_x*_num_y*_num_z*_num_groups;

  _array = new double[_num_rows];
  std::fill_n(_array, _num_rows, 0.0);

}


Vector::~Vector(){

  if (_array != NULL)
    delete [] _array;
}


void Vector::incrementValue(int row, double val){
  _array[row] += val;
}


void Vector::incrementValueByCoords(int x, int y, int z, int g, double val){
  int row = ((z*_num_y + y)*_num_x + x)*_num_groups + g;
  incrementValue(row, val);
}


void Vector::incrementValueByCell(int cell, int g, double val){
  int row = cell*_num_groups + g;
  incrementValue(row, val);
}


void Vector::setAll(double val){
  std::fill_n(_array, _num_rows, val);
}


void Vector::setValue(int row, double val){
  _array[row] = val;
}


void Vector::setValueByCoords(int x, int y, int z, int g, double val){
  int row = ((z*_num_y + y)*_num_x + x)*_num_groups + g;
  setValue(row, val);
}


void Vector::setValueByCell(int cell, int g, double val){
  int row = cell*_num_groups + g;
  setValue(row, val);
}


void Vector::clear(){
  setAll(0.0);
}


void Vector::scaleByValue(double val){

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


double Vector::getValue(int row){
  return _array[row];
}


double Vector::getValueByCoords(int x, int y, int z, int g){
  int row = ((z*_num_y + y)*_num_x + x)*_num_groups + g;
  return getValue(row);
}


double Vector::getValueByCell(int cell, int g){
  int row = cell*_num_groups + g;
  return getValue(row);
}


double* Vector::getArray(){
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
    _array[i] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}


double Vector::sum(){
  return pairwise_sum(_array, _num_rows);
}
