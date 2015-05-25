#include "Matrix.h"

Matrix::Matrix(int num_x, int num_y, int num_z, int num_groups){

  _num_x = num_x;
  _num_y = num_y;
  _num_z = num_z;
  _num_groups = num_groups;
  _num_rows = _num_x*_num_y*_num_z*_num_groups;
  
  /* Initialize variables */
  for (int i=0; i < _num_rows; i++){
    std::map<int, double> *values = new std::map<int, double>;
    _LIL.push_back(*values);
  }

  _A = NULL;
  _IA = NULL;
  _JA = NULL;
  _DIAG = NULL;
  _modified = true;
}


Matrix::~Matrix(){

  if (_A != NULL)
    delete [] _A;

  if (_IA != NULL)
    delete [] _IA;

  if (_JA != NULL)
    delete [] _JA;

  if (_DIAG != NULL)
    delete [] _DIAG;
  
  for (int i=0; i < _num_rows; i++)
    _LIL[i].clear();
}


void Matrix::incrementValue(int row, int col, double val){
  _LIL[row][col] += val;
  _modified = true;
}


void Matrix::incrementValueByCoords(int x_from, int y_from, int z_from, int g_from,
                                    int x_to, int y_to, int z_to, int g_to, double val){
  int row = ((z_to*_num_y + y_to)*_num_x + x_to)*_num_groups + g_to;
  int col = ((z_from*_num_y + y_from)*_num_x + x_from)*_num_groups + g_from;
  incrementValue(row, col, val);
}


void Matrix::incrementValueByCell(int cell_from, int g_from,
                                  int cell_to, int g_to, double val){
  int row = cell_to*_num_groups + g_to;
  int col = cell_from*_num_groups + g_from;
  incrementValue(row, col, val);
}


void Matrix::setValue(int row, int col, double val){
  _LIL[row][col] = val;
  _modified = true;
}


void Matrix::setValueByCoords(int x_from, int y_from, int z_from, int g_from,
                              int x_to, int y_to, int z_to, int g_to, double val){
  int row = ((z_to*_num_y + y_to)*_num_x + x_to)*_num_groups + g_to;
  int col = ((z_from*_num_y + y_from)*_num_x + x_from)*_num_groups + g_from;
  setValue(row, col, val);
}


void Matrix::setValueByCell(int cell_from, int g_from,
                            int cell_to, int g_to, double val){
  int row = cell_to*_num_groups + g_to;
  int col = cell_from*_num_groups + g_from;
  setValue(row, col, val);
}


void Matrix::clear(){
  for (int i=0; i < _num_rows; i++)
    _LIL[i].clear();

  _modified = true;
}


void Matrix::convertToCSR(){
  
  /* Get number of nonzero values */
  int NNZ = getNNZ();

  /* Allocate memory for arrays */
  if (_A != NULL)
    delete [] _A;

  if (_IA != NULL)
    delete [] _IA;

  if (_JA != NULL)
    delete [] _JA;

  if (_DIAG != NULL)
    delete [] _DIAG;

  _A = new double[NNZ];
  _IA = new int[_num_rows+1];
  _JA = new int[NNZ];
  _DIAG = new double[_num_rows];
  std::fill_n(_DIAG, _num_rows, 0.0);  
  
  /* Form arrays */

  int j = 0;
  std::map<int, double>::iterator iter;
  for (int row=0; row < _num_rows; row++){
    _IA[row] = j;
    for (iter = _LIL[row].begin(); iter != _LIL[row].end(); ++iter){
      if (iter->second != 0.0){
        _JA[j] = iter->first;
        _A[j] = iter->second;
        
        if (row == iter->first)
          _DIAG[row] = iter->second;
        
        j++;
      }
    }
  }

  _IA[_num_rows] = NNZ;

  _modified = false;
}



void Matrix::printString(){

  /* Convert to CSR form */
  convertToCSR();

  log_printf(NORMAL, "Matrix");

  std::cout << " Num rows: " << _num_rows << std::endl;
  std::cout << " NNZ     : " << getNNZ() << std::endl;

  for (int row=0; row < _num_rows; row++){
    for (int i = _IA[row]; i < _IA[row+1]; i++){
      std::cout << std::setprecision(6) << " ( " << row << ", " << _JA[i] << "): " << _A[i] << std::endl;
    }
  }

  log_printf(NORMAL, "End Matrix");
}


double Matrix::getValue(int row, int col){
  return _LIL[row][col];
}


double Matrix::getValueByCoords(int x_from, int y_from, int z_from, int g_from,
                                int x_to, int y_to, int z_to, int g_to){
  int row = ((z_to*_num_y + y_to)*_num_x + x_to)*_num_groups + g_to;
  int col = ((z_from*_num_y + y_from)*_num_x + x_from)*_num_groups + g_from;
  return getValue(row, col);
}


double Matrix::getValueByCell(int cell_from, int g_from,
                              int cell_to, int g_to){
  int row = cell_to*_num_groups + g_to;
  int col = cell_from*_num_groups + g_from;
  return getValue(row, col);
}


double* Matrix::getA(){

  if (_modified)
    convertToCSR();
  
  return _A;
}


int* Matrix::getIA(){

  if (_modified)
    convertToCSR();

  return _IA;
}


int* Matrix::getJA(){

  if (_modified)
    convertToCSR();
    
  return _JA;
}


double* Matrix::getDIAG(){

  if (_modified)
    convertToCSR();
    
  return _DIAG;
}


int Matrix::getNumX(){
  return _num_x;
}


int Matrix::getNumY(){
  return _num_y;
}


int Matrix::getNumZ(){
  return _num_z;
}


int Matrix::getNumGroups(){
  return _num_groups;
}


int Matrix::getNumRows(){
  return _num_rows;
}


void Matrix::random(){
  for (int i=0; i < _num_rows; i++){
    for (int j=0; j < _num_rows; j++){
      if (i == j)
        setValue(i, j, static_cast<double>(rand()) / static_cast<double>(RAND_MAX / (_num_rows*10)) + _num_rows);
      else
        setValue(i, j, -static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
    }
  }
}


int Matrix::getNNZ(){

  int NNZ = 0;
  std::map<int, double>::iterator iter;
  for (int row=0; row < _num_rows; row++){
    for (iter = _LIL[row].begin(); iter != _LIL[row].end(); ++iter){
      if (iter->second != 0.0)
        NNZ++;
    }
  }

  return NNZ;
}
