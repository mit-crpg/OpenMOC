#include "Matrix.h"

Matrix::Matrix(int num_x, int num_y, int num_groups){

  setNumX(num_x);
  setNumY(num_y);
  setNumGroups(num_groups);  
  _num_rows = _num_x*_num_y*_num_groups;
  
  /* Initialize variables */
  for (int i=0; i < _num_rows; i++){
    std::map<int, FP_PRECISION> *values = new std::map<int, FP_PRECISION>;
    _LIL.push_back(*values);
  }

  _A = NULL;
  _IA = NULL;
  _JA = NULL;
  _DIAG = NULL;
  _modified = true;

  /* Allocate memory for OpenMP locks for each Matrix cell */ 
  _cell_locks = new omp_lock_t[_num_x*_num_y];

  /* Loop over all Matrix cells to initialize OpenMP locks */
  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_x*_num_y; r++)
    omp_init_lock(&_cell_locks[r]);  
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

  if (_cell_locks != NULL)
    delete [] _cell_locks;
}


void Matrix::incrementValue(int cell_from, int group_from,
                            int cell_to, int group_to, FP_PRECISION val){

  if (cell_from >= _num_x*_num_y || cell_from < 0)
    log_printf(ERROR, "Unable to increment Matrix value for cell_from %i"
               " which is not between 0 and %i", cell_from, _num_x*_num_y-1);
  else if (cell_to >= _num_x*_num_y || cell_to < 0)
    log_printf(ERROR, "Unable to increment Matrix value for cell_to %i"
               " which is not between 0 and %i", cell_from, _num_x*_num_y-1);
  else if (group_from >= _num_groups || group_from < 0)
    log_printf(ERROR, "Unable to increment Matrix value for group_from %i"
               " which is not between 0 and %i", group_from, _num_groups-1);
  else if (group_to >= _num_groups || group_to < 0)
    log_printf(ERROR, "Unable to increment Matrix value for group_to %i"
               " which is not between 0 and %i", group_to, _num_groups-1);
  
  /* Atomically increment the Matrix value from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell_to]);

  int row = cell_to*_num_groups + group_to;
  int col = cell_from*_num_groups + group_from;
  _LIL[row][col] += val;
  
  /* Release Matrix cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell_to]);

  /* Set global modified flag to true */
  _modified = true;
}


void Matrix::setValue(int cell_from, int group_from,
                      int cell_to, int group_to, FP_PRECISION val){

  if (cell_from >= _num_x*_num_y || cell_from < 0)
    log_printf(ERROR, "Unable to set Matrix value for cell_from %i"
               " which is not between 0 and %i", cell_from, _num_x*_num_y-1);
  else if (cell_to >= _num_x*_num_y || cell_to < 0)
    log_printf(ERROR, "Unable to set Matrix value for cell_to %i"
               " which is not between 0 and %i", cell_from, _num_x*_num_y-1);
  else if (group_from >= _num_groups || group_from < 0)
    log_printf(ERROR, "Unable to set Matrix value for group_from %i"
               " which is not between 0 and %i", group_from, _num_groups-1);
  else if (group_to >= _num_groups || group_to < 0)
    log_printf(ERROR, "Unable to set Matrix value for group_to %i"
               " which is not between 0 and %i", group_to, _num_groups-1);
  
  /* Atomically set the Matrix value from the
   * temporary array using mutual exclusion locks */
  omp_set_lock(&_cell_locks[cell_to]);

  int row = cell_to*_num_groups + group_to;
  int col = cell_from*_num_groups + group_from;
  _LIL[row][col] = val;
  
  /* Release Matrix cell mutual exclusion lock */
  omp_unset_lock(&_cell_locks[cell_to]);

  /* Set global modified flag to true */
  _modified = true;  
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
  
  _A = new FP_PRECISION[NNZ];
  _IA = new int[_num_rows+1];
  _JA = new int[NNZ];
  _DIAG = new FP_PRECISION[_num_rows];
  std::fill_n(_DIAG, _num_rows, 0.0);  
  
  /* Form arrays */
  int j = 0;
  std::map<int, FP_PRECISION>::iterator iter;
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
      std::cout << std::setprecision(6) << " ( " << row << ", " << _JA[i]
                << "): " << _A[i] << std::endl;
    }
  }

  log_printf(NORMAL, "End Matrix");
}


FP_PRECISION Matrix::getValue(int cell_from, int group_from,
                              int cell_to, int group_to){
  int row = cell_to*_num_groups + group_to;
  int col = cell_from*_num_groups + group_from;
  return _LIL[row][col];
}


FP_PRECISION* Matrix::getA(){

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


FP_PRECISION* Matrix::getDiag(){

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


int Matrix::getNumGroups(){
  return _num_groups;
}


int Matrix::getNumRows(){
  return _num_rows;
}


int Matrix::getNNZ(){

  int NNZ = 0;
  std::map<int, FP_PRECISION>::iterator iter;
  for (int row=0; row < _num_rows; row++){
    for (iter = _LIL[row].begin(); iter != _LIL[row].end(); ++iter){
      if (iter->second != 0.0)
        NNZ++;
    }
  }

  return NNZ;
}


void Matrix::setNumX(int num_x) {

  if (num_x < 1)
    log_printf(ERROR, "Unable to set Matrix num x to non-positive value %i",
               num_x);

  _num_x = num_x;
}


void Matrix::setNumY(int num_y) {

  if (num_y < 1)
    log_printf(ERROR, "Unable to set Matrix num y to non-positive value %i",
               num_y);

  _num_y = num_y;
}


void Matrix::setNumGroups(int num_groups) {

  if (num_groups < 1)
    log_printf(ERROR, "Unable to set Matrix num groups to non-positive value"
               " %i", num_groups);

  _num_groups = num_groups;
}
