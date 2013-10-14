#include "Mesh.h"

/* Mesh constructor */
Mesh::Mesh(bool acceleration){

  _acceleration = acceleration;
  _num_groups = 0;
  _num_fsrs = 0;
  _cells_x = 0;
  _cells_y = 0;

  _boundaries = new boundaryType[4];
  _boundaries[0] = REFLECTIVE;
  _boundaries[1] = REFLECTIVE;
  _boundaries[2] = REFLECTIVE;
  _boundaries[3] = REFLECTIVE;

}


/* Mesh destructor */
Mesh::~Mesh(){
}


/* Make mesh cells and set their surface ids */
void Mesh::initialize(){

  /* allocate memory for volumes */
  _volumes = new double[_cells_x*_cells_y];
  
  /* allocate memory for fluxes */
  FP_PRECISION* new_flux;
  FP_PRECISION* old_flux;
  new_flux = new FP_PRECISION[_cells_x*_cells_y*_num_groups];
  old_flux = new FP_PRECISION[_cells_x*_cells_y*_num_groups];
  _fluxes.insert(std::pair<std::string, FP_PRECISION*>("new_flux", new_flux));
  _fluxes.insert(std::pair<std::string, FP_PRECISION*>("old_flux", old_flux));
  
  /* allocate memory for cell widths, heights, and bounds */
  _lengths_x = new double[_cells_x];
  _lengths_y = new double[_cells_y];
  _bounds_x  = new double[_cells_x+1];
  _bounds_y  = new double[_cells_y+1];
  
  /* set mesh surface cell id's */
  for (int y = 0; y < _cells_y; y++){
    for (int x = 0; x < _cells_x; x++){
      
        for (int g = 0; g < _num_groups; g++){
	    new_flux[(y*_cells_x+x)*_num_groups + g] = 1.0;
	    old_flux[(y*_cells_x+x)*_num_groups + g] = 1.0;
	}
      
	/* allocate memory for fsr vector */
	std::vector<int> *fsrs = new std::vector<int>;
	_cell_fsrs.push_back(*fsrs);      
    }
  }
  
}


/* Get mesh cell width
 * @return mesh cell width
 */
int Mesh::getCellsX(){
  return _cells_x;
}


/* Get mesh cell height
 * @return mesh cell height
 */
int Mesh::getCellsY(){
  return _cells_y;
}


/* Set the number of mesh cells in a row
 * @param mesh cell width
 */
void Mesh::setCellLengthX(int cell_num, double length_x){

  int x = cell_num % _cells_x;
  _lengths_x[x] = length_x;
}


/* Set the number of mesh cells in a column
 * @param mesh cell height
 */
void Mesh::setCellLengthY(int cell_num, double length_y){
  
  int y = cell_num / _cells_x;
  _lengths_y[y] = length_y;
}


/* Set the number of mesh cells in a row
 * @param mesh cell width
 */
void Mesh::setCellsX(int cells_x){
  _cells_x = cells_x;
}


/* Set the number of mesh cells in a column
 * @param mesh cell height
 */
void Mesh::setCellsY(int cells_y){
  _cells_y = cells_y;
}


/* given an x,y coordinate, find what mesh cell the point is in
 * @param x coordinate
 * @param y coordinate
 * @return mesh cell id
 */
int Mesh::findMeshCell(double pointX, double pointY){

  /* initialize variables */
  double left;
  double top = _length_y / 2.0;
  bool flag = false;
  int cell = 0;
  
  /* loop over mesh cells in y direction */
  for (int y = 0; y < _cells_y; y++){
    
    /* left bound and cell number */
    left = - _length_x / 2.0;
    cell = y * _cells_x;
    
    /* if point is in current row loop over cells in row */
    if (pointY <= top && pointY >= top - _lengths_x[y]){
      
      /* loop over cells in row */
      for (int x = 0; x < _cells_x; x++){
	
	/* set current cell */
	cell = y * _cells_x + x;
	
	/* check if point is in cell */
	if (pointX >= left && pointX <= left + _lengths_x[x]){
	  flag = true;
	  break;
	}
	
	/* reset left bound */
	left = left + _lengths_x[x];
      }
      
      /* break if cell found */
      if (flag == true)
	break;
    }
    
    /* reset top bound */
    top = top - _lengths_y[y];
  }
  
  return cell;
}


/* Set mesh width
 * @param mesh width
 */
void Mesh::setLengthX(double length_x){
  _length_x = length_x;
}


/* Set mesh height
 * @param mesh height
 */
void Mesh::setLengthY(double length_y){
  _length_y = length_y;
}


/* Get mesh width
 * @return mesh width
 */
double Mesh::getLengthX(){
  return _length_x;
}


/* Get mesh height
 * @return mesh height
 */
double Mesh::getLengthY(){
  return _length_y;
}


/* Set the fsr bounds and mesh cell boundary types for each mesh cell
 */
void Mesh::setFSRBounds(){
  
  /* initialize variables */
  int min;
  int max;
  int fsr;
  std::vector<int>::iterator iter;
  
  /* create arrays of fsr indices, cell bounds, and surfaces */
  _fsr_indices = new int[2 * _cells_x * _cells_y];
  
  /* loop over mesh cells */
  for (int i = 0; i < _cells_y * _cells_x; i++){
    
    /* set fsr min and max bounds */
    min = _cell_fsrs.at(i).front();
    max = _cell_fsrs.at(i).front();
    
    /* loop over fsrs and update min and max bounds */
    for (iter = _cell_fsrs.at(i).begin(); iter != _cell_fsrs.at(i).end(); ++iter) {
      fsr = *iter;
      min = std::min(fsr, min);
      max = std::max(fsr, max);
    }
    
    /* set fsr bounds */
    _fsr_indices[2*i] = min;
    _fsr_indices[2*i+1] = max;
  }
  
  setCellBounds();
}


/* Using an fsr_id and coordinate, find which surface a coordinate is on
 * @param fsr id
 * @param coordinate object
 * @return id representing surface type
 */
int Mesh::findMeshSurface(int fsr_id, LocalCoords* coord){

	/* initialize variables */
	int surface = -1;
	double x = coord->getX();
	double y = coord->getY();
	int i;
	bool break_cells = false;

	int _cell_width = _cells_x;
	int _cell_height = _cells_y;
	int* _FSR_bounds = _fsr_indices;

	std::vector<int>::iterator iter;

	/* loop over mesh cells */
	for (i = 0; i < _cells_x * _cells_y; i++){

		break_cells = false;

		if (fsr_id >= _FSR_bounds[2*i] && fsr_id <= _FSR_bounds[2*i+1]){

			for (iter = _cell_fsrs.at(i).begin(); iter < _cell_fsrs.at(i).end(); iter++){

				if (fsr_id == (int)*iter){

					/* check if coordinate is on left surface, left top, or left bottom corner */
					if (fabs(x -  _bounds_x[i % _cells_x]) < 1e-6){

						/* check if coordinate is on left surface */
						if ((y - _bounds_y[i / _cell_width + 1]) > 1e-6 && (y - _bounds_y[i/_cell_width]) < -1e-6){
							surface = i*8+0;
							break_cells = true;
							break;
						}
						/* check if coordinate is on left top corner */
						else if (fabs(y - _bounds_y[i/_cell_width]) < 1e-6){
							surface = i*8+7;
							break_cells = true;
							break;
						}
						/* check if coordinate is on left bottom corner */
						else{
							surface = i*8+4;
							break_cells = true;
							break;
						}
					}
					/* check if coordinate is on right surface, right top corner, or right bottom corner */
					else if (fabs(x - _bounds_x[i%_cell_width + 1]) < 1e-6){
						/* check if coordinate is on right surface */
						if ((y - _bounds_y[i/_cell_width+1]) > 1e-6 && (y - _bounds_y[i/_cell_width]) < -1e-6){
							surface = i*8+2;
							break_cells = true;
							break;
						}
						/* check if coordinate is on right top surface */
						else if (fabs(y - _bounds_y[i/_cell_width]) < 1e-6){
							surface = i*8+6;
							break_cells = true;
							break;
						}
						/* coordinate is on right bottom surface */
						else{
							surface = i*8+5;
							break_cells = true;
							break;
						}
					}
					/* check if coordinate is on top surface */
					else if (fabs(y - _bounds_y[i/_cell_width]) < 1e-6){
						surface = i*8+3;
						break_cells = true;
						break;
					}
					/* coordinate is on bottom surface */
					else if (fabs(y - _bounds_y[i/_cell_width + 1]) < 1e-6){
						surface = i*8+1;
						break_cells = true;
						break;
					}
				}
			}
		}

		if (break_cells)
			break;
	}


	return surface;
}


/* print surface currents */
void Mesh::printCurrents(){

  double current;
  
  /* loop over cells */
  for (int i = 0; i < _cells_x * _cells_y; i++){
    /* loop over surfaces */
    for (int s = 0; s < 8; s++){
      /* loop over groups */
      for (int g = 0; g < _num_groups; g++){
	current = _currents[i*_num_groups*8 + s*_num_groups + g];
	log_printf(NORMAL, "cell: %i, surface: %i, group: %i, current: %f", i, s, g, current);
      }
    }
  }
}


/* set the mesh boundary type for left surface
 * @param mesh boundary type
 */
void Mesh::setBoundary(int side, boundaryType boundary){
  _boundaries[side] = boundary;
}


/* split the currents of the mesh cell corners to the nearby surfaces
 * left bottom corner -> bottom surface and left surface of mesh cell below
 * right bottom corner -> bottom surface and right surface of mesh cell below
 * right top corner -> right surface and top surface of mesh cell to the right
 * left top corner -> left surface and top surface of mesh cell to the left
 */
void Mesh::splitCorners(){

  log_printf(INFO, "splitting corners...");
  
  for (int x = 0; x < _cells_x; x++){
    for (int y = 0; y < _cells_y; y++){
      
      /* split the LEFT BOTTOM CORNER */
      
      /* if cell is not on left or bottom geometry edge
       * give to bottom surface and left surface of mesh cell below */
      if (x > 0 && y < _cells_y - 1){
	
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + 1*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 4*_num_groups + e];
	  _currents[((y+1)*_cells_x+x)*_num_groups*8 + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 4*_num_groups + e];
	}
      }
      /* if cell is on left or bottom geometry edge
       * give to bottom surface and left surface */
      else{
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + 1*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 4*_num_groups + e];
	  _currents[(y*_cells_x+x)*_num_groups*8 + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 4*_num_groups + e];
	}
      }
      
      /* split the RIGHT BOTTOM CORNER */
      
      /* if cell is not on right or bottom geometry edge
       * give to bottom surface and right surface of mesh cell below */
      if (x < _cells_x - 1 && y < _cells_y - 1){
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + 1*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 5*_num_groups + e];
	  _currents[((y+1)*_cells_x+x)*_num_groups*8 + 2*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 5*_num_groups + e];
	}
      }
      /* if cell is on right or bottom geometry edge
       * give to bottom surface and right surface */
      else{
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + 1*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 5*_num_groups + e];
	  _currents[(y*_cells_x+x)*_num_groups*8 + 2*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 5*_num_groups + e];
	}
      }
      
      /* split the RIGHT TOP CORNER */
      
      /* if cell is not on right or top geometry edge
       * give to right surface and top surface of mesh cell to the right */
      if (x < _cells_x - 1 && y > 0){
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + 2*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 6*_num_groups + e];
	  _currents[(y*_cells_x+x+1)*_num_groups*8 + 3*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 6*_num_groups + e];
	}
      }
      /* if cell is on right or top geometry edge
       * give to right surface and top surface */
      else{
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + 2*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 6*_num_groups + e];
	  _currents[(y*_cells_x+x)*_num_groups*8 + 3*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 6*_num_groups + e];
	}
      }
      
      /* split the LEFT TOP CORNER */
      
      /* if cell is not on left or top geometry edge
       * give to left surface and top surface of mesh cell to the left */
      if (x > 0 && y > 0){
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 7*_num_groups + e];
	  _currents[(y*_cells_x+x-1)*_num_groups*8 + 3*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 7*_num_groups + e];
	}
      }
      /* if cell is on left or top geometry edge
       * give to top surface and left surface */
      else{
	for (int e = 0; e < _num_groups; e++){
	  _currents[(y*_cells_x+x)*_num_groups*8 + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 7*_num_groups + e];
	  _currents[(y*_cells_x+x)*_num_groups*8 + 3*_num_groups + e] += _currents[(y*_cells_x+x)*_num_groups*8 + 7*_num_groups + e];
	}
      }
    }
  }
}


int Mesh::getNumCells(){
	return _cells_x * _cells_y;
}

void Mesh::setNumGroups(int num_groups){
	_num_groups = num_groups;
}

int Mesh::getNumGroups(){
	return _num_groups;
}

void Mesh::setNumFSRs(int num_fsrs){
	_num_fsrs = num_fsrs;
}

int Mesh::getNumFSRs(){
	return _num_fsrs;
}

boundaryType Mesh::getBoundary(int side){
  return _boundaries[side];
}

FP_PRECISION Mesh::getFlux(std::string flux_name, int cell_id, int group){
  FP_PRECISION* fluxes = _fluxes.find(flux_name)->second;
  return fluxes[cell_id*_num_groups + group];
}

int Mesh::findCellId(LocalCoords* coord){

  /* initialize variables */
  double x = coord->getX();
  double y = coord->getY();
  
  /* loop over cells in y direction */
  for (int y = 0; y < _cells_y; y++){
    if (y - _bounds_y[y+1] >= -1.e-8 && y - _bounds_y[y] <= 1.e-8){
      break;
    }
  }
  
  /* loop over cells in y direction */
  for (int x = 0; x < _cells_x; x++){
    if (x - _bounds_x[x] >= -1.e-8 && x - _bounds_x[x+1] <= 1.e-8){
      break;
    }
  }
  
  return (y*_cells_x + x);
}

void Mesh::setSurfaceCurrents(FP_PRECISION* surface_currents){
  _currents = surface_currents;
}


bool Mesh::getAcceleration(){
  return _acceleration;
}


std::vector<std::vector<int> >* Mesh::getCellFSRs(){
  return &_cell_fsrs;
}

void Mesh::setCellBounds(){

  _bounds_x[0] = -_length_x / 2.0;
  _bounds_y[0] =  _length_y / 2.0;
  

  log_printf(DEBUG, "bounds x: %f", _bounds_x[0]);

  /* set x bounds */
  for (int x = 1; x < _cells_x+1; x++){
    _bounds_x[x] = _bounds_x[x-1] + _lengths_x[x-1];
    log_printf(DEBUG, "bounds x: %f", _bounds_x[x]);
  }  

  log_printf(DEBUG, "bounds y: %f", _bounds_y[0]);

  /* set y bounds */
  for (int y = 1; y < _cells_y+1; y++){
    _bounds_y[y] = _bounds_y[y-1] - _lengths_y[y-1];
    log_printf(DEBUG, "bounds y: %f", _bounds_y[y]);
  }

}


Material** Mesh::getMaterials(){
  return _materials;
}

double* Mesh::getVolumes(){
  return _volumes;
}

void Mesh::setVolume(double volume, int cell_num){
  _volumes[cell_num] = volume;
}

FP_PRECISION* Mesh::getFluxes(std::string flux_name){
  return _fluxes.at(flux_name);
}

double* Mesh::getLengthsX(){
  return _lengths_x;
}

double* Mesh::getLengthsY(){
  return _lengths_y;
}


int Mesh::getCellNext(int cell_num, int surface_id){

  int cell_next = -1;
  
  if (surface_id == 0){
    if (cell_num % _cells_x != 0)
      cell_next = cell_num - 1;
  }
  else if (surface_id == 1){
    if (cell_num / _cells_x != _cells_y - 1)
      cell_next = cell_num + _cells_x;
  }
  else if (surface_id == 2){
    if (cell_num % _cells_x != _cells_x - 1)
      cell_next = cell_num + 1;
  }
  else if (surface_id == 3){
    if (cell_num / _cells_x != 0)
      cell_next = cell_num - _cells_x;
  }
  
  return cell_next;
}

FP_PRECISION* Mesh::getCurrents(){
  return _currents;
}


void Mesh::initializeMaterials(std::map<int, Material*>* materials, int* fsrs_to_mats){

  log_printf(INFO, "initializing mesh materials...");

  _materials = new Material*[_cells_x*_cells_y];
  std::vector<int>::iterator iter;
  
  for (int y = 0; y < _cells_y; y++){
    for (int x = 0; x < _cells_x; x++){
      
      _materials[y*_cells_x+x] = materials->at(fsrs_to_mats[_cell_fsrs.at(y*_cells_x+x).at(0)])->clone();

    }
  }

  log_printf(INFO, "Done initializing mesh materials");
}
