#include "Mesh.h"



/**
 * @brief Constructor initializes boundaries and variables that describe
          the Mesh.
 * @details The construcor initializes the many variables to zero,
 *          initializes cmfd acceleration to off, and initializes
 *          the boundary conditions to REFLECTIVE.
 * @param solve_type solve method (MOC or DIFFUSION)
 * @param cmfd_on an optional boolean to turn on CMFD
 * @param relax_factor relaxation factor
 * @param mesh_level cmfd nested universe Mesh level
 */
Mesh::Mesh(solveType solve_type, bool cmfd_on,
           double relax_factor, int mesh_level){

  if (solve_type == DIFFUSION)
    cmfd_on = true;

  /* Initialize variables */
  _cmfd_on = cmfd_on;
  _acceleration = cmfd_on;
   _num_groups = 0;
  _num_fsrs = 0;
  _num_currents = 0;
  _num_x = 0;
  _num_y = 0;
  _mesh_level = mesh_level;
  _optically_thick = false;
  _relax_factor = relax_factor;
  _solve_method = solve_type;

  /* Initialize boundaries to be reflective */
  _boundaries = new boundaryType[4];
  _boundaries[0] = REFLECTIVE;
  _boundaries[1] = REFLECTIVE;
  _boundaries[2] = REFLECTIVE;
  _boundaries[3] = REFLECTIVE;

  _volumes = NULL;
  _bounds_x = NULL;
  _bounds_y = NULL;
  _lengths_x = NULL;
  _lengths_y = NULL;

}


/**
 * @brief Destructor deletes arrays of boundaries, volumes, lengths, and bounds.
 * @details Deallocates memory for all arrays allocated by the CMFD class
 *          including boundaries, volumes, lengths, and bounds.
 */
Mesh::~Mesh(){

  if (_boundaries != NULL)
    delete [] _boundaries;

  if (_volumes != NULL)
    delete [] _volumes;

  if (_bounds_x != NULL)
    delete [] _bounds_x;

  if (_bounds_y != NULL)
    delete [] _bounds_y;

  if (_lengths_x != NULL)
    delete [] _lengths_x;

  if (_lengths_y != NULL)
    delete [] _lengths_y;

}


/**
 * @brief Initializes the Mesh by allocating memory for various arrays.
 * @details This method is called by the geometry once the width of
 *          the Mesh has been determined. This method allocates memory
 *          for the volumes, old and new flux arrays, lengths, bounds,
 *          and cell FSR vectors.
 */
void Mesh::initialize(){

  /* Allocate memory for Mesh properties */
  try{
    _volumes = new double[_num_x*_num_y];
    _lengths_x = new double[_num_x];
    _lengths_y = new double[_num_y];
    _bounds_x  = new double[_num_x+1];
    _bounds_y  = new double[_num_y+1];
  }
  catch(std::exception &e){
        log_printf(ERROR, "Could not allocate memory for the Mesh properties"
                   ". Backtrace:%s", e.what());
  }

  /* Set number of Mesh surface currents */
  _num_currents = _num_x*_num_y*8;

  /* Set initial Mesh cell flux to 1.0 and allocate memory for FSR vectors */
  for (int y = 0; y < _num_y; y++){
    for (int x = 0; x < _num_x; x++){

        /* Allocate memory for FSR vector */
        std::vector<int> *fsrs = new std::vector<int>;
        _cell_fsrs.push_back(*fsrs);
    }
  }
}


/**
 * @brief Get Mesh cell width.
 * @return Mesh cell width
 */
int Mesh::getCellsX(){
  return _num_x;
}


/**
 * @brief Get Mesh cell height.
 * @return Mesh cell height
 */
int Mesh::getCellsY(){
  return _num_y;
}


/**
 * @brief Set the Mesh cell width for a particular cell.
 * @param cell_num Mesh cell number
 * @param length_x width of Mesh cell
 */
void Mesh::setCellLengthX(int cell_num, double length_x){
  int x = cell_num % _num_x;
  _lengths_x[x] = length_x;
}


/**
 * @brief Set the Mesh cell height for a particular cell.
 * @param cell_num Mesh cell number
 * @param length_y height of Mesh cell
 */
void Mesh::setCellLengthY(int cell_num, double length_y){
  int y = cell_num / _num_x;
  _lengths_y[y] = length_y;
}


/**
 * @brief Set the number of Mesh cells in a row.
 * @param cells_x number of Mesh cells in a row
 */
void Mesh::setCellsX(int cells_x){
  _num_x = cells_x;
}


/**
 * @brief Set the number of Mesh cells in a column
 * @param cells_y number of Mesh cells in a column
 */
void Mesh::setCellsY(int cells_y){
  _num_y = cells_y;
}


/**
 * @brief given an x,y coordinate, find what Mesh cell the point is in.
 * @param x_coord coordinate
 * @param y_coord coordinate
 * @return the Mesh cell id
 */
int Mesh::findMeshCell(double x_coord, double y_coord){

  int x = 0, y = 0;

  /* Loop over cells in y direction */
  for (y = 0; y < _num_y; y++){
    if (y_coord - _bounds_y[y+1] >= -1.e-8 && y_coord - _bounds_y[y] <= 1.e-8){
      break;
    }
  }

  /* Loop over cells in y direction */
  for (x = 0; x < _num_x; x++){
    if (x_coord - _bounds_x[x] >= -1.e-8 && x_coord - _bounds_x[x+1] <= 1.e-8){
      break;
    }
  }

  int cell = (y*_num_x + x);
  return cell;
}


/**
 * @brief Set Mesh width.
 * @param length_x physical width of Mesh
 */
void Mesh::setLengthX(double length_x){
  _length_x = length_x;
}


/**
 * @brief Set Mesh height.
 * @param length_y physical height of Mesh
 */
void Mesh::setLengthY(double length_y){
  _length_y = length_y;
}


/**
 * @brief Get Mesh width.
 * @return physical width of mesh
 */
double Mesh::getLengthX(){
  return _length_x;
}


/**
 * @brief Get mesh height
 * @return physical height of mesh
 */
double Mesh::getLengthY(){
  return _length_y;
}


/**
 * @brief Set the FSR bounds for each Mesh cell.
 */
void Mesh::setFSRBounds(){

  int min;
  int max;
  std::vector<int>::iterator iter;

  /* Create arrays of FSR indices, cell bounds, and surfaces */
  try{
    _fsr_indices = new int[2 * _num_x * _num_y];
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh fsr bounds. "
               "Backtrace:%s", e.what());
  }

  /* Loop over Mesh cells */
  for (int i = 0; i < _num_y * _num_x; i++){

    /* Set fsr min and max bounds */
    min = _cell_fsrs.at(i).front();
    max = _cell_fsrs.at(i).front();

    /* Loop over FSRs and update min and max bounds */
    for (iter = _cell_fsrs.at(i).begin();
         iter != _cell_fsrs.at(i).end(); ++iter) {

      min = std::min(*iter, min);
      max = std::max(*iter, max);
    }

    /* Set FSR bounds */
    _fsr_indices[2*i] = min;
    _fsr_indices[2*i+1] = max;
  }
}


/**
 * @brief Using an FSR ID and coordinate, find which surface a coordinate is on.
 * @param fsr_id the ID of the FSR of interest
 * @param coord coordinate of segment on Mesh surface
 * @return surface UID representing surface
 */
int Mesh::findMeshSurface(int fsr_id, LocalCoords* coord){

  int surface = -1;
  double x = coord->getX();
  double y = coord->getY();
  int i;
  bool break_cells = false;

  std::vector<int>::iterator iter;

  /* Loop over mesh cells */
  for (i = 0; i < _num_x * _num_y; i++){

    break_cells = false;

    /* Find if FSR is within bounds of cell */
    if (fsr_id >= _fsr_indices[2*i] && fsr_id <= _fsr_indices[2*i+1]){

      /* Loop over FSRs in cell to see if FSR is actually in the cell
       * since the FSR bounds can overlap */
      for (iter = _cell_fsrs.at(i).begin();
           iter < _cell_fsrs.at(i).end(); iter++){

        if (fsr_id == *iter){

          /* Check if coordinate is on left surface, left top, or
           * left bottom corner */
          if (fabs(x -  _bounds_x[i % _num_x]) < 1e-6){

            /* Check if coordinate is on left surface */
            if ((y - _bounds_y[i / _num_x+1]) > 1e-6 &&
                (y - _bounds_y[i/_num_x]) < -1e-6){

              surface = i*8+0;
              break_cells = true;
              break;
            }

            /* Check if coordinate is on left top corner */
            else if (fabs(y - _bounds_y[i/_num_x]) < 1e-6){
              surface = i*8+7;
              break_cells = true;
              break;
            }

            /* Check if coordinate is on left bottom corner */
            else{
              surface = i*8+4;
              break_cells = true;
              break;
            }
          }

          /* Check if coordinate is on right surface, right top
           * corner, or right bottom corner */
          else if (fabs(x - _bounds_x[i%_num_x + 1]) < 1e-6){

            /* Check if coordinate is on right surface */
            if ((y - _bounds_y[i/_num_x+1]) > 1e-6 &&
                (y - _bounds_y[i/_num_x]) < -1e-6){

              surface = i*8+2;
              break_cells = true;
              break;
            }

            /* Check if coordinate is on right top surface */
            else if (fabs(y - _bounds_y[i/_num_x]) < 1e-6){
              surface = i*8+6;
              break_cells = true;
              break;
            }

            /* Coordinate is on right bottom surface */
            else{
              surface = i*8+5;
              break_cells = true;
              break;
            }
          }

          /* Check if coordinate is on top surface */
          else if (fabs(y - _bounds_y[i/_num_x]) < 1e-6){
            surface = i*8+3;
            break_cells = true;
            break;
          }

          /* Coordinate is on bottom surface */
          else if (fabs(y - _bounds_y[i/_num_x+1]) < 1e-6){
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


/**
 * @brief Print the Mesh surface currents to the console.
 */
void Mesh::printCurrents(){

  double current;

  /* Loop over Mesh cells */
  for (int i = 0; i < _num_x * _num_y; i++){

    /* Loop over Mesh surfaces */
    for (int s = 0; s < 8; s++){

      /* Loop over energy groups */
      for (int g = 0; g < _num_groups; g++){
        current = _currents[i*_num_groups*8 + s*_num_groups + g];
        log_printf(NORMAL, "cell: %i, surface: %i, group: %i, "
                   "current: %f", i, s, g, current);
      }
    }
  }
}


/**
 * @brief Set the Mesh boundary type for left surface.
 * @param side the Mesh surface UID
 * @param boundary the boundaryType of the surface
 */
void Mesh::setBoundary(int side, boundaryType boundary){
  _boundaries[side] = boundary;
}


/** @brief Split the currents of the Mesh cell corners to the nearby surfaces.
 * @details left bottom corner -> bottom surface and left surface
 *          of mesh cell below; right bottom corner -> bottom surface
 *          and right surface of mesh cell below; right top corner ->
 *          right surface and top surface of mesh cell to the right;
 *          left top corner -> left surface and top surface of mesh
 *          cell to the left.
 */
void Mesh::splitCorners(){

  log_printf(INFO, "splitting corners...");

  for (int x = 0; x < _num_x; x++){
    for (int y = 0; y < _num_y; y++){

      /* Split the LEFT BOTTOM CORNER */

      /* If cell is not on left or bottom geometry edge, give to bottom
       * surface and left surface of Mesh cell below */
      if (x > 0 && y < _num_y - 1){

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f",
                     y*_num_x+x,e,
                     _currents[(y*_num_x+x)*_num_groups*8 + 4*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + 1*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 4*_num_groups + e];
          _currents[((y+1)*_num_x+x)*_num_groups*8 + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 4*_num_groups + e];

        }
      }
      /* If cell is on left or bottom geometry edge, give to bottom
       * surface and left surface */
      else{

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, LEFT BOTTOM current: %f", 
		     y*_num_x+x,e, 
		     _currents[(y*_num_x+x)*_num_groups*8 + 4*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + 1*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 4*_num_groups + e];
          _currents[(y*_num_x+x)*_num_groups*8 + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 4*_num_groups + e];

        }
      }

      /* Split the RIGHT BOTTOM CORNER */

      /* If cell is not on right or bottom geometry edge, give to bottom
       * surface and right surface of mesh cell below */
      if (x < _num_x - 1 && y < _num_y - 1){

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f",
                     y*_num_x+x,e,
                     _currents[(y*_num_x+x)*_num_groups*8 + 5*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + 1*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 5*_num_groups + e];
          _currents[((y+1)*_num_x+x)*_num_groups*8 + 2*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 5*_num_groups + e];
        }
      }

      /* If cell is on right or bottom geometry edgej, give to bottom
       * surface and right surface */
      else{

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, RIGHT BOTTOM current: %f",
                     y*_num_x+x,e,
                     _currents[(y*_num_x+x)*_num_groups*8 + 5*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + 1*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 5*_num_groups + e];
          _currents[(y*_num_x+x)*_num_groups*8 + 2*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 5*_num_groups + e];
        }
      }

      /* Split the RIGHT TOP CORNER */

      /* If cell is not on right or top geometry edge, give to right
       * surface and top surface of mesh cell to the right */
      if (x < _num_x - 1 && y > 0){

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f",
                     y*_num_x+x,e,
                     _currents[(y*_num_x+x)*_num_groups*8 + 6*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + 2*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 6*_num_groups + e];
          _currents[(y*_num_x+x+1)*_num_groups*8 + 3*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 6*_num_groups + e];
        }
      }

      /* If cell is on right or top geometry edge, give to right
       * surface and top surface */
      else{

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, RIGHT TOP current: %f",
                     y*_num_x+x,e,
                     _currents[(y*_num_x+x)*_num_groups*8 + 6*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + 2*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 6*_num_groups + e];
          _currents[(y*_num_x+x)*_num_groups*8 + 3*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 6*_num_groups + e];
        }
      }

      /* Split the LEFT TOP CORNER */

      /* If cell is not on left or top geometry edge, give to left
       * surface and top surface of mesh cell to the left */
      if (x > 0 && y > 0){

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f",
                     y*_num_x+x,e,
                     _currents[(y*_num_x+x)*_num_groups*8 + 7*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 7*_num_groups + e];
          _currents[(y*_num_x+x-1)*_num_groups*8 + 3*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 7*_num_groups + e];
        }
      }

      /* If cell is on left or top geometry edge, give to top
       * surface and left surface */
      else{

        for (int e = 0; e < _num_groups; e++){

          log_printf(DEBUG, "cell: %i, group: %i, LEFT TOP current: %f",
                     y*_num_x+x,e,
                     _currents[(y*_num_x+x)*_num_groups*8 + 7*_num_groups + e]);

          _currents[(y*_num_x+x)*_num_groups*8 + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 7*_num_groups + e];
          _currents[(y*_num_x+x)*_num_groups*8 + 3*_num_groups + e] +=
            _currents[(y*_num_x+x)*_num_groups*8 + 7*_num_groups + e];
        }
      }
    }
  }
}


/**
 * @brief Get the number of Mesh cells in the Geometry.
 * @return the number of Mesh cells
 */
int Mesh::getNumCells(){
  int num_cells = _num_x * _num_y;
  return num_cells;
}


/**
 * @brief Set the number of energy groups.
 * @param num_groups number of energy groups
 */
void Mesh::setNumGroups(int num_groups){
  _num_groups = num_groups;
}


/**
 * @brief Get the number of energy groups.
 * @return the number of energy groups
 */
int Mesh::getNumGroups(){
  return _num_groups;
}


/**
 * @brief Set the number of FSRs in the Geometry.
 * @param num_fsrs the number of FSRs
 */
void Mesh::setNumFSRs(int num_fsrs){
  _num_fsrs = num_fsrs;
}


/**
 * @brief Get the number of FSRs in the Geometry.
 * @return the number of FSRs
 */
int Mesh::getNumFSRs(){
  return _num_fsrs;
}


/**
 * @brief Get the boundaryType for one side of the Mesh.
 * @param side the Mesh surface ID
 * @return the boundaryType for the surface
 */
boundaryType Mesh::getBoundary(int side){
  return _boundaries[side];
}


/**
 * @brief Get the flux for a certain Mesh cell and energy group.
 * @param flux_type type of the flux
 * @param cell_id UID of Mesh cell
 * @param group energy group
 * @return the scalar flux
 */
double Mesh::getFlux(int cell_id, int group, fluxType flux_type){
  double* fluxes = _fluxes.find(flux_type)->second;
  return fluxes[cell_id*_num_groups + group];
}

/**
 * @brief Get the Mesh Cell ID given a LocalCoords object.
 * @param coord LocalCoords object
 * @return the number of Mesh cells
 */
int Mesh::findCellId(LocalCoords* coord){

  double x_coord = coord->getX();
  double y_coord = coord->getY();
  int x,y;

  /* Loop over cells in y direction */
  for (y = 0; y < _num_y; y++){
    if (y_coord - _bounds_y[y+1] >= -1.e-8 && y_coord - _bounds_y[y] <= 1.e-8){
      break;
    }
  }

  /* Loop over cells in y direction */
  for (x = 0; x < _num_x; x++){
    if (x_coord - _bounds_x[x] >= -1.e-8 && x_coord - _bounds_x[x+1] <= 1.e-8){
      break;
    }
  }

  return (y*_num_x + x);
}


/**
 * @brief Set the pointer to the Mesh surface currents array.
 * @param surface_currents pointer to Mesh surface currents array
 */
void Mesh::setSurfaceCurrents(double* surface_currents){
  _currents = surface_currents;
}


/**
 * @brief Return whether or not CMFD is in use.
 * @return true if CMFD is in use and false otherwise
 */
bool Mesh::getCmfdOn(){
  return _cmfd_on;
}


/**
 * @brief Return whether CMFD acceleration is in use.
 * @return true if acceleration is in use and false otherwise
 */
bool Mesh::getAcceleration(){
  return _acceleration;
}


/**
 * @brief Set the whether CMFD acceleration is being used.
 * @param accel the CMFD acceleration flag (true or false)
 */
void Mesh::setAcceleration(bool accel){
  _acceleration = accel;
}


/**
 * @brief Get pointer to a std::vector of Mesh cell FSRs.
 * @return point to nested std::vector of cell FSRs
 */
std::vector<std::vector<int> >* Mesh::getCellFSRs(){
  return &_cell_fsrs;
}


/**
 * @brief Compute the physical bounds of each Mesh cell.
 */
void Mesh::setCellBounds(){

  _bounds_x[0] = -_length_x / 2.0;
  _bounds_y[0] = _length_y / 2.0;

  log_printf(DEBUG, "bounds x: %f", _bounds_x[0]);

  /* Set x bounds */
  for (int x = 1; x < _num_x+1; x++){
    _bounds_x[x] = _bounds_x[x-1] + _lengths_x[x-1];
    log_printf(DEBUG, "bounds x: %f", _bounds_x[x]);
  }

  log_printf(DEBUG, "bounds y: %f", _bounds_y[0]);

  /* Set y bounds */
  for (int y = 1; y < _num_y+1; y++){
    _bounds_y[y] = _bounds_y[y-1] - _lengths_y[y-1];
    log_printf(DEBUG, "bounds y: %f", _bounds_y[y]);
  }

  for (int x = 0; x < _num_x; x++){
    for (int y = 0; y < _num_y; y++){
      _volumes[y*_num_x+x] = _lengths_x[x] * _lengths_y[x];
    }
  }
}


/**
 * @brief Get pointer to the FSR Materials array.
 * @return pointer to array of Materials
 */
Material** Mesh::getMaterials(){
  return _materials;
}


/**
 * @brief Get pointer to FSR volumes array.
 * @return pointer to array of FSR volumes
 */
double* Mesh::getVolumes(){
  return _volumes;
}


/**
 * @brief Set the volume of a Mesh cell.
 * @param volume volume of Mesh cell
 * @param cell_num Mesh cell id
 */
void Mesh::setVolume(double volume, int cell_num){
  _volumes[cell_num] = volume;
}


/**
 * @brief Get a Mesh cell scalar flux array.
 * @param flux_type type of flux array (PRIMAL, PRIMAL_UPDATE, ADJOINT)
 * @return array of Mesh cell scalar fluxes
 */
double* Mesh::getFluxes(fluxType flux_type){
  return _fluxes.at(flux_type);
}


/**
 * @brief Get array of Mesh cell lengths in x direction.
 * @return array of Mesh cell lengths in x direction
 */
double* Mesh::getLengthsX(){
  return _lengths_x;
}


/**
 * @brief Get array of Mesh cell lengths in y direction.
 * @return array of Mesh cell lengths in y direction
 */
double* Mesh::getLengthsY(){
  return _lengths_y;
}


/**
 * @brief Get the ID of the Mesh cell next to given Mesh cell.
 * @param cell_num current Mesh cell ID
 * @param surface_id Mesh cell surface ID to look across for neighboring cell
 * @return neighboring Mesh cell ID
 */
int Mesh::getCellNext(int cell_num, int surface_id){

  int cell_next = -1;

  if (surface_id == 0){
    if (cell_num % _num_x != 0)
      cell_next = cell_num - 1;
  }

  else if (surface_id == 1){
    if (cell_num / _num_x != _num_y - 1)
      cell_next = cell_num + _num_x;
  }

  else if (surface_id == 2){
    if (cell_num % _num_x != _num_x - 1)
      cell_next = cell_num + 1;
  }

  else if (surface_id == 3){
    if (cell_num / _num_x != 0)
      cell_next = cell_num - _num_x;
  }

  return cell_next;
}


/**
 * @brief Get array of Mesh cell surface currents.
 * @return array of Mesh cell surface currents
 */
double* Mesh::getCurrents(){
  return _currents;
}


/**
 * @brief Initialize the Mesh cell Materials.
 */
void Mesh::initializeMaterialsMOC(){

  Material* material;

  try{
    _materials = new Material*[_num_x*_num_y];

    for (int y = 0; y < _num_y; y++){
      for (int x = 0; x < _num_x; x++){
        material = new Material(y*_num_x+x);
        material->setNumEnergyGroups(_num_groups);
        _materials[y*_num_x+x] = material;
      }
    }
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh cell materials. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initialize the Mesh cell Materials.
 * @param materials std::map of FSR Materials
 * @param fsrs_to_mats array of Material ID indexed by FSR ID
 */
void Mesh::initializeMaterialsDiffusion(std::map<int, Material*>* materials,
                               int* fsrs_to_mats){

  try{
    _materials = new Material*[_num_x*_num_y];
    std::vector<int>::iterator iter;

    for (int y = 0; y < _num_y; y++){
      for (int x = 0; x < _num_x; x++)
        _materials[y*_num_x+x] = materials->at
          (fsrs_to_mats[_cell_fsrs.at(y*_num_x+x).at(0)])->clone();
    }
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh cell materials. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Get the number of Mesh surface currents.
 * @return number of Mesh surface currents
 */
int Mesh::getNumCurrents(){
  return _num_currents;
}


/**
 * @brief Initializes the Mesh surface currents.
 */
void Mesh::initializeSurfaceCurrents(){

  try{
    _currents = new double[8*_num_x*_num_y*_num_groups];
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh currents. "
               "Backtrace:%s", e.what());
  }

  for (int i = 0; i < 8*_num_x*_num_y*_num_groups; i++)
    _currents[i] = 0.0;
}


/**
 * @brief Gets the Mesh nested universe level.
 * @return Mesh nested universe level
 */
int Mesh::getMeshLevel(){
  return _mesh_level;
}


/**
 * @brief Sets the Mesh nested universe level.
 * @param mesh_level Mesh nested universe level
 */
void Mesh::setMeshLevel(int mesh_level){
  _mesh_level = mesh_level;
}


/**
 * @brief Inform whether to use optically thick diffusion correction factor.
 * @param thick flag to turn on/off (true/false) correction factor.
 */
void Mesh::setOpticallyThick(bool thick){
  _optically_thick = thick;
}


/**
 * @brief Return whether optically thick diffusion correction factor is in use.
 * @return whether optically thick diffusion correction factor is in use
 */
bool Mesh::getOpticallyThick(){
  return _optically_thick;
}


/**
 * @brief Set the relaxation factor.
 * @param relax_factor the relaxation factor
 */
void Mesh::setRelaxFactor(double relax_factor){
  _relax_factor = relax_factor;
}


/**
 * @brief Get the relaxation factor.
 * @return the relaxation factor
 */
double Mesh::getRelaxFactor(){
  return _relax_factor;
}


/**
 * @brief Get the solution method (DIFFUSION or MOC).
 * @return the solution method
 */
solveType Mesh::getSolveType(){
  return _solve_method;
}


/**
 * @brief Initializes the Mesh cell PRIMAL, PRIMAL_UPDATE and ADJOINT fluxes.
 */
void Mesh::initializeFlux(){

  /* Create pointers to flux arrays */
  double* new_flux;
  double* old_flux;
  double* adj_flux;

  /* Allocate memory for fluxes and volumes */
  try{
    new_flux = new double[_num_x*_num_y*_num_groups];
    old_flux = new double[_num_x*_num_y*_num_groups];
    adj_flux = new double[_num_x*_num_y*_num_groups];
  }
  catch(std::exception &e){
    log_printf(ERROR, "Could not allocate memory for the Mesh cell fluxes, "
               "lengths, and volumes. Backtrace:%s", e.what());
  }

  /* Insert flux arrays */
  _fluxes.insert(std::pair<fluxType, double*>(PRIMAL, new_flux));
  _fluxes.insert(std::pair<fluxType, double*>(PRIMAL_UPDATE, old_flux));
  _fluxes.insert(std::pair<fluxType, double*>(ADJOINT, adj_flux));

  /* Set number of Mesh surface currents */
  _num_currents = _num_x*_num_y*8;

  /* Set initial Mesh cell flux to 1.0 and allocate memory for FSR vectors */
  for (int y = 0; y < _num_y; y++){
    for (int x = 0; x < _num_x; x++){
      for (int g = 0; g < _num_groups; g++){
        new_flux[(y*_num_x+x)*_num_groups + g] = 1.0;
        old_flux[(y*_num_x+x)*_num_groups + g] = 1.0;
        adj_flux[(y*_num_x+x)*_num_groups + g] = 1.0;
      }
    }
  }
}
