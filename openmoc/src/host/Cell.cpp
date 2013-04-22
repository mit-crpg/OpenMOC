#include "Cell.h"


short int Cell::_n = 0;


/**
 * @brief Default constructor used in rings/sectors subdivision of cells.
 */
Cell::Cell() { }


/**
 * @brief Constructor sets the unique and user-specifed IDs for this cell.
 * @param id the user-specified cell ID
 * @param universe the ID of the universe within which this cell resides
 */
Cell::Cell(short int id, short int universe) {
    _uid = _n;
    _id = id;
    _n++;
    _universe = universe;       
}


/**
 * @brief Destructor frees all surfaces making up cell.
 */
Cell::~Cell() {
    _surfaces.clear();
}


/**
 * @brief Return the cell's unique ID.
 * @return the cell's unique ID
 */
short int Cell::getUid() const {
    return _uid;
}


/**
 * @brief Return the cell's user-specified ID.
 * @return the cell's user-specified ID
 */
short int Cell::getId() const {
    return _id;
}


/**
 * @brief Return the cell type (FILL or MATERIAL).
 * @return the cell type
 */
cellType Cell::getType() const {
    return _cell_type;
}


/**
 * @brief Return the ID of the universe within which this cell resides.
 * @return the universe ID
 */
short int Cell::getUniverse() const {
    return _universe;
}


/**
 * @brief Return the number of surfaces in the cell.
 * @return the number of surfaces
 */
short int Cell::getNumSurfaces() const {
    return this->_surfaces.size();
}


/**
 * @brief Return the hashtable of surfaces IDs and surface pointers for all 
 *        surfaces making up the cell.
 * @return map of surface IDs and surface pointers
 */
std::map<short int,Surface*> Cell::getSurfaces() const {
    return _surfaces;
}


/**
 * @brief Set the ID for teh universe within whic this cell resides.
 * @param universe the universe's user-specified ID
 */
void Cell::setUniverse(short int universe) {
    _universe = universe;
}


/**
 * @brief Insert the a surface into this cell's container.
 * @param halfspace the surface halfspace (positive/negative for surface side)
 * @param surface a pointer to the surface
 */
void Cell::addSurface(short int halfspace, Surface* surface) {
    if (halfspace != -1 && halfspace != +1)
        log_printf(ERROR, "Unable to add surface %d to cell %d since the "
                    "halfspace %d is not -1 or 1", surface->getId(), 
                                                    _id, halfspace);
    _surfaces.insert(std::pair<short int, Surface*>(halfspace*surface->getId(),
                                                                     surface));
}


/**
 * @brief Registers a surface pointer with the Cell's surfaces map. 
 * @details This method is used by the geometry whenever a new cell is 
 *          added to it.
 * @param surface the surface pointer
 */
void Cell::setSurfacePointer(Surface* surface) {

    /* if _surfaces does not contain this surface id throw an error */
    if (_surfaces.find(surface->getId()) == _surfaces.end() &&
          _surfaces.find(-surface->getId()) == _surfaces.end())

      log_printf(WARNING, "Unable to set surface pointer for cell id = %d "
                 "for surface id = %d since cell does not contain this surface",
                 _id, surface->getId());

    try{
        /* If the cell contains the positive side of the surface */
        if (_surfaces.find(surface->getId()) != _surfaces.end())
            _surfaces[surface->getId()] = surface;

        /* If the cell contains the negative side of the surface */
        else
            _surfaces[-1*surface->getId()] = surface;

        log_printf(INFO, "Set the surface pointer for cell id = %d for "
                   "surface id = %d", _id, surface->getId());
    }
    catch (std::exception &e) {
        log_printf(ERROR, "Unable to add surface with id = %d to cell with "
		 "id = %d. Backtrace:\n%s", surface, _id, e.what());
    }
}


/**
 * @brief Determines whether a point is contained inside a cell. 
 * @details Queries each surface inside the cell to determine if the point
 *          particle is on the same side of the surface. This point is 
 *          only inside the cell if it is on the same side of every surface 
 *          in the cell.
 * @param point a pointer to a point
 */
bool Cell::cellContains(Point* point) {

    /* Loop over all surfaces inside the cell */
    std::map<short int, Surface*>::iterator iter;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
      
        /* If the surface evaluated point is not the same sign as the surface
         * or within a threshold, return false */
        if (iter->second->evaluate(point) * iter->first < -ON_SURFACE_THRESH)
            return false;
    }

    return true;
}


/**
 * @brief Determines whether a point is contained inside a cell. 
 * @details Queries each surface inside the cell to determine if the particle 
 *          is on the same side of the surface. This particle is only inside 
 *          the cell if it is on the same side of every surface in the cell.
 * @param coords a pointer to a localcoord
 */
bool Cell::cellContains(LocalCoords* coords) {
    return this->cellContains(coords->getPoint());
}


/**
 * @brief Computes the minimum distance to a surface from a point with a given
 *        trajectory at a certain angle. 
 * @details If the trajectory will not intersect any of the surfaces in the cell*           returns INFINITY.
 * @param point the point of interest
 * @param angle the angle of the trajectory (in radians from 0 to 2*PI)
 * @param min_intersection a pointer to the intersection point that is found
 */
double Cell::minSurfaceDist(Point* point, double angle,
                            Point* min_intersection) {

    double min_dist = INFINITY; 
    double d;
    Point intersection;

    std::map<short int, Surface*>::iterator iter;
    
    /* Loop over all of the cell's surfaces */
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

        /* Find the minimum distance from this surface to this point */
        d = iter->second->getMinDistance(point, angle, &intersection);

        /* If the distance to cell is less than current min distance, update */
        if (d < min_dist) {
            min_dist = d;
            min_intersection->setX(intersection.getX());
            min_intersection->setY(intersection.getY());
        }
    }

    return min_dist;
}


/**
 * Constructor sets the user-specified and unique IDs for this cellbasic.
 * @param id the user-specified cell ID
 * @param universe the ID for the universe within which this cellbasic resides
 * @param material the ID for the material which fills this cellbasic
 */
CellBasic::CellBasic(short int id, short int universe, short int material,
		     int num_rings, int num_sectors): Cell(id, universe) {
    _cell_type = MATERIAL;
    _material = material;
    _num_rings = num_rings;
    _num_sectors = num_sectors;
}



/**
 * @brief Return the ID of the material filling the cellbasic.
 * @return the material's ID
 */
short int CellBasic::getMaterial() const {
    return _material;
}


/**
 * @brief Return the number of rings in the cell.
 * @return the number of rings
 */
short int CellBasic::getNumRings() {
    return this->_num_rings;
}


/**
 * @brief Return the number of sectors in the cell.
 * @return the number of sectors
 */
short int CellBasic::getNumSectors() {
    return this->_num_sectors;
}


/**
 * @brief Return the number of flat source regions in this cellbasic. 
 * @details This method is used when the geometry recursively constructs flat
 *           source regions. By definition, cellbasic's are the lowest level
 *           entity in the geometry and thus only have one flat source region
 *           within them, so this method always returns 1.
 * @return the number of FSRs in this cell
 */
int CellBasic::getNumFSRs() {
    return 1;
}


/**
 * @brief Set the cell's number of rings.
 * @param num_rings the number of rings in this cell
 */
void CellBasic::setNumRings(short int num_rings) {
    _num_rings = num_rings;
}


/**
 * @brief Set the cell's number of sectors.
 * @param num_sectors the number of sectors in this cell
 */
void CellBasic::setNumSectors(short int num_sectors) {
    _num_sectors = num_sectors;
}


/**
 * @brief Create a duplicate of the cellbasic.
 * @param new_id the clone's user-specified cell ID
 * @param num_rings the number of rings to put in the clone
 * @param num_sectors the number of sectors to put in the clone
 * @return a pointer to the clone
 */
CellBasic* CellBasic::clone(short int new_id, short int num_rings,
                            short int num_sectors) {

    /* Construct new cell */
    CellBasic* new_cell = new CellBasic(new_id, _universe, _material);
    new_cell->setNumSectors(_num_sectors);
    new_cell->setNumRings(_num_rings);

    /* Loop over all of this cell's surfaces and add them to the clone */
    std::map<short int, Surface*>::iterator iter;
    int halfspace;
    Surface* surface;
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
        halfspace = iter->first / abs(iter->first);
        surface = iter->second;
        new_cell->addSurface(halfspace, surface);
    }
    return new_cell;
}


/**
 * @brief Convert this cellbasic's attributes to a string format.
 * @return a character array of this cellbasic's attributes
 */
std::string CellBasic::toString() {

    std::stringstream string;

    string << "Cell id = " << _id 
           << ", type = MATERIAL, material id = " << _material 
           << ", universe = " << _universe  
           << ", num_surfaces = " << getNumSurfaces() 
           << ", num of rings = " << _num_rings 
           << ", num of sectors = " << _num_sectors;
    
    /* Append each of the surface ids to the string */
    std::map<short int, Surface*>::iterator iter;
    string << ", surface ids = ";
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
        string << iter->first << ", ";

    return string.str();
}


/**
 * @brief Prints a string representation of all of the cellbasic's objects to
 *        the console.
 */
void CellBasic::printString() {
    log_printf(RESULT, toString().c_str());
}


/**
 *  @brief CellFill constructor
 *  @param id the user-specified cell ID
 *  @param universe the ID of the universe within which this cell resides
 *  @param universe_fill the ID of the universe filling this cell
 */
CellFill::CellFill(short int id, short int universe, short int universe_fill):
    Cell(id, universe) {
        _cell_type = FILL;
        _universe_fill.first = universe_fill;
}


/**
 * @brief Return the ID of the universe filling this cell.
 * @return the universe's id
 */
short int CellFill::getUniverseFillId() const {
    return _universe_fill.first;
}


/**
 * @brief Return a pointer to the universe filling this cell.
 * @return the universe pointer
 */
Universe* CellFill::getUniverseFill() const {
    return _universe_fill.second;
}


/**
 * @brief Return the number of flat source regions in this cellfill. 
 * @details This method is used when the geometry recursively constructs flat
 *          source regions.
 * @return the number of FSRs in this cellfill
 */
int CellFill::getNumFSRs() {
    Universe *univ = getUniverseFill();
    if (univ->getType() == SIMPLE)
        return univ->computeFSRMaps();
    else
        return static_cast<Lattice*>(univ)->computeFSRMaps();
}


/**
 * @brief Set the ID of the universe filling this cellfill.
 * @param universe_fill the universe's ID
 */
void CellFill::setUniverseFill(short int universe_fill) {
    _universe_fill.first = universe_fill;
}


/**
 * @brief Set a pointer to the universe filling this cellfill.
 * @param universe the universe's pointer
 */
void CellFill::setUniverseFillPointer(Universe* universe) {
    _universe_fill.second = universe;
}


/**
 * @brief Convert this cellfill's attributes to a string format.
 * @return a character array of this cell's attributes
 */
std::string CellFill::toString() {

    std::stringstream string;

    string << "Cell id = " << _id << ", type = FILL, universe_fill = " <<
        _universe_fill.first << ", universe = " << _universe <<
      ", num_surfaces = " << getNumSurfaces();

    /** Add the IDs for the surfaces in this cell */
    std::map<short int, Surface*>::iterator iter;
    string << ", surface ids = ";
    for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
        string << iter->first << ", ";

    return string.str();
}


/**
 * @brief Prints a string representation of all of the cellfill's objects to
 *        the console.
 */
void CellFill::printString() {
    log_printf(RESULT, toString().c_str());
}
