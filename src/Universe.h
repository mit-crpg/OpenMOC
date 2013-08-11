/**
 * @file Universe.h
 * @brief The Universe class.
 * @date January 9, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#ifdef __cplusplus
#include <map>
#include <vector>
#include "Cell.h"
#include "LocalCoords.h"
#endif

/** Error threshold for determining how close to the boundary of a lattice cell
 * a point needs to be to be considered on it */
#define ON_LATTICE_CELL_THRESH 1E-12

/** Distance a point is moved to cross over a surface into a new cell during
 * track segmentation */
#define TINY_MOVE 1E-10

class LocalCoords;
class Cell;


/**
 * @enum universeType
 * @brief The type of universe
 */
enum universeType{

    /** A simple non-repeating universe */
    SIMPLE,

    /** A collection of universes in a repeating lattice */
    LATTICE
};


/**
 * @class Universe Universe.h "openmoc/src/host/Universe.h"
 * @brief A universe represents an unbounded space in the 2D xy-plane.
 * @details A universe contains cell which are bounded subspaces in the 2D
 *          xy-plane and which together form the universe. Universes allow
 *          for complex, repeating (ie lattices) geometries to be simply
 *          represented with as few data structures as possible.
 */
class Universe {

protected:

    /** A static counter for the number of universes in a simulation */
    static short int _n;

    /** A monotonically increasing unique ID for each universe created */
    short int _uid;

    /** A user-defined id for each universe created */
    short int _id;

    /** The type of universe (ie, SIMLE or LATTICE) */
    universeType _type;

    /** A hash table of cell IDs and cell pointers */
    std::map<short int, Cell*> _cells;

    /** The coordinates of the origin for the universe */
    Point _origin;

    /** A hash table of cell IDs and their corresponding flat source region
     *  IDs. This helps for computing FSR IDs and building FSR maps for 
     *  plotting FSR-based quantities such as the scalar flux and pin powers. */
    std::map<int, int> _region_map;

    /** A boolean representing whether or not this universe contains a material
     *  with a non-zero fission cross-section and is fissionable */
    bool _fissionable;

public:
    Universe(const short int id);
    virtual ~Universe();

    void addCell(Cell* cell);
    std::map<short int, Cell*> getCells() const;
    short int getUid() const;
    short int getId() const;
    universeType getType();
    short int getNumCells() const;
    int getFSR(short int cell_id);
    Point* getOrigin();
    std::vector<short int> getMaterialIds();
    std::vector<short int> getNestedUniverseIds();
    bool isFissionable();

    void setType(universeType type);
    void setOrigin(Point* origin);
    void setFissionability(bool fissionable);
    Cell* findCell(LocalCoords* coords,
		   std::map<short int, Universe*> universes);
    int computeFSRMaps();
    void subdivideCells();
    std::string toString();
    void printString();
};


/** 
 * @class Lattice Universe.h "openmoc/src/host/Universe.h"
 * @brief Represents a repeating 2D lattice array of universes.
 */
class Lattice: public Universe {

private:

    /** The number of lattice cells along the x-axis */
    short int _num_x;

    /** The number of lattice cells along the y-axis */
    short int _num_y;

    /** The width of each lattice cell (cm) along the x-axis */
    double _width_x;

    /** The width of each lattice cell (cm) along the y-axis */
    double _width_y;

    /** A container of universes ? */
    std::vector< std::vector< std::pair<short int, Universe*> > > _universes;

    /** A container of the number of FSRs in each lattice cell */
    std::vector< std::vector< std::pair<int, int> > > _region_map;

public:

    Lattice(const short int id, const double width_x, const double width_y);
    virtual ~Lattice();

    short int getNumX() const;
    short int getNumY() const;
    Point* getOrigin();
    std::vector< std::vector< std::pair<short int, Universe*> > >
                                                         getUniverses() const;
    Universe* getUniverse(short int lattice_x, short int lattice_y) const;
    double getWidthX() const;
    double getWidthY() const;
    int getFSR(short int lat_x, short int lat_y);
    std::vector<short int> getNestedUniverseIds();

    void setLatticeCells(int num_x, int num_y, short* universes);
    void setUniversePointer(Universe* universe);

    bool withinBounds(Point* point);
    Cell* findCell(LocalCoords* coords, 
		   std::map<short int, Universe*> universes);
    Cell* findNextLatticeCell(LocalCoords* coords, double angle,
			      std::map<short int, Universe*> universes);	
    int computeFSRMaps();

    std::string toString();
    void printString();
};

#endif /* UNIVERSE_H_ */

