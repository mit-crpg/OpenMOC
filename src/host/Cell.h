/**
 * Cell.h
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#ifndef CELL_H_
#define CELL_H_

#include "Surface.h"
#include "Point.h"
#include "LocalCoords.h"

class Universe;
class Surface;
class LocalCoords;

/* Represents cell type */
enum cellType {
	MATERIAL,
	FILL
};

/**
 * Represents a cell inside of a universe
 * @param _n static counter of the number of cells instantiated
 * @param _uid monotonically increasing unique id
 * @param _id cell id from geometry input file
 * @param _type Cell type (MATERIAL or FILL)
 * @param _universe id for the universe this cell is in
 * @param _surfaces map of surface ids (keys) to surface pointers (values)
 */
class Cell {

protected:

	static short int _n;
	short int _uid;
	short int _id;
	cellType _type;
	short int _universe;

	/* keys are +/- depending on side of surface */
	std::map<short int, Surface*> _surfaces;

public:

	Cell();
	Cell(short int id, cellType type, short int universe,
			short int num_surfaces, short int *surfaces);
	virtual ~Cell();
	void addSurface(short int surface_id, Surface* surface);
	void setSurfacePointer(Surface* surface);
	short int getUid() const;
	short int getId() const;
	cellType getType() const;
	short int getUniverse() const;
	short int getNumSurfaces() const;
	std::map<short int, Surface*> getSurfaces() const;
	void setUniverse(short int universe);
	bool cellContains(Point* point);
	bool cellContains(LocalCoords* coords);
	double minSurfaceDist(Point* point, double angle, 
			      Point* min_intersection);
	virtual std::string toString() =0;
	virtual int getNumFSRs() =0;

};


/**
 * Represents a cell filled with a material as a Cell subclass
 * @param _material id for the material filling this cell
 * @param _num_rings the number of rings within this cell
 * @param _num_sectors the number of sectors within this cell
 */
class CellBasic: public Cell {

private: 

	short int _material;
	short int _num_rings;
	short int _num_sectors;

public:

	CellBasic(short int id, short int universe, short int num_surfaces,
				short int *surfaces, short int material, short int num_rings,
				short int num_sectors);
	CellBasic(short int id, short int universe, short int material,
			  short int num_rings, short int num_sectors);
	CellBasic(short int id, short int universe, short int material);
	short int getMaterial() const;
	void addSurface(short int surface_id, Surface* surface);
	std::string toString();
	int getNumFSRs();
	CellBasic* clone(short int new_id, short int num_rings,
						short int num_sectors);
	short int getNumRings();
	short int getNumSectors();
	void setNumSectors(short int num);

};


/**
 * Represents a cell filled with a universe as a Cell subclass
 * @param _universe_fill id for the universe filling this cell
 */
class CellFill: public Cell {
private:
	std::pair<short int, Universe*> _universe_fill;
public:
	CellFill(short int id, short int universe, short int num_surfaces,
		 short int *surfaces, short int universe_fill);
	short int getUniverseFillId() const;
	Universe* getUniverseFill() const;
	void setUniverseFill(short int universe_Fill);
	void setUniverseFillPointer(Universe* universe_fill);
	std::string toString();
	int getNumFSRs();
};

#endif /* CELL_H_ */
