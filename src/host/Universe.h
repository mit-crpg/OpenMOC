/*
 * Universe.h
 *
 *  Created on: Jan 9, 2012
 *      Author: wbinventor
 */


#ifndef UNIVERSE_H_
#define UNIVERSE_H_


#include <map>
#include "Cell.h"
#include "LocalCoords.h"
#include "silo.h"

class LocalCoords;
class Cell;

enum universeType{
	SIMPLE,
	LATTICE
};


class Universe {

protected:

	static short int _n;		/* Counts the number of universes */
	short int _uid;		/* monotonically increasing id based on n */
	short int _id;
	universeType _type;
	std::map<short int, Cell*> _cells;
	Point _origin;
	std::map<int, int> _region_map;

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
	void setId(const short int id);
	void setType(universeType type);
	void setNumCells(const short int num_cells);
	void setOrigin(Point* origin);
	virtual Cell* findCell(LocalCoords* coords,
			       std::map<short int, Universe*> universes);
	std::string toString();
	virtual int computeFSRMaps();

};

#endif /* UNIVERSE_H_ */

