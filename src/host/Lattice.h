/*
 * Lattice.h
 *
 *  Created on: Jan 18, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#ifndef LATTICE_H_
#define LATTICE_H_

#include <vector>
#include "Universe.h"
#include "silo.h"


class Lattice: public Universe {

private:

	short int _num_x, _num_y;
	double _width_x, _width_y;
	std::vector< std::vector< std::pair<short int, Universe*> > > _universes;
	std::vector< std::vector< std::pair<int, int> > > _region_map;
	friend class Universe;

public:

	Lattice(const short int id, const short int num_x, const short int num_y,
		const double width_x, const double width_y,
		short int universes_count, short int *universes);
	virtual ~Lattice();
	void setUniversePointer(Universe* universe);
	short int getNumX() const;
	short int getNumY() const;
	Point* getOrigin();
	std::vector< std::vector< std::pair<short int, Universe*> > >
		getUniverses() const;
	Universe* getUniverse(short int lattice_x, short int lattice_y) const;
	double getWidthX() const;
	double getWidthY() const;
	int getFSR(short int lat_x, short int lat_y);

	bool withinBounds(Point* point);
	Cell* findCell(LocalCoords* coords, std::map<short int, Universe*> universes);
	Cell* findNextLatticeCell(LocalCoords* coords, 
				  double angle,
				  std::map<short int, Universe*> universes);
	std::string toString();
	
	int virtual computeFSRMaps();

};
#endif /* LATTICE_H_ */
