/*
 * Geometry.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <limits.h>
#include <limits>
#include <sys/types.h>
#include <sys/stat.h>
#include "Parser.h"
#include "LocalCoords.h"
#include "Track.h"
#include "silo.h"


class Geometry {

private:

	double _x_min, _y_min, _x_max, _y_max; 		/* the corners */
	bool _top_bc, _bottom_bc;			/* 0 for vacuum, 1 for reflective */
	bool _left_bc, _right_bc;
	int _num_FSRs;
	int* _FSRs_to_cells;
	int* _FSRs_to_materials;
	double _max_seg_length;
	double _min_seg_length;
	std::map<short int, Material*> _materials;
	std::map<short int, Surface*> _surfaces;
	std::map<short int, Cell*> _cells;
	std::map<short int, Universe*> _universes;
	std::map<short int, Lattice*> _lattices;

public:

	Geometry(Parser* parser);
	virtual ~Geometry();
	double getWidth() const;
	double getHeight() const;
	bool getBCTop() const;
	bool getBCBottom() const;
	bool getBCLeft() const;
	bool getBCRight() const;
	short int getNumRings() const;
	short int getNumSectors() const;
	double getSectorOffset() const;
	int getNumFSRs() const;
	int getNumMaterials() const;
	double getMaxSegmentLength() const;
	double getMinSegmentLength() const;
	int* getFSRtoCellMap() const;
	int* getFSRtoMaterialMap() const;

	void addMaterial(Material* material);
	Material* getMaterial(short int id);
	std::map<short int, Material*> getMaterials();
	void addSurface(Surface* surface);
	Surface* getSurface(short int id);
	void addCell(Cell *cell);
	void initializeCellFillPointers();
	Cell* getCell(short int id);
	void addUniverse(Universe* universe);
	Universe* getUniverse(short int id);
	void addLattice(Lattice* lattice);
	Lattice* getLattice(short int id);
	std::string toString();
	void printString();

	Cell* findCell(LocalCoords* coords);
	Cell* findFirstCell(LocalCoords* coords, double angle);
	Cell* findCell(int fsr_id);
	Cell* findCell(Universe* univ, int fsr_id);
	Cell* findNextCell(LocalCoords* coords, double angle);
	int findFSRId(LocalCoords* coords);
	void segmentize(Track* track);

	void computePinPowers(FP_PRECISION* FSRs_to_powers,
									FP_PRECISION* FSRs_to_pin_powers);
	FP_PRECISION computePinPowers(Universe* univ, char* output_file_prefix,
								int FSR_id, FP_PRECISION* FSRs_to_powers,
										FP_PRECISION* FSRs_to_pin_powers);

	template <class K, class V>
	bool mapContainsKey(std::map<K, V> map, K key);

};

#endif /* GEOMETRY_H_ */
