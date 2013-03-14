/*
 * TrackGenerator.h
 *
 *  Created on: Jan 23, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#ifndef TRACKGENERATOR_H_
#define TRACKGENERATOR_H_

#define _USE_MATH_DEFINES
#include "Track.h"
#include "Geometry.h"
#include "Plotter.h"

// writing on a text file
#include <iostream>
#include <fstream>

class TrackGenerator {

private:

	int _num_azim;			/* number of azimuthal angles */
	double _spacing;		/* track spacing */
	int* _num_tracks;		/* number of tracks at each angle */
	int* _num_x;			/* number of tracks starting on x-axis */
	int* _num_y;			/* number of tracks starting on y-axis */
	FP_PRECISION* _azim_weights;	/* azimuthal weights */
	Track** _tracks;
	Geometry* _geom;
	Plotter* _plotter;
	bool _use_inputfile;
	std::string _tracks_filename;

public:

	TrackGenerator(Geometry* geom, Plotter* plotter,
					const int num_azim,const double spacing);
	virtual ~TrackGenerator();
    FP_PRECISION *getAzimWeights() const;
    int getNumAzim() const;
    int *getNumTracks() const;
    double getSpacing() const;
    Track **getTracks() const;
    void generateTracks(bool use_inputfile, std::string tracks_filename);
	void computeEndPoint(Point* start, Point* end,  const double phi,
			const double width, const double height);
	void initializeBoundaryConditions();
	void segmentize();
	void printTrackingTimers();
	void dumpTracksToFile();
	void readTracksFromFile();

};

#endif /* TRACKGENERATOR_H_ */
