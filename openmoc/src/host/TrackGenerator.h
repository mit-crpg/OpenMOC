/**
 * @file TrackGenerator.h
 * @brief The TrackGenerator class.
 * @date January 23, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef TRACKGENERATOR_H_
#define TRACKGENERATOR_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omp.h>
#include "Track.h"
#include "Geometry.h"
#endif


/**
 * @class TrackGenerator TrackGenerator.h "openmoc/src/host/TrackGenerator.h"
 * @brief The track generator is dedicated to generating tracks which cyclically
 *        wrap across the geometry.
 * @details The track generator creates track and initializes boundary 
 *          conditions (vacuum or reflective) for each track.
 */
class TrackGenerator {

private:
    /** Number of azimuthal angles in \f$ [0, \pi] \f$ */
    int _num_azim;
    /** The user-specified track spacing (cm) */
    double _spacing;
    /** An array of the number of tracks for each azimuthal angle */
    int* _num_tracks;
    /** The total number of tracks */
    int _tot_num_tracks;
    /** An array of the number of track segments per track */
    int* _num_segments;
    /** The total number of track segments */
    int _tot_num_segments;
    /** An array of the number of tracks starting on the x-axis for each
     *  azimuthal angle */
    int* _num_x;
    /** An array of the number of tracks starting on the y-axis for each
     *  azimuthal angle */
    int* _num_y;
    /** An array of the weights for each azimuthal angle */
    double* _azim_weights;
    /** A 2D ragged array of tracks */
    Track** _tracks;
    /** Pointer to the geometry */
    Geometry* _geometry;
    /** Boolean or whether to use track input file (true) or not (false) */
    bool _use_input_file;
    /** Filename for the *.tracks input / output file */
    std::string _tracks_filename;
    /** Boolean whether the tracks have been generated (true) or not (false) */
    bool _contains_tracks;

    void computeEndPoint(Point* start, Point* end,  const double phi,
                        const double width, const double height);

    void initializeBoundaryConditions();
    void segmentize();
    void dumpTracksToFile();
    bool readTracksFromFile();

public:
    TrackGenerator();
    TrackGenerator(Geometry* geometry, int num_azim, double spacing);
    virtual ~TrackGenerator();

    int getNumAzim();
    double getTrackSpacing();
    Geometry* getGeometry();
    int getNumTracks();
    int* getNumTracksArray();
    int getNumSegments();
    int* getNumSegmentsArray();
    Track** getTracks();
    double* getAzimWeights();
    bool containsTracks();
    void retrieveTrackCoords(double* coords, int num_tracks);
    void retrieveSegmentCoords(double* coords, int num_segments);

    void setNumAzim(int num_azim);
    void setTrackSpacing(double spacing);
    void setGeometry(Geometry* geometry);

    void generateTracks();
};

#endif /* TRACKGENERATOR_H_ */
