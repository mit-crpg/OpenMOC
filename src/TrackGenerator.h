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
#ifdef SWIG
#include "Python.h"
#endif
#include "Track.h"
#include "Geometry.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omp.h>
#endif


/**
 * @class TrackGenerator TrackGenerator.h "src/TrackGenerator.h"
 * @brief The TrackGenerator is dedicated to generating and storing Tracks
 *        which cyclically wrap across the Geometry.
 * @details The TrackGenerator creates Track and initializes boundary
 *          conditions (vacuum or reflective) for each Track.
 */
class TrackGenerator {

private:

  /** The number of shared memory OpenMP threads */
  int _num_threads;

  /** Number of azimuthal angles in \f$ [0, \pi] \f$ */
  int _num_azim;

  /** The track spacing (cm) */
  double _spacing;

  /** An integer array of the number of Tracks for each azimuthal angle */
  int* _num_tracks;

  /** An integer array with the Track uid separating the azimuthal and periodic
   * halfspaces */
  int* _num_tracks_by_halfspace;

  /** An integer array of the number of Tracks starting on the x-axis for each
   *  azimuthal angle */
  int* _num_x;

  /** An integer array of the number of Tracks starting on the y-axis for each
   *  azimuthal angle */
  int* _num_y;

  /** An array of the azimuthal angle quadrature weights */
  FP_PRECISION* _azim_weights;

  /** A 2D ragged array of Tracks */
  Track** _tracks;

  /** Pointer to the Geometry */
  Geometry* _geometry;

  /** Boolean for whether to use Track input file (true) or not (false) */
  bool _use_input_file;

  /** Filename for the *.tracks input / output file */
  std::string _tracks_filename;

  /** Boolean whether the Tracks have been generated (true) or not (false) */
  bool _contains_tracks;

  void computeEndPoint(Point* start, Point* end,  const double phi,
                       const double width, const double height);

  void initializeTrackFileDirectory();
  void initializeTracks();
  void recalibrateTracksToOrigin();
  void initializeTrackUIDs();
  void initializeBoundaryConditions();
  void segmentize();
  void dumpTracksToFile();
  bool readTracksFromFile();

public:

  TrackGenerator(Geometry* geometry, int num_azim, double spacing);
  virtual ~TrackGenerator();

  /* Get parameters */
  int getNumAzim();
  double getTrackSpacing();
  Geometry* getGeometry();
  int getNumTracks();
  int getNumX(int azim);
  int getNumY(int azim);
  int* getNumTracksArray();
  int* getNumTracksByHalfspaceArray();
  int getNumSegments();
  Track** getTracks();
  FP_PRECISION* getAzimWeights();
  int getNumThreads();
  FP_PRECISION* getFSRVolumes();
  FP_PRECISION getFSRVolume(int fsr_id);
  FP_PRECISION getMaxOpticalLength();

  /* Set parameters */
  void setNumAzim(int num_azim);
  void setTrackSpacing(double spacing);
  void setGeometry(Geometry* geometry);
  void setNumThreads(int num_threads);

  /* Worker functions */
  bool containsTracks();
  void retrieveTrackCoords(double* coords, int num_tracks);
  void retrieveSegmentCoords(double* coords, int num_segments);
  void generateTracks();
  void correctFSRVolume(int fsr_id, FP_PRECISION fsr_volume);
  void generateFSRCentroids();
  void splitSegments(FP_PRECISION max_optical_length);
};

#endif /* TRACKGENERATOR_H_ */
