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
#include "Quadrature.h"
#include "Timer.h"
#include "segmentation_type.h"
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
 *          conditions (vacuum, reflective, or periodic) for each Track.
 */
class TrackGenerator {

protected:

  /** The number of shared memory OpenMP threads */
  int _num_threads;

  /** Half the number of azimuthal angles in \f$ [0, 2\pi] \f$ */
  int _num_azim_2;

  /** The requested track azimuthal spacing (cm) */
  double _azim_spacing;

  /** An integer array of the number of Tracks for each azimuthal angle */
  int* _num_tracks;

  /** An integer array of the number of Tracks starting on the x-axis for each
   *  azimuthal angle */
  int* _num_x;

  /** An integer array of the number of Tracks starting on the y-axis for each
   *  azimuthal angle */
  int* _num_y;

  /** A boolean indicating if a user-defined Quadrature was assigned */
  bool _user_quadrature;

  /** The associated Quadrature object */
  Quadrature* _quadrature;

  /** A 2D ragged array of Tracks */
  Track** _tracks;

  /** A 1D array of Track pointers arranged by UID */
  Track** _tracks_array;

  /** Pointer to the Geometry */
  Geometry* _geometry;

  /** Boolean for whether to use Track input file (true) or not (false) */
  bool _use_input_file;

  /** Filename for the *.tracks input / output file */
  std::string _tracks_filename;

  /** OpenMP mutual exclusion locks for atomic FSR operations */
  omp_lock_t* _FSR_locks;

  /** Boolean whether the Tracks have been generated (true) or not (false) */
  bool _contains_tracks;

  /** Boolean whether segments have been generated for every Track (true) or
      not (false) */
  bool _contains_segments;

  /** Boolean whether to print segments to a tracking file (true) or not
   *  (false) */
  bool _dump_segments;

  /** The z-coord where the 2D Tracks should be created */
  double _z_coord;

  /** Determines the type of track segmentation to use */
  segmentationType _segment_formation;

  /** Max optical path length for segments before splitting */
  FP_PRECISION _max_optical_length;

  /** A timer to record timing data for ray tracing */
  Timer* _timer;

  /** A buffer holding the computed FSR volumes */
  FP_PRECISION* _FSR_volumes;

  void computeEndPoint(Point* start, Point* end,  const double phi,
                       const double width_x, const double width_y);

  void initializeTrackFileDirectory();
  void initializeTracks();
  void recalibrateTracksToOrigin();
  void initializeTrackUids();
  void initializeBoundaryConditions();
  void initializeTrackCycleIndices(boundaryType bc);
  void initializeVolumes();
  void initializeFSRLocks();
  void segmentize();
  void dumpTracksToFile();
  bool readTracksFromFile();
  void clearTimerSplits();
  void calculateFSRVolumes();
  void resetStatus();

public:

  TrackGenerator(Geometry* geometry, int num_azim, double azim_spacing);
  virtual ~TrackGenerator();

  /* Get parameters */
  int getNumAzim();
  double getDesiredAzimSpacing();
  double getPhi(int azim);
  Geometry* getGeometry();
  Quadrature* getQuadrature();
  int getNumTracks();
  int getNumX(int azim);
  int getNumY(int azim);
  int getNumSegments();
  Track** getTracks();
  Track** getTracksArray();
  FP_PRECISION retrieveMaxOpticalLength();
  int getNumThreads();
  FP_PRECISION* getFSRVolumes();
  FP_PRECISION getFSRVolume(int fsr_id);
  FP_PRECISION getMaxOpticalLength();
  double getZCoord();
  omp_lock_t* getFSRLocks();
  segmentationType getSegmentFormation();

  /* Set parameters */
  void setNumAzim(int num_azim);
  void setDesiredAzimSpacing(double azim_spacing);
  void setGeometry(Geometry* geometry);
  void setQuadrature(Quadrature* quadrature);
  void setNumThreads(int num_threads);
  void setZCoord(double z_coord);
  void setMaxOpticalLength(FP_PRECISION tau);

  /* Worker functions */
  bool containsTracks();
  bool containsSegments();
  void retrieveTrackCoords(double* coords, int num_tracks);
  void retrieveSegmentCoords(double* coords, int num_segments);
  void generateTracks(bool store=true, bool neighbor_cells=false);
  void correctFSRVolume(int fsr_id, FP_PRECISION fsr_volume);
  void generateFSRCentroids();
  void splitSegments(FP_PRECISION max_optical_length);
  void initializeSegments();
  void printTimerReport();
  void resetFSRVolumes();
};

#endif /* TRACKGENERATOR_H_ */
