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
#include "Track2D.h"
#include "Geometry.h"
#include "MOCKernel.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omp.h>
#include <tuple>
#include "segmentation_type.h"
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

  /** Number of azimuthal angles in \f$ [0, 2 \pi] \f$ */
  int _num_azim;

  /** Number of polar angles in \f$ [0, \pi] \f$ */
  int _num_polar;

  /** The requested track azimuthal spacing (cm) */
  double _azim_spacing;

  /** The actual track azimuthal spacing for each azimuthal angle (cm) */
  double* _azim_spacings;

  /** An integer array of the number of Tracks in a cycle for each azim angle */
  int* _tracks_per_cycle;

  /** An array of the number of cycles for each azimuthal angle */
  int* _cycles_per_azim;

  /** An array of the cycle length of each cycle for each azimuthal angle */
  double* _cycle_length;

  /** The total number of Tracks for all azimuthal angles */
  int _num_2D_tracks;

  /** An integer array with the Track uid separating the azimuthal, polar, and
   * periodic halfspaces */
  int* _num_tracks_by_parallel_group;

  /** The number of parallel groups of tracks */
  int _num_parallel_track_groups;

  /** Boolen to indicate whether a periodic BC exists */
  bool _periodic;

  /** An integer array of the number of Tracks starting on each axis */
  int* _num_x;
  int* _num_y;
  double* _dx_eff;
  double* _dy_eff;

  /** The quadrature set */
  Quadrature* _quadrature;

  /** A 2D ragged array of 2D tracks (azim, track index) */
  Track2D** _tracks_2D;

  /** A 2D ragged array of 2D tracks (azim, cycle, train index) */
  Track2D**** _tracks_2D_cycle;

  /** An array of track pointers used in the Solver */
  Track** _tracks_2D_array;

  /** Pointer to the Geometry */
  Geometry* _geometry;

  /** Boolean for whether to use Track input file (true) or not (false) */
  bool _use_input_file;

  /** Determines the type of track segmentation to use */
  segmentationType _segment_formation;

  /** The z coord where the 2D tracks should be generated */
  double _z_coord;

  /** Filename for the *.tracks input / output file */
  std::string _tracks_filename;

  /** Max optical segment length for Tracks before splitting */
  //FIXME: make this true
  FP_PRECISION _max_optical_length;

  /** Maximum number of track segmenets in a single Track */
  int _max_num_segments;

  /** Boolean to indicate whether the segments should be dumped to file */
  bool _dump_segments;

  /** OpenMP mutual exclusion locks for atomic FSR operations */
  omp_lock_t* _FSR_locks;

  /** A buffer holding the computed FSR volumes */
  FP_PRECISION* _FSR_volumes;

  /** Booleans to indicate whether the Tracks and segments have been generated
   *  (true) or not (false) */
  bool _contains_2D_tracks;
  bool _contains_2D_segments;

  /** Private class methods */
  void initialize2DTracks();
  void initialize2DTrackReflections();
  void initialize2DTrackCycles();
  void recalibrate2DTracksToOrigin();
  virtual void segmentize();

public:

  TrackGenerator(Geometry* geometry, int num_azim, int num_polar,
                 double azim_spacing);
  virtual ~TrackGenerator();

  /* Get parameters */
  int getNumAzim();
  int getNumPolar();
  double getDesiredAzimSpacing();
  Geometry* getGeometry();
  int getNum2DTracks();
  int getNum2DSegments();
  void countSegments();
  int* getNumTracksByParallelGroupArray();
  int getNumParallelTrackGroups();
  bool getPeriodic();
  Track** get2DTracksArray();
  virtual Track** getTracksArray();
  Track2D** get2DTracks();
  double* getAzimSpacings();
  double getAzimSpacing(int azim);
  FP_PRECISION getMaxOpticalLength();
  int getMaxNumSegments();
  int getMaxNumTracksPerStack();
  int getNumThreads();
  int* getTracksPerCycle();
  int*** getTracksPerStack();
  int* getCyclesPerAzim();
  double getCycleLength(int azim);
  int getNumX(int azim);
  int getNumY(int azim);
  double getDxEff(int azim);
  double getDyEff(int azim);
  void exportFSRVolumes(double* out_volumes, int num_fsrs);
  FP_PRECISION* getFSRVolumesBuffer();
  FP_PRECISION* getFSRVolumes();
  FP_PRECISION getFSRVolume(int fsr_id);
  double getZCoord();
  Quadrature* getQuadrature();
  bool getCycleDirection(int azim, int cycle, int track_index);
  FP_PRECISION retrieveMaxOpticalLength();
  omp_lock_t* getFSRLocks();
  bool contains2DTracksArray();
  bool contains2DTracks();
  bool contains2DSegments();

  /* Set parameters */
  void setNumThreads(int num_threads);
  void setNumAzim(int num_azim);
  void setNumPolar(int num_polar);
  void setDesiredAzimSpacing(double spacing);
  void setGeometry(Geometry* geometry);
  void setZCoord(double z_coord);
  void setQuadrature(Quadrature* quadrature);
  void setMaxOpticalLength(FP_PRECISION tau);
  void setMaxNumSegments(int max_num_segments);
  void setDumpSegments(bool dump_segments);

  /* Worker functions */
  virtual void createDefaultQuadrature();
  void retrieve2DTrackCoords(double* coords, int num_tracks);
  void retrieve2DPeriodicCycleCoords(double* coords, int num_tracks);
  void retrieve2DReflectiveCycleCoords(double* coords, int num_tracks);
  void retrieve2DSegmentCoords(double* coords, int num_segments);
  void generateFSRCentroids(FP_PRECISION* FSR_volumes);
  void generateTracks();
  void splitSegments(FP_PRECISION max_optical_length);
  double leastCommonMultiple(double a, double b);
  void dumpSegmentsToFile();
  bool readSegmentsFromFile();
  void initializeTrackFileDirectory();
  void initialize2DTrackPeriodicIndices();
  void initialize2DTracksArray();
  void checkBoundaryConditions();
  void initialize2DTrackCycleIds();
};



#endif /* TRACKGENERATOR_H_ */
