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
#include "Track3D.h"
#include "Geometry.h"
#include "MOCKernel.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omp.h>
#include <tuple>
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

  /** The requested track polar spacing (cm) */
  double _polar_spacing;

  /** The actual track azimuthal spacing for each azimuthal angle (cm) */
  double* _azim_spacings;

  /** The actual track polar spacing (cm) by (azim, polar) */
  double** _polar_spacings;

  /** An integer array of the number of Tracks in a cycle for each azim angle */
  int* _tracks_per_cycle;

  /** An array of the number of cycles for each azimuthal angle */
  int* _cycles_per_azim;

  /** An array of the # of 3D tracks in each z-stack (azim, 2D track, polar) */
  int*** _tracks_per_stack;

  /** An array of the # of 3D tracks in each train
   *  (azim, cycle, polar, lz track) */
  int**** _tracks_per_train;

  /** An array of the cycle length of each cycle for each azimuthal angle */
  double* _cycle_length;

  /** The total number of Tracks for all azimuthal and polar angles */
  int _num_2D_tracks;
  int _num_3D_tracks;

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
  int** _num_z;
  int** _num_l;
  double* _dx_eff;
  double* _dy_eff;
  double** _dz_eff;
  double** _dl_eff;

  /** The quadrature set */
  Quadrature* _quadrature;

  /** A 2D ragged array of 2D tracks (azim, track index) */
  Track2D** _tracks_2D;

  /** A 2D ragged array of 2D tracks (azim, cycle, train index) */
  Track2D**** _tracks_2D_cycle;

  /** An array of 3D tracks (azim, 2D track, polar, z-stack) */
  Track3D**** _tracks_3D;

  /** An array of 3D tracks (azim, cycle, polar, lz track, train index) */
  Track3D****** _tracks_3D_cycle;

  /** An array of track pointers used in the Solver */
  Track** _tracks;
  Track** _flattened_tracks;

  /** Pointer to the Geometry */
  Geometry* _geometry;

  /** Boolean for whether to use Track input file (true) or not (false) */
  bool _use_input_file;

  /** Boolean for whether to solve 3D (true) or 2D (false) problem */
  bool _solve_3D;

  /** Boolean for whether to ray trace on-the-fly (true) or explicitly generate
   *  segments (false) */
  bool _OTF;

  /** Bool for whether to trace whole axial track stacks on-the-fly (true) or
      not (false) */
  bool _OTF_stacks;

  /** The z coord where the 2D tracks should be generated */
  double _z_coord;

  /** The mesh defining axial heights of radial planes segmented in on-the-fly
      calculations */
  std::vector<FP_PRECISION> _segmentation_heights;

  /** The global axial mesh to use in on-the-fly calculations */
  std::vector<FP_PRECISION> _global_z_mesh;

  /** Filename for the *.tracks input / output file */
  std::string _tracks_filename;

  /** The method to use for generating 3D tracks */
  int _track_generation_method;

  /** Max optical segment length for 3D tracks to be split during on-the-fly
   *  computation */
  FP_PRECISION _max_optical_length;

  /** Maximum number of track segmenets in a single 3D track for on-the-fly
    computation */
  int _max_num_segments;

  /** Maximum number of tracks a single 3D track stack for on-the-fly
    computation */
  int _max_num_tracks_per_stack;

  /** Boolean to indicate whether the segments should be dumped */
  bool _dump_segments;

  /** OpenMP mutual exclusion locks for atomic FSR operations */
  omp_lock_t* _FSR_locks;

  /** Booleans to indicate whether the Tracks and segments have been generated
   *  (true) or not (false) */
  bool _contains_2D_tracks;
  bool _contains_flattened_tracks;
  bool _contains_3D_tracks;
  bool _contains_2D_segments;
  bool _contains_3D_segments;
  bool _contains_global_z_mesh;
  bool _contains_segmentation_heights;

  /** Private class methods */
  void initialize2DTracks();
  void initialize3DTracks();
  void initialize2DTrackReflections();
  void initialize3DTrackReflections();
  void initialize2DTrackCycles();
  void recalibrate2DTracksToOrigin();
  void recalibrate3DTracksToOrigin();
  void segmentize2D();
  void segmentize3D();
  void segmentizeExtruded();
  void decomposeLZTrack(Track3D* track, double l_start, double l_end,
                        int azim, int cycle, int polar, int lz_index,
                        bool create_tracks);
  double findTrackEndPoint(Track2D* track, double phi, int azim_index);
  double convertLtoX(double l, int azim, int cycle);
  double convertLtoY(double l, int azim, int cycle);

public:

  TrackGenerator(Geometry* geometry, int num_azim, int num_polar,
                 double azim_spacing, double polar_spacing);
  virtual ~TrackGenerator();

  /* Get parameters */
  int getNumAzim();
  int getNumPolar();
  double getDesiredAzimSpacing();
  double getDesiredPolarSpacing();
  Geometry* getGeometry();
  int getNum2DTracks();
  int getNum3DTracks();
  int getNum2DSegments();
  int getNum3DSegments();
  void countSegments();
  int* getNumTracksByParallelGroupArray();
  int getNumParallelTrackGroups();
  bool getPeriodic();
  Track** getTracksArray();
  Track** getFlattenedTracksArray();
  Track2D** get2DTracks();
  Track3D**** get3DTracks();
  double* getAzimSpacings();
  double getAzimSpacing(int azim);
  double** getPolarSpacings();
  double getPolarSpacing(int azim, int polar);
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
  int getNumZ(int azim, int polar);
  int getNumL(int azim, int polar);
  double getDxEff(int azim);
  double getDyEff(int azim);
  FP_PRECISION* get2DFSRVolumes();
  FP_PRECISION get2DFSRVolume(int fsr_id);
  void export3DFSRVolumes(double* out_volumes, int num_fsrs);
  FP_PRECISION* get3DFSRVolumes();
  FP_PRECISION* get3DFSRVolumesOTF();
  FP_PRECISION get3DFSRVolume(int fsr_id);
  double getZCoord();
  Quadrature* getQuadrature();
  int getTrackGenerationMethod();
  bool getCycleDirection(int azim, int cycle, int track_index);
  FP_PRECISION retrieveMaxOpticalLength();
  omp_lock_t* getFSRLocks();

  /* Set parameters */
  void setNumThreads(int num_threads);
  void setNumAzim(int num_azim);
  void setNumPolar(int num_polar);
  void setDesiredAzimSpacing(double spacing);
  void setDesiredPolarSpacing(double spacing);
  void setGeometry(Geometry* geometry);
  void setSolve2D();
  void setSolve3D();
  void setOTF();
  void setOTFStacks();
  void setZCoord(double z_coord);
  void setQuadrature(Quadrature* quadrature);
  void setTrackGenerationMethod(int method);
  void setSegmentationHeights(std::vector<double> z_mesh);
  void setGlobalZMesh();
  void setMaxOpticalLength(FP_PRECISION tau);
  void setDumpSegments(bool dump_segments);

  /* Worker functions */
  bool contains2DTracks();
  bool containsFlattenedTracks();
  bool contains3DTracks();
  bool contains2DSegments();
  bool contains3DSegments();
  void retrieve2DTrackCoords(double* coords, int num_tracks);
  void retrieve2DPeriodicCycleCoords(double* coords, int num_tracks);
  void retrieve2DReflectiveCycleCoords(double* coords, int num_tracks);
  void retrieve3DPeriodicCycleCoords(double* coords, int num_tracks);
  void retrieve3DReflectiveCycleCoords(double* coords, int num_tracks);
  void retrieve3DTrackCoords(double* coords, int num_tracks);
  void retrieveSingle3DTrackCoords(double coords[6], int track_id);
  void retrieve2DSegmentCoords(double* coords, int num_segments);
  void retrieve3DSegmentCoords(double* coords, int num_segments);
  void generateFSRCentroids(FP_PRECISION* FSR_volumes);
  void generateTracks();
  void splitSegments(FP_PRECISION max_optical_length);
  double leastCommonMultiple(double a, double b);
  int binarySearch(FP_PRECISION* values, int size, FP_PRECISION val, int sign);
  bool isSolve2D();
  bool isSolve3D();
  bool isOTF();
  bool isOTFStacks();
  void dump2DSegmentsToFile();
  void dump3DSegmentsToFile();
  bool read2DSegmentsFromFile();
  bool read3DSegmentsFromFile();
  void initializeTrackFileDirectory();
  void traceSegmentsOTF(Track* flattened_track, Point* start,
                        double theta, MOCKernel* kernel);
  void traceStackOTF(Track* flattened_track, int polar_index, 
                     MOCKernel** kernels);
  void initialize2DTrackPeriodicIndices();
  void initialize3DTrackPeriodicIndices();
  void initializeTracksArray();
  void checkBoundaryConditions();
  void initialize2DTrackCycleIds();
  void initialize3DTrackCycleIds();
  void create3DTracksArrays();
};



#endif /* TRACKGENERATOR_H_ */
