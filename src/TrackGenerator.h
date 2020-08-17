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
#include "MOCKernel.h"
#include "segmentation_type.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <omp.h>
#endif

/* Modified from so/33745364/sched-getcpu-equivalent-for-os-x */
#ifndef __linux__
#include <cpuid.h>

inline int sched_getcpu() {
  int CPU;
  uint32_t CPUInfo[4];
  __cpuid_count(1, 0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
  /* CPUInfo[1] is EBX, bits 24-31 are APIC ID */
  if ( (CPUInfo[3] & (1 << 9)) == 0)
    CPU = -1;  /* no APIC on chip */
  else
    CPU = (unsigned)CPUInfo[1] >> 24;
  if (CPU < 0) CPU = 0;
  return CPU;
}
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

  /** Number of azimuthal angles in \f$ [0, 2 \pi] \f$ */
  int _num_azim;

  /** The requested track azimuthal spacing (cm) */
  double _azim_spacing;

  /** The total number of Tracks for all azimuthal angles */
  int _num_2D_tracks;

  /** An integer array of the number of Tracks starting on the x-axis for each
   *  azimuthal angle */
  int* _num_x;

  /** An integer array of the number of Tracks starting on the y-axis for each
   *  azimuthal angle */
  int* _num_y;

  /** A long integer array of the number of Tracks for each azimuthal angle */
  long* _tracks_per_azim;

  /** A 2D ragged array of 2D tracks (azim, track index) */
  Track** _tracks_2D;

  /** A 1D array of Track pointers arranged by UID */
  Track** _tracks_2D_array;

  /** Pointer to the Geometry */
  Geometry* _geometry;

  /** Boolean for whether to use Track input file (true) or not (false) */
  bool _use_input_file;

  /** Filename for the *.tracks input / output file */
  std::string _tracks_filename;

  /** OpenMP mutual exclusion locks for atomic FSR operations */
  omp_lock_t* _FSR_locks;

  /** Boolean indicating whether the Tracks have been generated (true) or not
   * (false) */
  bool _contains_2D_tracks;

  /** Boolean indicating whether 2D segments have been generated (true) or not
    * (false) */
  bool _contains_2D_segments;

  /** The quadrature set */
  Quadrature* _quadrature;

  /** The z-coord where the 2D Tracks should be created */
  double _z_coord;

  /** Boolen to indicate whether a periodic BC exists */
  bool _periodic;

  /** Determines the type of track segmentation to use */
  segmentationType _segment_formation;

  /** Max optical segment length for Tracks before splitting */
  FP_PRECISION _max_optical_length;

  /** Maximum number of track segmenets in a single Track */
  int _max_num_segments;

  /** Boolean to indicate whether the segments should be dumped to file */
  bool _dump_segments;

  /** Boolean to indicate whether the segments have been centered around their
   * centroid or not */
  bool _segments_centered;

  /** A buffer holding the computed FSR volumes */
  FP_PRECISION* _FSR_volumes;

  /** A timer to record timing data for track generation */
  Timer* _timer;

  /** Geometry boundaries for this domain */
  double _x_min;
  double _y_min;
  double _z_min;
  double _x_max;
  double _y_max;
  double _z_max;

  /** Private class methods */
  virtual void initializeTracks();
  void initializeTrackReflections();
  virtual void segmentize();
  virtual void setContainsSegments(bool contains_segments);
  virtual void allocateTemporarySegments();
  virtual void resetStatus();
  virtual void initializeDefaultQuadrature();
  virtual void writeExtrudedFSRInfo(FILE* out);
  virtual void readExtrudedFSRInfo(FILE* in);
  virtual std::string getTestFilename(std::string directory);

public:

  TrackGenerator(Geometry* geometry, int num_azim, double azim_spacing);
  virtual ~TrackGenerator();

  /* Get parameters */
  int getNumAzim();
  double getDesiredAzimSpacing();
  Geometry* getGeometry();
  virtual long getNumTracks();
  virtual long getNumSegments();
  long getNum2DTracks();
  long getNum2DSegments();
  void countSegments();
  bool getPeriodic();
  Track** get2DTracksArray();
  Track** getTracksArray();
  Track** get2DTracks();
  FP_PRECISION getMaxOpticalLength();
  int getMaxNumSegments();
  int getNumThreads();
  int getNumX(int azim);
  int getNumY(int azim);
  void exportFSRVolumes(double* out_volumes, int num_fsrs);
  void initializeVolumes();
  void initializeFSRVolumesBuffer();
  FP_PRECISION* getFSRVolumesBuffer();
  FP_PRECISION* getFSRVolumes();
  FP_PRECISION getFSRVolume(long fsr_id);
  double getZCoord();
  Quadrature* getQuadrature();
  FP_PRECISION retrieveMaxOpticalLength();
  omp_lock_t* getFSRLocks();
  segmentationType getSegmentFormation();
  virtual bool containsTracks();
  virtual bool containsSegments();
  int get2DTrackID(int a, int x);
  long* getTracksPerAzim();

  /* Set parameters */
  void setNumThreads(int num_threads);
  void setNumAzim(int num_azim);
  void setDesiredAzimSpacing(double spacing);
  void setGeometry(Geometry* geometry);
  void setZCoord(double z_coord);
  void setQuadrature(Quadrature* quadrature);
  void setMaxOpticalLength(FP_PRECISION tau);
  void setMaxNumSegments(int max_num_segments);
  void setDumpSegments(bool dump_segments);

  /* Worker functions */
  virtual void retrieveTrackCoords(double* coords, long num_tracks);
  void retrieve2DTrackCoords(double* coords, long num_tracks);
  virtual void retrieveSegmentCoords(double* coords, long num_segments);
  void retrieve2DSegmentCoords(double* coords, long num_segments);
  void generateFSRCentroids(FP_PRECISION* FSR_volumes);
  void generateTracks();
  void splitSegments(FP_PRECISION max_optical_length);
  double leastCommonMultiple(double a, double b);
  void dumpSegmentsToFile();
  bool readSegmentsFromFile();
  void initializeTrackFileDirectory();
  void initializeTracksArray();
  virtual void checkBoundaryConditions();

  /* Log functions */
  void printTimerReport(bool mpi_reduce);
  void printMemoryReport();
};


#endif /* TRACKGENERATOR_H_ */
