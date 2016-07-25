/**
 * @file TrackGenerator3D.h
 * @brief The TrackGenerator3D class.
 * @date March 13, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef TRACKGENERATOR3D_H_
#define TRACKGENERATOR3D_H_

#include "TrackGenerator.h"
#include "Track3D.h"


/**
 * @struct CycleTrackIndexes
 * @brief A cycle track represents a track in the l-z plane and this struct
 *        contains the azim, cycle, polar, lz, and train indices of a track.
 */
struct CycleTrackIndexes {

  /** The azimuthal index (in 0 to _num_azim / 4) */
  int _azim;

  /** The x index (in 0 to _num_x[_azim]) */
  int _x;

  /** The polar index (in 0 to _num_polar) */
  int _polar;

  /** The lz index (in 0 to _num_l[_azim][_polar] + _num_z[_azim][_polar]) */
  int _lz;

  /** The train index */
  int _train;

  /** Constructor initializes each attribute to -1 */
  CycleTrackIndexes() {
    _azim  = -1;
    _x = -1;
    _polar = -1;
    _lz    = -1;
    _train = -1;
  }
};


/**
 * @struct StackTrackIndexes
 * @brief A stack track represents a track in a z-stack of tracks and this
 *        struct contains the azim, xy, polar, and z indices of a track.
 */
struct StackTrackIndexes {

  /** The azimuthal index (in 0 to _num_azim / 2) */
  int _azim;

  /** The xy index (in 0 to _num_x[_azim] + _num_y[_azim]) */
  int _xy;

  /** The polar index (in 0 to _num_polar) */
  int _polar;

  /** The z index in the z-stack (in 0 to _tracks_per_stack[_azim][_xy][_polar]) */
  int _z;

  /** Constructor initializes each attribute to -1 */
  StackTrackIndexes() {
    _azim  = -1;
    _xy    = -1;
    _polar = -1;
    _z     = -1;
  }
};


/**
 * @class TrackGenerator3D TrackGenerator3D.h "src/TrackGenerator3D.h"
 * @brief The TrackGenerator3D is dedicated to generating and storing Tracks
 *        which cyclically wrap across the Geometry in three dimensions.
 * @details The TrackGenerator creates Track and initializes boundary
 *          conditions (vacuum, reflective, or periodic) for each Track in
 *          three dimensions.
 */
class TrackGenerator3D : public TrackGenerator {

private:

  /** The requested track polar spacing (cm) */
  double _polar_spacing;

  /** A 2D ragged array of 2D tracks (azim, x index, stack index) */
  Track**** _tracks_2D_cycle;

  /** An array of the # of 3D tracks in each z-stack (azim, 2D track, polar) */
  int*** _tracks_per_stack;

  // FIXME
  long*** _cum_tracks_per_stack;
  long** _cum_tracks_per_xy;
  int*** _first_lz_of_stack;

  /** The total number of Tracks for all azimuthal and polar angles */
  int _num_3D_tracks;

  /** An integer array of the number of Tracks starting on each axis */
  int** _num_z;
  int** _num_l;
  double** _dz_eff;
  double** _dl_eff;

  /** An array of 3D tracks (azim, 2D track, polar, z-stack) */
  Track3D**** _tracks_3D;

  /** The mesh defining axial heights of radial planes segmented in on-the-fly
      calculations */
  std::vector<FP_PRECISION> _segmentation_heights;

  /** The global axial mesh to use in on-the-fly calculations */
  std::vector<FP_PRECISION> _global_z_mesh;

  /** The method to use for generating 3D tracks */
  int _track_generation_method;

  /** Dimensions of temporary segments storage matrix */
  int _num_seg_matrix_rows;
  int _num_seg_matrix_columns;

  /** A matrix of temporary segments are created for on-the-fly methods to
    * improve efficiency */
  std::vector<segment*> _temporary_segments;
  bool _contains_temporary_segments;

  /** An array temporary Tracks are created for on-the-fly methods to
    * improve efficiency for every thread */
  std::vector<Track3D*> _temporary_3D_tracks;
  std::vector<Track**> _temporary_tracks_array;
  bool _contains_temporary_tracks;

  /** Maximum number of tracks a single 3D track stack for on-the-fly
   *  computation */
  int _max_num_tracks_per_stack;

  /** Booleans to indicate whether the Tracks and segments have been generated
   *  (true) or not (false) */
  bool _contains_3D_tracks;
  bool _contains_3D_segments;
  bool _contains_global_z_mesh;
  bool _contains_segmentation_heights;

  /** Private class methods */
  void initializeTracks();
  void initializeTrackCycles();
  void recalibrateTracksToOrigin();
  void segmentize();
  void setContainsSegments(bool contains_segments);
  void allocateTemporarySegments();
  void allocateTemporaryTracks();
  void resetStatus();
  void initializeDefaultQuadrature();
  std::string getTestFilename(std::string directory);
  void getCycleTrackData(CycleTrackIndexes* ctis, int num_cycles,
                         bool save_tracks);

  void segmentizeExtruded();
  void get3DTrack(CycleTrackIndexes* cti, Track3D* track,
                  bool create_arrays, bool save_tracks);
  long get3DTrackID(StackTrackIndexes* sti);
  double getLStart(CycleTrackIndexes* cti);
  int getFirstStack(CycleTrackIndexes* cti, Track3D* track_3D);

public:

  TrackGenerator3D(Geometry* geometry, int num_azim, int num_polar,
                   double azim_spacing, double polar_spacing);
  virtual ~TrackGenerator3D();

  /* Get parameters */
  double getDesiredPolarSpacing();
  int getNumTracks();
  int getNumSegments();
  int getNum3DTracks();
  long getNum3DSegments();
  Track3D**** get3DTracks();
  double getZSpacing(int azim, int polar);
  int getNumRows();
  int getNumColumns();
  segment* getTemporarySegments(int thread_id);
  Track3D* getTemporary3DTracks(int thread_id);
  Track** getTemporaryTracksArray(int thread_id);
  int getNumZ(int azim, int polar);
  int getNumL(int azim, int polar);
  int getTrackGenerationMethod();
  int*** getTracksPerStack();
  int getMaxNumTracksPerStack();
  bool containsTracks();
  bool containsSegments();
  bool containsTemporaryTracks();
  void convertCTItoSTI(CycleTrackIndexes* cti, StackTrackIndexes* sti);
  void convertSTItoCTI(StackTrackIndexes* sti, CycleTrackIndexes* cti);
  void getTrainIndex(CycleTrackIndexes* cti, StackTrackIndexes* sti);
  int getStackIndex(CycleTrackIndexes* cti);
  int getNumTracksPerLZ(CycleTrackIndexes* cti);
  void get3DTrackData(StackTrackIndexes* sti, CycleTrackIndexes* cti,
                      bool outgoing, Track3D* track);
  void getSTIByIndex(long int id, StackTrackIndexes* sti);

  void getTrackOTF(Track3D* track, StackTrackIndexes* sti);

  /* Set parameters */
  void setDesiredPolarSpacing(double spacing);
  void setSegmentFormation(segmentationType segmentation_type);
  void setTrackGenerationMethod(int method);
  void setSegmentationHeights(std::vector<FP_PRECISION> z_mesh);
  void useGlobalZMesh();

  /* Worker functions */
  void retrieveTrackCoords(double* coords, int num_tracks);
  void retrieve3DTrackCoords(double* coords, int num_tracks);
  void retrieveGlobalZMesh(FP_PRECISION*& z_mesh, int& num_fsrs);
  void retrieveSingle3DTrackCoords(double coords[6], int track_id);
  void retrieveSegmentCoords(double* coords, int num_segments);
  void retrieve3DSegmentCoords(double* coords, int num_segments);
  void create3DTracksArrays();
  void checkBoundaryConditions();
};

#endif /* TRACKGENERATOR3D_H_ */
