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


  /** An array of the # of 3D tracks in each z-stack (azim, 2D track, polar) */
  int*** _tracks_per_stack;

  /** An array of the # of 3D tracks in each train
   *  (azim, cycle, polar, lz track) */
  int**** _tracks_per_train;

  /** The total number of Tracks for all azimuthal and polar angles */
  int _num_3D_tracks;

  /** An integer array of the number of Tracks starting on each axis */
  int** _num_z;
  int** _num_l;
  double** _dz_eff;
  double** _dl_eff;

  /** An array of 3D tracks (azim, 2D track, polar, z-stack) */
  Track3D**** _tracks_3D;

  /** An array of 3D tracks (azim, cycle, polar, lz track, train index) */
  Track3D****** _tracks_3D_cycle;

  /** An array of track pointers used in the Solver */
  Track** _tracks_3D_array;

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
  void initializeTracksArray();
  void initializeTrackReflections();
  void initializeTrackPeriodicIndices();
  void recalibrateTracksToOrigin();
  void segmentize();
  void setContainsSegments(bool contains_segments);
  void allocateTemporarySegments();
  void resetStatus();
  void initializeDefaultQuadrature();
  std::string getTestFilename(std::string directory);

  void segmentizeExtruded();
  void decomposeLZTrack(Track3D* track, double l_start, double l_end,
                        int azim, int cycle, int polar, int lz_index,
                        bool create_tracks);
  double convertLtoX(double l, int azim, int cycle);
  double convertLtoY(double l, int azim, int cycle);

public:

  TrackGenerator3D(Geometry* geometry, int num_azim, int num_polar,
                   double azim_spacing, double polar_spacing);
  virtual ~TrackGenerator3D();


  /* Get parameters */
  double getDesiredPolarSpacing();
  int getNumTracks();
  int getNumSegments();
  int getNum3DTracks();
  int getNum3DSegments();
  Track** getTracksArray();
  Track3D**** get3DTracks();
  double getZSpacing(int azim, int polar);
  int getNumRows();
  int getNumColumns();
  segment* getTemporarySegments(int thread_id);
  int getNumZ(int azim, int polar);
  int getNumL(int azim, int polar);
  int getTrackGenerationMethod();
  int*** getTracksPerStack();
  int getMaxNumTracksPerStack();
  bool containsTracks();
  bool containsSegments();

  /* Set parameters */
  void setDesiredPolarSpacing(double spacing);
  void setSegmentFormation(segmentationType segmentation_type);
  void setTrackGenerationMethod(int method);
  void setSegmentationHeights(std::vector<FP_PRECISION> z_mesh);
  void useGlobalZMesh();

  /* Worker functions */
  void retrieve3DPeriodicCycleCoords(double* coords, int num_tracks);
  void retrieve3DReflectiveCycleCoords(double* coords, int num_tracks);
  void retrieveTrackCoords(double* coords, int num_tracks);
  void retrieve3DTrackCoords(double* coords, int num_tracks);
  void retrieveGlobalZMesh(FP_PRECISION*& z_mesh, int& num_fsrs);
  void retrieveSingle3DTrackCoords(double coords[6], int track_id);
  void retrieveSegmentCoords(double* coords, int num_segments);
  void retrieve3DSegmentCoords(double* coords, int num_segments);
  void create3DTracksArrays();
  void initializeTrackCycleIds();
  void checkBoundaryConditions();
  void deleteTemporarySegments();
};

#endif /* TRACKGENERATOR3D_H_ */
