/**
 * @file TrackGenerator.h
 * @brief The TrackGenerator class.
 * @date January 23, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
FIXME
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
 *          conditions (vacuum or reflective) for each Track.
 */
class TrackGenerator3D {

private:

  /** The requested track polar spacing (cm) */
  double _polar_spacing;

  /** The actual track polar spacing (cm) by (azim, polar) */
  double** _polar_spacings;

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
  int _num_rows;
  int _num_columns;

  /** A matrix of temporary segments are created for on-the-fly methods to
    * improve efficiency */
  std::vector<segment**> _temporary_segments;
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
  void initialize3DTracks();
  void initialize3DTrackReflections();
  void recalibrate3DTracksToOrigin();
  void segmentize3D();
  void segmentizeExtruded();
  void decomposeLZTrack(Track3D* track, double l_start, double l_end,
                        int azim, int cycle, int polar, int lz_index,
                        bool create_tracks);
  double convertLtoX(double l, int azim, int cycle);
  double convertLtoY(double l, int azim, int cycle);

public:

  TrackGenerator(Geometry* geometry, int num_azim, int num_polar,
                 double azim_spacing, double polar_spacing);
  virtual ~TrackGenerator();


  /* Get parameters */
  double getDesiredPolarSpacing();
  int getNum3DTracks();
  int getNum3DSegments();
  Track** get3DTracksArray();
  Track3D**** get3DTracks();
  double** getPolarSpacings();
  double getPolarSpacing(int azim, int polar);
  double getZSpacing(int azim, int polar);
  int getNumRows();
  int getNumColumns();
  segment* getTemporarySegments(int thread_id, int row_num);
  int getNumZ(int azim, int polar);
  int getNumL(int azim, int polar);
  int getTrackGenerationMethod();
  segmentationType getSegmentFormation();
  bool contains3DTracksArray();
  bool contains3DTracks();
  bool contains3DSegments();

  /* Set parameters */
  void setDesiredPolarSpacing(double spacing);
  void setSegmentFormation(segmentationType segmentation_type);
  void setTrackGenerationMethod(int method);
  void setSegmentationHeights(std::vector<double> z_mesh);
  void setGlobalZMesh();

  /* Worker functions */
  void retrieve3DPeriodicCycleCoords(double* coords, int num_tracks);
  void retrieve3DReflectiveCycleCoords(double* coords, int num_tracks);
  void retrieve3DTrackCoords(double* coords, int num_tracks);
  void retrieveGlobalZMesh(double*& z_mesh, int& num_fsrs);
  void retrieveSingle3DTrackCoords(double coords[6], int track_id);
  void retrieve3DSegmentCoords(double* coords, int num_segments);
  void dump2DSegmentsToFile();
  void initialize3DTrackPeriodicIndices();
  void initialize3DTrackCycleIds();
  void create3DTracksArrays();
};

#endif /* TRACKGENERATOR3D_H_ */
