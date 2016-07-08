/**
 * @file Track.h
 * @brief The generic Track class.
 * @date May 16, 2015
 * @author Samuel Shaner, MIT Course, 22 (shaner@mit.edu)
 */

#ifndef TRACK_H_
#define TRACK_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "Point.h"
#include "Material.h"
#include "boundary_type.h"
#include <vector>
#include <algorithm>
#endif

/**
 * @struct segment
 * @brief A segment represents a line segment within a single flat source
 *        region along a track.
 */
struct segment {

  /** The length of the segment (cm) */
  FP_PRECISION _length;

  /** A pointer to the material in which this segment resides */
  Material* _material;

  /** The ID for flat source region in which this segment resides */
  int _region_id;

  // FIXME
  int _track_idx;

  /** The ID for the mesh surface crossed by the Track end point */
  int _cmfd_surface_fwd;

  /** The ID for the mesh surface crossed by the Track start point */
  int _cmfd_surface_bwd;

  /** Constructor initializes CMFD surfaces */
  segment() {
    _track_idx = 0;
    _cmfd_surface_fwd = -1;
    _cmfd_surface_bwd = -1;
  }
};


/**
 * @class Track Track.h "src/Track.h"
 * @brief A Track represents a characteristic line across the geometry.
 * @details A Track has particular starting and ending points on the
 *          boundaries of the geometry and an azimuthal and polar angle.
 */
class Track {

protected:

  /** A monotonically increasing unique ID for each Track created */
  int _uid;

  /** The Track's start point */
  Point _start;

  /** The Track's end point */
  Point _end;

  /** The azimuthal angle for the Track */
  double _phi;

  /** A dynamically sized vector of segments making up this Track */
  std::vector<segment> _segments;

  /** Number of segments recorded during volume calculation */
  int _num_segments;

  /** An enum to indicate whether the outgoing angular flux along this
   *  Track's "forward" direction should be zeroed out for vacuum boundary
   *  conditions or sent to a periodic or reflective track. */
  boundaryType _bc_fwd;

  /** An enum to indicate whether the outgoing angular flux along this
   *  Track's "reverse" direction should be zeroed out for vacuum boundary
   *  conditions or sent to a periodic or reflective track. */
  boundaryType _bc_bwd;

  /* Indices that are used to locate the track in the various track arrays */
  int _azim_index;
  int _xy_index;
  int _cycle_index;
  int _stack_index;

  /** Pointers to reflective and periodic Tracks in the forward and reverse
   *  directions */
  long _track_next_fwd;
  long _track_next_bwd;
  long _track_refl_fwd;
  long _track_refl_bwd;

  /** Booleans to indicate wheter the reflective Tracks in the forward and
   *  and backward direction enter into Tracks pointed in the forward
   *  direction. */
  bool _next_fwd_fwd;
  bool _next_bwd_fwd;
  bool _refl_fwd_fwd;
  bool _refl_bwd_fwd;

  /** Boolean indicating whether the track is pointing fwd (True) or bwd
   *  (False) in the cycle of tracks */
  bool _direction_in_cycle;

public:
  Track();
  virtual ~Track();

  /* Setter methods */
  virtual void setValues(const double start_x, const double start_y,
                         const double end_x, const double end_y,
                         const double phi);
  virtual void setCoords(double x0, double y0,
                         double x1, double y1);
  void setUid(int uid);
  void setPhi(const double phi);
  void setBCFwd(const boundaryType bc_fwd);
  void setBCBwd(const boundaryType bc_bwd);
  void setTrackNextFwd(long track_id);
  void setTrackNextBwd(long track_id);
  void setTrackReflFwd(long track_id);
  void setTrackReflBwd(long track_id);
  void setNextFwdFwd(bool fwd);
  void setNextBwdFwd(bool fwd);
  void setReflFwdFwd(bool fwd);
  void setReflBwdFwd(bool fwd);
  void setXYIndex(int index);
  void setAzimIndex(int index);
  void setCycleIndex(int index);
  void setStackIndex(int index);
  void setDirectionInCycle(bool fwd);

  /* Getter methods */
  int getUid();
  Point* getEnd();
  Point* getStart();
  double getPhi() const;
  double getLength();
  long getTrackNextFwd();
  long getTrackNextBwd();
  long getTrackReflFwd();
  long getTrackReflBwd();
  bool getNextFwdFwd();
  bool getNextBwdFwd();
  bool getReflFwdFwd();
  bool getReflBwdFwd();
  int getXYIndex();
  int getAzimIndex();
  int getCycleIndex();
  int getStackIndex();
  boundaryType getBCFwd() const;
  boundaryType getBCBwd() const;
  segment* getSegment(int s);
  segment* getSegments();
  int getNumSegments();
  bool getDirectionInCycle();

  /* Worker methods */
  void addSegment(segment* segment);
  void removeSegment(int index);
  void insertSegment(int index, segment* segment);
  void clearSegments();
  void setNumSegments(int num_segments);
  virtual std::string toString();
};


/**
 * @brief Return the Track's unique ID
 * @return the Track's unique ID
 */
inline int Track::getUid() {
  return _uid;
}


/**
 * @brief Returns a pointer to a segment with a given index.
 * @details Returns a pointer to the segment or ends program if Track does
 *          not have the requested segment.
 * @param segment index into the Track's segments container
 * @return a pointer to the requested segment
 */
inline segment* Track::getSegment(int segment) {

  /* If Track doesn't contain this segment, exits program */
  if (segment >= (int)_segments.size())
    log_printf(ERROR, "Attempted to retrieve segment s = %d but Track only "
               "has %d segments", segment, _segments.size());

  return &_segments[segment];
}



/**
 * @brief Returns a vector of pointers to the Track's segments.
 * @return vector of segment pointers
 */
inline segment* Track::getSegments() {
  return &_segments[0];
}


/**
 * @brief Return the number of segments along this Track.
 * @return the number of segments
 */
inline int Track::getNumSegments() {
  if (_num_segments == 0)
    return _segments.size();
  else
    return _num_segments;
}

/**
 * @brief Sets the number of segments in a track
 * @details This function sets the number of segments in a track. It's purpose
 *          is to be used for 3D tracks with on-the-fly ray tracing where
 *          segments are not explicitly created, but there is a need to know
 *          how many segments exist.
 */
inline void Track::setNumSegments(int num_segments) {
    _num_segments = num_segments;
}

#endif /* TRACK_H_ */
