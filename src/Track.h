/**
 * @file Track.h
 * @brief The Track class.
 * @date January 19, 2012
 * @author William Boyd, MIT Course, 22 (wboyd@mit.edu)
 */

#ifndef TRACK_H_
#define TRACK_H_

#ifdef __cplusplus
#include "Python.h"
#include "Point.h"
#include "Material.h"
#include <vector>
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

  /** The ID for the mesh surface crossed by the Track end point */
  int _cmfd_surface_fwd;

  /** The ID for the mesh surface crossed by the Track start point */
  int _cmfd_surface_bwd;
};


/**
 * @class Track Track.h "src/Track.h"
 * @brief A Track represents a characteristic line across the geometry.
 * @details A Track has particular starting and ending points on the
 *          boundaries of the geometry and an azimuthal angle.
 */
class Track {

private:
  /** A monotonically increasing unique ID for each Track created */
  int _uid;

  /** The Track's start point */
  Point _start;

  /** The Track's end point */
  Point _end;

  /** The azimuthal angle for the Track */
  double _phi;

  /** The azimuthal angle index into the global 2D ragged array of Tracks */
  int _azim_angle_index;

  /** A dynamically sized vector of segments making up this Track */
  std::vector<segment> _segments;

  /** The Track which reflects out of this Track along its "forward"
   * direction for reflective boundary conditions. */
  Track* _track_in;

  /** The Track which reflects out of this Track along its "reverse"
   * direction for reflective boundary conditions. */
  Track* _track_out;

  /** The first index into the global 2D ragged array of Tracks for the Track
   *  that reflects out of this Track along its "forward" direction for
   *  reflective boundary conditions. */
  int _track_in_i;

  /** The second index into the global 2D ragged array of Tracks for the Track
   *  that reflects out of this Track along its "forward" direction for
   *  reflective boundary conditions. */
  int  _track_in_j;

  /** The first index into the global 2D ragged array of Tracks for the Track
   *  that reflects out of this Track along its "reverse" direction for
   *  reflective boundary conditions. */
  int _track_out_i;

  /** The second index into the global 2D ragged array of Tracks for the Track
   *  that reflects out of this Track along its "reverse" direction for
   *  reflective boundary conditions */
  int _track_out_j;

  /** A boolean to indicate whether to give the flux to the "forward" (false)
   *  or "reverse" (true) direction of the Track reflecting out of this one
   *  along its "forward" direction for reflective boundary conditions.*/
  bool _refl_in;

  /** A boolean to indicate whether to give the flux to the "forward" (false)
   *  or "reverse" (true) direction of the Track reflecting out of this one
   *  along its "forward" direction for reflective boundary conditions. */
  bool _refl_out;

  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "forward" direction should be zeroed out for vacuum boundary
   *  conditions. */
  bool _bc_in;

  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "reverse" direction should be zeroed out for vacuum boundary
   *  conditions. */
  bool  _bc_out;

public:
  Track();
  virtual ~Track();
  void setValues(const double start_x, const double start_y,
                 const double end_x, const double end_y, const double phi);
  void setUid(int uid);
  void setPhi(const double phi);
  void setAzimAngleIndex(const int index);
  void setReflIn(const bool refl_in);
  void setReflOut(const bool refl_out);
  void setBCIn(const bool bc_in);
  void setBCOut(const bool bc_out);
  void setTrackIn(Track *track_in);
  void setTrackOut(Track *track_out);
  void setTrackInI(int i);
  void setTrackInJ(int j);
  void setTrackOutI(int i);
  void setTrackOutJ(int j);

  int getUid();
  Point* getEnd();
  Point* getStart();
  double getPhi() const;
  int getAzimAngleIndex() const;
  segment* getSegment(int s);
  segment* getSegments();
  int getNumSegments();
  Track *getTrackIn() const;
  Track *getTrackOut() const;
  int getTrackInI() const;
  int getTrackInJ() const;
  int getTrackOutI() const;
  int getTrackOutJ() const;
  bool isReflIn() const;
  bool isReflOut() const;
  bool getBCIn() const;
  bool getBCOut() const;

  bool contains(Point* point);
  void addSegment(segment* segment);
  void clearSegments();
  std::string toString();
};


/**
 * @brief Return the Track's unique ID
 * @return the Track's unique ID
 */
inline int Track::getUid() {
  return _uid;
}


/**
 * @brief Returns the incoming Track.
 * @return a pointer to the incoming Track
 */
inline Track* Track::getTrackIn() const {
  return _track_in;
}


/**
 * @brief Returns the outgoing Track
 * @return a pointer to the outgoing Track
 */
inline Track* Track::getTrackOut() const {
  return _track_out;
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
    log_printf(ERROR, "Attempted to retrieve segment s = %d but Track only"
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
  return _segments.size();
}


#endif /* TRACK_H_ */
