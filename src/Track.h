/**
 * @file Track.h
 * @brief The generic Track class.
 * @date May 16, 2015
 * @author Samuel Shaner, MIT Course, 22 (shaner@mit.edu)
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
  
  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "forward" direction should be zeroed out for vacuum boundary
   *  conditions. */
  bool _bc_in;

  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "reverse" direction should be zeroed out for vacuum boundary
   *  conditions. */
  bool  _bc_out;

  Track* _track_out;
  Track* _track_in;

  int _azim_index;
  
public:
  Track();
  virtual ~Track();

  void setUid(int uid);
  void setPhi(const double phi);

  void setBCIn(const bool bc_in);
  void setBCOut(const bool bc_out);

  void setTrackIn(Track* track_in);
  void setTrackOut(Track* track_out);

  void setAzimIndex(int azim_index);
  
  int getUid();
  Point* getEnd();
  Point* getStart();
  double getPhi() const;
  double getLength();
  
  Track* getTrackIn();
  Track* getTrackOut();

  int getAzimIndex();
  
  bool getBCIn() const;
  bool getBCOut() const;

  segment* getSegment(int s);
  segment* getSegments();
  int getNumSegments();

  /* Worker functions */
  void addSegment(segment* segment);
  void clearSegments();
  virtual std::string toString()=0;
};


/**
 * @brief Returns the incoming Track.
 * @return a pointer to the incoming Track
 */
inline Track* Track::getTrackIn() {
  return _track_in;
}


/**
 * @brief Returns the outgoing Track
 * @return a pointer to the outgoing Track
 */
inline Track* Track::getTrackOut() {
  return _track_out;
}


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
