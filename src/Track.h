/**
 * @file Track.h
 * @brief The Track class.
 * @date January 19, 2012
 * @author William Boyd, MIT Course, 22 (wboyd@mit.edu)
 */

#ifndef TRACK_H_
#define TRACK_H_

#ifdef __cplusplus
#include <vector>
#include "Point.h"
#include "Material.h"
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

    /** The id for flat source region in which this segment resides */
    int _region_id;
};


/**
 * @class Track Track.h "openmoc/src/host/Track.h"
 * @brief A track represents a characteristic line across the geometry.
 * @details A track has particular starting and ending points on the 
 *          boundaries of the geometry and an azimuthal angle.
 */
class Track {

private:
    /** A monotonically increasing unique ID for each track created */
    int _uid;

    /** The track's start point */
    Point _start;

    /** The track's end point */
    Point _end;

    /** The azimuthal angle for the track */
    double _phi;

    /** The azimuthal angle index into the global 2D ragged array of tracks */
    int _azim_angle_index;

    /** A dynamically sized vector of segments making up this track */
    std::vector<segment> _segments;

    /** The track which reflects out of this track along its "forward"
     * direction for reflective boundary conditions. */
    Track* _track_in;

    /** The track which reflects out of this track along its "reverse"
     * direction for reflective boundary conditions. */
    Track* _track_out;

    /** The first index into the global 2D ragged array of tracks for the track
     *  that reflects out of this track along its "forward" direction for
     *  reflective boundary conditions. */
    int _track_in_i;

    /** The second index into the global 2D ragged array of tracks for the track
     *  that reflects out of this track along its "forward" direction for 
     *  reflective boundary conditions. */
    int  _track_in_j;

    /** The first index into the global 2D ragged array of tracks for the track
     *  that reflects out of this track along its "reverse" direction for
     *  reflective boundary conditions. */ 
    int _track_out_i;

    /** The second index into the global 2D ragged array of tracks for the track
     *  that reflects out of this track along its "reverse" direction for
     *  reflective boundary conditions */
    int _track_out_j;

    /** A boolean to indicate whether to give the flux to the "forward" 
     *  (false) or "reverse" (true) direction of the track reflecting out of 
     *  this one along its "forward" direction for reflective boundary 
     *  conditions.*/
    bool _refl_in;

    /** A boolean to indicate whether to give the flux to the "forward" 
     *  (false) or "reverse" (true) direction of the track reflecting out of 
     *  this one along its "forward" direction for reflective boundary 
     *  conditions. */
    bool _refl_out;

    /** A boolean to indicate whether the outgoing angular flux along this 
     *  track's "forward" direction should be zeroed out for vacuum boundary
     *  conditions. */
    bool _bc_in;

    /** A boolean to indicate whether the outgoing angular flux along this 
     *  track's "reverse" direction should be zeroed out for vacuum boundary
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
 * @brief Return the track's unique ID
 * @return the track's unique ID
 */
inline int Track::getUid() {
    return _uid;
}


/**
 * @brief Returns the incoming track.
 * @return a pointer to the incoming track
 */
inline Track* Track::getTrackIn() const {
    return _track_in;
}


/**
 * @brief Returns the outgoing track
 * @return a pointer to the outgoing track
 */
inline Track* Track::getTrackOut() const {
    return _track_out;
}


/**
 * @brief Returns a pointer to a segment with a given index.
 * @details Returns a pointer to the segment or ends program if track does 
 *          not have the requested segment.
 * @param segment index into the track's segments container
 * @return a pointer to the requested segment
 */
inline segment* Track::getSegment(int segment) {

    /* If track doesn't contain this segment, exits program */
    if (segment >= (int)_segments.size())
        log_printf(ERROR, "Attempted to retrieve segment s = %d but track only"
                   "has %d segments", segment, _segments.size());

    return &_segments[segment];
}



/**
 * @brief Returns a vector of pointers to the track's segments.
 * @return vector of segment pointers
 */
inline segment* Track::getSegments() {
    return &_segments[0];
}


/**
 * @brief Return the number of segments along this track.
 * @return the number of segments
 */
inline int Track::getNumSegments() {
    return _segments.size();
}


#endif /* TRACK_H_ */
