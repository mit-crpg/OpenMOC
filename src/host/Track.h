/*
 * Track.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#ifndef TRACK_H_
#define TRACK_H_

#include <vector>
#include "Point.h"
#include "Material.h"


/* Represent a segment along a given track */
struct segment {
	FP_PRECISION _length;
	Material* _material;
	int _region_id;
#if STORE_PREFACTORS
	FP_PRECISION _prefactors[NUM_ENERGY_GROUPS][NUM_POLAR_ANGLES];
#endif
};


class Track {

private:

	Point _start;
	Point _end;
	double _phi;
	int _azim_angle_index;
	FP_PRECISION _azim_weight;
	FP_PRECISION _polar_weights[NUM_POLAR_ANGLES];
	FP_PRECISION _polar_fluxes[2 * GRP_TIMES_ANG];
	std::vector<segment*> _segments;
	Track *_track_in, *_track_out;
	int _track_in_i, _track_in_j;
	int _track_out_i, _track_out_j;
	bool _refl_in, _refl_out;
	bool _bc_in, _bc_out;

public:

	Track();
	virtual ~Track();
	void setValues(const double start_x, const double start_y,
			const double end_x, const double end_y, const double phi);
    void setAzimuthalWeight(const FP_PRECISION azim_weight);
    void setPolarWeight(const int angle, FP_PRECISION polar_weight);
    void setPolarFluxes(bool direction, int start_index, FP_PRECISION* polar_fluxes);
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

    Point* getEnd();
    Point* getStart();
    double getPhi() const;
    int getAzimAngleIndex() const;
    FP_PRECISION getAzimuthalWeight() const;
    FP_PRECISION* getPolarWeights();
    FP_PRECISION* getPolarFluxes();
	segment* getSegment(int s);
	std::vector<segment*> getSegments();
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

    void normalizeFluxes(FP_PRECISION factor);
    bool contains(Point* point);
    void addSegment(segment* segment);
    void clearSegments();
    std::string toString();

    void cloneOnDevice();

};

/**
 * Set this track's polar fluxes for a particular direction (0 or 1)
 * @param direction incoming/outgoing (0/1) flux for forward/reverse directions
 * @param polar_fluxes pointer to an array of fluxes
 */
inline void Track::setPolarFluxes(bool direction, int start_index,
												FP_PRECISION* polar_fluxes) {

	int start = direction * GRP_TIMES_ANG;

	if (direction != true && direction != false)
		log_printf(ERROR, "Tried to set this track's polar flux in a direction"
				"which does not exist: direction = %b", direction);

	if (!direction) {
		for (int i = 0; i < GRP_TIMES_ANG; i++)
			_polar_fluxes[start + i] = polar_fluxes[i+start_index] * _bc_in;
	}
	else {
		for (int i = 0; i < GRP_TIMES_ANG; i++)
			_polar_fluxes[start + i] = polar_fluxes[i+start_index] * _bc_out;
	}

	return;
}


/**
 * Return a pointer to this track's polar flux array
 * @return a pointer to the polar flux array
 */
inline FP_PRECISION* Track::getPolarFluxes() {
	return _polar_fluxes;
}


/**
 * Returns the incoming track
 * @return a pointer to the incoming track
 */
inline Track *Track::getTrackIn() const {
    return _track_in;
}


/**
 * Return the track's azimuthal weight
 * @param the track's azimuthal weight
 */
inline FP_PRECISION Track::getAzimuthalWeight() const {
    return _azim_weight;
}


/**
 * Returns the outgoing track
 * @return a pointer to the outgoing track
 */
inline Track *Track::getTrackOut() const {
    return _track_out;
}


/**
 * Returns a pointer to a segment with a given index or ends program if
 * track does not have the segment
 * @param segment index into the track's segments container
 * @return a pointer to the requested segment
 */
inline segment* Track::getSegment(int segment) {

	/* Checks to see if segments container contains this segment index */
	if (segment < (int)_segments.size())
		return _segments.at(segment);

	/* If track doesn't contain this segment, exits program */
	else
		log_printf(ERROR, "Attempted to retrieve segment s = %d but track only"
				"has %d segments", segment, _segments.size());
	exit(1);
}



/**
 * Returns a vector of this track's segments
 * @return vector of segment pointer
 */
inline std::vector<segment*> Track::getSegments() {
	return _segments;
}


/**
 * Return the number of segments along this track
 * @return the number of segments
 */
inline int Track::getNumSegments() {
	return _segments.size();
}


/**
 * Normalizes all of the polar flux values by multiplying by a factor
 * @param factor the factor to scale the flux by
 */
inline void Track::normalizeFluxes(FP_PRECISION factor)  {

	/* Loop over all polar fluxes */
	for (int i = 0; i < 2*GRP_TIMES_ANG; i++)
		_polar_fluxes[i] *= factor;

	return;
}


#endif /* TRACK_H_ */
