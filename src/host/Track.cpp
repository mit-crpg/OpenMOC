/*
 * Track.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Track.h"


/*
 * Default track constructor
 */
Track::Track() { }



/**
 * Track destructor
 */
Track::~Track() {
	clearSegments();
}


/**
 * Set the values for the track' start and end point and angle
 * @param start_x the x-coordinate at the starting point
 * @param start_y the y-coordinate at the starting point
 * @param end_x the x-coordinate at the ending point
 * @param end_y the y-coordinate at the ending point
 * @param phi the track's azimuthal angle
  */
 void Track::setValues(const double start_x, const double start_y,
		 const double end_x, const double end_y, const double phi) {

	 _start.setCoords(start_x, start_y);
	 _end.setCoords(end_x, end_y);
	 _phi = phi;
 }


 /**
  * Set the track azimuthal weight
  * @param weight the azimuthal weight
  */
 void Track::setAzimuthalWeight(const FP_PRECISION azim_weight) {
     _azim_weight = azim_weight;
 }


 /**
  * Sets the weight of this track at one of the quadrature polar angles
  * @param angle polar angle
  * @param polar_weight the weight of that angle
  */
 void Track::setPolarWeight(const int angle, FP_PRECISION polar_weight) {
	 _polar_weights[angle] = polar_weight;
 }


 /*
  * Set the track azimuthal angle
  * @param phi the azimuthal angle
  */
 void Track::setPhi(const double phi) {
	_phi = phi;
}


void Track::setAzimAngleIndex(const int index) {
	_azim_angle_index = index;
}


/**
 * Adds a segment pointer to this Track's list of segments
 * IMPORTANT: assumes that segments are added in order of their starting
 * location from the track's start point
 * @param segment a pointer to the segment
 */
void Track::addSegment(segment* segment) {

	try {
		_segments.push_back(segment);
	}
	catch (std::exception &e) {
		log_printf(ERROR, "Unable to add a segment to track. Backtrace:"
				"\n%s", e.what());
	}
}


/**
 * Sets whether the incoming flux is at the beginning (false) or
 * end (true) of this Track
 * @param relf_in - beginning (false)/end (true)
 */
void Track::setReflIn(const bool refl_in) {
    _refl_in = refl_in;
}


/**
 * Sets whether the outgoing flux is at the beginning (false) or
 * end (true) of the outgoing Track
 * @param relf_out - beginning (false)/end (true)
 */
void Track::setReflOut(const bool refl_out) {
    _refl_out = refl_out;
}


/**
 * Sets the boundary condition for the incoming flux
 * (vacuum (0) or reflective (1))
 * @param bc_in - boundary condition for the incoming flux
 */
void Track::setBCIn(const bool bc_in) {
	_bc_in = bc_in;
}



/**
 * Sets the boundary condition for the outgoing flux
 * (vacuum (0) or reflective (1))
 * @param bc_out - boundary condition for the outgoing flux
 */
void Track::setBCOut(const bool bc_out) {
	_bc_out = bc_out;
}


/**
 * Sets the incoming track for boundary conditions
 * @param track_in pointer to the incoming track
 */
void Track::setTrackIn(Track *track_in) {
    _track_in = track_in;
}


/**
 * Sets the outgoing track for boundary conditions
 * @param track_out pointer to the outgoing track
 */
void Track::setTrackOut(Track *track_out) {
    _track_out = track_out;
}


/**
 * Sets the ith index in the 2D jagged array of all tracks
 * for this track's incoming track
 * @param i the ith index of the incoming track
 */
void Track::setTrackInI(int i) {
  _track_in_i = i;
}


/**
 * Sets the jth index in the 2D jagged array of all tracks
 * for this track's incoming track
 * @param j the jth index of the incoming track
 */
void Track::setTrackInJ(int j) {
  _track_in_j = j;
}


/**
 * Sets the ith index in the 2D jagged array of all tracks
 * for this track's outgoing track
 * @param i the ith index of the outgoing track
 */
void Track::setTrackOutI(int i) {
  _track_out_i = i;
}


/**
 * Sets the jth index in the 2D jagged array of all tracks
 * for this track's incoming track
 * @param j the jth index of the outgoing track
 */
void Track::setTrackOutJ(int j) {
  _track_out_j = j;
}


/**
 * Returns whether the incoming flux is at the beginning (false) or
 * end (true) of this Track
 * @return beginning (false)/end (true)
 */
bool Track::isReflIn() const {
  return _refl_in;
}


/**
 * Returns whether the outgoing flux is at the beginning (false) or
 * end (true) of the outgoing Track
 * @return beginning (false)/end (true)
 */
bool Track::isReflOut() const {
  return _refl_out;
}


/**
 * Returns the boundary condition for the incoming flux
 * (vacuum (0) or reflective (1))
 * @return boundary condition for incoming flux
 */
bool Track::getBCIn() const {
	return _bc_in;
}


/**
 * Returns the boundary condition for the outgoing flux
 * (vacuum (0) or reflective (1))
 * @return boundary condition for outgoing flux
 */
bool Track::getBCOut() const {
	return _bc_out;
}


/**
 * Returns the track's end point
 * @return a pointer to the track's end point
 */
Point* Track::getEnd() {
    return &_end;
}


/**
 * Returns the track's start point
 * @return a pointer to the track's start point
 */
Point* Track::getStart() {
    return &_start;
}


/**
 * Return the track's azimuthal angle (with respect to the x-axis)
 * @return the aximuthal angle
 */
double Track::getPhi() const {
    return _phi;
}


int Track::getAzimAngleIndex() const {
	return _azim_angle_index;
}

/**
 * Return an array pointer to the track's polar weights
 * @return pointer to the tracks' polar weights
 */
FP_PRECISION* Track::getPolarWeights() {
	return _polar_weights;
}


/**
 * Returns the ith index of this track's incoming track for boundary
 * conditions in the 2D jagged array of all tracks
 * @return the ith index
 */
int Track::getTrackInI() const {
  return _track_in_i;
}


/**
 * Returns the jth index of this track's incoming track for boundary
 * conditions in the 2D jagged array of all tracks
 * @return the jth index
 */
int Track::getTrackInJ() const {
  return _track_in_j;
}


/**
 * Returns the ith index of this track's outgoing track for boundary
 * conditions in the 2D jagged array of all tracks
 * @return the ith index
 */
int Track::getTrackOutI() const {
  return _track_out_i;
}


/**
 * Returns the jth index of this track's outgoing track for boundary
 * conditions in the 2D jagged array of all tracks
 * @return the jth index
 */
int Track::getTrackOutJ() const {
  return _track_out_j;
}


/**
 * Checks whether a point is contained along this track
 * @param a pointer to the point of interest
 * @return true if the point is on the track, false otherwise
 */
bool Track::contains(Point* point) {

	double m; 		// the slope of the track

	/* If the point is outside of the bounds of the start and end points of the
	 * track it does not lie on the track */
	if (!(((point->getX() <= _start.getX()+1.0E-2 &&
			point->getX() >= _end.getX()-1.0E-2)
		|| (point->getX() >= _start.getX()-1.0E-2 &&
				point->getX() <= _end.getX()+1.0E-2))
		&&
		((point->getY() <= _start.getY()+1.0E-2 &&
				point->getY() >= _end.getY()-1.0E-2)
		|| (point->getY() >= _start.getY()-1.0E-2 &&
				point->getY() <= _end.getY()+1.0E-2)))) {

		return false;
	}


	/* If the track is vertical */
	if (fabs(_phi - M_PI / 2) < 1E-10) {
		if (fabs(point->getX() - _start.getX()) < 1E-10)
			return true;
		else
			return false;
	}
	/* If the track is not vertical */
	else {
		m = sin(_phi) / cos(_phi);

		/* Use point-slope formula */
		if (fabs(point->getY() - (_start.getY() +
				m * (point->getX() - _start.getX()))) < 1e-10) {

			return true;
		}
		else
			return false;
	}
}


/**
 * Deletes each of this track's segments
 */
void Track::clearSegments() {

	for (int i=0; i < (int)_segments.size(); i++)
		delete _segments.at(i);

	_segments.clear();
}


/**
 * Convert this track's attributes to a character array
 * @return a character array of this track's attributes
 */
std::string Track::toString() {

	std::stringstream string;
	string << "Track: start, x = " << _start.getX() << ", y = " <<
			_start.getY() << " end, x = " << _end.getX() << ", y = "
			<< _end.getY() << ", phi = " << _phi << " azim_weight = " <<
			_azim_weight;

	return string.str();
}
