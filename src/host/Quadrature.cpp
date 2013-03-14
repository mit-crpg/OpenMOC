/*
 * Quadrature.cpp
 *
 *  Created on: Jan 20, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Quadrature.h"


/**
 * Quadrature constructor
 * @param type the quadrature type (Leonard or Tabuchi)
 */
Quadrature::Quadrature(quadratureType type) {

	/* If TabuchiYamomoto */
	if (type == TABUCHI) {

		_type = TABUCHI;

		if (NUM_POLAR_ANGLES == 1) {
			_sinthetas[0] = 0.798184;
			_weights[0] = 1.0;
			_multiples[0] = _sinthetas[0] * _weights[0];
		}

		else if (NUM_POLAR_ANGLES == 2) {
			_sinthetas[0] = 0.363900;
			_sinthetas[1] = 0.899900;
			_weights[0] = 0.212854;
			_weights[1] = 0.787146;
			_multiples[0] = _sinthetas[0] * _weights[0];
			_multiples[1] = _sinthetas[1] * _weights[1];
		}

		else if (NUM_POLAR_ANGLES == 3) {
			_sinthetas[0] = 0.166648;
			_sinthetas[1] = 0.537707;
			_sinthetas[2] = 0.932954;
			_weights[0] = 0.046233;
			_weights[1] = 0.283619;
			_weights[2] = 0.670148;
			_multiples[0] = _sinthetas[0] * _weights[0];
			_multiples[1] = _sinthetas[1] * _weights[1];
			_multiples[2] = _sinthetas[2] * _weights[2];
		}

		else
			log_printf(ERROR, "Tabuchi type quadrature supports 1, 2, or 3 polar"
					" angles but %d are defined", NUM_POLAR_ANGLES);
	}

	/* If Leonard */
	else if (type == LEONARD) {

		_type = LEONARD;

		if (NUM_POLAR_ANGLES == 2) {
			_sinthetas[0] = 0.273658;
			_sinthetas[1] = 0.865714;
			_weights[0] = 0.139473;
			_weights[1] = 0.860527;
			_multiples[0] = _sinthetas[0] * _weights[0];
			_multiples[1] = _sinthetas[1] * _weights[1];
		}

		else if (NUM_POLAR_ANGLES == 3) {
			_sinthetas[0] = 0.099812;
			_sinthetas[1] = 0.395534;
			_sinthetas[2] = 0.891439;
			_weights[0] = 0.017620;
			_weights[1] = 0.188561;
			_weights[2] = 0.793819;
			_multiples[0] = _sinthetas[0] * _weights[0];
			_multiples[1] = _sinthetas[1] * _weights[1];
			_multiples[2] = _sinthetas[2] * _weights[2];
		}

		else {
			log_printf(ERROR, "Leonard type quadrature supports 2, or 3 polar"
					" angles but %d are defined", NUM_POLAR_ANGLES);
		}

	}

	else
		log_printf(ERROR, "Leonard and Tabuchi quadrature types supported, but unknown"
				" type given");
}


/**
 * Quadrature destructor
 */
Quadrature::~Quadrature() { }



/**
 * Returns the quadrature type (Leonard or Tabuchi)
 * @return the quadrature type
 */
quadratureType Quadrature::getType() const {
	return _type;
}


/**
 * Returns the sintheta value for a particular polar angle
 * @param n the polar angle of interest
 * @return the sintheta value
 */
FP_PRECISION Quadrature::getSinTheta(const int n) const {

	if (n > -1 && n < NUM_POLAR_ANGLES)
		return _sinthetas[n];

	else {
		log_printf(ERROR, "Attempted to retrieve sintheta for polar angle = %d"
				" but only %d polar angles are defined", NUM_POLAR_ANGLES);
	}

	exit(0);
}


/**
 * Returns the polar weight value for a particular polar angle
 * @param n the polar angle of interest
 * @return the polar weight value
 */
FP_PRECISION Quadrature::getWeight(const int n) const {

	if (n > -1 && n < NUM_POLAR_ANGLES)
		return _weights[n];

	else
		log_printf(ERROR, "Attempted to retrieve the weight for polar angle = %d"
				" but only %d polar angles are defined", NUM_POLAR_ANGLES);
	exit(0);
}


/**
 * Returns the multiple value for a particular polar angle
 * @param n the polar angle of interest
 * @return the multiple value
 */
FP_PRECISION Quadrature::getMultiple(const int n) const {

	if (n > -1 && n < NUM_POLAR_ANGLES)
		return _multiples[n];

	else
		log_printf(ERROR, "Attempted to retrieve the multiple for polar angle = %d"
				" but only %d polar angles are defined", NUM_POLAR_ANGLES);
	exit(0);
}


/**
 * Returns a pointer to the quadrature's sintheta array
 * @return a pointer to the sintheta array
 */
FP_PRECISION* Quadrature::getSinThetas() {
	return _sinthetas;
}


/**
 * Returns a pointer to the quadrature's polar weights array
 * @return a pointer to the polar weights array
 */
FP_PRECISION* Quadrature::getWeights() {
	return _weights;
}



/**
 * Returns a pointer to the quadrature's multiples array
 * @return a pointer to the multiples array
 */
FP_PRECISION* Quadrature::getMultiples() {
	return _multiples;
}


/**
 * Converts this quadrature to a character array of its attributes
 * @param a character array of this quadrature's attributes
 */
std::string Quadrature::toString() {
	std::stringstream string;

	string << "Quadrature type = ";

	if (_type == LEONARD)
		string << "LEONARD";
	else if (_type == TABUCHI)
		string << "TABUCHI";

	string << "\nsinthetas = ";
	for (int p = 0; p < NUM_POLAR_ANGLES; p++)
		string << _sinthetas[p] << ", ";

	string << "\nweights = ";
	for (int p = 0; p < NUM_POLAR_ANGLES; p++)
		string << _weights[p] << ", ";

	string << "\nmultiples= ";
	for (int p = 0; p < NUM_POLAR_ANGLES; p++)
		string << _multiples[p] << ", ";

	return string.str();
}
