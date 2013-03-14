/*
 * Quadrature.h
 *
 *  Created on: Jan 20, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <sstream>
#include "configurations.h"
#include "log.h"

enum quadratureType {
	LEONARD,
	TABUCHI
};


class Quadrature {

private:

	quadratureType _type;
	FP_PRECISION _sinthetas[NUM_POLAR_ANGLES];
	FP_PRECISION _weights[NUM_POLAR_ANGLES];
	FP_PRECISION _multiples[NUM_POLAR_ANGLES];

public:

	Quadrature(quadratureType type);
	virtual ~Quadrature();
	quadratureType getType() const;
	FP_PRECISION getSinTheta(const int n) const;
	FP_PRECISION getWeight(const int n) const;
	FP_PRECISION getMultiple(const int n) const;
	FP_PRECISION* getSinThetas();
    FP_PRECISION* getWeights();
    FP_PRECISION* getMultiples();
    std::string toString();

};

#endif /* QUADRATURE_H_ */
