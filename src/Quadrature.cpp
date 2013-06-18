#include "Quadrature.h"


/**
 * @brief Constructor initializes arrays with polar angles and weights.
 * @param type the quadrature type (LEONARD or TABUCHI)
 * @param num_polar_angles the number of polar angles (1, 2, or 3)
 */
Quadrature::Quadrature(quadratureType type, int num_polar_angles) {

    if (num_polar_angles != 1 && num_polar_angles != 2 &&
        num_polar_angles != 3)
        log_printf(ERROR, "Unable to set the number of polar angles to %d. Only"
               "1, 2 and 3 polar angles are supported by OpenMOC.");

    _num_polar_angles = num_polar_angles;

    /** Initialize memory */
    _sinthetas = new FP_PRECISION[_num_polar_angles];
    _weights = new FP_PRECISION[_num_polar_angles];
    _multiples = new FP_PRECISION[_num_polar_angles];

    /* If quadrature type is TABUCHI */
    if (type == TABUCHI) {

        _type = TABUCHI;

        if (_num_polar_angles == 1) {
            _sinthetas[0] = 0.798184;
            _weights[0] = 1.0;
            _multiples[0] = _sinthetas[0] * _weights[0];
        }

        else if (_num_polar_angles == 2) {
            _sinthetas[0] = 0.363900;
            _sinthetas[1] = 0.899900;
            _weights[0] = 0.212854;
            _weights[1] = 0.787146;
            _multiples[0] = _sinthetas[0] * _weights[0];
            _multiples[1] = _sinthetas[1] * _weights[1];
        }

        else if (_num_polar_angles == 3) {
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
            log_printf(ERROR, "TABUCHI type polar quadrature supports 1, 2, or "
                       "3 polar angles but %d are defined", _num_polar_angles);
    }

    /* If quadrature type is LEONARD */
    else if (type == LEONARD) {

        _type = LEONARD;
                
        if (_num_polar_angles == 2) {
            _sinthetas[0] = 0.273658;
            _sinthetas[1] = 0.865714;
            _weights[0] = 0.139473;
            _weights[1] = 0.860527;
            _multiples[0] = _sinthetas[0] * _weights[0];
            _multiples[1] = _sinthetas[1] * _weights[1];
        }

        else if (_num_polar_angles == 3) {
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
            log_printf(ERROR, "LEONARD type polar quadrature supports 2, or 3"
                       "polar angles but %d are defined", _num_polar_angles);
        }

    }

    else
        log_printf(ERROR, "LEONARD and TABUCHI polar quadrature types "
                   "supported, but unknown type given");
}


/**
 * @brief Destructor deletes arrray of sines of the polar angles, the weights
 *        of the polar angles and the products of the sines and weights.
 */
Quadrature::~Quadrature() { 
    delete [] _sinthetas;
    delete [] _weights;
    delete [] _multiples;
}


/**
 * @brief Returns the polar quadrature type (LEONARD or TABUCHI)
 * @return the polar quadrature type
 */
quadratureType Quadrature::getType() const {
    return _type;
}


/**
 * @brief Returns the number of polar angles (1, 2 or 3).
 * @return the number of polar angles
 */
int Quadrature::getNumPolarAngles() const {
    return _num_polar_angles;
}


/**
 * @brief Returns the \f$ sin(\theta)\f$ value for a particular polar angle.
 * @param n index of the polar angle of interest
 * @return the value of \f$ \sin(\theta) \f$ for this polar angle
 */
FP_PRECISION Quadrature::getSinTheta(const int n) const {

  if (n < 0 || n >= _num_polar_angles)
        log_printf(ERROR, "Attempted to retrieve sintheta for polar angle = %d"
                   " but only %d polar angles are defined", _num_polar_angles);

    return _sinthetas[n];
}


/**
 * @brief Returns the weight value for a particular polar angle.
 * @param n index of the polar angle of interest
 * @return the weight for a polar angle 
 */
FP_PRECISION Quadrature::getWeight(const int n) const {

    if (n < 0 || n >= _num_polar_angles)
        log_printf(ERROR, "Attempted to retrieve the weight for polar "
                   "angle = %d but only %d polar angles are defined", 
                   _num_polar_angles);

    return _weights[n];
}


/**
 * @brief Returns the multiple value for a particular polar angle. (the 
 *        sine of a polar angle multiplied by its weight).
 * @param n index of the polar angle of interest
 * @return the value of the sine of the polar angle multiplied with its weight
 */
FP_PRECISION Quadrature::getMultiple(const int n) const {

    if (n < 0 || n >= _num_polar_angles)
        log_printf(ERROR, "Attempted to retrieve the multiple for polar "
                   "angle = %d but only %d polar angles are defined", 
                   _num_polar_angles);

    return _multiples[n];
}


/**
 * @brief Returns a pointer to the quadrature's sintheta array.
 * @return a pointer to the sintheta array
 */
FP_PRECISION* Quadrature::getSinThetas() {
    return _sinthetas;
}


/**
 * @brief Returns a pointer to the quadrature's polar weights array.
 * @return a pointer to the polar weights array
 */
FP_PRECISION* Quadrature::getWeights() {
    return _weights;
}



/**
 * @brief Returns a pointer to the quadrature's multiples array.
 * @return a pointer to the multiples array
 */
FP_PRECISION* Quadrature::getMultiples() {
    return _multiples;
}


/**
 * @brief Converts this quadrature to a character array of its attributes.
 * @details The character array includes the polar quadrature type, the number
 *          of polar angles, and the values of the sine of each polar angle,
 *          the weight of each polar angle, and the product of the sine
 *          and weight of each polar angle.
 * @return a character array of this quadrature's attributes
 */
std::string Quadrature::toString() {
    std::stringstream string;

    string << "Quadrature type = ";

    if (_type == LEONARD)
        string << "LEONARD";
    else if (_type == TABUCHI)
        string << "TABUCHI";

    string << "\nsinthetas = ";
    for (int p = 0; p < _num_polar_angles; p++)
        string << _sinthetas[p] << ", ";
    
    string << "\nweights = ";
    for (int p = 0; p < _num_polar_angles; p++)
        string << _weights[p] << ", ";

    string << "\nmultiples= ";
    for (int p = 0; p < _num_polar_angles; p++)
        string << _multiples[p] << ", ";

    return string.str();
}
