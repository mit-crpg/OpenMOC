/** 
 * @file    Neutron.cpp
 * @brief   contains functions for the Neutron class
 * @author  Luke Eure
 * @date    January 8 2016
 */

#include "Neutron.h"

/**
 * @brief   constructor for Neutron class
 * @param   neutron_num to be saved as the neutron's id number
 */
Neutron::Neutron(int neutron_num) {
  _neutron_alive = true;
  _neutron_direction.resize(3);
  _id = neutron_num;
  const int global_seed = 12;
  _seed = _id + global_seed;
  rand_r(&_seed);
}


/**
 * @brief   deconstructor for Neutron class
 */
Neutron::~Neutron() {}


/**
 * @brief   tells if the neutron is alive
 * @return  a bool (true=neutron alive. false=neutron dead)
 */
bool Neutron::alive(){
  return _neutron_alive;
}


/**
 * @brief   returns the energy group of the neutron
 * @return  an int: the energy group of the neutron
 */
int Neutron::getGroup() {
  return _neutron_group;
}


/**
 * @brief   moves the neutron a given distance
 * @param   distance the distance the neutron should be moved
 */
void Neutron::move(double distance) {
  _xyz.setX(_xyz.getX() + _neutron_direction[0] * distance);
  _xyz.setY(_xyz.getY() + _neutron_direction[1] * distance);
  _xyz.setZ(_xyz.getZ() + _neutron_direction[2] * distance);
}


/**
 * @brief   reverses a neutron's direction of motion along an axis
 * @param   axis the axis along which the nuetron should be reflected
 */
void Neutron::reflect(int axis) {
  _neutron_direction[axis] *= -1;
}


/**
 * @brief   sets the neutron's position along an axis
 * @param   axis the axis along which the position will be set
 * @param   value the value to which the position will be set
 */ 
void Neutron::setPosition(int axis, double value) {
  if (axis==0)
    _xyz.setX(value);
  if (axis==1)
    _xyz.setY(value);
  if (axis==2)
    _xyz.setZ(value);
}


/**
 * @brief   kills the neutron
 */
void Neutron::kill() {
  _neutron_alive = false;
}


/**
 * @brief   gets the position of the neutron along a certain axis
 * @param   axis an int containing the axis along which the position will be
 *          returned
 * @return  a double denoting the neutron's position along axis
 */
double Neutron::getPosition(int axis) {
  double pos;
  if (axis==0)
    pos = _xyz.getX();
  if (axis==1)
    pos = _xyz.getY();
  if (axis==2)
    pos = _xyz.getZ();
  return pos;
}


/**
 * @brief   gets the position vector of the neutron
 * @param   a Point* to be pointed at a Point containing the neutron's position
 */
void Neutron::getPositionVector(Point* &position) {
  position = &_xyz;
}


/**
 * @brief   gets the direction of the neutron along a certain axis
 * @param   axis an int containing the axis along which the direction will be
 *          returned
 * @return  a double denoting the neutron's direction along axis
 */
double Neutron::getDirection(int axis) {
  return _neutron_direction[axis];
}


/**
 * @brief   sets the direction of the neutron along a certain axis
 * @param   axis the axis along which the direction will be set
 * @param   magnitude the magnitude of the nuetron's motion along axis
 */
void Neutron::setDirection(int axis, double magnitude) {
  _neutron_direction[axis] = magnitude;
}


/**
 * @brief   gets the neutron's distance from a given point
 * @param   coord a vector denoting the point to find the neutron's distance from
 * @return  the neutron's distance from that point
 */
double Neutron::getDistance(Point* coord) {
  double dis;
  dis = _xyz.distanceToPoint(coord);
  return dis;
}


/**
 * @brief   set the neutron's group
 * @param   new_group the new energy group of the neutron
 */
void Neutron::setGroup(int new_group) {
  _neutron_group = new_group;
}


/**
 * @brief   sets the neutron's direction to a random direction based on
 *          the neutron's random number seed
 */
void Neutron::sampleDirection() {

  // sample azimuthal angle
  double phi = 2 * M_PI * arand();

  // sample cosine of the polar angle
  double mu = 2 * arand() - 1.0;

  // set diection
  _neutron_direction[0] = sqrt(1 - mu*mu) * cos(phi);
  _neutron_direction[1] = sqrt(1 - mu*mu) * sin(phi);
  _neutron_direction[2] = mu;
}


/**
 * @brief   returns a pseudo-random number using the seed between 0 and 1
 * @return  a psuedo-random number between 0 and 1
 */
double Neutron::arand() {
  double r = rand_r(&_seed);
  return (double) r / (double) RAND_MAX;
}


/**
 * @brief   returns a pseudo-random number using the seed
 * @return  a psuedo-random number
 */
int Neutron::rand() {
  return rand_r(&_seed);
}


/**
 * @brief   samples the neutron energy group after a scattering event
 * @param   scattering_matrix the scattering cross section matrix
 * @return  the neutron group after scattering
 */
int Neutron::sampleScatteredGroup(std::vector <double> &scattering_matrix) {

  // get the total scattering cross-section from this group
  int num_groups = scattering_matrix.size();
  double scattering_total = 0;
  for (int g=0; g < num_groups; ++g)
    scattering_total += scattering_matrix[g];

  // sample the outgoing scattered energy group
  double r = arand() * scattering_total;
  double scatter_sum = 0.0;
  for (int g=0; g<num_groups; ++g) {
    scatter_sum += scattering_matrix[g];
    if (r<scatter_sum)
      return g;
  }

  // throw error if no group has been found
  log_printf(ERROR,
      "no group found to scatter in to %s");
}


/**
 * @brief   samples an initial neutron energy group after fission
 * @param   chi the neutron emission spectrum from fission
 * @return  the group number of the emitted neutron
 */
int Neutron::sampleEnergyGroup(std::vector <double> chi) {
  double r = arand();
  double chi_sum = 0.0;
  for (int g=0; g<chi.size(); ++g) {
    chi_sum += chi[g];
    if (r<chi_sum) {
      return g;
    }
  }
  return chi.size() - 1;
}
