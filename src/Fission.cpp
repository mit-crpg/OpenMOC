/* 
 @file    Fission.cpp
 @brief   contains functions for the Fission class
 @author  Luke Eure
 @date    February 17 2016
*/

#include "Fission.h"

/*
 @brief   constructor for Fission class
*/
Fission::Fission(int num_sites) {

  _old_fission_bank = new std::vector <Point*>;
  _new_fission_bank = new std::vector <Point*>;

  // reserve space and allocate each point
  _old_fission_bank->reserve(num_sites);
  _new_fission_bank->reserve(num_sites);
  _num_sites = num_sites;
  for (int i=0; i<_num_sites; ++i) {
    Point* dummy_site_new = new Point();
    Point* dummy_site_old = new Point();
    _new_fission_bank->push_back(dummy_site_new);
    _old_fission_bank->push_back(dummy_site_old);
  }
}


/*
 @brief   desconstructor for Fission Class
*/
Fission::~Fission() {
  delete _old_fission_bank;
  delete _new_fission_bank;
}


/*
 @brief   sets the old fission bank to contain the locations from the most
      recent batch, clears the new fission bank for the next round
*/
void Fission::newBatch() {

  _temp_fission_bank = _old_fission_bank;
  _old_fission_bank = _new_fission_bank;
  _new_fission_bank = _temp_fission_bank;
  _old_actual_num_sites = _current_num_sites;
  _site_sampled = 0;
  _current_num_sites = 0;
}


/*
 @brief   assigns a given neutron a position sampled from the fission bank
 @param   neutron a pointer to a Neutron object
*/
void Fission::sampleSite(Neutron *neutron) {
  int index; 
  
  // pick an index
  if (_old_actual_num_sites == 0)
    log_printf(ERROR, "tried to sample fission site, but none exist %s");
  if (_site_sampled < _old_actual_num_sites and _site_sampled < _num_sites)
    index = _site_sampled;
  else
    index = neutron->rand() % _old_actual_num_sites;
  
  // set neutron's position from the fission bank
  Point* site = _old_fission_bank->at(index);
  neutron->setPosition(0, site->getX());
  neutron->setPosition(1, site->getY());
  neutron->setPosition(2, site->getZ());

  _site_sampled++;
}


/*
 @brief   add a location to new_fission_bank
 @param   position a position to be added
*/
void Fission::add(Point* position, Neutron* neutron) {

  if (_current_num_sites < 6 and _current_num_sites > 0) {

  }
  //Point* dummy_site = new Point();
  double x = double(position->getX());
  double y = double(position->getY());
  double z = double(position->getZ());

  int index;

  // add the location to the end of the fission bank if there is room
  if (_current_num_sites < _num_sites)
  {
    //_new_fission_bank->push_back(dummy_site);
   _new_fission_bank->at(_current_num_sites)->setX(x);
    _new_fission_bank->at(_current_num_sites)->setY(y);
    _new_fission_bank->at(_current_num_sites)->setZ(z);
   index = _current_num_sites;
  }

  // otherwise replace a random site with the new site
  else {
    index = neutron->rand() % _num_sites;
    _new_fission_bank->at(index)->setX(x);
    _new_fission_bank->at(index)->setY(y);
    _new_fission_bank->at(index)->setZ(z);
  }

  _current_num_sites++;
}
