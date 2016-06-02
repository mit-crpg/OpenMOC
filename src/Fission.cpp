/* 
 @file      Fission.cpp
 @brief     contains functions for the Fission class
 @author    Luke Eure
 @date      February 17 2016
*/

#include "Fission.h"

/*
 @brief     constructor for Fission class
*/
Fission::Fission() {
    _old_fission_bank = new std::vector <Point*>;
    _new_fission_bank = new std::vector <Point*>;
}

/*
 @brief     desconstructor for Fission Class
*/
Fission::~Fission() {
    delete _old_fission_bank;
    delete _new_fission_bank;
}

/*
 @brief     sets the old fission bank to contain the locations from the most
            recent batch, clears the new fission bank for the next round
*/
void Fission::newBatch() {
    _temp_fission_bank = _old_fission_bank;
    _old_fission_bank = _new_fission_bank;
    _new_fission_bank = _temp_fission_bank;
    for (int i; i=0; i < _new_fission_bank->size())
        delete _new_fission_bank->at(i);
    _new_fission_bank->clear();
}

/*
 @brief     assigns a given neutron a position sampled from the fission bank
 @param     neutron a pointer to a Neutron object
*/
void Fission::sampleSite(Neutron *neutron) {
    int index = neutron->rand() % _old_fission_bank->size();
    Point* site;
    site = _old_fission_bank->at(index);
    neutron->setPosition(0, site->getX());
    neutron->setPosition(1, site->getY());
    neutron->setPosition(2, site->getZ());
}

/*
 @brief     add a location to new_fission_bank
 @param     position a position to be added
*/
void Fission::add(Point* position) {
    Point* fission_site = new Point();
    fission_site->setX(position->getX());
    fission_site->setY(position->getY());
    fission_site->setZ(position->getZ());
    _new_fission_bank->push_back(fission_site);
}
