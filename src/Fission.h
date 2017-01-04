/* 
 @file    Fission.h
 @brief   contains Fission class
 @author  Luke Eure
 @date    February 17 2016
*/

#ifndef FISSION_H
#define FISSION_H

#include <vector>
#include <iostream>

#include "Neutron.h"
#include "Point.h"

class Fission {
public:
  Fission(int num_sites);
  virtual ~Fission();

  void newBatch();
  void sampleSite(Neutron *neutron);
  void add(Point* position, Neutron* neutron);

private:

  /** the fission locations from the last batch */
  std::vector <Point*>* _old_fission_bank;

  /** the fission locations from this batch */
  std::vector <Point*>* _new_fission_bank;
  
  /** a placeholder used when switching banks */
  std::vector <Point*>* _temp_fission_bank;

  /** the number of locations, the size of the fission banks */
  int _num_sites;

  /** the number of fission sites currently populating the vector */
  int _current_num_sites;

  /** the actual number of fission sites in the old fission bank */
  int _old_actual_num_sites;

  /** the site being sampled by a neutron */
  int _site_sampled;

};
#endif
