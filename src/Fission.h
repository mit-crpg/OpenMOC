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
  Fission();
  virtual ~Fission();

  void newBatch();
  void sampleSite(Neutron *neutron);
  void add(Point* position);

private:

  /** the fission locations from the last batch */
  std::vector <Point*>* _old_fission_bank;

  /** the fission locations from this batch */
  std::vector <Point*>* _new_fission_bank;
  
  /** a placeholder used when switching banks */
  std::vector <Point*>* _temp_fission_bank;

};
#endif
