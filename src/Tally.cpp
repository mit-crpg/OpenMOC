/* 
 @file      Tally.cpp
 @brief     contains functions for the Tally class
 @author    Luke Eure
 @date      January 8 2016
*/

#include "Tally.h"

/*
 @brief     constructor for Tally class
*/
Tally::Tally() {}


/*
 @brief     deconstructor
*/
Tally::~Tally() {}


/*
 @brief     sets the tally and _tally_squared equal to zero
*/
void Tally::clear() {
    _tally_count = 0;
    _tally_squared = 0;
    _n = 0;
}


/*
  @brief    gets the current tally count
  @return   a double: the number stored in the tally
*/
double Tally::getCount() {
    return _tally_count;
}


/*
  @brief    returns the standard deviation from the amount held in the tally
  @return   the standard deviation from the amount held in the tally
*/
double Tally::getStandardDeviation() {
    double standardDev = 
        _tally_squared/_n - (_tally_count / _n) * (_tally_count / _n);
    return standardDev;
}


/*
 @brief     overload for += that allows an amount to be added to the tally. Its
            square is added to _tally_squared
 @param     tally_addition an amount to be added
*/
Tally Tally::operator+=(double tally_addition) {
    _tally_count += tally_addition;
    _tally_squared += tally_addition * tally_addition;
    _n++;
    return *this;
}


/*
  @brief    returns the mean of all the entries into the tally
  @return   the mean of all the entries into the tallyl
*/
double Tally::getMean() {
    return _tally_count/_n;
}
