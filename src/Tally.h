/* 
 @file      Tally.h
 @brief     contains Tally class
 @author    Luke Eure
 @date      January 8 2016
*/

#ifndef TALLY_H
#define TALLY_H

class Tally {

public:    
    Tally();
    virtual ~Tally();
    
    void clear();
    double getCount();
    double getStandardDeviation();
    double getMean();
    Tally operator+=(double tally_addition);

private:

    /** count stored in the tally */
    double _tally_count;

    /** sum of squares of tallied amounts */
    double _tally_squared;

    /** number of items added to the tally */
    int _n;
};

#endif
