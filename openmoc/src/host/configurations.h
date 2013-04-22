/**
 * @file configurations.h
 * @date January 16, 2012
 * @author William Boyd, Course 22, MIT (wboyd@mit.edu)
 */

#ifndef CONFIGURATIONS_H_
#define CONFIGURATIONS_H_


/******************************************************************************
 ****************************** USER DEFINED **********************************
 *****************************************************************************/

/** Floating point precision (DOUBLE or SINGLE) */
#define SINGLE

/** The floating point precision (float or double) to use */
#ifdef DOUBLE
#define FP_PRECISION double
#else
#define FP_PRECISION float
#endif

#endif
