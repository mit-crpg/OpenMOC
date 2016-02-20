/**
 * @file ringify_type.h
 * @details The ringifyType enum.
 * @date February 16, 2016
 * @author Derek Gaston, MIT, Course 22 (gastdr@mit.edu)
 */

#ifndef RINGIFY_TYPE_H_
#define RINGIFY_TYPE_H_

/**
 * @enum ringifyType
 * @brief The method to use when subdividing cells into rings
 */
enum ringifyType {
  /** Max ring radius chosen such that it has the same area as the
   * bounding box of the universe it lies in */
  EQUIVALENT_AREA,

  /** Max ring radius chosen to be the distance from the center to
   * a corner of the bounding box for the universe it lies in */
  MAX_DISTANCE
};

#endif /* RINGIFY_TYPE_H_ */
