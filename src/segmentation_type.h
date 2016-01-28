/**
 * @file segmentation_type.h
 * @details The segmentationType enum.
 * @date January 27, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef SEGMENTATION_TYPE_H_
#define SEGMENTATION_TYPE_H_

/**
 * @enum segmentationType
 * @brief The types of Track segmentation supported by OpenMOC.
 */
enum segmentationType {
  /** Explicit 3D segments */
  EXPLICIT,

  /** Axial on-the-fly 3D segment formation by 3D track */
  OTF_TRACKS,

  /** Axial on-the-fly 3D segment formation by z-stack */
  OTF_STACKS

};

#endif /* SEGMENTATION_TYPE_H_ */
