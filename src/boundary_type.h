/**
 * @file boundary_type.h
 * @details The boundaryType enum.
 * @date January 10, 2015
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef BOUNDARY_TYPE_H_
#define BOUNDARY_TYPE_H_

/**
 * @enum boundaryType
 * @brief The types of boundary conditions supported by OpenMOC for Surfaces.
 */
enum boundaryType {
  /** A vacuum boundary condition */
  VACUUM,

  /** A reflective boundary condition */
  REFLECTIVE,

  /** A periodic boundary condition */
  PERIODIC,

  /* Boundary between two domains (only in domain-decomposed geometry) */
  INTERFACE,

  /** No boundary type (typically an interface between flat source regions) */
  BOUNDARY_NONE
};

#endif /* BOUNDARY_TYPE_H_ */
