/**
 * @file Track2D.h
 * @brief The 2D Track class.
 * @date May 9, 2015
 * @author Samuel Shaner, MIT Course, 22 (shaner@mit.edu)
 */

#ifndef TRACK2D_H_
#define TRACK2D_H_

#ifdef __cplusplus
#include "Python.h"
#include "Point.h"
#include "Material.h"
#include "Track.h"
#include <vector>
#endif

/**
 * @class Track2D Track2D.h "src/Track2D.h"
 * @brief A 2D Track represents a characteristic line across the geometry.
 * @details A 2D Track has particular starting and ending points on the
 *          boundaries of the geometry and an azimuthal and polar angle.
 */
class Track2D : public Track {

protected:

public:
  Track2D();
  virtual ~Track2D();
  void setValues(const double start_x, const double start_y,
                 const double end_x, const double end_y,
                 const double phi);

  void setCoords(double x0, double y0, double x1, double y1);

  std::string toString();
};



#endif /* TRACK2D_H_ */
