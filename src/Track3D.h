/**
 * @file Track3D.h
 * @brief The 3D Track class.
 * @date May 9, 2015
 * @author Samuel Shaner, MIT Course, 22 (shaner@mit.edu)
 */

#ifndef TRACK3D_H_
#define TRACK3D_H_

#ifdef __cplusplus
#include "Python.h"
#include "Point.h"
#include "Material.h"
#include "Track.h"
#include <vector>
#endif



/**
 * @class Track3D Track3D.h "src/Track3D.h"
 * @brief A 3D Track represents a characteristic line across the geometry.
 * @details A 3D Track has particular starting and ending points on the
 *          boundaries of the geometry and an azimuthal and polar angle.
 */
class Track3D : public Track {

protected:

  /** The polar angle for the Track */
  double _theta;

  /** The polar angle index into the the TrackCycle _track_stacks 2D ragged array */
  int _polar_index;
  
public:
  Track3D();
  virtual ~Track3D();

  /* Setters */
  void setValues(const double start_x, const double start_y,
                 const double start_z, const double end_x,
                 const double end_y, const double end_z,
                 const double phi, const double theta);
  void setTheta(const double theta);
  void setPolarIndex(const int index);
  void setCoords(double x0, double y0, double z0, double x1, double y1, double z1);
  
  /* Getters */
  double getTheta() const;
  int getPolarIndex() const;

  std::string toString();
};


#endif /* TRACK3D_H_ */
