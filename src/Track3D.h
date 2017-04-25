/**
 * @file Track3D.h
 * @brief The 3D Track class.
 * @date May 9, 2015
 * @author Samuel Shaner, MIT Course, 22 (shaner@mit.edu)
 */

#ifndef TRACK3D_H_
#define TRACK3D_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "Point.h"
#include "Material.h"
#include "Track.h"
#include "Track.h"
#include "LocalCoords.h"
#include <vector>
#endif



/** Forward declaration of ExtrudedFSR struct */
struct ExtrudedFSR;


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

  /* Indices that are used to locate the track in the various track arrays */
  int _polar_index;
  int _z_index;
  int _lz_index;
  int _cycle_index;
  int _cycle_track_index;
  int _train_index;

  /** Boolean to indicate whether track is in the same direction as the
   *  reflective track cycle. */
  bool _cycle_fwd;

public:
  Track3D();
  virtual ~Track3D();

  /* Setter methods */
  void setValues(const double start_x, const double start_y,
                 const double start_z, const double end_x,
                 const double end_y, const double end_z,
                 const double phi, const double theta);
  void setTheta(const double theta);
  void setCoords(double x0, double y0, double z0, double x1, double y1,
                 double z1);
  void setPolarIndex(int index);
  void setZIndex(int index);
  void setLZIndex(int index);
  void setCycleIndex(int index);
  void setCycleTrackIndex(int index);
  void setTrainIndex(int index);
  void setCycleFwd(bool fwd);

  /* Getter methods */
  double getTheta() const;
  int getPolarIndex();
  int getZIndex();
  int getLZIndex();
  int getCycleIndex();
  int getCycleTrackIndex();
  int getTrainIndex();
  bool getCycleFwd();

  /* Worker methods */
  std::string toString();
};


#endif /* TRACK3D_H_ */
