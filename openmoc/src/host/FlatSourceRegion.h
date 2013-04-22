/**
 * @file FlatSourceRegion.h
 * @brief The FlatSourceRegion class.
 * @date February 3, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef FLATSOURCEREGION_H_
#define FLATSOURCEREGION_H_

#ifdef __cplusplus
#include <omp.h>
#include "Material.h"
#endif

/**
 * @class FlatSourceRegion FlatSourceRegion.h
 *        "openmoc/src/host/FlatSourceRegion.h"
 * @brief This class represents a unique discretized region within the 
 *        geometry in which OpenMOC approximates the scalar flux and source
 *        as constant.
 */
class FlatSourceRegion {

private:
    /** A static counter for the number of flat source regions 
     *  in a simulation */    
    static int _n;
    /** A monotonically increasing unique ID for each flat source 
     *  region created */
    int _uid;
    /** A pointer to this flat source region's material */
    Material* _material;
    /** The flat source region's volume approximated by the sum of track
     *  segment lengths within the FSR */
    FP_PRECISION _volume;

public:
    FlatSourceRegion();
    virtual ~FlatSourceRegion();
    int getUid() const;
    Material* getMaterial();
    FP_PRECISION getVolume() const;

    void setMaterial(Material* material);
    void setVolume(FP_PRECISION volume);
    void incrementVolume(FP_PRECISION volume);
};


#endif /* FLATSOURCEREGION_H_ */
