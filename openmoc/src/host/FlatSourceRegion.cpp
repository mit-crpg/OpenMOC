#include "FlatSourceRegion.h"

int FlatSourceRegion::_n = 0;

/**
 * @brief Constructor initializes a unique ID and NULL material pointer for
 *        the flat source region.
 */
FlatSourceRegion::FlatSourceRegion() {
    _material = NULL;
    _volume = 0.0;
    _uid = _n;
    _n++;
}


/**
 * @brief Destructor.
 */
FlatSourceRegion::~FlatSourceRegion() { }


/**
 * @brief Return the flat source regoin's unique ID.
 * @return the flat source region's unique ID
 */
int FlatSourceRegion::getUid() const {
    return _uid;
}

/**
 * @brief Returns a pointer to this flat source region's  material.
 * @return a pointer to the material
 */
Material* FlatSourceRegion::getMaterial() {
    return _material;
}


/**
 * @brief Returns the region's volume as computed by the sum of all 
 *        track segment lenths in the flat source region.
 * @return the flat source region's volume
 */
FP_PRECISION FlatSourceRegion::getVolume() const {
    return _volume;
}


/**
 * @brief Set the pointer to the flat source region's material.
 * @param material pointer to the material
 */
void FlatSourceRegion::setMaterial(Material* material) {
    _material = material;
}


/**
 * @brief Sets the flat source region's volume.
 * @param volume the flat source region's volume
 */
void FlatSourceRegion::setVolume(FP_PRECISION volume) {
    _volume = volume;
}


/**
 * @brief Increment this flat source regions's volume by some amount.
 * @details This method is called during segmentation as the volume for the flat
 *          source region is successively approximated by the sum of track
 *          segment lengths within the region. 
 * @param volume the amount to increment by (ie, a track segment length (cm))
 */
void FlatSourceRegion::incrementVolume(FP_PRECISION volume) {
    _volume += volume;
}
