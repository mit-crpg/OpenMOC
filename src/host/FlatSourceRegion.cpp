/*
 * FlatSourceRegion.cu
 *
 *  Created on: Feb 3, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#include "FlatSourceRegion.h"

/* _n keeps track of the number FSRs of instantiated */
int FlatSourceRegion::_n = 0;

/**
 * FlatSourceRegion constructor
 */
FlatSourceRegion::FlatSourceRegion() {

	_material = NULL;
	_volume = 0.0;
	_uid = _n;
	_n++;

	/* Zero the region's flux and source */
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
		_flux[e] = 0.0;
		_old_flux[e] = 0.0;
		_source[e] = 0.0;
		_old_source[e] = 0.0;
	}

#if USE_OPENMP
	omp_init_lock(&_flux_lock);
#endif

}


/**
 * Default destructor
 */
FlatSourceRegion::~FlatSourceRegion() {
#if USE_OPENMP
	omp_destroy_lock(&_flux_lock);
#endif
}


/**
 * Return the FSR's uid
 * @return the FSR's uid
 */
int FlatSourceRegion::getUid() const {
	return _uid;
}

/**
 * Returns a pointer to this region's material
 * @return a pointer to a material
 */
Material* FlatSourceRegion::getMaterial() {
	return _material;
}


/**
 * Returns the region's volume as computed by the number of track segments
 * and their cumulative lengths in this region
 * @preturn the region's volume
 */
FP_PRECISION FlatSourceRegion::getVolume() const {
    return _volume;
}


/**
 * Returns an array of the multi energy group fluxes tallied in this region
 * @return a flux array
 */
FP_PRECISION* FlatSourceRegion::getFlux() {
    return _flux;
}


/**
 * Returns an array of the multi energy group fluxes tallied in this region
 * from the previous iteration
 * @return a flux array
 */
FP_PRECISION* FlatSourceRegion::getOldFlux() {
    return _old_flux;
}


/**
 * Returns a array of the multi-energy group source term in this region
 * computed by the solver
 * @return the old source array
 */
FP_PRECISION* FlatSourceRegion::getSource() {
    return _source;
}


/**
 * Returns a array of the multi-energy group source term in this region from
 * the previous iteration computed by the solver
 * @return the old source array
 */
FP_PRECISION* FlatSourceRegion::getOldSource() {
    return _old_source;
}


/**
 * Return an array of ratios of source / sigma_t for this flat source region
 * for each energy group
 * @return array of source / sigma_t ratio
 */
FP_PRECISION* FlatSourceRegion::getRatios() {
	return _ratios;
}


/**
 * Set this region's material
 * @param material pointer to a material
 */
void FlatSourceRegion::setMaterial(Material* material) {
	_material = material;
}


/**
 * Sets the reigon's volume
 * @param volume the region's volume
 */
void FlatSourceRegion::setVolume(FP_PRECISION volume) {
    _volume = volume;
}


/**
 * Increment this FSR's volume by some amount corresponding to a segment length
 * @param volume the amount to increment by
 */
void FlatSourceRegion::incrementVolume(FP_PRECISION volume) {
	_volume += volume;
}


/**
 * Set the scalar flux for one energy group inside this FSR
 * @param energy the energy group index
 * @param flux the scalar flux
 */
void FlatSourceRegion::setFlux(int energy, FP_PRECISION flux) {

	if (energy < -1 || energy >= NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Attempted to set the scalar flux for FSR uid = %d "
				"in an energy group which does not exist: %d", _uid, energy);

#if USE_OPENMP
	omp_set_lock(&_flux_lock);
#endif

	_flux[energy] = flux;

#if USE_OPENMP
	omp_unset_lock(&_flux_lock);
#endif

	return;
}



/**
 * Set the scalar flux for one energy group inside this FSR
 * @param energy the energy group index
 * @param flux the scalar flux from the previous iteration
 */
void FlatSourceRegion::setOldFlux(int energy, FP_PRECISION old_flux) {

	if (energy < -1 || energy >= NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Attempted to set the old scalar flux for FSR "
				"uid = %d in an energy group which does not exist: %d",
															_uid, energy);

	_old_flux[energy] = old_flux;

	return;
}



/**
 * Increment the scalar flux for one energy group inside this FSR
 * @param energy the energy group index
 * @param flux the scalar flux
 */
void FlatSourceRegion::incrementFlux(int energy, FP_PRECISION flux) {

	if (energy < -1 || energy >= NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Attempted to increment the scalar flux for FSR uid = "
				"%d in an energy group which does not exist: %d", _uid, energy);

#if USE_OPENMP
	omp_set_lock(&_flux_lock);
#endif

	_flux[energy] += flux;

#if USE_OPENMP
	omp_unset_lock(&_flux_lock);
#endif

	return;
}


/**
 * Increment the scalar flux for all energy groups inside this FSR
 * @param flux the scalar flux for all energy groups
 */
void FlatSourceRegion::incrementFlux(FP_PRECISION* flux) {

#if USE_OPENMP
	omp_set_lock(&_flux_lock);
#endif

	for (int e=0; e < NUM_ENERGY_GROUPS; e++)
		_flux[e] += flux[e];

#if USE_OPENMP
	omp_unset_lock(&_flux_lock);
#endif

	return;
}



/**
 * Set the source for one energy group for this FSR
 * @param energy the energy group index
 * @param source the source
 */
void FlatSourceRegion::setSource(int energy, FP_PRECISION source) {
	if (energy < -1 || energy >= NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Attempted to set the source for FSR uid = %d "
				"in an energy group which does not exist: %d", _uid, energy);

	#pragma omp critical
	{
		_old_source[energy] = source;
	}

	return;
}


/**
 * Set the old source from the previous iteration for one energy group
 * inside this FSR
 * @param energy the energy group index
 * @param old_source the old source
 */
void FlatSourceRegion::setOldSource(int energy, FP_PRECISION old_source) {
	if (energy < -1 || energy >= NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Attempted to set the old source for FSR uid = %d "
				"in an energy group which does not exist: %d", _uid, energy);

	#pragma omp critical
	{
		_old_source[energy] = old_source;
	}

	return;
}


/**
 * Normalizes all of the scalar flux values by multiplying by a factor
 * @param factor the factor to scale the flux by
 */
void FlatSourceRegion::normalizeFluxes(FP_PRECISION factor) {

	/* Loop over all energy groups */
	#pragma omp critical
	{
		for (int e=0; e < NUM_ENERGY_GROUPS; e++)
			_flux[e] *= factor;
	}

	return;
}


/**
 * Compute source / sigma_t for each energy group
 */
void FlatSourceRegion::computeRatios() {

	FP_PRECISION* sigma_t = _material->getSigmaT();

	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		_ratios[e] = _source[e]/sigma_t[e];

	return;
}


/**
 * Compute the volumetric fission rate in this flat source region by adding
 * up the fission rates in each energy group. This mehtod assumes that fixed
 * source iteration has already been run since it uses the flux stored in
 * this region
 */
FP_PRECISION FlatSourceRegion::computeFissionRate() {

	FP_PRECISION power = 0.0;
	FP_PRECISION* sigma_f = _material->getSigmaF();

	/* Add the fission rates from each energy group */
	for (int e=0; e < NUM_ENERGY_GROUPS; e++)
		power += sigma_f[e] * _flux[e];

	/* Multiply by volume of FSR */
	power *= _volume;

	return power;
}
