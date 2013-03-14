/*
 * FlatSourceRegion.h
 *
 *  Created on: Feb 3, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef FLATSOURCEREGION_H_
#define FLATSOURCEREGION_H_


#include "Material.h"

#if USE_OPENMP
	#include <omp.h>
#endif


class FlatSourceRegion {
private:
	static int _n;
	int _uid;
	Material* _material;
	FP_PRECISION _volume;
	FP_PRECISION _flux[NUM_ENERGY_GROUPS];
	FP_PRECISION _old_flux[NUM_ENERGY_GROUPS];
	FP_PRECISION _source[NUM_ENERGY_GROUPS];
	FP_PRECISION _old_source[NUM_ENERGY_GROUPS];

	/* Pre-computed Ratio of source / sigma_t */
	FP_PRECISION _ratios[NUM_ENERGY_GROUPS];

	#if USE_OPENMP
		omp_lock_t _flux_lock;
	#endif

public:
	FlatSourceRegion();
	virtual ~FlatSourceRegion();
	int getUid() const;
    Material* getMaterial();
    FP_PRECISION getVolume() const;
    FP_PRECISION* getFlux();
    FP_PRECISION* getOldFlux();
    FP_PRECISION* getOldSource();
    FP_PRECISION* getSource();
    FP_PRECISION* getRatios();
    void setMaterial(Material* material);
    void setVolume(FP_PRECISION volume);
    void incrementVolume(FP_PRECISION volume);
    void setFlux(int energy, FP_PRECISION flux);
    void incrementFlux(int energy, FP_PRECISION flux);
    void incrementFlux(FP_PRECISION* flux);
    void setOldFlux(int energy, FP_PRECISION old_flux);
    void setSource(int energy, FP_PRECISION source);
    void setOldSource(int energy, FP_PRECISION old_source);
    void normalizeFluxes(FP_PRECISION factor);
    void computeRatios();
    FP_PRECISION computeFissionRate();
};


#endif /* FLATSOURCEREGION_H_ */
