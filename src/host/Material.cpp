/*
 * Material.cpp
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Material.h"

/* _n keeps track of the number of materials instantiated */
int Material::_n = 0;

/**
 * Material constructor
 * @param id the material's id
 */
Material::Material(short int id,
				   FP_PRECISION *sigma_a, int sigma_a_cnt,
				   FP_PRECISION *sigma_t, int sigma_t_cnt,
				   FP_PRECISION *sigma_f, int sigma_f_cnt,
				   FP_PRECISION *nu_sigma_f, int nu_sigma_f_cnt,
				   FP_PRECISION *chi, int chi_cnt,
				   FP_PRECISION *sigma_s, int sigma_s_cnt) {
	_uid = _n;
	_id = id;
	_n++;

	if (sigma_a_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_a");
	memcpy(_sigma_a, sigma_a, NUM_ENERGY_GROUPS*sizeof(*_sigma_a));

	if (sigma_t_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_t");
	memcpy(_sigma_t, sigma_t, NUM_ENERGY_GROUPS*sizeof(*_sigma_t));

	if (nu_sigma_f_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_f");
	memcpy(_sigma_f, sigma_f, NUM_ENERGY_GROUPS*sizeof(*_sigma_f));

	if (nu_sigma_f_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of nu_sigma_f");
	memcpy(_nu_sigma_f, nu_sigma_f, NUM_ENERGY_GROUPS*sizeof(*_nu_sigma_f));

	if (chi_cnt != NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of chi");
	memcpy(_chi, chi, NUM_ENERGY_GROUPS*sizeof(*_chi));

	if (sigma_s_cnt != NUM_ENERGY_GROUPS*NUM_ENERGY_GROUPS)
		log_printf(ERROR, "Wrong number of sigma_s");
	
	/* Set the material's scattering matrix. This assumes that the scattering
	 * matrix passed in has the standard notation: the ij element is for scattering
	 * from group i to j. For efficient caching of the elements of this matrix
	 * during fixed source iteration, the matrix transpose is what is actually
	 * stored in the material
     */
	for (int i=0; i<NUM_ENERGY_GROUPS; i++) {
		for (int j=0; j<NUM_ENERGY_GROUPS; j++) {
			_sigma_s[i][j] = sigma_s[j*NUM_ENERGY_GROUPS+i];
		}
	}

	return;
}

/**
 * Destructor
 */
Material::~Material() { }



/**
 * Get the material's chi array
 * @return the material's chi array
 */
FP_PRECISION* Material::getChi() {
    return _chi;
}



/**
 * Return the material's uid
 * @return the material's uid
 */
int Material::getUid() const {
	return _uid;
}

/**
 * Return the material's id
 * @return the material's id
 */
short int Material::getId() const {
    return _id;
}



/**
 * Return the material's fission cross-section array
 * @return the material's fission cross-section array
 */
FP_PRECISION* Material::getSigmaF() {
    return _sigma_f;
}


/**
 * Return the material's nu*sigma_f array
 * @return the material's nu*sigma_f array
 */
FP_PRECISION* Material::getNuSigmaF() {
    return _nu_sigma_f;
}


/**
 * Return the material's scattering matrix
 * @return the material's scattering matrix
 */
FP_PRECISION* Material::getSigmaS() {
    return *_sigma_s;
}


/**
 * Return the material's total cross-section array
 * @return the material's total cross-section array
 */
FP_PRECISION* Material::getSigmaT() {
    return _sigma_t;
}


/**
 * Return the material's absorption cross section array
 * @return the material's absorption cross-section array
 */
FP_PRECISION* Material::getSigmaA() {
	return _sigma_a;
}


/**
 * Set the material's chi array
 * @param chi the chi array
 */
void Material::setChi(FP_PRECISION chi[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_chi[i] = chi[i];
}


/**
 * Set the material's fission cross-section array
 * @param nu_sigma_f the fission cross-section array
 */
void Material::setSigmaF(FP_PRECISION sigma_f[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_sigma_f[i] = sigma_f[i];
}


/**
 * Set the material's nu*sigma_f array
 * @param nu_sigma_f the nu*sigma_f array
 */
void Material::setNuSigmaF(FP_PRECISION nu_sigma_f[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_nu_sigma_f[i] = nu_sigma_f[i];
}


/**
 * Set the material's scattering matrix. This assumes that the scattering
 * matrix passed in has the standard notation: the ij element is for scattering
 * from group i to j. For efficient caching of the elements of this matrix
 * during fixed source iteration, the matrix transpose is what is actually
 * stored in the material
 * @param sigma_s the material's scattering matrix
 */
void Material::setSigmaS(FP_PRECISION sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
			for (int j=0; j < NUM_ENERGY_GROUPS; j++)
			_sigma_s[i][j] = sigma_s[j][i];
	}
}


/**
 * Set the material's total scattering cross-section array
 * @param sigma_t the material's total scattering cross-section
 */
void Material::setSigmaT(FP_PRECISION sigma_t[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_sigma_t[i] = sigma_t[i];
}


/**
 * Set the material's absorption scattering cross-section array
 * @param sigma_a the material's absorption scattering cross-section
 */
void Material::setSigmaA(FP_PRECISION sigma_a[NUM_ENERGY_GROUPS]) {
	for (int i=0; i < NUM_ENERGY_GROUPS; i++)
		_sigma_a[i] = sigma_a[i];
}


/**
 * Checks if the total cross-section for this material is equal to the
 * absorption plus scattering cross-sections for all energy groups
 */
void Material::checkSigmaT() {

	FP_PRECISION calc_sigma_t;

	/* Loop over all energy groups */
	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {

		/* Initialize the calculated total xs to the absorption xs */
		calc_sigma_t = _sigma_a[i];

		/* Increment calculated total xs by scatter xs for each energy group */
		for (int j=0; j < NUM_ENERGY_GROUPS; j++)
			calc_sigma_t += _sigma_s[j][i];

		/* Check if the calculated and total match up to certain threshold */
		if (fabs(calc_sigma_t - _sigma_t[i]) > SIGMA_T_THRESH) {
			log_printf(ERROR, "Material id = %d has a different total "
					"cross-section than the sum of its scattering and "
					"absorption cross-sections for group %d: "
					"sigma_t = %f, calc_sigma_t = %f", _id, i, _sigma_t[i],
					   calc_sigma_t);
		}
	}

	return;
}


/**
 * Converts this material's attributes to a character array representation
 * @param a character array of this member's attributes
 */
std::string Material::toString() {
	std::stringstream string;

	string << "Material id = " << _id;

	string << "\n\t\tSigma_a = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _sigma_a[e] << ", ";

	string << "\n\t\tSigma_t = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _sigma_t[e] << ", ";

	string << "\n\t\tnu_sigma_f = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _nu_sigma_f[e] << ", ";

	string << "\n\t\tSigma_s = \n\t\t";
	for (int G = 0; G < NUM_ENERGY_GROUPS; G++) {
		for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
			string << _sigma_s[G][g] << "\t\t ";
		string << "\n\t\t";
	}

	string << "Chi = ";
	for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		string << _chi[e] << ", ";

	return string.str();
}
