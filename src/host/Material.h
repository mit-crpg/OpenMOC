/*
 * Material.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "configurations.h"
#include "log.h"


class Material {
private:
	static int _n; /* Counts the number of materials */
	int _uid;      /* monotonically increasing id based on n */
	short int _id;
	FP_PRECISION _sigma_t[NUM_ENERGY_GROUPS];
	FP_PRECISION _sigma_a[NUM_ENERGY_GROUPS];
	FP_PRECISION _sigma_f[NUM_ENERGY_GROUPS];
	FP_PRECISION _nu_sigma_f[NUM_ENERGY_GROUPS];
	FP_PRECISION _chi[NUM_ENERGY_GROUPS];
	/* first index is row number; second index is column number */
	FP_PRECISION _sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]; 

public:
	Material(short int id,
			 FP_PRECISION *sigma_a, int sigma_a_cnt,
			 FP_PRECISION *sigma_t, int sigma_t_cnt,
			 FP_PRECISION *sigma_f, int sigma_f_cnt,
			 FP_PRECISION *nu_sigma_f, int nu_sigma_f_cnt,
			 FP_PRECISION *chi, int chi_cnt,
			 FP_PRECISION *sigma_s, int sigma_s_cnt);
	virtual ~Material();

	int getUid() const;
	short int getId() const;
	FP_PRECISION* getSigmaF();
	FP_PRECISION* getNuSigmaF();
	FP_PRECISION* getSigmaS();
	FP_PRECISION* getSigmaT();
	FP_PRECISION* getSigmaA();
	FP_PRECISION* getChi();

	void setChi(FP_PRECISION chi[NUM_ENERGY_GROUPS]);
	void setSigmaF(FP_PRECISION sigma_f[NUM_ENERGY_GROUPS]);
	void setNuSigmaF(FP_PRECISION nu_sigma_f[NUM_ENERGY_GROUPS]);
	void setSigmaS(FP_PRECISION sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS]);
	void setSigmaT(FP_PRECISION sigma_t[NUM_ENERGY_GROUPS]);
	void setSigmaA(FP_PRECISION sigma_a[NUM_ENERGY_GROUPS]);

	void checkSigmaT();
	std::string toString();
};

#endif /* MATERIAL_H_ */
