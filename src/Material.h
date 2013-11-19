/**
 * @file Material.h
 * @brief
 * @date January 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#ifdef __cplusplus
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "log.h"
#endif

#ifndef ICPC
#define _mm_free(array) free(array)
#define _mm_malloc(size,alignment) malloc(size)
#endif

/** Error threshold for determining how close the sum of \f$ \Sigma_a \f$ 
 *  and \f$ \Sigma_s \f$ must match that of \f$ \Sigma_t \f$ for each energy 
 *  group 
 */
#define SIGMA_T_THRESH 1E-3


int material_id();


/**
 * @class Material Material.h "openmoc/src/host/Material.h"
 * @brief The material class represents a unique material and its relevant
 *        nuclear data (ie, multigroup cross-sections) for neutron transport.
 */
class Material {

private:

    /** A static counter for the number of materials in a simulation */
    static int _n;

    /** A monotonically increasing unique ID for each material created */
    int _uid;

    /** A user-defined ID for each material created */
    short int _id;

    /** The number of energy groups */
    int _num_groups;

    /** An array of the total cross-sections for each energy group */
    FP_PRECISION* _sigma_t;

    /** An array of the absorption cross-sections for each energy group */
    FP_PRECISION* _sigma_a;

    /** A 2D array of the scattering cross-section matrix. The first index is 
     *  row number and second index is column number */
    FP_PRECISION* _sigma_s; 

    /** An array of the fission cross-sections for each energy group */
    FP_PRECISION* _sigma_f;

    /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$ 
     *  for each energy group */
    FP_PRECISION* _nu_sigma_f;

    /** An array of the chi \f$ \chi \f$ values for each energy group */
    FP_PRECISION* _chi;

    /** An array of the diffusion coefficients values for each energy group */
    FP_PRECISION* _dif_coef;

    /** An array of the diffusion coefficients values for each energy group */
    FP_PRECISION* _buckling;

    /** An array of the diffusion coefficient values
     *  for each energy group for each surface */
    FP_PRECISION* _dif_hat;

    /** An array of the CMFD correction to the diffusion coefficient values
     *  for each energy group for each surface */
    FP_PRECISION* _dif_tilde;

    /** A boolean representing whether or not this material contains a non-zero
     *  fission cross-section and is fissionable */
    bool _fissionable;

    /** A boolean to indicate whether or not the data has been 
     * allocated to be vector aligned for SIMD instructions */
    bool _data_aligned;

    /** The number of vector widths needed to fit all energy groups */
    int _num_vector_groups;

public:
    Material(short int id);
    virtual ~Material();

    int getUid() const;
    short int getId() const;
    int getNumEnergyGroups() const;
    FP_PRECISION* getSigmaT();
    FP_PRECISION* getSigmaA(); 
    FP_PRECISION* getSigmaS();
    FP_PRECISION* getSigmaF();
    FP_PRECISION* getNuSigmaF();
    FP_PRECISION* getChi();
    FP_PRECISION* getDifCoef();
    FP_PRECISION* getBuckling();
    FP_PRECISION* getDifHat();
    FP_PRECISION* getDifTilde();
    bool isFissionable();
    bool isDataAligned();
    int getNumVectorGroups();

    void setNumEnergyGroups(const int num_groups);

    void setSigmaT(double* xs, int num_groups);
    void setSigmaA(double* xs, int num_groups);
    void setSigmaS(double* xs, int num_groups);
    void setSigmaF(double* xs, int num_groups);
    void setNuSigmaF(double* xs, int num_groups);
    void setChi(double* xs, int num_groups);
    void setBuckling(double* xs, int num_groups);
    void setDifCoef(double* xs, int num_groups);
    void setDifHat(double* xs, int num_groups);
    void setDifTilde(double* xs, int num_groups);
    
    void setSigmaTByGroup(double xs, int group);
    void setSigmaAByGroup(double xs, int group);
    void setSigmaFByGroup(double xs, int group);
    void setNuSigmaFByGroup(double xs, int group);
    void setSigmaSByGroup(double xs, int group1, int group2);
    void setChiByGroup(double xs, int group);
    void setBucklingByGroup(double xs, int group);
    void setDifCoefByGroup(double xs, int group);
    void setDifHatByGroup(double xs, int group, int surface);
    void setDifTildeByGroup(double xs, int group, int surface);


    void checkSigmaT();
    std::string toString();
    void printString();

    void alignData();
    
    Material* clone();
};


#endif /* MATERIAL_H_ */
