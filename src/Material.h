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
    double* _sigma_t;

    /** An array of the absorption cross-sections for each energy group */
    double* _sigma_a;

    /** A 2D array of the scattering cross-section matrix. The first index is 
     *  row number and second index is column number */
    double* _sigma_s; 

    /** An array of the fission cross-sections for each energy group */
    double* _sigma_f;

    /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$ 
     *  for each energy group */
    double* _nu_sigma_f;

    /** An array of the chi \f$ \chi \f$ values for each energy group */
    double* _chi;

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
    double* getSigmaT();
    double* getSigmaA(); 
    double* getSigmaS();
    double* getSigmaF();
    double* getNuSigmaF();
    double* getChi();
    bool isDataAligned();
    int getNumVectorGroups();

    void setNumEnergyGroups(const int num_groups);
    void setSigmaT(double* xs, int num_groups);
    void setSigmaA(double* xs, int num_groups);
    void setSigmaS(double* xs, int num_groups);
    void setSigmaF(double* xs, int num_groups);
    void setNuSigmaF(double* xs, int num_groups);
    void setChi(double* xs, int num_groups);
    
    void setSigmaTByGroup(double xs, int group);
    void setSigmaAByGroup(double xs, int group);
    void setSigmaFByGroup(double xs, int group);
    void setNuSigmaFByGroup(double xs, int group);
    void setSigmaSByGroup(double xs, int group1, int group2);
    void setChiByGroup(double xs, int group);


    void checkSigmaT();
    std::string toString();
    void printString();

    void alignData();
};


#endif /* MATERIAL_H_ */
