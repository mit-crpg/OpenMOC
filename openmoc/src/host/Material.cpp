#include "Material.h"

int Material::_n = 0;


static int auto_id = 10000;


/**
 * @brief Returns an auto-generated unique material ID.
 * @details This method is intended as a utility method for user's writing
 *          OpenMOC input files. The method makes use of a static material
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined material IDs greater
 *          than or equal to 10000 is prohibited.
 */
int material_id() {
    int id = auto_id;
    auto_id++;
    return id;
}


/**
 * @brief Constructor sets the ID and unique ID for the material.
 * @param id the user-defined id for the material
 */
Material::Material(short int id) {

    _uid = _n;
    _id = id;
    _n++;

    _sigma_t = NULL;
    _sigma_a = NULL;
    _sigma_s = NULL;
    _sigma_f = NULL;
    _nu_sigma_f = NULL;
    _chi = NULL;

    return;
}


/**
 * @brief Destructor deletes all cross-section data.
 */
Material::~Material() { 

   /* This assumes that if the number of energy groups has been set, all 
    * cross-sections have been set */
    if (_sigma_t != NULL)
        delete [] _sigma_t;
    if (_sigma_a != NULL)
        delete [] _sigma_a;
    if (_sigma_s != NULL)
        delete [] _sigma_s;
    if (_sigma_f != NULL)
        delete [] _sigma_f;
    if (_nu_sigma_f != NULL)
        delete [] _nu_sigma_f;
    if (_chi != NULL)
        delete [] _chi;
}


/**
 * @brief Return the material's unique ID.
 * @return the material's unique ID
 */
int Material::getUid() const {
    return _uid;
}

/**
 * @brief Return the material's user-defined ID
 * @return the material's user-defined ID
 */
short int Material::getId() const {
    return _id;
}


/**
 * @brief Returns the number of energy groups for this material's nuclear data.
 * @return the number of energy groups
 */
int Material::getNumEnergyGroups() const {
    return _num_energy_groups;
}


/**
 * @brief Return the array of the material's total cross-sections.
 * @return the pointer to the material's array of total cross-sections
 */
double* Material::getSigmaT() {
    if (_sigma_t == NULL)
        log_printf(ERROR, "Unable to return material %d's total "
                        "cross-section since it has not yet been set", _id);
    return _sigma_t;
}



/**
 * @brief Return the array of the material's absorption cross-sections.
 * @return the pointer to the material's array of absorption cross-sections
 */
double* Material::getSigmaA() {
    if (_sigma_a == NULL)
        log_printf(ERROR, "Unable to return material %d's absorption "
                        "cross-section since it has not yet been set", _id);

    return _sigma_a;
}


/**
 * @brief Return the array of the material's scattering cross-section matrix.
 * @return the pointer to the material's array of scattering cross-sections
 */
double* Material::getSigmaS() {
    if (_sigma_s == NULL)
        log_printf(ERROR, "Unable to return material %d's scattering "
                        "cross-section since it has not yet been set", _id);

    return _sigma_s;
}


/**
 * @brief Return the array of the material's fission cross-sections.
 * @return the pointer to the material's array of fission cross-sections
 */
double* Material::getSigmaF() {
    if (_sigma_f == NULL)
        log_printf(ERROR, "Unable to return material %d's fission "
                        "cross-section since it has not yet been set", _id);

    return _sigma_f;
}


/**
 * @brief Return the array of the material's fission cross-sections
 *        multiplied by nu \f$ \nu \f$.
 * @return the pointer to the material's array of fission cross-sections 
 *         multiplied by nu \f$ \nu \f$
 */
double* Material::getNuSigmaF() {
    if (_nu_sigma_f == NULL)
        log_printf(ERROR, "Unable to return material %d's nu times fission "
                        "cross-section since it has not yet been set", _id);

    return _nu_sigma_f;
}


/**
 * @brief Return the array of the material's chi \f$ \chi \f$.
 * @return the pointer to the material's array of chi \f$ \chi \f$ values
 */
double* Material::getChi() {
    if (_chi == NULL)
        log_printf(ERROR, "Unable to return material %d's chi spectrum "
                        "since it has not yet been set", _id);

    return _chi;
}


/**
 * @brief Set the number of energy groups for this material.
 * @param num_energy_groups the number of energy groups.
 */
void Material::setNumEnergyGroups(const int num_energy_groups) {
    if (num_energy_groups < 0)
        log_printf(ERROR, "Unable to set the number of energy groups for "
                   "material %d to %d", num_energy_groups);

    _num_energy_groups = num_energy_groups;
}

/**
 * @brief Set the material's array of total cross-sections.
 * @param sigma_t the array of total cross-sections
 * @param num_energy_groups the number of energy groups
 */
void Material::setSigmaT(double* sigma_t, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set sigma_t with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _sigma_t = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _sigma_t[i] = sigma_t[i];
}


/**
 * @brief Set the material's array of absorption scattering cross-sections.
 * @details This method is intended to be called from 
 * @param sigma_a the array of absorption scattering cross-sections
 * @param num_energy_groups the number of energy groups
 */
void Material::setSigmaA(double* sigma_a, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set sigma_a with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _sigma_a = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _sigma_a[i] = sigma_a[i];
}



/**
 * @brief Set the material's 2D array of scattering cross-sections. 
 * @details This assumes that the scattering matrix passed in has the standard 
 *          notation: the ij element is for scattering from group i to j. For 
 *          efficient caching of the elements of this matrix during fixed 
 *          source iteration, the matrix transpose is what is actually stored 
 *          in the material
 * @param sigma_s the array of scattering cross-sections
 * @param num_energy_groups_squared the number of energy groups squared
 */
void Material::setSigmaS(double* sigma_s, int num_energy_groups_squared) {
 
    if (_num_energy_groups*_num_energy_groups != num_energy_groups_squared)
        log_printf(ERROR, "Unable to set sigma_s with %f groups for material "
		   "%d which contains %d energy groups", 
		   float(sqrt(num_energy_groups_squared)), _num_energy_groups);

    _sigma_s = new double[_num_energy_groups*_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++) {
        for (int j=0; j < _num_energy_groups; j++)
            _sigma_s[j*_num_energy_groups+i] = sigma_s[i*_num_energy_groups+j];
    }
    
}


/**
 * @brief Set the material's array of fission cross-sections.
 * @param sigma_f the array of fission cross-sections
 * @param num_energy_groups the number of energy groups
 */
void Material::setSigmaF(double* sigma_f, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _sigma_f = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _sigma_f[i] = sigma_f[i];
}


/**
 * @brief Set the material's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param nu_sigma_f the array of fission cross-sections multiplied by nu 
 *        \f$ \nu \f$
 * @param num_energy_groups the number of energy groups 
*/
void Material::setNuSigmaF(double* nu_sigma_f, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _nu_sigma_f = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _nu_sigma_f[i] = nu_sigma_f[i];
}


/**
 * @brief Set the material's array of \f$ \chi \f$ values.
 * @param chi the array of chi \f$ \chi \f$ values
 * @param num_energy_groups the number of energy groups 
 */
void Material::setChi(double* chi, int num_energy_groups) {

    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set chi with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _chi = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _chi[i] = chi[i];
}


/**
 * @brief Set the material's array of total cross-sections.
 * @param sigma_t the array of total cross-sections
 * @param num_energy_groups the number of energy groups
 */
void Material::setSigmaT(float* sigma_t, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set sigma_t with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _sigma_t = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _sigma_t[i] = sigma_t[i];
}


/**
 * @brief Set the material's array of absorption scattering cross-sections.
 * @details This method is intended to be called from 
 * @param sigma_a the array of absorption scattering cross-sections
 * @param num_energy_groups the number of energy groups
 */
void Material::setSigmaA(float* sigma_a, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set sigma_a with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _sigma_a = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _sigma_a[i] = sigma_a[i];
}



/**
 * @brief Set the material's 2D array of scattering cross-sections. 
 * @details This assumes that the scattering matrix passed in has the standard 
 *          notation: the ij element is for scattering from group i to j. For 
 *          efficient caching of the elements of this matrix during fixed 
 *          source iteration, the matrix transpose is what is actually stored 
 *          in the material
 * @param sigma_s the array of scattering cross-sections
 * @param num_energy_groups_squared the number of energy groups squared
 */
void Material::setSigmaS(float* sigma_s, int num_energy_groups_squared) {
 
    if (_num_energy_groups * _num_energy_groups != num_energy_groups_squared)
        log_printf(ERROR, "Unable to set sigma_s with %f groups for material "
		   "%d which contains %d energy groups", 
		   float(sqrt(num_energy_groups_squared)), _num_energy_groups);

    _sigma_s = new double[_num_energy_groups*_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++) {
        for (int j=0; j < _num_energy_groups; j++)
            _sigma_s[i*_num_energy_groups+j] = sigma_s[i*_num_energy_groups+j];    
    }
}


/**
 * @brief Set the material's array of fission cross-sections.
 * @param sigma_f the array of fission cross-sections
 * @param num_energy_groups the number of energy groups
 */
void Material::setSigmaF(float* sigma_f, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _sigma_f = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _sigma_f[i] = sigma_f[i];
}


/**
 * @brief Set the material's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param nu_sigma_f the array of fission cross-sections multiplied by nu 
 *        \f$ \nu \f$
 * @param num_energy_groups the number of energy groups 
*/
void Material::setNuSigmaF(float* nu_sigma_f, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _nu_sigma_f = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _nu_sigma_f[i] = nu_sigma_f[i];
}


/**
 * @brief Set the material's array of \f$ \chi \f$ values.
 * @param chi the array of chi \f$ \chi \f$ values
 * @param num_energy_groups the number of energy groups 
 */
void Material::setChi(float* chi, int num_energy_groups) {
    if (_num_energy_groups != num_energy_groups)
      log_printf(ERROR, "Unable to set chi with %d groups for material "
                 "%d which contains %d energy groups", num_energy_groups,
                 _num_energy_groups);

    _chi = new double[_num_energy_groups];

    for (int i=0; i < _num_energy_groups; i++)
        _chi[i] = chi[i];
}


/**
 * @brief Checks if the total cross-section for this material is equal to the
 *        absorption plus scattering cross-sections for all energy groups.
 * @details If the total cross-section does not equal the absorption plus 
 *          scattering cross-section within SIGMA_T_THRESH (defined in
 *          openmoc/src/host/configurations.h) then this method exits OpenMOC.
 */
void Material::checkSigmaT() {

    if (_num_energy_groups == 0)
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since the number of energy groups has not been set", _id);
    if (_sigma_t == NULL) 
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since its total cross-section has not been set", _id);
    if (_sigma_a == NULL) 
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since its absorption cross-section has not been set", _id);
    if (_sigma_s == NULL)
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since its scattering cross-section has not been set", _id);

    double calc_sigma_t;

    /* Loop over all energy groups */
    for (int i=0; i < _num_energy_groups; i++) {

        /* Initialize the calculated total xs to the absorption xs */
        calc_sigma_t = _sigma_a[i];

        /* Increment calculated total xs by scatter xs for each energy group */
        for (int j=0; j < _num_energy_groups; j++)
            calc_sigma_t += _sigma_s[i+j*_num_energy_groups];

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
 * @brief Converts this material's attributes to a character array 
 *        representation.
 * @details The character array returned includes the user-defined ID,
 *          and each of the absorption, total, fission, nu multiplied by
 *          fission and scattering cross-sections and chi for all energy
 *          groups.
 * @return character array of this member's attributes
 */
std::string Material::toString() {

    std::stringstream string;

    string << "Material id = " << _id;

    if (_sigma_a != NULL) {
        string << "\n\t\tSigma_a = ";
        for (int e = 0; e < _num_energy_groups; e++)
            string << _sigma_a[e] << ", ";
    }

    if (_sigma_t != NULL) {
        string << "\n\t\tSigma_t = ";
        for (int e = 0; e < _num_energy_groups; e++)
            string << _sigma_t[e] << ", ";
    }

    if (_sigma_f != NULL) {
        string << "\n\t\tSigma_f = ";
        for (int e = 0; e < _num_energy_groups; e++)
            string << _sigma_f[e] << ", ";
    }


    if (_nu_sigma_f != NULL) {
        string << "\n\t\tnu_sigma_f = ";
        for (int e = 0; e < _num_energy_groups; e++)
            string << _nu_sigma_f[e] << ", ";
    }

    if (_sigma_s != NULL) {
        string << "\n\t\tSigma_s = \n\t\t";
        for (int G = 0; G < _num_energy_groups; G++) {
  	    for (int g = 0; g < _num_energy_groups; g++)
                string << _sigma_s[G+g*_num_energy_groups] << "\t\t ";
            string << "\n\t\t";
        }
    }

    if (_chi != NULL) {
        string << "Chi = ";
        for (int e = 0; e < _num_energy_groups; e++)
            string << _chi[e] << ", ";
    }

    return string.str();
}


/**
 * @brief Prints a string representation of all of the material's objects to
 *        the console.
 */
void Material::printString() {
    printf("%s", toString().c_str());
}
