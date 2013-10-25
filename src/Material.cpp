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
    _dif_coef = NULL;
    _dif_hat = NULL;
    _dif_tilde = NULL;
    _buckling = NULL;

    _fissionable = false;

    _data_aligned = false;

    return;
}


/**
 * @brief Destructor deletes all cross-section data.
 */
Material::~Material() { 

    /* If data is vector aligned */
    if (_data_aligned) {
        if (_sigma_t != NULL)
	    _mm_free(_sigma_t);

	if (_sigma_a != NULL)
	    _mm_free(_sigma_a);

	if (_sigma_s != NULL)
	    _mm_free(_sigma_s);
	
	if (_sigma_f != NULL)
	    _mm_free(_sigma_f);
	
	if (_nu_sigma_f != NULL)
	    _mm_free(_nu_sigma_f);
	
	if (_chi != NULL)
	    _mm_free(_chi);

	if (_dif_coef != NULL)
	    delete [] _dif_coef;

	if (_dif_hat != NULL)
	    delete [] _dif_hat;

	if (_dif_tilde != NULL)
	    delete [] _dif_tilde;

	if (_buckling != NULL)
	    delete [] _buckling;

    }

    /* Data is not vector aligned */
    else {
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

	if (_dif_coef != NULL)
	    delete [] _dif_coef;

	if (_dif_hat != NULL)
	    delete [] _dif_hat;

	if (_dif_tilde != NULL)
	    delete [] _dif_tilde;

	if (_buckling != NULL)
	    delete [] _buckling;
    }
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
    return _num_groups;
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
 * @brief Return the array of the material's diffusion coefficients.
 * @return the pointer to the material's array of diffusion coefficients
 */
double* Material::getDifCoef() {

  if (_dif_coef == NULL){
    _dif_coef = new double[_num_groups];
    for (int e = 0; e < _num_groups; e++)
      _dif_coef[e] = 0.0;
  }

    return _dif_coef;
}


/**
 * @brief Return the array of the material's surface diffusion coefficients.
 * @return the pointer to the material's array of surface diffusion coefficients
 */
double* Material::getDifHat() {

  if (_dif_hat == NULL){
    _dif_hat = new double[4*_num_groups];
    for (int e = 0; e < 4*_num_groups; e++)
      _dif_hat[e] = 0.0;
  }


    return _dif_hat;
}

/**
 * @brief Return the array of the material's CMFD correction to the
 * surface diffusion coefficients.
 * @return the pointer to the material's array of CMFD correction to
 * the surface diffusion coefficients
 */
double* Material::getDifTilde() {

  if (_dif_tilde == NULL){
    _dif_tilde = new double[4*_num_groups];
    for (int e = 0; e < 4*_num_groups; e++)
      _dif_tilde[e] = 0.0;
  }


    return _dif_tilde;
}


/**
 * @brief Return the array of the material's CMFD correction to the
 * surface diffusion coefficients.
 * @return the pointer to the material's array of CMFD correction to
 * the surface diffusion coefficients
 */
double* Material::getBuckling() {

  if (_buckling == NULL){
    _buckling = new double[_num_groups];
    for (int e = 0; e < _num_groups; e++)
      _buckling[e] = 0.0;
  }

    return _buckling;
}


/**
 * @brief Returns whether or not the material contains a fissionable (non-zero)
 *        fission cross-section.
 * @return true if fissionable, false otherwise
 */
bool Material::isFissionable() {
    return _fissionable;
}


/**
 * @brief Returns true if the data is vector aligned, false otherwise (default).
 * @return Whether or not the materials data is vector aligned
 */
bool Material::isDataAligned() {
    return _data_aligned;
}


/**
 * @brief Returns the rounded up number of energy groups to fill an integral
 *        number of vector lengths
 * @return The number of vector-aligned energy groups
 */
int Material::getNumVectorGroups() {
    return _num_vector_groups;
}



/**
 * @brief Set the number of energy groups for this material.
 * @param num_groups the number of energy groups.
 */
void Material::setNumEnergyGroups(const int num_groups) {

    if (num_groups < 0)
        log_printf(ERROR, "Unable to set the number of energy groups for "
                   "material %d to %d", _num_groups);

    _num_groups = num_groups;

    /* Free old memory arrays if they have already been allocated 
     * for a previous simulation */
    /* If data is vector aligned */
    if (_data_aligned) {
        if (_sigma_t != NULL)
	    _mm_free(_sigma_t);

	if (_sigma_a != NULL)
	    _mm_free(_sigma_a);

	if (_sigma_s != NULL)
	    _mm_free(_sigma_s);
	
	if (_sigma_f != NULL)
	    _mm_free(_sigma_f);
	
	if (_nu_sigma_f != NULL)
	    _mm_free(_nu_sigma_f);
	
	if (_chi != NULL)
	    _mm_free(_chi);
    }

    /* Data is not vector aligned */
    else {
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

	if (_dif_coef != NULL)
	    delete [] _dif_coef;

	if (_dif_hat != NULL)
	    delete [] _dif_hat;

	if (_dif_tilde != NULL)
	    delete [] _dif_tilde;

	if (_buckling != NULL)
	    delete [] _buckling;
    }


    _sigma_t = new double[_num_groups];
    _sigma_a = new double[_num_groups];
    _sigma_f = new double[_num_groups];
    _nu_sigma_f = new double[_num_groups];
    _chi = new double[_num_groups];
    _sigma_s = new double[_num_groups*_num_groups];
}


/**
 * @brief Set the material's array of total cross-sections.
 * @param xs the array of total cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaT(double* xs, int num_groups) {
    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set sigma_t with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++)
        _sigma_t[i] = xs[i];
}


/**
 * @brief Set the material's total cross-section for some energy group.
 * @param xs the total cross-section (\f$ \Sigma_t [cm^1] \f$)
 * @param group the energy group
 */
void Material::setSigmaTByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_t for group %d for material "
                 "%d which contains %d energy groups", group, _uid, _num_groups);

    _sigma_t[group] = xs;
}


/**
 * @brief Set the material's array of absorption scattering cross-sections.
 * @details This method is intended to be called from 
 * @param xs the array of absorption scattering cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaA(double* xs, int num_groups) {
    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set sigma_a with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++){
        _sigma_a[i] = xs[i];

	if (_buckling != NULL & _dif_coef != NULL)
	  _sigma_a[i] += _buckling[i] * _dif_coef[i];
    }
}


/**
 * @brief Set the material's absorption cross-section for some energy group.
 * @param xs the absorption cross-section (\f$ \Sigma_a [cm^1] \f$)
 * @param group the energy group
 */
void Material::setSigmaAByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_a for group %d for material "
                 "%d which contains %d energy groups", group, _uid, _num_groups);

    _sigma_a[group] = xs;
}




/**
 * @brief Set the material's 2D array of scattering cross-sections. 
 * @details This assumes that the scattering matrix passed in has the standard 
 *          notation: the ij element is for scattering from group i to j. For 
 *          efficient caching of the elements of this matrix during fixed 
 *          source iteration, the matrix transpose is what is actually stored 
 *          in the material
 * @param xs the array of scattering cross-sections
 * @param num_groups_squared the number of energy groups squared
 */
void Material::setSigmaS(double* xs, int num_groups_squared) {
 
    if (_num_groups*_num_groups != num_groups_squared)
        log_printf(ERROR, "Unable to set sigma_s with %f groups for material "
		   "%d which contains %d energy groups", 
		   float(sqrt(num_groups_squared)), _num_groups);

    for (int i=0; i < _num_groups; i++) {
        for (int j=0; j < _num_groups; j++)
            _sigma_s[j*_num_groups+i] = xs[i*_num_groups+j];
    }
    
}


/**
 * @brief Set the material's scattering cross-section for some energy group.
 * @param xs the scattering cross-section (\f$ \Sigma_s [cm^1] \f$)
 * @param group1 the row index in the scattering matrix
 * @param group2 the column index in the scattering matrix
 */
void Material::setSigmaSByGroup(double xs, int group1, int group2) {

    if (group1 < 0 || group2 < 0 || group1 >= _num_groups || group2 >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_s for group %d,%d for material "
                 "%d which contains %d energy groups", group1, group2, _uid,
                 _num_groups);

    _sigma_s[_num_groups*group1 + group2] = xs;
}



/**
 * @brief Set the material's array of fission cross-sections.
 * @param xs the array of fission cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaF(double* xs, int num_groups) {
    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++)
        _sigma_f[i] = xs[i];

    /* Determine whether or not this material is fissionable */
    _fissionable = false;

    for (int i=0; i < _num_groups; i++) {
        if (_sigma_f[i] > 0.0) {
	    _fissionable = true;
	    return;
	}
    }
}


/**
 * @brief Set the material's fission cross-section for some energy group.
 * @param xs the fission cross-section (\f$ \Sigma_f [cm^1] \f$)
 * @param group the energy group
 */
void Material::setSigmaFByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_f for group %d for material "
                 "%d which contains %d energy groups", group, _uid, _num_groups);

    _sigma_f[group] = xs;

    /* Determine whether or not this material is fissionable */
    _fissionable = false;

    for (int i=0; i < _num_groups; i++) {
        if (_sigma_f[i] > 0.0) {
	    _fissionable = true;
	    return;
	}
    }
}


/**
 * @brief Set the material's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param xs the array of fission cross-sections multiplied by nu 
 *        \f$ \nu \f$
 * @param num_groups the number of energy groups 
*/
void Material::setNuSigmaF(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_groups, _uid, _num_groups);

    for (int i=0; i < _num_groups; i++)
        _nu_sigma_f[i] = xs[i];
}


/**
 * @brief Set the material's fission cross-section multiplied by \f$ \nu \f$
 *        for some energy group.
 * @param xs the fission cross-section (\f$ \nu\Sigma_f [cm^1] \f$)
 * @param group the energy group
 */
void Material::setNuSigmaFByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set nu_sigma_f for group %d for material "
                 "%d which contains %d energy groups", group, _uid);

    _nu_sigma_f[group] = xs;
}



/**
 * @brief Set the material's array of \f$ \chi \f$ values.
 * @param xs the array of chi \f$ \chi \f$ values
 * @param num_groups the number of energy groups 
 */
void Material::setChi(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set chi with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++)
        _chi[i] = xs[i];
}


/**
 * @brief Set the material's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
void Material::setChiByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set chi for group %d for material "
		   "%d which contains %d energy groups", 
		   group, _num_groups, _uid);

    _chi[group] = xs;
}


/**
 * @brief Set the material's array of diffusion coefficients.
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setDifCoef(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
    		  "for material %d which contains %d energy groups", num_groups,
                 _num_groups);

    if (_dif_coef == NULL)
      _dif_coef = new double[_num_groups];

    for (int i=0; i < _num_groups; i++)
        _dif_coef[i] = xs[i];
}



/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDifCoefByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    if (_dif_coef == NULL){
      _dif_coef = new double[_num_groups];
      for (int i=0; i < _num_groups; i++)
        _dif_coef[i] = 0.0;
    }

    _dif_coef[group] = xs;
}


/**
 * @brief Set the material's array of diffusion coefficients.
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setBuckling(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
    		  "for material %d which contains %d energy groups", num_groups,
                 _num_groups);

    if (_buckling == NULL){
      _buckling = new double[_num_groups];
    }

    for (int i=0; i < _num_groups; i++)
        _buckling[i] = xs[i];
}


/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setBucklingByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    if (_buckling == NULL){
      _buckling = new double[_num_groups];
      for (int i=0; i < _num_groups; i++)
        _buckling[i] = 0.0;
    }

    _buckling[group] = xs;
}


/**
 * @brief Set the material's array of diffusion coefficients.
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setDifHat(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
    		  "for material %d which contains %d energy groups", num_groups,
                 _num_groups);

    if (_dif_hat == NULL)
      _dif_hat = new double[4*_num_groups];

    for (int i=0; i < _num_groups*4; i++)
        _dif_hat[i] = xs[i];
}


/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDifHatByGroup(double xs, int group, int surface) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    if (_dif_hat == NULL){
      _dif_hat = new double[4*_num_groups];
      for (int i=0; i < _num_groups*4; i++)
        _dif_hat[i] = 0.0;
    }

    _dif_hat[surface*_num_groups + group] = xs;
}


/**
 * @brief Set the material's array of diffusion coefficients.
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setDifTilde(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
    		  "for material %d which contains %d energy groups", num_groups,
                 _num_groups);

    if (_dif_tilde == NULL){
      _dif_tilde = new double[4*_num_groups];
    }

    for (int i=0; i < _num_groups*4; i++)
        _dif_tilde[i] = xs[i];
}


/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDifTildeByGroup(double xs, int group, int surface) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient correction for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    if (_dif_tilde == NULL){
      _dif_tilde = new double[4*_num_groups];
      for (int i=0; i < _num_groups*4; i++)
        _dif_tilde[i] = 0.0;
    }

    _dif_tilde[surface*_num_groups + group] = xs;
}


/**
 * @brief Checks if the total cross-section for this material is equal to the
 *        absorption plus scattering cross-sections for all energy groups.
 * @details If the total cross-section does not equal the absorption plus 
 *          scattering cross-section within SIGMA_T_THRESH (defined in
 *          openmoc/src/host/configurations.h) then this method exits OpenMOC.
 */
void Material::checkSigmaT() {

    if (_num_groups == 0)
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
    for (int i=0; i < _num_groups; i++) {

        /* Initialize the calculated total xs to the absorption xs */
        calc_sigma_t = _sigma_a[i];

        /* Increment calculated total xs by scatter xs for each energy group */
        for (int j=0; j < _num_groups; j++)
            calc_sigma_t += _sigma_s[i+j*_num_groups];

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
        for (int e = 0; e < _num_groups; e++)
            string << _sigma_a[e] << ", ";
    }

    if (_sigma_t != NULL) {
        string << "\n\t\tSigma_t = ";
        for (int e = 0; e < _num_groups; e++)
            string << _sigma_t[e] << ", ";
    }

    if (_sigma_f != NULL) {
        string << "\n\t\tSigma_f = ";
        for (int e = 0; e < _num_groups; e++)
            string << _sigma_f[e] << ", ";
    }


    if (_nu_sigma_f != NULL) {
        string << "\n\t\tnu_sigma_f = ";
        for (int e = 0; e < _num_groups; e++)
            string << _nu_sigma_f[e] << ", ";
    }

    if (_sigma_s != NULL) {
        string << "\n\t\tSigma_s = \n\t\t";
        for (int G = 0; G < _num_groups; G++) {
  	    for (int g = 0; g < _num_groups; g++)
                string << _sigma_s[G+g*_num_groups] << "\t\t ";
            string << "\n\t\t";
        }
    }

    if (_chi != NULL) {
        string << "Chi = ";
        for (int e = 0; e < _num_groups; e++)
            string << _chi[e] << ", ";
    }

    if (_dif_coef != NULL) {
        string << "Diffusion Coefficient = ";
        for (int e = 0; e < _num_groups; e++)
            string << _dif_coef[e] << ", ";
    }

    if (_buckling != NULL) {
        string << "Buckling = ";
        for (int e = 0; e < _num_groups; e++)
            string << _buckling[e] << ", ";
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


/**
 * @brief Aligns the cross-section data structures 
 */
void Material::alignData() {

    if (_data_aligned)
        return;

    if (_num_groups <= 0)
        log_printf(ERROR, "Unable to align material %d data since the "
		   "cross-sections have not yet been set\n", _uid);

    _num_vector_groups = (_num_groups / VEC_LENGTH) + 1;

    /* Allocate memory for the new aligned xs data */
    int size = _num_vector_groups * VEC_LENGTH * sizeof(double);

    double* new_sigma_t = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_sigma_a = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_sigma_f = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_nu_sigma_f = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_chi = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    
    /* The scattering matrix will be the number of vector groups 
     * wide (SIMD) and the actual number of groups long since 
     * instructions are not SIMD in this dimension */
    size = _num_vector_groups * VEC_LENGTH * _num_vector_groups * VEC_LENGTH * sizeof(double);
    double* new_sigma_s = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    
    /* Initialize data structures to ones for sigma_t since it is used to
     * divide the source in the solver, and zeroes for everything else */
    size = _num_vector_groups * VEC_LENGTH * sizeof(double);
    for (int i=0; i < _num_vector_groups * VEC_LENGTH; i++) {
        new_sigma_t[i] = 1.0;
	new_sigma_a[i] = 0.0;
	new_sigma_f[i] = 0.0;
	new_nu_sigma_f[i] = 0.0;
	new_chi[i] = 0.0;
    }
    
    size *= _num_vector_groups * VEC_LENGTH;
    memset(new_sigma_s, 0.0, size);
    
    /* Copy materials data from unaligned arrays into new aligned arrays */
    size = _num_groups * sizeof(double);
    memcpy(new_sigma_t, _sigma_t, size);
    memcpy(new_sigma_a, _sigma_a, size);
    memcpy(new_sigma_f, _sigma_f, size);
    memcpy(new_nu_sigma_f, _nu_sigma_f, size);
    memcpy(new_chi, _chi, size);

    for (int e=0; e < _num_groups; e++) {
        memcpy(new_sigma_s, _sigma_s, size);
	new_sigma_s += _num_vector_groups * VEC_LENGTH;
	_sigma_s += _num_groups;
    }

    _sigma_s -= _num_groups * _num_groups;
    
    /* Reset the new scattering cross section array pointer */
    new_sigma_s -= _num_vector_groups * VEC_LENGTH * _num_groups;

    /* Delete the old unaligned arrays */
    delete [] _sigma_t;
    delete [] _sigma_a;
    delete [] _sigma_f;
    delete [] _nu_sigma_f;
    delete [] _chi;
    delete [] _sigma_s;
    	
    /* Set the material's array pointers to the new aligned arrays */
    _sigma_t = new_sigma_t;
    _sigma_a = new_sigma_a;
    _sigma_f = new_sigma_f;
    _nu_sigma_f = new_nu_sigma_f;
    _chi = new_chi;
    _sigma_s = new_sigma_s;
    
    _data_aligned = true;

    return;
}


Material* Material::clone(){

  Material* material_clone = new Material(getId());
  
  material_clone->setNumEnergyGroups(_num_groups);
  material_clone->setSigmaT(_sigma_t, _num_groups);
  material_clone->setSigmaF(_sigma_f, _num_groups);
  material_clone->setNuSigmaF(_nu_sigma_f, _num_groups);
  material_clone->setChi(_chi, _num_groups);
  
  if (_dif_coef != NULL)
    material_clone->setDifCoef(_dif_coef, _num_groups);  
  
  if (_buckling != NULL)
    material_clone->setBuckling(_buckling, _num_groups);  

  if (_dif_hat != NULL)
    material_clone->setDifHat(_dif_hat, _num_groups);  

  if (_dif_tilde != NULL)
    material_clone->setDifTilde(_dif_tilde, _num_groups);  

  material_clone->setSigmaA(_sigma_a, _num_groups);

  copySigmaS(material_clone);

  return material_clone;

}


void Material::copySigmaS(Material* material){

  for (int i=0; i < _num_groups; i++){
    for (int j=0; j < _num_groups; j++){
      material->setSigmaSByGroup(_sigma_s[i*_num_groups+j], i, j);
    }
  }
}

