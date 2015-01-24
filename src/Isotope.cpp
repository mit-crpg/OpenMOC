#include "Isotope.h"

int Isotope::_n = 0;


/**
 * @brief Constructor sets the ID, unique ID, and name for the Isotope.
 * @param id the user-specified optional Isotope ID
 * @param name the user-specified name for the Isotope
 */
Isotope::Isotope(int id, const char* name) {
  _uid = _n;
  _id = id;
  _n++;
  
  _name = NULL;
  setName(name);
  
  _sigma_t = NULL;
  _sigma_a = NULL;
  _sigma_s = NULL;
  _sigma_f = NULL;
  _nu_sigma_f = NULL;
  _chi = NULL;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data structures from memory.
 */
Isotope::~Isotope() {
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
 * @brief Return the isotope's unique ID.
 * @return the isotope's unique ID
 */
int Isotope::getUid() const {
  return _uid;
}


/**
 * @brief Return the isotope's user-defined ID
 * @return the isotope's user-defined ID
 */
int Isotope::getId() const {
  return _id;
}


/**
 * @brief Return the user-defined name of the Isotope
 * @return the Isotope name
 */
char* Isotope::getName() const {
  return _name;
}


/**
 * @brief Sets the name of the Isotope
 * @param name the Isotope name string
 */
void Isotope::setName(const char* name) {
  int length = strlen(name);

  if (_name != NULL)
    delete [] _name;

  /* Initialize a character array for the Isotope's name */
  _name = new char[length+1];

  /* Copy the input character array Isotope name to the class attribute name */
  for (int i=0; i <= length; i++)
    _name[i] = name[i];
}


/**
 * @brief Returns the number of energy groups for this isotope's nuclear data.
 * @return the number of energy groups
 */
int Isotope::getNumEnergyGroups() const {
  return _num_groups;
}


/**
 * @brief Return the array of the isotope's total cross-sections.
 * @return the pointer to the isotope's array of total cross-sections
 */
FP_PRECISION* Isotope::getSigmaT() {
  if (_sigma_t == NULL)
    log_printf(ERROR, "Unable to return isotope %d's total "
                      "cross-section since it has not yet been set", _id);

  return _sigma_t;
}

/**
 * @brief Return the value of the isotope's total cross section in given group
 * @return The value of the total cross section
 * @param group The index of the group 
 */
FP_PRECISION Isotope::getSigmaTByGroup(int group) {
  if (_sigma_t == NULL)
    log_printf(ERROR, "Unable to return isotope %d's total "
                      "cross-section since it has not yet been set", _uid);
                      
  if (group < 1 || group > _num_groups)
    log_printf(ERROR, "Unable to return isotope %d's total cross section "
                      " for group %d, since it only has %d groups.",
                      _uid, group, _num_groups);
  
  return _sigma_t[group-1];
}

/**
 * @brief Return the array of the isotope's absorption cross-sections.
 * @return the pointer to the isotope's array of absorption cross-sections
 */
FP_PRECISION* Isotope::getSigmaA() {
  if (_sigma_a == NULL)
    log_printf(ERROR, "Unable to return isotope %d's absorption "
                      "cross-section since it has not yet been set", _uid);

  return _sigma_a;
}

/**
 * @brief Return the value of the isotope's absorption cross section in given group
 * @return The value of the absorption cross section
 * @param group The index of the group 
 */
FP_PRECISION Isotope::getSigmaAByGroup(int group) {    
  if (_sigma_a == NULL)
    log_printf(ERROR, "Unable to return isotope %d's absorption "
                      "cross-section since it has not yet been set", _uid);
                      
  if (group < 1 || group > _num_groups)
    log_printf(ERROR, "Unable to return isotope %d's absorption cross section "
                      " for group %d, since it only has %d groups.",
                      _uid, group, _num_groups);
  
  return _sigma_a[group-1];
}

/**
 * @brief Return the array of the isotope's scattering cross-sections.
 * @return the pointer to the isotope's array of scattering cross-sections
 */
FP_PRECISION* Isotope::getSigmaS() {
  if (_sigma_s == NULL)
    log_printf(ERROR, "Unable to return isotope %d's scattering "
                      "cross-section since it has not yet been set", _uid);

  return _sigma_s;
}


/**
 * @brief Return the value of the isotope's scattering cross section in given group
 * @return The value of the scattering cross section
 * @param group1 The index of the origin group
 * @param group2 The index of the destination group
 */
FP_PRECISION Isotope::getSigmaSByGroup(int origin, int destination) {    
  if (_sigma_t == NULL)
    log_printf(ERROR, "Unable to return isotope %d's scattering "
                      "cross-section since it has not yet been set", _uid);
                      
  if (origin < 1 || origin > _num_groups || destination < 1 
                 || destination > _num_groups)
    log_printf(ERROR, "Unable to return isotope %d's scattering cross section "
                      " for group %d-%d, since it only has %d groups.",
                      _uid, origin, destination, _num_groups);  
  
  return _sigma_s[_num_groups*(origin-1) + (destination-1)];
}



/**
 * @brief Return the array of the isotope's fission cross-sections.
 * @return the pointer to the isotope's array of fission cross-sections
 */
FP_PRECISION* Isotope::getSigmaF() {
  if (_sigma_f == NULL)
    log_printf(ERROR, "Unable to return isotope %d's fission "
                      "cross-section since it has not yet been set", _uid);

  return _sigma_f;
}

/**
 * @brief Return the value of the isotope's fission cross section in given group
 * @return The value of the fission cross section
 * @param group The index of the group 
 */
FP_PRECISION Isotope::getSigmaFByGroup(int group) {    
  if (_sigma_f == NULL)
    log_printf(ERROR, "Unable to return isotope %d's fission "
                      "cross-section since it has not yet been set", _uid);
                      
  if (group < 1 || group > _num_groups)
    log_printf(ERROR, "Unable to return isotope %d's fission cross section "
                      " for group %d, since it only has %d groups.",
                      _uid, group, _num_groups);
  
  return _sigma_f[group-1];
}    

/**
 * @brief Return the array of the isotope's nu-fission cross-sections.
 * @return the pointer to the isotope's array of nu-fission cross-sections
 */
FP_PRECISION* Isotope::getNuSigmaF() {
  if (_nu_sigma_f == NULL)
    log_printf(ERROR, "Unable to return isotope %d's nu-fission "
                      "cross-section since it has not yet been set", _uid);

  return _nu_sigma_f;
}

/**
 * @brief Return the value of the isotope's nu-fission cross section in given group
 * @return The value of the nu-fission cross section
 * @param group The index of the group 
 */
FP_PRECISION Isotope::getNuSigmaFByGroup(int group) { 
  if (_nu_sigma_f == NULL)
    log_printf(ERROR, "Unable to return isotope %d's nu-fission "
                      "cross-section since it has not yet been set", _uid);
                      
  if (group < 1 || group > _num_groups)
    log_printf(ERROR, "Unable to return isotope %d's nu-fission cross section "
                      " for group %d, since it only has %d groups.",
                      _uid, group, _num_groups);
     
  return _nu_sigma_f[group-1];
}    

/**
 * @brief Return the array of the isotope's fission spectrum.
 * @return the pointer to the isotope's array of fission spectrum
 */
FP_PRECISION* Isotope::getChi() {
  if (_chi == NULL)
    log_printf(ERROR, "Unable to return isotope %d's fission "
                      "spectrum since it has not yet been set", _uid);

  return _chi;
}

/**
 * @brief Return the value of the isotope's fission spectrum in given group
 * @return The value of the fission spectrum
 * @param group The index of the group 
 */
FP_PRECISION Isotope::getChiByGroup(int group) {   
  if (_chi == NULL)
    log_printf(ERROR, "Unable to return isotope %d's fission "
                      "spectrum since it has not yet been set", _uid);
                      
  if (group < 1 || group > _num_groups)
    log_printf(ERROR, "Unable to return isotope %d's fission spectrum "
                      " for group %d, since it only has %d groups.",
                      _uid, group, _num_groups);
       
  return _chi[group-1];
}  


/**
 * @brief Set the number of energy groups for this isotope.
 * @param num_groups the number of energy groups.
 */
void Isotope::setNumEnergyGroups(const int num_groups) {

  if (num_groups < 0)
    log_printf(ERROR, "Unable to set the number of energy groups for "
                 "isotope %d to %d", _id, _num_groups);

  _num_groups = num_groups;

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


  _sigma_t = new FP_PRECISION[_num_groups];
  _sigma_a = new FP_PRECISION[_num_groups];
  _sigma_s = new FP_PRECISION[_num_groups*_num_groups];
  _sigma_f = new FP_PRECISION[_num_groups];
  _nu_sigma_f = new FP_PRECISION[_num_groups];
  _chi = new FP_PRECISION[_num_groups];
  

  for (int g = 0; g < _num_groups; g++){
    _chi[g] = 0.0;
    _nu_sigma_f[g] = 0.0;
  }
}


/**
 * @brief Set the isotope's array of total cross-sections.
 * @param xs the array of total cross-sections
 * @param num_groups the number of energy groups
 */
void Isotope::setSigmaT(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_t with %d groups for isotope "
       "%d which contains %d energy groups", num_groups, _uid,
       _num_groups);

  

  for (int i=0; i < _num_groups; i++) 
    _sigma_t[i] = FP_PRECISION(xs[i]);
}


/**
 * @brief Set the isotope's total cross-section for some energy group.
 * @param xs the total cross-section (\f$ \Sigma_t [cm^1] \f$)
 * @param group the energy group
 */
void Isotope::setSigmaTByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set sigma_t for group %d for isotope "
              "%d which contains %d energy groups", group, _uid, _num_groups);

  _sigma_t[group-1] = xs;
}

/**
 * @brief Set the isotope's array of absorption scattering cross-sections.
 * @details This method is intended to be called from Python
 * @param xs the array of absorption scattering cross-sections
 * @param num_groups the number of energy groups
 */
void Isotope::setSigmaA(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_a with %d groups for isotope "
               "%d which contains %d energy groups", num_groups, _uid,
               _num_groups);

  for (int i=0; i < _num_groups; i++){
    _sigma_a[i] = FP_PRECISION(xs[i]);
  }
}


/**
 * @brief Set the isotope's absorption cross-section for some energy group.
 * @param xs the absorption cross-section (\f$ \Sigma_a [cm^1] \f$)
 * @param group the energy group
 */
void Isotope::setSigmaAByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set sigma_a for group %d for isotope "
              "%d which contains %d energy groups", group, _uid, _num_groups);

  _sigma_a[group-1] = xs;
}




/**
 * @brief Set the isotope's 2D array of scattering cross-sections. 
 * @details This assumes that the scattering matrix passed in has the standard 
 *          notation: the ij element is for scattering from group i to j. For 
 *          efficient caching of the elements of this matrix during fixed 
 *          source iteration, the matrix transpose is what is actually stored 
 *          in the isotope
 * @param xs the array of scattering cross-sections
 * @param num_groups_squared the number of energy groups squared
 */
void Isotope::setSigmaS(double* xs, int num_groups_squared) {
 
  if (_num_groups*_num_groups != num_groups_squared)
    log_printf(ERROR, "Unable to set sigma_s with %f groups squared for "
         "matrix isotope %d which contains %d energy groups", 
         float(sqrt(num_groups_squared)), _uid, _num_groups);

  for (int i=0; i < _num_groups; i++) {
    for (int j=0; j < _num_groups; j++)
          _sigma_s[j*_num_groups+i] = xs[i*_num_groups+j];
  }
    
}



/**
 * @brief Set the isotope's scattering cross-section for some energy group.
 * @param xs the scattering cross-section (\f$ \Sigma_s [cm^1] \f$)
 * @param origin the column index in the scattering matrix
 * @param destination the row index in the scattering matrix
 */
void Isotope::setSigmaSByGroup(double xs, int origin, int destination) {

  if (origin <= 0 || destination <= 0 || origin > _num_groups 
                  || destination > _num_groups)
    log_printf(ERROR, "Unable to set sigma_s for group %d,%d for isotope "
                 "%d which contains %d energy groups", origin, destination, 
                 _uid, _num_groups);

  _sigma_s[_num_groups*(destination-1) + (origin-1)] = xs;
}



/**
 * @brief Set the isotope's array of fission cross-sections.
 * @param xs the array of fission cross-sections
 * @param num_groups the number of energy groups
 */
void Isotope::setSigmaF(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_f with %d groups for isotope "
               "%d which contains %d energy groups", num_groups, _uid,
               _num_groups);

  for (int i=0; i < _num_groups; i++)
    _sigma_f[i] = xs[i];

  /* Determine whether or not this isotope is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_groups; i++) {
    if (_sigma_f[i] > 0.0) {
      _fissionable = true;
      return;
    }
  }
}


/**
 * @brief Set the isotope's fission cross-section for some energy group.
 * @param xs the fission cross-section (\f$ \Sigma_f [cm^1] \f$)
 * @param group the energy group
 */
void Isotope::setSigmaFByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set sigma_f for group %d for isotope "
              "%d which contains %d energy groups", group, _uid, _num_groups);

  _sigma_f[group-1] = xs;

  /* Determine whether or not this isotope is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_groups; i++) {
    if (_sigma_f[i] > 0.0) {
      _fissionable = true;
      return;
    }
  }
}


/**
 * @brief Set the isotope's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param xs the array of fission cross-sections multiplied by nu 
 *        \f$ \nu \f$
 * @param num_groups the number of energy groups 
*/
void Isotope::setNuSigmaF(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for isotope "
               "%d which contains %d energy groups", num_groups, _uid, _num_groups);

  for (int i=0; i < _num_groups; i++)
    _nu_sigma_f[i] = xs[i];
}


/**
 * @brief Set the isotope's fission cross-section multiplied by \f$ \nu \f$
 *        for some energy group.
 * @param xs the fission cross-section (\f$ \nu\Sigma_f [cm^1] \f$)
 * @param group the energy group
 */
void Isotope::setNuSigmaFByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set nu_sigma_f for group %d for isotope "
               "%d which contains %d energy groups", group, _uid, _num_groups);

  _nu_sigma_f[group-1] = xs;
}



/**
 * @brief Set the isotope's array of \f$ \chi \f$ values.
 * @param xs the array of chi \f$ \chi \f$ values
 * @param num_groups the number of energy groups 
 */
void Isotope::setChi(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set chi with %d groups for isotope "
               "%d which contains %d energy groups", num_groups, _uid,
               _num_groups);

  for (int i=0; i < _num_groups; i++)
    _chi[i] = xs[i];
}


/**
 * @brief Set the isotope's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
void Isotope::setChiByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set chi for group %d for isotope "
         "%d which contains %d energy groups", 
         group, _uid, _num_groups);

  _chi[group-1] = xs;
}


/**
 * @brief Get the scatter source for an isotope
 * @param group The energy group
 * @param flux Pointer to the array of flux values
 * @return The scatter source
 */
FP_PRECISION Isotope::getScatterSource(int group, FP_PRECISION* flux) {
  
  FP_PRECISION scatter_source = 0;
  
  /* loop over origin groups */
  /* Note: could switch to pairwise sum if needed */
  for (int g=0; g < _num_groups; g++)
    scatter_source += getSigmaSByGroup(g+1, group+1) * flux[g];
                                   
  return scatter_source;  
  
}


/**
 * @brief Returns whether or not the Isotope contains a fissionable (non-zero)
 *        fission cross-section.
 * @return true if fissionable, false otherwise
 */
bool Isotope::isFissionable() {
  return _fissionable;
}
