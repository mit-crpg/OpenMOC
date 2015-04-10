#include "Material.h"

static int auto_id = 10000;


/**
 * @brief Returns an auto-generated unique Material ID.
 * @details This method is intended as a utility method for user's writing
 *          OpenMOC input files. The method makes use of a static Material
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
 * @brief Resets the auto-generated unique Material ID counter to 10000.
 */
void reset_material_id() {
  auto_id = 10000;
}


/**
 * @brief Constructor sets the ID and unique ID for the Material.
 * @param id the user-specified optional Material ID
 * @param name the user-specified optional Material name
 */
Material::Material(int id, const char* name) {

  /* If the user did not define an optional ID, create one */
  if (id == 0)
    _id = material_id();

  /* Use the user-defined ID */
  else
    _id = id;

  _name = NULL;
  setName(name);

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
 * @brief Destructor deletes all cross-section data structures from memory.
 */
Material::~Material() {

  if (_name != NULL)
    delete [] _name;

  /* If data is vector aligned */
  if (_data_aligned) {
    if (_sigma_t != NULL)
      MM_FREE(_sigma_t);

    if (_sigma_a != NULL)
      MM_FREE(_sigma_a);

    if (_sigma_s != NULL)
      MM_FREE(_sigma_s);

    if (_sigma_f != NULL)
      MM_FREE(_sigma_f);

    if (_nu_sigma_f != NULL)
      MM_FREE(_nu_sigma_f);

    if (_chi != NULL)
      MM_FREE(_chi);

    /* Whomever SIMD'izes OpenMOC will need to properly free these
     * using mm_free if they are vector aligned */
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
 * @brief Return the Material's user-defined ID
 * @return the Material's user-defined ID
 */
int Material::getId() const {
  return _id;
}


/**
 * @brief Return the user-defined name of the Material
 * @return the Material name
 */
char* Material::getName() const {
  return _name;
}


/**
 * @brief Returns the number of energy groups for this Material's nuclear data.
 * @return the number of energy groups
 */
int Material::getNumEnergyGroups() const {
  return _num_groups;
}


/**
 * @brief Return the array of the Material's total cross-sections.
 * @return the pointer to the Material's array of total cross-sections
 */
FP_PRECISION* Material::getSigmaT() {
  if (_sigma_t == NULL)
    log_printf(ERROR, "Unable to return Material %d's total "
               "cross-section since it has not yet been set", _id);

  return _sigma_t;
}



/**
 * @brief Return the array of the Material's absorption cross-sections.
 * @return the pointer to the Material's array of absorption cross-sections
 */
FP_PRECISION* Material::getSigmaA() {
  if (_sigma_a == NULL)
      log_printf(ERROR, "Unable to return Material %d's absorption "
                 "cross-section since it has not yet been set", _id);

  return _sigma_a;
}


/**
 * @brief Return the array of the Material's scattering cross-section matrix.
 * @return the pointer to the Material's array of scattering cross-sections
 */
FP_PRECISION* Material::getSigmaS() {
  if (_sigma_s == NULL)
    log_printf(ERROR, "Unable to return Material %d's scattering "
               "cross-section since it has not yet been set", _id);

  return _sigma_s;
}


/**
 * @brief Return the array of the Material's fission cross-sections.
 * @return the pointer to the Material's array of fission cross-sections
 */
FP_PRECISION* Material::getSigmaF() {
  if (_sigma_f == NULL)
    log_printf(ERROR, "Unable to return material %d's fission "
               "cross-section since it has not yet been set", _id);

  return _sigma_f;
}


/**
 * @brief Return the array of the Material's fission cross-sections
 *        multiplied by nu \f$ \nu \f$.
 * @return the pointer to the Material's array of fission cross-sections
 *         multiplied by nu \f$ \nu \f$
 */
FP_PRECISION* Material::getNuSigmaF() {
  if (_nu_sigma_f == NULL)
    log_printf(ERROR, "Unable to return Material %d's nu times fission "
               "cross-section since it has not yet been set", _id);

  return _nu_sigma_f;
}


/**
 * @brief Return the array of the Material's chi \f$ \chi \f$.
 * @return the pointer to the Material's array of chi \f$ \chi \f$ values
 */
FP_PRECISION* Material::getChi() {
  if (_chi == NULL)
    log_printf(ERROR, "Unable to return Material %d's chi spectrum "
               "since it has not yet been set", _id);

  return _chi;
}


/**
 * @brief Return the array of the Material's diffusion coefficients.
 * @return the pointer to the Material's array of diffusion coefficients
 */
FP_PRECISION* Material::getDifCoef() {

  if (_dif_coef == NULL){

    _dif_coef = new FP_PRECISION[_num_groups];

    for (int e = 0; e < _num_groups; e++)
      _dif_coef[e] = 0.0;
  }

  return _dif_coef;
}


/**
 * @brief Return the array of the Material's surface diffusion coefficients.
 * @details Returns the diffusion coefficients \f$ \hat{D} \f$ for this Material
 *           on each of the four surfaces of a CMFD mesh cell.
 * @return the pointer to the Material's array of surface diffusion coefficients
 */
FP_PRECISION* Material::getDifHat() {

  if (_dif_hat == NULL){

    _dif_hat = new FP_PRECISION[4*_num_groups];

    for (int e = 0; e < 4*_num_groups; e++)
      _dif_hat[e] = 0.0;
  }

  return _dif_hat;
}


/**
 * @brief Return the array of the Material's CMFD correction to the surface
 *        diffusion coefficients.
 * @return the pointer to the Material's array of CMFD correction to the
 *         surface diffusion coefficients
 */
FP_PRECISION* Material::getDifTilde() {

  if (_dif_tilde == NULL){

    _dif_tilde = new FP_PRECISION[4*_num_groups];

    for (int e = 0; e < 4*_num_groups; e++)
      _dif_tilde[e] = 0.0;
  }

  return _dif_tilde;
}


/**
 * @brief Return the array of the Material's CMFD correction to the surface
 *        diffusion coefficients.
 * @return the pointer to the Material's array of CMFD correction to
 *         the surface diffusion coefficients
 */
FP_PRECISION* Material::getBuckling() {

  if (_buckling == NULL){

    _buckling = new FP_PRECISION[_num_groups];

    for (int e = 0; e < _num_groups; e++)
      _buckling[e] = 0.0;
  }

  return _buckling;
}

/**
 * @brief Get the Material's total cross section for some energy group.
 * @param group the energy group
 * @return the total cross section
 */
FP_PRECISION Material::getSigmaTByGroup(int group) {    
  if (_sigma_t == NULL)
    log_printf(ERROR, "Unable to return Material %d's total "
               "cross section since it has not yet been set", _id);
               
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get sigma_t for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
               
  return _sigma_t[group-1];
}


/**
 * @brief Get the Material's absorption cross section for some energy group.
 * @param group the energy group
 * @return the absorption cross section
 */
FP_PRECISION Material::getSigmaAByGroup(int group) {    
  if (_sigma_a == NULL)
    log_printf(ERROR, "Unable to return Material %d's absorption "
               "cross section since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get sigma_a for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
  
  return _sigma_a[group-1];
}


/**
 * @brief Get the Material's scattering cross section for some energy group.
 * @param origin the incoming energy group
 * @param destination the outgoing energy group
 * @return the scattering cross section
 */
FP_PRECISION Material::getSigmaSByGroup(int origin, int destination) {   
  if (_sigma_s == NULL)
    log_printf(ERROR, "Unable to return Material %d's scattering "
               "cross section since it has not yet been set", _id);
  
  if (origin <= 0 || destination <= 0 || origin > _num_groups || destination > _num_groups)
    log_printf(ERROR, "Unable to get sigma_s for group %d,%d for Material %d "
               "which contains %d energy groups",
               origin, destination, _id, _num_groups);
   
  return getSigmaSByGroupInline(origin-1,destination-1);
  
}


/**
 * @brief Get the Material's fission cross section for some energy group.
 * @param group the energy group
 * @return the fission cross section
 */
FP_PRECISION Material::getSigmaFByGroup(int group) {    
  if (_sigma_f == NULL)
    log_printf(ERROR, "Unable to return material %d's fission "
               "cross section since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);  
  
  return _sigma_f[group-1];
}    


/**
 * @brief Get the Material's nu-fission cross section for some energy group.
 * @param group the energy group
 * @return the nu-fission cross section
 */
FP_PRECISION Material::getNuSigmaFByGroup(int group) {    
  if (_nu_sigma_f == NULL)
    log_printf(ERROR, "Unable to return Material %d's nu-fission "
               "cross section since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get nu_sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
  
  return _nu_sigma_f[group-1];
}    


/**
 * @brief Get the Material's fission spectrum for some energy group.
 * @param group the energy group
 * @return the fission spectrum
 */
FP_PRECISION Material::getChiByGroup(int group) {        
  if (_chi == NULL)
    log_printf(ERROR, "Unable to return Material %d's chi spectrum "
               "since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get chi for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
  
  return _chi[group-1];
}  


/**
 * @brief Get the Material's diffustion coefficient for some energy group.
 * @param group the energy group
 * @return the diffusion coefficient
 */
FP_PRECISION Material::getDifCoefByGroup(int group) {
  if (_dif_coef == NULL)
    log_printf(ERROR, "Unable to return Material %d's diffusion coefficient "
               "since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get dif_coef for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
  
  return _dif_coef[group-1];
}


/**
 * @brief Get the Material's dif_hat for some energy group.
 * @param group the energy group
 * @param surface the index of the surface
 * @return the dif_hat value
 */
FP_PRECISION Material::getDifHatByGroup(int group, int surface) {
  if (_dif_hat == NULL)
    log_printf(ERROR, "Unable to return Material %d's dif_hat "
               "since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get dif_hat for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
               
  return _dif_hat[surface*_num_groups + (group-1)];
}


/**
 * @brief Get the Material's dif_tilde for some energy group.
 * @param group the energy group
 * @return the dif_tilde value
 */
FP_PRECISION Material::getDifTildeByGroup(int group) {
  if (_dif_tilde == NULL)
    log_printf(ERROR, "Unable to return Material %d's dif_tilde "
               "since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get dif_tilde for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
  
  return _dif_tilde[group-1];
}


/**
 * @brief Get the Material's buckling for some energy group.
 * @param group the energy group
 * @return the buckling value
 */
FP_PRECISION Material::getBucklingByGroup(int group) {
  if (_buckling == NULL)
    log_printf(ERROR, "Unable to return Material %d's buckling "
               "since it has not yet been set", _id);
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get buckling for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);
  
  return _buckling[group-1];
}



/**
 * @brief Returns whether or not the Material contains a fissionable (non-zero)
 *        fission cross-section.
 * @return true if fissionable, false otherwise
 */
bool Material::isFissionable() {
  return _fissionable;
}


/**
 * @brief Returns true if the data is vector aligned, false otherwise (default).
 * @return Whether or not the Material's data is vector aligned
 */
bool Material::isDataAligned() {
  return _data_aligned;
}


/**
 * @brief Returns the rounded up number of energy groups to fill an integral
 *        number of vector lengths.
 * @return The number of vector-aligned energy groups
 */
int Material::getNumVectorGroups() {
  return _num_vector_groups;
}


/**
 * @brief Sets the name of the Material
 * @param name the Material name string
 */
void Material::setName(const char* name) {
  int length = strlen(name);

  if (_name != NULL)
    delete [] _name;

  /* Initialize a character array for the Material's name */
  _name = new char[length+1];

  /* Copy the input character array Material name to the class attribute name */
  for (int i=0; i <= length; i++)
    _name[i] = name[i];
}


/**
 * @brief Set the number of energy groups for this Material.
 * @param num_groups the number of energy groups.
 */
void Material::setNumEnergyGroups(const int num_groups) {

  if (num_groups < 0)
    log_printf(ERROR, "Unable to set the number of energy groups for "
               "material %d to %d", _id, num_groups);

  _num_groups = num_groups;

  /* Free old data arrays if they were allocated for a previous simulation */

  /* If data is vector aligned */
  if (_data_aligned) {
    if (_sigma_t != NULL)
      MM_FREE(_sigma_t);

    if (_sigma_a != NULL)
      MM_FREE(_sigma_a);

    if (_sigma_s != NULL)
      MM_FREE(_sigma_s);

    if (_sigma_f != NULL)
      MM_FREE(_sigma_f);

    if (_nu_sigma_f != NULL)
      MM_FREE(_nu_sigma_f);

    if (_chi != NULL)
      MM_FREE(_chi);
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

  /* Allocate memory for data arrays */
  _sigma_t = new FP_PRECISION[_num_groups];
  _sigma_a = new FP_PRECISION[_num_groups];
  _sigma_f = new FP_PRECISION[_num_groups];
  _nu_sigma_f = new FP_PRECISION[_num_groups];
  _chi = new FP_PRECISION[_num_groups];
  _sigma_s = new FP_PRECISION[_num_groups*_num_groups];


  /* Assign the null vector to each data array */
  memset(_sigma_t, 0.0, sizeof(FP_PRECISION) * _num_groups);
  memset(_sigma_a, 0.0, sizeof(FP_PRECISION) * _num_groups);
  memset(_sigma_f, 0.0, sizeof(FP_PRECISION) * _num_groups);
  memset(_nu_sigma_f, 0.0, sizeof(FP_PRECISION) * _num_groups);
  memset(_chi, 0.0, sizeof(FP_PRECISION) * _num_groups);
  memset(_sigma_s, 0.0, sizeof(FP_PRECISION) * _num_groups * _num_groups);
}


/**
 * @brief Set the Material's array of total cross-sections.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's total cross-sections. An example of how this function
 *          might be called in Python is as follows:
 *
 * @code
 *          sigma_t = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaT(sigma_t)
 * @endcode
 *
 *          NOTE: This routine will override an zero-valued cross-sections 
 *          (e.g., in void or gap regions) with a minimum value of 1E-10 to
 *          void numerical issues in the MOC solver.
 *
 * @param xs the array of total cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaT(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_t with %d groups for Material "
               "%d which contains %d energy groups", num_groups, 
               _id, _num_groups);

  for (int i=0; i < _num_groups; i++) {

    /* If the cross-section is near zero (e.g., within (-1E-10, 1E-10)) */
    if (fabs(xs[i]) < ZERO_SIGMA_T) {
      log_printf(WARNING, "Overriding zero cross-section in "
                 "group %d for Material %d with 1E-10", i, _id);
      _sigma_t[i] = FP_PRECISION(ZERO_SIGMA_T);
    }
    else
      _sigma_t[i] = FP_PRECISION(xs[i]);
  }
}


/**
 * @brief Set the Material's total cross-section for some energy group.
 * @param xs the total cross-section
 * @param group the energy group
 */
void Material::setSigmaTByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set sigma_t for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);

    /* If the cross-section is near zero (e.g., within (-1E-10, 1E-10)) */
    if (fabs(xs) < ZERO_SIGMA_T) {
      log_printf(WARNING, "Overriding zero cross-section in "
                 "group %d for Material %d with 1E-10", group, _id);
      _sigma_t[group-1] = FP_PRECISION(ZERO_SIGMA_T);
    }
    else
      _sigma_t[group-1] = FP_PRECISION(xs);
}


/**
 * @brief Set the Material's array of absorption scattering cross-sections.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's absorption cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          sigma_a = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaA(sigma_a)
 * @endcode
 *
 * @param xs the array of absorption scattering cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaA(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_a with %d groups for Material "
               "%d which contains %d energy groups", num_groups, _id, _num_groups);

  for (int i=0; i < _num_groups; i++){
    _sigma_a[i] = FP_PRECISION(xs[i]);

    if (_buckling != NULL & _dif_coef != NULL)
      _sigma_a[i] += FP_PRECISION(_buckling[i] * _dif_coef[i]);
  }
}


/**
 * @brief Set the Material's absorption cross-section for some energy group.
 * @param xs the absorption cross-section
 * @param group the energy group
 */
void Material::setSigmaAByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set sigma_a for group %d for material "
               "%d which contains %d energy groups", group, _id, _num_groups);

  _sigma_a[group-1] = xs;
}




/**
 * @brief Set the Material's 2D array of scattering cross-sections.
 * @details The array should be passed to OpenMOC as a 1D array in 
 *          column-major order.  This assumes the standard convention, 
 *          where column index is the origin group and the row index is
 *          the destination group.  That is, the array should be ordered
 *          as follows:
 *              1 -> 1
 *              1 -> 2
 *              1 -> 3
 *                ...
 *              2 -> 1
 *              2 -> 2
 *                ...         
 *
 *          Note that if the scattering matrix is defined in NumPy by
 *          the standard convention, "flat" will put the matrix into row
 *          major order.  Thus, one should transpose the matrix before
 *          flattening. 
 * 
 *          For cache efficiency, the transpose of the input is actually
 *          stored in OpenMOC.
 * 
 *          This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the square of the number of energy groups) as input to this 
 *          function. This function then fills the NumPy array with the data 
 *          values for the Material's scattering cross-sections. An example 
 *          of how this function might be called in Python is as follows:
 *
 * @code
 *          sigma_s = numpy.array([[0.05,    0,    0,     ... ],
 *                                 [0.10, 0.08,   ...     ... ],
 *                                             ...
 *                                 [...        ...        ... ]])
 *          sigma_s = numpy.transpose(sigma_s)
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaS(sigma_s.flat)
 * @endcode
 *
 * @param xs the array of scattering cross-sections
 * @param num_groups_squared the number of energy groups squared
 */
void Material::setSigmaS(double* xs, int num_groups_squared) {

  if (_num_groups*_num_groups != num_groups_squared)
    log_printf(ERROR, "Unable to set sigma_s with %f groups for Material %d "
               "which contains %d energy groups",
                float(sqrt(num_groups_squared)), _id, _num_groups);

  for (int dest=0; dest < _num_groups; dest++) {
    for (int orig=0; orig < _num_groups; orig++)
      _sigma_s[dest*_num_groups+orig] = xs[orig*_num_groups+dest];
  }
}


/**
 * @brief Set the Material's scattering cross-section for some energy group.
 * @param xs the scattering cross-section
 * @param origin the column index in the scattering matrix
 * @param destination the row index in the scattering matrix
 */
void Material::setSigmaSByGroup(double xs, int origin, int destination) {

  if (origin <= 0 || destination <= 0 || origin > _num_groups || destination > _num_groups)
    log_printf(ERROR, "Unable to set sigma_s for group %d -> %d for Material %d "
               "which contains %d energy groups",
               origin, destination, _id, _num_groups);

  _sigma_s[_num_groups*(destination-1) + (origin-1)] = xs;
}


/**
 * @brief Set the Material's array of fission cross-sections.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's fission cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          sigma_f = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaF(sigma_f)
 * @endcode
 *
 * @param xs the array of fission cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaF(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_f with %d groups for Material "
               "%d which contains %d energy groups", num_groups, _id, _num_groups);

  for (int i=0; i < _num_groups; i++)
    _sigma_f[i] = xs[i];

  /* Determine whether or not this Material is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_groups; i++) {
    if (_sigma_f[i] > 0.0 || _nu_sigma_f[i] > 0.0) {
      _fissionable = true;
      return;
    }
  }
}


/**
 * @brief Set the Material's fission cross-section for some energy group.
 * @param xs the fission cross-section
 * @param group the energy group
 */
void Material::setSigmaFByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);

  _sigma_f[group-1] = xs;

  /* Determine whether or not this Material is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_groups; i++) {
    if (_sigma_f[i] > 0.0 || _nu_sigma_f[i] > 0.0) {
      _fissionable = true;
      return;
    }
  }
}


/**
 * @brief Set the Material's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param xs the array of fission cross-sections multiplied by nu
 *        \f$ \nu \f$
 * @param num_groups the number of energy groups
*/
void Material::setNuSigmaF(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for Material %d "
              "which contains %d energy groups", num_groups, _id, _num_groups);

  for (int i=0; i < _num_groups; i++)
    _nu_sigma_f[i] = xs[i];
    
  /* Determine whether or not this Material is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_groups; i++) {
    if (_sigma_f[i] > 0.0 || _nu_sigma_f[i] > 0.0) {
      _fissionable = true;
      return;
    }
  }
}


/**
 * @brief Set the Material's fission cross-section multiplied by \f$ \nu \f$
 *        for some energy group.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's nu*fission cross-sections. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          nu_sigma_f = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setNuSigmaF(nu_sigma_f)
 * @endcode
 *
 * @param xs the fission cross-section multiplied by nu \f$ \nu \f$
 * @param group the energy group
 */
void Material::setNuSigmaFByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set nu_sigma_f for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);

  _nu_sigma_f[group-1] = xs;

  /* Determine whether or not this Material is fissionable */
  _fissionable = false;

  for (int i=0; i < _num_groups; i++) {
    if (_sigma_f[i] > 0.0 || _nu_sigma_f[i] > 0.0) {
      _fissionable = true;
      return;
    }
  }
}



/**
 * @brief Set the Material's array of chi \f$ \chi \f$ values.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's chi distribution. An example of how this function might
 *          be called in Python is as follows:
 *
 * @code
 *          chi = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setChi(chi)
 * @endcode
 *
 * @param xs the array of chi \f$ \chi \f$ values
 * @param num_groups the number of energy groups
 */
void Material::setChi(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set chi with %d groups for Material "
               "%d which contains %d energy groups", num_groups, _id, _num_groups);

  double chi_sum = 0.0;
  for (int i=0; i < _num_groups; i++)
    chi_sum += xs[i];

  for (int i=0; i < _num_groups; i++){
    if (chi_sum == 0)
      _chi[i] = xs[i];
    else
      _chi[i] = xs[i] / chi_sum;
  }
}


/**
 * @brief Set the Material's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
void Material::setChiByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set chi for group %d for Material "
              "%d which contains %d energy groups", group, _id, _num_groups);

  _chi[group-1] = xs;
}


/**
 * @brief Set the Material's array of diffusion coefficients.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's diffusion coefficients. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          dif_coef = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setDifCoef(dif_coef)
 * @endcode
 *
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setDifCoef(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient with %d groups for "
               "Material %d which contains %d energy groups", num_groups,
               _id, _num_groups);

  if (_dif_coef == NULL)
    _dif_coef = new FP_PRECISION[_num_groups];

  for (int i=0; i < _num_groups; i++)
    _dif_coef[i] = xs[i];
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDifCoefByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_groups);

  if (_dif_coef == NULL){
    _dif_coef = new FP_PRECISION[_num_groups];

    for (int i=0; i < _num_groups; i++)
      _dif_coef[i] = 0.0;
  }

  _dif_coef[group-1] = xs;
}


/**
 * @brief Set the Material's array of diffusion coefficients.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's buckling coefficients. An example of how this function
 *          might be called in Python is as follows:
 *
 * @code
 *          buckling = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setBuckling(buckling)
 * @endcode
 *
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setBuckling(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient with %d groups for "
               "Material %d which contains %d energy groups", num_groups,
               _id, _num_groups);

  if (_buckling == NULL)
    _buckling = new FP_PRECISION[_num_groups];

  for (int i=0; i < _num_groups; i++)
    _buckling[i] = xs[i];
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setBucklingByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient for group %d for "
               "Material %d which contains %d energy groups",
               group, _id, _num_groups);

  if (_buckling == NULL){
    _buckling = new FP_PRECISION[_num_groups];
    for (int i=0; i < _num_groups; i++)
      _buckling[i] = 0.0;
  }

  _buckling[group-1] = xs;
}


/**
 * @brief Set the Material's array of diffusion coefficients.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's surface diffusion coefficients. An example of how this
 *          function might be called in Python is as follows:
 *
 * @code
 *          dif_hat = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setDifHat(dif_hat)
 * @endcode
 *
 * @param xs the array of surface diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setDifHat(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
               "for Material %d which contains %d energy groups", num_groups,
               _id, _num_groups);

  if (_dif_hat == NULL)
    _dif_hat = new FP_PRECISION[4*_num_groups];

  for (int i=0; i < _num_groups*4; i++)
    _dif_hat[i] = xs[i];
}


/**
 * @brief Set the Material's diffusion coefficient for some energy group along
 *        some CMFD mesh cell.
 * @param xs the diffusion coefficient along some mesh cell surface
 * @param group the energy group
 * @param surface the Surface corresponding to this coefficient
 */
void Material::setDifHatByGroup(double xs, int group, int surface) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient for group %d for "
              "Material %d which contains %d energy groups",
              group, _id, _num_groups);

  if (_dif_hat == NULL){

    _dif_hat = new FP_PRECISION[4*_num_groups];

    for (int i=0; i < _num_groups*4; i++)
      _dif_hat[i] = 0.0;
  }

  _dif_hat[surface*_num_groups + (group-1)] = xs;
}


/**
 * @brief Set the Material's array of diffusion coefficients.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Material's nuclear data in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of energy groups) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Material's CMFD corrected diffusion coefficients. An example of how
 *          this function might be called in Python is as follows:
 *
 * @code
 *          dif_tilde = numpy.array([0.05, 0.1, 0.15, ... ])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setDifTilde(dif_tilde)
 * @endcode
 *
 * @param xs the array of CMFD corrected diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setDifTilde(double* xs, int num_groups) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
              "for Material %d which contains %d energy groups", num_groups,
              _id, _num_groups);

  if (_dif_tilde == NULL){
    _dif_tilde = new FP_PRECISION[4*_num_groups];
  }

  for (int i=0; i < _num_groups*4; i++)
    _dif_tilde[i] = xs[i];
}


/**
 * @brief Set the Material's CMFD corrected diffusion coefficient for some
 *        energy group.
 * @param xs the CMFD corrected diffusion coefficient
 * @param group the energy group
 * @param surface the Surface corresponding to this coefficient
 */
void Material::setDifTildeByGroup(double xs, int group, int surface) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set diffusion coefficient correction for "
              "group %d for Material %d which contains %d energy groups",
               group, _id, _num_groups);

  if (_dif_tilde == NULL){
    _dif_tilde = new FP_PRECISION[4*_num_groups];

    for (int i=0; i < _num_groups*4; i++)
      _dif_tilde[i] = 0.0;
  }

  _dif_tilde[surface*_num_groups + (group-1)] = xs;
}


/**
 * @brief Checks if the total cross-section for this Material is equal to the
 *        absorption plus scattering cross-sections for all energy groups.
 * @details If the total cross-section does not equal the absorption plus
 *          scattering cross-section within SIGMA_T_THRESH then this method
 *          exits OpenMOC.
 */
void Material::checkSigmaT() {

  if (_num_groups == 0)
    log_printf(ERROR, "Unable to verify Material %d's total cross-section "
              "since the number of energy groups has not been set", _id);
  if (_sigma_t == NULL)
    log_printf(ERROR, "Unable to verify Material %d's total cross-section "
              "since its total cross-section has not been set", _id);
  if (_sigma_a == NULL)
    log_printf(ERROR, "Unable to verify Material %d's total cross-section "
               "since its absorption cross-section has not been set", _id);
  if (_sigma_s == NULL)
    log_printf(ERROR, "Unable to verify Material %d's total cross-section "
              "since its scattering cross-section has not been set", _id);

  FP_PRECISION calc_sigma_t;

  /* Loop over all energy groups */
  for (int i=0; i < _num_groups; i++) {

    /* Initialize the calculated total xs to the absorption xs */
    calc_sigma_t = _sigma_a[i];

    /* Increment calculated total xs by scatter xs for each energy group */
    for (int j=0; j < _num_groups; j++)
      calc_sigma_t += _sigma_s[i+j*_num_groups];

    /* Check if the calculated and total match up to certain threshold */
    if (fabs(calc_sigma_t - _sigma_t[i]) > SIGMA_T_THRESH) {
      log_printf(ERROR, "Material id = %d has a different total cross-section "
                 "than the sum of its scattering and absorption cross-sections "
                 "for group %d: sigma_t = %f, calc_sigma_t = %f",
                 _id, i+1, _sigma_t[i], calc_sigma_t);
    }
  }

  return;
}


/**
 * @brief Converts this Material's attributes to a character array
 *        representation.
 * @details The character array returned includes the user-defined ID, and each
 *          of the absorption, total, fission, nu multiplied by fission and
 *          scattering cross-sections and chi for all energy groups.
 * @return character array of this Material's attributes
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
 * @brief Prints a string representation of all of the Material's attributes to
 *        the console.
 */
void Material::printString() {
  log_printf(NORMAL, toString().c_str());
}


/**
 * @brief Reallocates the Material's cross-section data structures along
 *        word-aligned boundaries
 * @details This method is used to assist with SIMD auto-vectorization of the
 *          MOC routines in the Solver classes. Rather than using the assigned
 *          number of energy groups, this method adds "dummy" energy groups
 *          such that the total number of groups is some multiple of VEC_LENGTH
 *          (typically 4, 8, or 16). As a result, the SIMD-vectorized Solver
 *          subclasses can break up loops over energy groups in such a way
 *          to "expose" the SIMD nature of the algorithm.
 */
void Material::alignData() {

  /* If the data has already been aligned, do nothing */
  if (_data_aligned)
    return;

  if (_num_groups <= 0)
    log_printf(ERROR, "Unable to align Material %d data since the "
               "cross-sections have not yet been set\n", _id);

  _num_vector_groups = (_num_groups / VEC_LENGTH) + 1;

  /* Allocate memory for the new aligned xs data */
  int size = _num_vector_groups * VEC_LENGTH * sizeof(FP_PRECISION);

  /* Allocate word-aligned memory for cross-section data arrays */
  FP_PRECISION* new_sigma_t = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_sigma_a = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_sigma_f = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_nu_sigma_f=(FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_chi = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

  /* The scattering matrix will be the number of vector groups
   * wide (SIMD) and the actual number of groups long since
   * instructions are not SIMD in this dimension */

  size = _num_vector_groups * VEC_LENGTH * _num_vector_groups;
  size *= VEC_LENGTH * sizeof(FP_PRECISION);
  FP_PRECISION* new_sigma_s = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

  /* Initialize data structures to ones for sigma_t since it is used to
   * divide the source in the solver, and zeroes for everything else */
  size = _num_vector_groups * VEC_LENGTH * sizeof(FP_PRECISION);
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
  size = _num_groups * sizeof(FP_PRECISION);
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


/**
 * @brief Create a duplicate of the Material.
 * @return a pointer to the clone
 */
Material* Material::clone(){

  Material* clone = new Material(getId());

  clone->setNumEnergyGroups(_num_groups);

  for (int i=0; i < _num_groups; i++) {
    clone->setSigmaTByGroup((double)_sigma_t[i], i+1);
    clone->setSigmaAByGroup((double)_sigma_a[i], i+1);
    clone->setSigmaFByGroup((double)_sigma_f[i], i+1);
    clone->setNuSigmaFByGroup((double)_nu_sigma_f[i], i+1);
    clone->setChiByGroup((double)_chi[i], i+1);

    for (int j=0; j < _num_groups; j++)
      clone->setSigmaSByGroup(
        (double)getSigmaSByGroupInline(i,j), i+1, j+1);

    if (_dif_coef != NULL)
      clone->setDifCoefByGroup((double)_dif_coef[i], i+1);

    if (_buckling != NULL)
      clone->setBucklingByGroup((double)_buckling[i], i+1);

    for (int j=0; j < 4; j++) {

      if (_dif_hat != NULL)
        clone->setDifHatByGroup((double)_dif_hat[i*4+j], i+1, j);

      if (_dif_tilde != NULL)
        clone->setDifTildeByGroup((double)_dif_tilde[i*4+j], i+1, j);
    }
  }

  return clone;
}
