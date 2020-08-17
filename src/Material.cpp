#include "Material.h"

static int auto_id = DEFAULT_INIT_ID;


/**
 * @brief Returns an auto-generated unique Material ID.
 * @details This method is intended as a utility method for user's writing
 *          OpenMOC input files. The method makes use of a static Material
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 1,000,000. Hence, user-defined material IDs greater
 *          than or equal to 1,000,000 is prohibited.
 */
int material_id() {
  int id = auto_id;
  auto_id++;
  return id;
}


/**
 * @brief Resets the auto-generated unique Material ID counter to 1,000,000.
 */
void reset_material_id() {
  auto_id = DEFAULT_INIT_ID;
}


/**
 * @brief Maximize the auto-generated unique Material ID counter.
 * @details This method updates the auto-generated unique Material ID
 *          counter if the input parameter is greater than the present
 *          value. This is useful for the OpenMC compatibility module
 *          to ensure that the auto-generated Material IDs do not
 *          collide with those created in OpenMC.
 * @param material_id the id assigned to the auto-generated counter
 */
void maximize_material_id(int material_id) {
  if (material_id > auto_id)
    auto_id = material_id;
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
  if (name != NULL)
    setName(name);

  _volume = 0.;
  _num_instances = 0;

  /* Initialize a dummy number groups */
  _num_groups = -1;

  _sigma_t = NULL;
  _max_sigma_t = -1;
  _sigma_s = NULL;
  _sigma_a = NULL;
  _sigma_f = NULL;
  _nu_sigma_f = NULL;
  _chi = NULL;
  _fiss_matrix = NULL;

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

    if (_sigma_s != NULL)
      MM_FREE(_sigma_s);

    if (_sigma_f != NULL)
      MM_FREE(_sigma_f);

    if (_nu_sigma_f != NULL)
      MM_FREE(_nu_sigma_f);

    if (_chi != NULL)
      MM_FREE(_chi);

    if (_fiss_matrix != NULL)
      MM_FREE(_fiss_matrix);
  }

  /* Data is not vector aligned */
  else {
    if (_sigma_t != NULL)
      free(_sigma_t);

    if (_sigma_s != NULL)
      delete [] _sigma_s;

    if (_sigma_a != NULL)
      delete [] _sigma_a;

    if (_sigma_f != NULL)
      delete [] _sigma_f;

    if (_nu_sigma_f != NULL)
      delete [] _nu_sigma_f;

    if (_chi != NULL)
      delete [] _chi;

    if (_fiss_matrix != NULL)
      delete [] _fiss_matrix;
  }
}


/**
 * @brief Return the Material's user-defined ID.
 * @return the Material's user-defined ID
 */
int Material::getId() const {
  return _id;
}


/**
 * @brief Return the user-defined name of the Material.
 * @return the Material name
 */
char* Material::getName() const {
  return _name;
}


/**
 * @brief Return the aggregate volume/area of all instances of this Material.
 * @details The volume/area of the Material is computed from track segments
 *          which overlap this Material during track generation.
 * @return the volume/area of the Material
 */
double Material::getVolume() {
  return _volume;
}


/**
 * @brief Return the number of instances of this Material in the Geometry.
 * @details The number of instances of this Material in the Geometry is
 *          determined during track generation.
 * @return the number of material instances
 */
int Material::getNumInstances() {
  return _num_instances;
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
 * @brief Return the maximum of the Material's total cross-sections.
 * @return the maximum of the Material's total cross-sections.
 */
FP_PRECISION Material::getMaxSigmaT() {
  return _max_sigma_t;
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
 * @brief Return the array of the Material's absorption group cross-sections.
 * @return the pointer to the Material's array of absorption cross-sections
 */
FP_PRECISION* Material::getSigmaA() {

  if (_num_groups <= 0)
    log_printf(ERROR, "Unable to return Material %d's absorption "
               "cross-section since it has %d groups", _id, _num_groups);

  if (_sigma_s == NULL || _sigma_t == NULL)
    log_printf(ERROR, "Unable to return Material %d's absorption "
               "cross-section since scattering and total have not been set",
               _id);

  /* If not initialized, compute _sigma_a the absorption cross section */
  if (_sigma_a == NULL) {
    _sigma_a = new FP_PRECISION[_num_groups];
    for (int g=0; g < _num_groups; g++) {
      _sigma_a[g] = _sigma_t[g];
      for (int gp=0; gp < _num_groups; gp++)
        _sigma_a[g] -= _sigma_s[gp*_num_groups + g];
    }
  }

  return _sigma_a;
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
 * @brief Return the array of the Material's fission matrix.
 * @return the pointer to the Material's fission matrix array
 */
FP_PRECISION* Material::getFissionMatrix() {
  if (_fiss_matrix == NULL)
    buildFissionMatrix();

  return _fiss_matrix;
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

  else if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get sigma_t for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);

  return _sigma_t[group-1];
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

  else if (origin <= 0 || destination <= 0 ||
            origin > _num_groups || destination > _num_groups)
    log_printf(ERROR, "Unable to get sigma_s for group %d,%d for Material %d "
               "which contains %d energy groups",
               origin, destination, _id, _num_groups);

  return _sigma_s[(destination-1)*_num_groups + (origin-1)];
}


/**
 * @brief Get the Material's absorption cross section for some energy group.
 * @param group the energy group
 * @return the absorption cross section
 */
FP_PRECISION Material::getSigmaAByGroup(int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get sigma_a for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);

  if (_sigma_s == NULL || _sigma_t == NULL)
    log_printf(ERROR, "Unable to return Material %d's absorption "
               "cross-section since scattering and total have not been set",
               _id);

  /* If not initialized, compute _sigma_a the absorption cross section */
  if (_sigma_a == NULL) {
    _sigma_a = new FP_PRECISION[_num_groups];
    for (int g=0; g < _num_groups; g++) {
      _sigma_a[g] = _sigma_t[g];
      for (int gp=0; gp < _num_groups; gp++)
        _sigma_a[g] -= _sigma_s[gp*_num_groups + g];
    }
  }

  return _sigma_a[group-1];
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

  else if (group <= 0 || group > _num_groups)
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

  else if (group <= 0 || group > _num_groups)
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

  else if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to get chi for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);

  return _chi[group-1];
}


/**
 * @brief Get the Material's fission matrix for some energy group.
 * @param origin the incoming energy group \f$ E_{0} \f$
 * @param destination the outgoing energy group \f$ E_{1} \f$
 * @return the fission matrix entry \f$ \nu\Sigma_{f}(E_{0}) * \chi(E_{1})\f$
 */
FP_PRECISION Material::getFissionMatrixByGroup(int origin, int destination) {
  if (_fiss_matrix == NULL)
    buildFissionMatrix();

  else if (origin <= 0 || destination <= 0 ||
           origin > _num_groups || destination > _num_groups)
    log_printf(ERROR, "Unable to get fission matrix for group %d,%d for "
               "Material %d which contains %d energy groups",
               origin, destination, _id, _num_groups);

  return _fiss_matrix[(destination-1)*_num_groups + (origin-1)];
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
 * @brief Set the volume/area of the Material.
 * @param volume the volume/area of the Material
 */
void Material::setVolume(double volume) {
  _volume = volume;
}


/**
 * @brief Increment the volume/area of the Material by some amount.
 * @details This routine is called by the TrackGenerator during track
 *          generation and segmentation.
 * @param volume the amount to increment the current volume by
 */
void Material::incrementVolume(double volume) {
  _volume += volume;
}


/**
 * @brief Set the number of instances of this Material.
 * @param num_instances the number of instances of this Material in the Geometry
 */
void Material::setNumInstances(int num_instances) {
  _num_instances = num_instances;
}


/**
 * @brief Increment the number of instances of this Material.
 * @details This routine is called by the TrackGenerator during track
 *          generation and segmentation.
 */
void Material::incrementNumInstances() {
  _num_instances++;
}


/**
 * @brief Set the number of energy groups for this Material.
 * @param num_groups the number of energy groups.
 */
void Material::setNumEnergyGroups(const int num_groups) {

  if (num_groups < 1)
    log_printf(ERROR, "Unable to set the number of energy groups for "
               "material %d to %d", _id, num_groups);

  _num_groups = num_groups;

  /* Free old data arrays if they were allocated for a previous simulation */

  /* If data is vector aligned */
  if (_data_aligned) {
    if (_sigma_t != NULL)
      MM_FREE(_sigma_t);

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
      free(_sigma_t);

    if (_sigma_s != NULL)
      delete [] _sigma_s;

    if (_sigma_f != NULL)
      delete [] _sigma_f;

    if (_nu_sigma_f != NULL)
      delete [] _nu_sigma_f;

    if (_chi != NULL)
      delete [] _chi;
  }

  /* Allocate memory for data arrays */
  _sigma_t = (FP_PRECISION*) MM_MALLOC(_num_groups*sizeof(FP_PRECISION),
                                       VEC_ALIGNMENT);
  _sigma_f = new FP_PRECISION[_num_groups];
  _nu_sigma_f = new FP_PRECISION[_num_groups];
  _chi = new FP_PRECISION[_num_groups];
  _sigma_s = new FP_PRECISION[_num_groups*_num_groups];

  /* Assign the null vector to each data array */
  memset(_sigma_t, 0, sizeof(FP_PRECISION) * _num_groups);
  memset(_sigma_f, 0, sizeof(FP_PRECISION) * _num_groups);
  memset(_nu_sigma_f, 0, sizeof(FP_PRECISION) * _num_groups);
  memset(_chi, 0, sizeof(FP_PRECISION) * _num_groups);
  memset(_sigma_s, 0, sizeof(FP_PRECISION) * _num_groups * _num_groups);
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
 *          sigma_t = numpy.array([0.05, 0.1, 0.15, ...])
 *          material = openmoc.Material(openmoc.material_id())
 *          material.setSigmaT(sigma_t)
 * @endcode
 *
 *          NOTE: This routine will override any zero-valued cross-sections
 *          (e.g., in void or gap regions) with a minimum value of 1E-10 to
 *          avoid numerical issues in the MOC solver.
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
      log_printf(INFO_ONCE, "Overriding zero cross-section in "
                 "group %d for Material %d with 1E-10", i, _id);
      _sigma_t[i] = FP_PRECISION(ZERO_SIGMA_T);
    }
    else
      _sigma_t[i] = FP_PRECISION(xs[i]);

    /* Update the maximum total cross-section */
    _max_sigma_t = std::max(_max_sigma_t, _sigma_t[i]);
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
      log_printf(INFO_ONCE, "Overriding zero cross-section in "
                 "group %d for Material %d with 1E-10", group, _id);
      _sigma_t[group-1] = FP_PRECISION(ZERO_SIGMA_T);
    }
    else
      _sigma_t[group-1] = FP_PRECISION(xs);

  /* Update the maximum total cross-section */
  _max_sigma_t = std::max(_max_sigma_t, _sigma_t[group-1]);
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
    log_printf(ERROR, "Unable to set sigma_s with %.4e groups for Material %d "
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

  if (origin <= 0 || destination <= 0 || origin > _num_groups ||
      destination > _num_groups)
    log_printf(ERROR, "Unable to set sigma_s for group %d -> %d for Material %d"
               " which contains %d energy groups",
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
               "%d which contains %d energy groups", num_groups, _id,
               _num_groups);

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

  /* If the SigmaF input is positive, this material is fissionable */
  if (xs > 0)
    _fissionable = true;
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

  /* If the SigmaF input is positive, this material is fissionable */
  if (xs > 0)
    _fissionable = true;
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

  for (int i=0; i < _num_groups; i++) {
    if (fabs(chi_sum) < FLT_EPSILON)
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
 * @brief Set the Material's absorption cross-section for some energy group.
 * @details The absorption cross section is computed by OpenMOC when needed,
 *          but the user can set it manually, to replace absorption by another
 *          cross section as absorption is not needed when determining Keff from
 *          fission rates.
 * @param xs the absorption cross-section
 * @param group the energy group
 */
void Material::setSigmaAByGroup(double xs, int group) {

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR, "Unable to set sigma_a for group %d for Material "
               "%d which contains %d energy groups", group, _id, _num_groups);

  if (_sigma_a == NULL)
    _sigma_a = new FP_PRECISION[_num_groups];

  _sigma_a[group-1] = FP_PRECISION(xs);
}


/**
 * @brief Builds the fission matrix from chi and the fission cross-section.
 * @details The fission matrix is constructed as the outer product of the
 *          chi and fission cross-section vectors. This routine is intended
 *          for internal use and is called by the Solver at runtime.
 */
void Material::buildFissionMatrix() {

  if (_num_groups == 0)
    log_printf(ERROR, "Unable to build Material %d's fission matrix "
           "since the number of energy groups has not been set", _id);
  else if (_nu_sigma_f == NULL)
    log_printf(ERROR, "Unable to build Material %d's fission matrix "
              "since its nu-fission cross-section has not been set", _id);
  else if (_chi == NULL)
    log_printf(ERROR, "Unable to build Material %d's fission matrix "
               "since its chi spectrum has not been set", _id);

  /* Deallocate memory for old fission matrix if needed */
  if (_fiss_matrix != NULL)
    delete [] _fiss_matrix;

  _fiss_matrix = new FP_PRECISION[_num_groups*_num_groups];

  /* Compute vector outer product of chi and the fission cross-section */
  for (int G=0; G < _num_groups; G++) {
    for (int g=0; g < _num_groups; g++)
      _fiss_matrix[G*_num_groups+g] = _chi[G] * _nu_sigma_f[g];
  }
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

  if (_fiss_matrix == NULL)
    buildFissionMatrix();

  _num_vector_groups = _num_groups / VEC_LENGTH + (_num_groups % VEC_LENGTH != 0);

  /* Allocate memory for the new aligned xs data */
  int size = _num_vector_groups * VEC_LENGTH * sizeof(FP_PRECISION);

  /* Allocate word-aligned memory for cross-section data arrays */
  FP_PRECISION* new_sigma_t = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_sigma_f = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_nu_sigma_f=(FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_chi = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  //FIXME Sigma_a is not copied

  /* The fission and scattering matrices will be the number of vector
   * groups wide (SIMD) and the actual number of groups long since
   * instructions are not SIMD in this dimension */

  size = _num_vector_groups * VEC_LENGTH * _num_vector_groups;
  size *= VEC_LENGTH * sizeof(FP_PRECISION);
  FP_PRECISION* new_fiss_matrix = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  FP_PRECISION* new_sigma_s = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

  /* Initialize data structures to ones for sigma_t since it is used to
   * divide the source in the solver, and zeroes for everything else */
  size = _num_vector_groups * VEC_LENGTH * sizeof(FP_PRECISION);
  for (int i=0; i < _num_vector_groups * VEC_LENGTH; i++) {
    new_sigma_t[i] = 1.0;
    new_sigma_f[i] = 0.0;
    new_nu_sigma_f[i] = 0.0;
    new_chi[i] = 0.0;
  }

  size *= _num_vector_groups * VEC_LENGTH;
  memset(new_fiss_matrix, 0, size);
  memset(new_sigma_s, 0, size);

  /* Copy materials data from unaligned arrays into new aligned arrays */
  size = _num_groups * sizeof(FP_PRECISION);
  memcpy(new_sigma_t, _sigma_t, size);
  memcpy(new_sigma_f, _sigma_f, size);
  memcpy(new_nu_sigma_f, _nu_sigma_f, size);
  memcpy(new_chi, _chi, size);

  for (int e=0; e < _num_groups; e++) {
    memcpy(new_fiss_matrix, _fiss_matrix, size);
    memcpy(new_sigma_s, _sigma_s, size);
    new_fiss_matrix += _num_vector_groups * VEC_LENGTH;
    new_sigma_s += _num_vector_groups * VEC_LENGTH;
    _fiss_matrix += _num_groups;
    _sigma_s += _num_groups;
  }

  _fiss_matrix -= _num_groups * _num_groups;
  _sigma_s -= _num_groups * _num_groups;

  /* Reset the new fission / scattering matrix array pointers */
  new_fiss_matrix -= _num_vector_groups * VEC_LENGTH * _num_groups;
  new_sigma_s -= _num_vector_groups * VEC_LENGTH * _num_groups;

  /* Delete the old unaligned arrays */
  delete [] _sigma_t;
  delete [] _sigma_f;
  delete [] _nu_sigma_f;
  delete [] _chi;
  delete [] _fiss_matrix;
  delete [] _sigma_s;

  /* Set the material's array pointers to the new aligned arrays */
  _sigma_t = new_sigma_t;
  _sigma_f = new_sigma_f;
  _nu_sigma_f = new_nu_sigma_f;
  _chi = new_chi;
  _fiss_matrix = new_fiss_matrix;
  _sigma_s = new_sigma_s;

  _data_aligned = true;

  return;
}


/**
 * @brief Transposes the scattering and fission matrices.
 * @details This routine is used by the Solver when performing
 *          adjoint flux caclulations.
 */
void Material::transposeProductionMatrices() {

  int num_groups;
  if (_data_aligned)
    num_groups = _num_vector_groups * VEC_LENGTH;
  else
    num_groups = _num_groups;

  /* Perform matrix transpose on each matrix that has been allocated */
  if (_fiss_matrix != NULL)
    matrix_transpose<FP_PRECISION>(_fiss_matrix, num_groups, num_groups);
  if (_sigma_s != NULL)
    matrix_transpose<FP_PRECISION>(_sigma_s, num_groups, num_groups);
}


/**
 * @brief Create a duplicate of the Material.
 * @return a pointer to the clone
 */
Material* Material::clone() {

  Material* clone = new Material(material_id(), _name);

  /* Set the number of groups if this Material's groups have been set */
  if (_num_groups > 0)
    clone->setNumEnergyGroups(_num_groups);

  for (int i=0; i < _num_groups; i++) {
    clone->setSigmaTByGroup((double)_sigma_t[i], i+1);
    clone->setSigmaFByGroup((double)_sigma_f[i], i+1);
    clone->setNuSigmaFByGroup((double)_nu_sigma_f[i], i+1);
    clone->setChiByGroup((double)_chi[i], i+1);
    if (_sigma_a != NULL)
      clone->setSigmaAByGroup((double)_sigma_a[i], i+1);

    for (int j=0; j < _num_groups; j++)
      clone->setSigmaSByGroup((double)getSigmaSByGroup(i+1,j+1), i+1, j+1);
  }

  if (_fiss_matrix != NULL)
    clone->buildFissionMatrix();

  return clone;
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

  if (_fiss_matrix != NULL) {
    string << "\n\t\tFiss. Matrix = \n\t\t";
    for (int G = 0; G < _num_groups; G++) {
      for (int g = 0; g < _num_groups; g++)
        string << _fiss_matrix[G+g*_num_groups] << "\t\t ";
      string << "\n\t\t";
    }
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
