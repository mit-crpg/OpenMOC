##
# @file materialize.py
# @package openmoc.materialize
# @brief The materialize module provides utility functions to read and write
#        multi-group materials cross section data from HDF5 binary file format.
# @author William Boyd (wboyd@mit.edu)
# @date October 9, 2015

import sys
import os

import numpy as np

import openmoc

# For Python 2.X.X
if sys.version_info[0] == 2:
    from log import py_printf
# For Python 3.X.X
else:
    from openmoc.log import py_printf


def _get_domain(domains, domain_spec):
    """A helper routine to find materials/cells for load_from_hdf5(...)"""

    # If domain_spec is an integer, it must be a domain ID
    if isinstance(domain_spec, int) and domain_spec in domains:
        return domains[domain_spec]

    # If domain_spec is a string, it must be a domain name
    elif isinstance(domain_spec, str):
        for domain_id, domain in domains.items():
            if domain_spec == domain.getName():
                return domain

    # If domain could not be found
    return None


def _get_numpy_array(hdf5_group, key, suffix):
    """A helper routine to ensure that MGXS data is a proper NumPy array"""

    sigma = np.array(hdf5_group['{}/'.format(key) + suffix][...])
    sigma = np.atleast_1d(sigma)
    sigma = sigma.flatten()
    return sigma


##
# @brief This routine loads an HDF5 file of multi-group cross section data.
# @details The routine instantiates Material objects with multi-group cross
#          section data and returns a dictionary of each Material object keyed
#          by its ID. An OpenMOC Geometry may optionally be given and the
#          routine will directly insert the multi-group cross sections into each
#          Material in the Geometry. If a Geometry is passed in, Material
#          objects from the Geometry will be used in place of those instantiated
#          by this routine. An optional parameter for the domain types may be
#          used to define whether multi-group cross sections are tabulated by
#          material or cell in the HDF5 binary file.
# @param filename filename for cross sections HDF5 file (default is 'mgxs.h5')
# @param directory directory for cross sections HDF5 file (default is 'mgxs')
# @param geometry an optional geometry populated with materials, cells, etc.
# @param domain_type the domain type ('material' or 'cell') upon which the
#        cross sections are defined (default is 'material')
# @param suffix an optional string suffix to index the HDF5 file beyond the
#        assumed domain_type/domain_id/mgxs_type group sequence (default is '')
# @return a dictionary of Material objects keyed by ID
def load_from_hdf5(filename='mgxs.h5', directory='mgxs',
                   geometry=None, domain_type='material', suffix=''):

    # Create a h5py file handle for the file
    import h5py
    filename = os.path.join(directory, filename)
    f = h5py.File(filename, 'r')

    # Check that the file has an 'energy groups' attribute
    if '# groups' not in f.attrs:
        py_printf('ERROR', 'Unable to load HDF5 file "%s" since it does '
                           'not contain an \'# groups\' attribute', filename)

    if domain_type not in f.keys():
        py_printf('ERROR', 'Unable to load HDF5 file "%s" since it does '
                           'not contain domain type "%s"', filename, domain_type)

    if geometry and 'openmoc.Geometry' not in str(type(geometry)):
        py_printf('ERROR', 'Unable to load HDF5 file "%s" for "%s" which is not '
                           'an OpenMOC Geometry', filename, str(type(geometry)))

    # Instantiate dictionary to hold Materials to return to user
    materials = {}
    num_groups = int(f.attrs['# groups'])

    # If a Geometry was passed in, extract all cells or materials from it
    if geometry:
        if domain_type == 'material':
            domains = geometry.getAllMaterials()
        elif domain_type == 'cell':
            domains = geometry.getAllMaterialCells()
        else:
            py_printf('ERROR', 'Domain type "%s" is not supported', domain_type)

    # Iterate over all domains (e.g., materials or cells) in the HDF5 file
    for domain_spec in f[domain_type]:

        py_printf('INFO', 'Importing cross sections for %s "%s"',
                          domain_type, str(domain_spec))

        # Create shortcut to HDF5 group for this domain
        domain_group = f[domain_type][domain_spec]

        # If domain_spec is an integer, it is an ID; otherwise a string name
        if domain_spec.isdigit():
            domain_spec = int(domain_spec)
        else:
            domain_spec = str(domain_spec)

        # If using an OpenMOC Geometry, extract a Material from it
        if geometry:

            if domain_type == 'material':
                material = _get_domain(domains, domain_spec)

            elif domain_type == 'cell':
                cell = _get_domain(domains, domain_spec)
                material = cell.getFillMaterial()

                # If the user filled the Cell with a Material, clone it
                if material != None:
                    material = material.clone()

                # If the Cell does not contain a Material, creat one for it
                else:
                    if isinstance(domain_spec, int):
                        material = openmoc.Material(id=domain_spec)
                    else:
                        material = openmoc.Material(name=domain_spec)

                # Fill the Cell with the new Material
                cell.setFill(material)

        # If not Geometry, instantiate a new Material with the ID/name
        else:
            if isinstance(domain_spec, int):
                material = openmoc.Material(id=domain_spec)
            else:
                material = openmoc.Material(name=domain_spec)

        # Add material to the collection
        materials[domain_spec] = material
        material.setNumEnergyGroups(num_groups)

        # Search for the total/transport cross section
        if 'transport' in domain_group:
            sigma = _get_numpy_array(domain_group, 'transport', suffix)
            material.setSigmaT(sigma)
        elif 'total' in domain_group:
            sigma = _get_numpy_array(domain_group, 'total', suffix)
            material.setSigmaT(sigma)
        else:
            py_printf('WARNING', 'No "total" or "transport" MGXS found for'
                                 '"%s %s"', domain_type, domain_spec)

        # Search for the fission production cross section
        if 'nu-fission' in domain_group:
            sigma = _get_numpy_array(domain_group, 'nu-fission', suffix)
            material.setNuSigmaF(sigma)
        else:
            py_printf('WARNING', 'No "nu-fission" MGXS found for'
                                 '"%s %s"', domain_type, domain_spec)

        # Search for the scattering matrix cross section
        if 'nu-scatter matrix' in domain_group:
            sigma = _get_numpy_array(domain_group, 'nu-scatter matrix', suffix)
            material.setSigmaS(sigma)
        elif 'scatter matrix' in domain_group:
            sigma = _get_numpy_array(domain_group, 'scatter matrix', suffix)
            material.setSigmaS(sigma)
        else:
            py_printf('WARNING', 'No "scatter matrix" or "nu-scatter matrix" '
                                 'found for "%s %s"', domain_type, domain_spec)

        # Search for chi (fission spectrum)
        if 'chi' in domain_group:
            sigma = _get_numpy_array(domain_group, 'chi', suffix)
            material.setChi(sigma)
        else:
            py_printf('WARNING', 'No "chi" MGXS found for "%s %s"',
                                 domain_type, domain_spec)

        # Search for optional cross sections
        if 'absorption' in domain_group:
            sigma = _get_numpy_array(domain_group, 'absorption', suffix)
            material.setSigmaA(sigma)
        if 'fission' in domain_group:
            sigma = _get_numpy_array(domain_group, 'fission', suffix)
            material.setSigmaF(sigma)
        if 'diffusion' in domain_group:
            sigma = _get_numpy_array(domain_group, 'diffusion', suffix)
            material.setDifCoef(sigma)
        if 'buckling' in domain_group:
            sigma = _get_numpy_array(domain_group, 'buckling', suffix)
            material.setBuckling(sigma)

    # Return collection of materials
    return materials



##
# @brief This routine loads an OpenMC Library of multi-group cross section data.
# @details The routine instantiates Material objects with multi-group cross
#          section data and returns a dictionary of each Material object keyed
#          by its ID. An OpenMOC Geometry may optionally be given and the
#          routine will directly insert the multi-group cross sections into each
#          Material in the Geometry. If a Geometry is passed in, Material
#          objects from the Geometry will be used in place of those instantiated
#          by this routine.
# @param mgxs_lib an openmc.mgxs.Library object with cross section data
# @param geometry an optional geometry populated with materials, cells, etc.
# @return a dictionary of Material objects keyed by ID
def load_openmc_mgxs_lib(mgxs_lib, geometry=None):

    # Attempt to import openmc
    try:
        import openmc
    except ImportError:
        py_printf('ERROR', 'The OpenMC code must be installed on your system')

    if not isinstance(mgxs_lib, openmc.mgxs.Library):
        py_printf('ERROR', 'Unable to load cross sections from %s which is not '
                           'an openmc.mgxs.Library object', str(type(mgxs_lib)))

    if geometry and 'openmoc.Geometry' not in str(type(geometry)):
        py_printf('ERROR', 'Unable to load cross sections for "%s" which '
                           'is not an OpenMOC Geometry', str(type(geometry)))

    # Instantiate dictionary to hold Materials to return to user
    materials = {}
    num_groups = mgxs_lib.num_groups
    domain_type = mgxs_lib.domain_type

    # If a Geometry was passed in, extract all cells or materials from it
    if geometry:
        if domain_type == 'material':
            domains = geometry.getAllMaterials()
        elif domain_type == 'cell':
            domains = geometry.getAllMaterialCells()
        else:
            py_printf('ERROR', 'Unable to load a cross sections library with '
                               'domain type %s', mgxs_lib.domain_type)

    # Iterate over all domains (e.g., materials or cells) in the HDF5 file
    for domain in mgxs_lib.domains:

        py_printf('INFO', 'Importing cross sections for %s "%d"',
                          domain_type, domain.id)

        # If using an OpenMOC Geometry, extract a Material from it
        if geometry:

            if domain_type == 'material':
                material = _get_domain(domains, domain.id)

            elif domain_type == 'cell':
                print(domain, domains)
                cell = _get_domain(domains, domain.id)
                print(cell)
                material = cell.getFillMaterial()

                # If the user filled the Cell with a Material, clone it
                if material != None:
                    material = material.clone()

                # If the Cell does not contain a Material, create one for it
                else:
                    material = openmoc.Material(id=domain)

                # Fill the Cell with the new Material
                cell.setFill(material)

        # If not Geometry, instantiate a new Material with the ID/name
        else:
            material = openmoc.Material(id=domain.id)

        # Add material to the collection
        materials[domain.id] = material
        material.setNumEnergyGroups(num_groups)

        # Search for the total/transport cross section
        if 'transport' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'transport')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaT(sigma)
        elif 'total' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'total')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaT(sigma)
        else:
            py_printf('WARNING', 'No "total" or "transport" MGXS found for'
                                 '"%s %d"', domain_type, domain.id)

        # Search for the fission production cross section
        if 'nu-fission' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'nu-fission')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setNuSigmaF(sigma)
        else:
            py_printf('WARNING', 'No "nu-fission" MGXS found for'
                                 '"%s %d"', domain_type, domain.id)

        # Search for the scattering matrix cross section
        if 'nu-scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'nu-scatter matrix')
            sigma = mgxs.get_xs(nuclides='sum').flatten()
            material.setSigmaS(sigma)
        elif 'scatter matrix' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'scatter matrix').flatten()
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaS(sigma)
        else:
            py_printf('WARNING', 'No "scatter matrix" or "nu-scatter matrix" '
                                 'found for "%s %d"', domain_type, domain.id)

        # Search for chi (fission spectrum)
        if 'chi' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'chi')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setChi(sigma)
        else:
            py_printf('WARNING', 'No "chi" MGXS found for "%s %d"',
                                 domain_type, domain.id)

        # Search for optional cross sections
        if 'absorption' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'absorption')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaA(sigma)
        if 'fission' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'fission')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setSigmaF(sigma)
        if 'diffusion' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'diffusion')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setDifCoef(sigma)
        if 'buckling' in mgxs_lib.mgxs_types:
            mgxs = mgxs_lib.get_mgxs(domain, 'buckling')
            sigma = mgxs.get_xs(nuclides='sum')
            material.setBuckling(sigma)

    # Return collection of materials
    return materials
