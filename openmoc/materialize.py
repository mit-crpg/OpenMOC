##
# @file materialize.py
# @package openmoc.materialize
# @brief The materialize module provides utility functions to read and write
#        multi-group materials cross-section data from HDF5 binary file format.
# @author William Boyd (wboyd@mit.edu)
# @date October 9, 2015

import sys
import os

import openmoc

# For Python 2.X.X
if sys.version_info[0] == 2:
    from log import py_printf
# For Python 3.X.X
else:
    from openmoc.log import py_printf



def _get_domain(domains, domain_spec):

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


##
# @brief This routine processes an HDF5 file of multi-group cross section data.
# @details The routine instantiates Material objects with that data and
#          returns a dictionary of each Material object keyed by its ID.
#          An OpenMOC Geometry may optionally be given and the routine will
#          directly insert the multi-group cross sections into each Material in
#          the Geometry. If a Geometry is passed in, Material objects from the
#          Geometry will be used in place of those instantiated by this routine.
#          A second optional parameter for the domain types may be used to
#          define whether multi-group cross sections are tabulated by material
#          or cell in the HDF5 binary file.
# @param filename the file of nuclear cross-sections for each Material
# @param geometry an optional geometry populated with materials, cells, etc.
# @param domain_type the domain type ('material' or 'cell') upon which the
#        cross sections are defined (default is 'material')
# @param suffice
# @return a dictionary of Material objects keyed by ID
def load_from_hdf5(filename='mgxs.h5', directory='mgxs',
                   geometry=None, domain_type='material', suffix=''):

    # Create a h5py file handle for the file
    import h5py
    filename = os.path.join(directory, filename)
    f = h5py.File(filename, 'r')

    # Check that the file has an 'energy groups' attribute
    if '# groups' not in f.attrs:
        py_printf('ERROR', 'Unable to materialize file "%s" since it does ' +
                           'not contain an \'# groups\' attribute', filename)

    if domain_type not in f.keys():
        py_printf('ERROR', 'Unable to materialize file "%s" since it does ' +
                           'not contain domain type "%s"', filename, domain_type)

    if geometry and 'openmoc.Geometry' not in str(type(geometry)):
        py_printf('ERROR', 'Unable to materialize file "%s"  with "%s" which '
                           'is not an OpenMOC Geometry', filename, str(geometry))

    materials = {}
    num_groups = int(f.attrs['# groups'])

    if geometry:
        if domain_type == 'material':
            domains = geometry.getAllMaterials()
        elif domain_type == 'cell':
            domains = geometry.getAllMaterialCells()
        else:
            py_printf('ERROR', 'Domain type "%s" is not supported', domain_type)

    for domain_spec in f[domain_type]:

        py_printf('INFO', 'Importing material for %s "%s"',
                          domain_type, str(domain_spec))

        domain_group = f[domain_type][domain_spec]

        # If the domain spec is an integer, it is an ID
        if domain_spec.isdigit():
            domain_spec = int(domain_spec)

        # If using an OpenMOC Geometry, extract a Material from it
        if geometry:
            if domain_type == 'material':
                material = _get_domain(domains, domain_spec)
            elif domain_type == 'cell':
                material = _get_domain(domains, domain_spec).getFillMaterial()

        # Instantiate a new Material with an appropriate ID or name
        else:
            if isinstance(domain_spec, int):
                material = openmoc.Material(id=domain_spec)
            else:
                material = openmoc.Material(name=domain_spec)

        # Add material to the collection
        materials[domain_spec] = material
        material.setNumEnergyGroups(num_groups)

        # Total cross section
        if 'total' in domain_group:
            material.setSigmaT(domain_group['total/' + suffix][...])
        elif 'transport' in domain_group:
            material.setSigmaT(domain_group['transport/' + suffix][...])
        else:
            py_printf('WARNING', 'No "total" or "transport" MGXS found for'
                                 '"%s %s"', domain_type, domain_spec)

        # Nu-Fission cross section
        if 'nu-fission' in domain_group:
            material.setNuSigmaF(domain_group['nu-fission/' + suffix][...])
        else:
            py_printf('WARNING', 'No "nu-fission" MGXS found for'
                                 '"%s %s"', domain_type, domain_spec)

        # Scattering matrix cross section
        if 'scatter matrix' in domain_group:
            material.setSigmaS(domain_group['scatter matrix/' + suffix][...])
        elif 'nu-scatter matrix' in domain_group:
            material.setSigmaS(domain_group['nu-scatter matrix/' + suffix][...])
        else:
            py_printf('WARNING', 'No "scatter matrix" or "nu-scatter matrix" '
                                 'found for "%s %s"', domain_type, domain_spec)

        # Scattering matrix cross section
        if 'chi' in domain_group:
            material.setChi(domain_group['chi/' + suffix][...])
        else:
            py_printf('WARNING', 'No "chi" MGXS found for "%s %s"',
                                 domain_type, domain_spec)

    # Return collection of materials
    return materials