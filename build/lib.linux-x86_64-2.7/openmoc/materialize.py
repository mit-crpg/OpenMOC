##
# @file materialize.py
# @package openmoc.materialize
# @brief The materialize module provides utility functions to read and write
#        multi-group materials cross-section data from HDF5 binary file format.
# @author William Boyd (wboyd@mit.edu)
# @date April 23, 2013

from openmoc import *
from log import *


##
# @brief
# @param filename
# @return a list of materials for OpenMOC
def materialize(filename):

    xs_types = ['Total XS', 'Absorption XS', 'Scattering XS', \
                    'Fission XS', 'Nu Fission XS', 'Chi']
    materials = {}

    # Check that the filename is a string
    if not isinstance(filename, str):
        py_printf('ERROR', 'Unable to materialize using filename %s ' + \
                      'since it is not a string', str(filename))


    ############################################################################
    #                               HDF5 DATA FILES
    ############################################################################
    if filename.endswith('.hdf5'):

        import h5py
        import numpy as np

        # Create a h5py file handle for the file
        f = h5py.File(filename)

        # Check that the file has an 'energy groups' attribute
        if not 'Energy Groups' in f.attrs:
            py_printf('ERROR', 'Unable to materialize file %s since it ' + \
                          'does not contain an \'Energy Groups\' attribute', \
                          filename)
    
        num_groups = f.attrs['Energy Groups']

        # Check that the number of energy groups is an integer
        if not isinstance(num_groups, int):
            py_printf('ERROR', 'Unable to materialize file %s since the ' + \
                          'number of energy groups %s is not an integer', \
                          filename, str(num_groups))

        material_names = list(f)

        # Loop over each material and 
        for name in material_names:

            py_printf('INFO', 'Importing material %s', str(name))

            new_material = Material(material_id())
            new_material.setNumEnergyGroups(int(num_groups))

            # Retrieve xs data from file for this material
            sigma_t = f[name]['Total XS'][...]
            sigma_a = f[name]['Absorption XS'][...]
            sigma_s = f[name]['Scattering XS'][...]
            sigma_f = f[name]['Fission XS'][...]
            nu_sigma_f = f[name]['Nu Fission XS'][...]
            chi = f[name]['Chi'][...]

            # Load the cross-section data into the material object
            new_material.setSigmaT(sigma_t)
            new_material.setSigmaA(sigma_a)
            new_material.setSigmaS(sigma_s)
            new_material.setSigmaF(sigma_f)
            new_material.setNuSigmaF(nu_sigma_f)
            new_material.setChi(chi)

            # Add this material to the list
            materials[name] = new_material


    ############################################################################
    #                      PYTHON DICTIONARY DATA FILES
    ############################################################################
    elif filename.endswith('.py'):

        import imp
        data = imp.load_source(filename, filename).dataset

        # Check that the file has an 'energy groups' attribute
        if not 'Energy Groups' in data.keys():
            py_printf('ERROR', 'Unable to materialize file %s since it ' + \
                      'does not contain an \'Energy Groups\' attribute', \
                      filename)
    
        num_groups = data['Energy Groups']

        # Check that the number of energy groups is an integer
        if not isinstance(num_groups, int):
            py_printf('ERROR', 'Unable to materialize file %s since the ' + \
                          'number of energy groups %s is not an integer', \
                          filename, str(num_groups))

        data = data['Materials']
        material_names = data.keys()

        # Loop over each material and 
        for name in material_names:

            py_printf('INFO', 'Importing material %s', str(name))

            new_material = Material(material_id())
            new_material.setNumEnergyGroups(int(num_groups))

            # Retrieve xs data from file for this material
            sigma_t = data[name]['Total XS']
            sigma_a = data[name]['Absorption XS']
            sigma_s = data[name]['Scattering XS']
            sigma_f = data[name]['Fission XS']
            nu_sigma_f = data[name]['Nu Fission XS']
            chi = data[name]['Chi']

            # Setting the cross-sections for this material
            new_material.setSigmaT(sigma_t)
            new_material.setSigmaA(sigma_a)
            new_material.setSigmaS(sigma_s)
            new_material.setSigmaF(sigma_f)
            new_material.setNuSigmaF(nu_sigma_f)
            new_material.setChi(chi)

            # Add this material to the list
            materials[name] = new_material



    ############################################################################
    #                      UNSUPPORTED DATA FILE TYPES
    ############################################################################
    else:
        py_printf('ERROR', 'Unable to materialize using filename %s ' + \
                      'since it has an unkown extension. Supported ' + \
                      'extension types are .hdf5 and .py', filename)



    # Return the list of materials
    return materials
