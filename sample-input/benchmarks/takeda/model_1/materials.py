"""
This file writes all of the materials data (multi-group nuclear
cross-sections) for the Takeda
benchmark problem to an HDF5 file. The script uses the h5py Python package
to interact with the HDF5 file format. This may be a good example for those
wishing to write their nuclear data to an HDF5 file to import using the
OpenMOC 'materialize' Python module.
"""


# Create a Python dictionary to store Takeda multi-group cross-sections
dataset = {}

dataset['Energy Groups'] = 2
dataset['Materials'] = {}

takeda_materials = dataset['Materials']

###############################################################################
##################################   core   ###################################
###############################################################################

# Create a subdictionary for region 1 material data
takeda_materials['Core'] = {}

takeda_materials['Core']['Absorption XS'] = [0.00852709, 0.158196]

takeda_materials['Core']['Total XS'] = [0.223775, 1.03864]

takeda_materials['Core']['Scattering XS'] = [0.192423, 0.0228253, 0.00, 0.880439]

takeda_materials['Core']['Fission XS'] = [0.0004, 0.1]

takeda_materials['Core']['Nu Fission XS'] = [0.00909319, 0.290183]

takeda_materials['Core']['Chi'] = [1.0, 0.0]

###############################################################################
##################################   reflector  ###############################
###############################################################################

# Create a subdictionary for region 1 material data
takeda_materials['Reflector'] = {}

takeda_materials['Reflector']['Absorption XS'] = [0.000416392, 0.0202999]

takeda_materials['Reflector']['Total XS'] = [0.250367, 1.64482]

takeda_materials['Reflector']['Scattering XS'] = [0.193446, 0.0565042, 0.00, 1.62452]

takeda_materials['Reflector']['Fission XS'] = [0.0, 0.0]

takeda_materials['Reflector']['Nu Fission XS'] = [0.0, 0.0]

takeda_materials['Reflector']['Chi'] = [1.0, 0.0]

###############################################################################
################################   control rod    #############################
###############################################################################

# Create a subdictionary for region 1 material data
takeda_materials['Control Rod'] = {}

takeda_materials['Control Rod']['Absorption XS'] = [0.0174439, 0.182224]

takeda_materials['Control Rod']['Total XS'] = [0.0852325, 0.217460]

takeda_materials['Control Rod']['Scattering XS'] = [0.0677241, 0.0000645461, 0.00, 0.0352358]

takeda_materials['Control Rod']['Fission XS'] = [0.0, 0.0]

takeda_materials['Control Rod']['Nu Fission XS'] = [0.0, 0.0]

takeda_materials['Control Rod']['Chi'] = [1.0, 0.0]

###############################################################################
################################       Void       #############################
###############################################################################

# Create a subdictionary for region 1 material data
takeda_materials['Void'] = {}

takeda_materials['Void']['Absorption XS'] = [0.0000465132, 0.00132890]

takeda_materials['Void']['Total XS'] = [0.0128407, 0.0120676]

takeda_materials['Void']['Scattering XS'] = [0.01277, 0.0000240997, 0.00, 0.0107387]

takeda_materials['Void']['Fission XS'] = [0.0, 0.0]

takeda_materials['Void']['Nu Fission XS'] = [0.0, 0.0]

takeda_materials['Void']['Chi'] = [1.0, 0.0]
