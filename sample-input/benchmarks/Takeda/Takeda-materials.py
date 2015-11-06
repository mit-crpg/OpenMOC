"""
This file writes all of the materials data (multi-group nuclear 
cross-sections) for the Takeda
benchmark problem to an HDF5 file. The script uses the h5py Python package
to interact with the HDF5 file format. This may be a good example for those
wishing ot write their nuclear data to an HDF5 file to import using the
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
takeda_materials['core'] = {}

takeda_materials['core']['Absorption XS'] = [0.00852709, 0.158196]

takeda_materials['core']['Total XS'] = [0.223775, 1.03864]

takeda_materials['core']['Scattering XS'] = [0.192423, 0.0228253, 0.00, 0.880439]

takeda_materials['core']['Fission XS'] = [0.0004, 0.1]

takeda_materials['core']['Nu Fission XS'] = [0.00909319, 0.290183]

takeda_materials['core']['Chi'] = [1.0, 0.0]

###############################################################################
##################################   reflector  ###############################
###############################################################################

# Create a subdictionary for region 1 material data
takeda_materials['reflector'] = {}

takeda_materials['reflector']['Absorption XS'] = [0.000416392, 0.0202999]

takeda_materials['reflector']['Total XS'] = [0.250367, 1.64482]

takeda_materials['reflector']['Scattering XS'] = [0.193446, 0.0565042, 0.00, 1.62452]

takeda_materials['reflector']['Fission XS'] = [0.0, 0.0]

takeda_materials['reflector']['Nu Fission XS'] = [0.0, 0.0]

takeda_materials['reflector']['Chi'] = [1.0, 0.0]

###############################################################################
################################   control rod    #############################
###############################################################################

# Create a subdictionary for region 1 material data
takeda_materials['control_rod'] = {}

takeda_materials['control_rod']['Absorption XS'] = [0.0174439, 0.182224]

takeda_materials['control_rod']['Total XS'] = [0.0852325, 0.217460]

takeda_materials['control_rod']['Scattering XS'] = [0.0677241, 0.0000645461, 0.00, 0.0352358]

takeda_materials['control_rod']['Fission XS'] = [0.0, 0.0]

takeda_materials['control_rod']['Nu Fission XS'] = [0.0, 0.0]

takeda_materials['control_rod']['Chi'] = [1.0, 0.0]

###############################################################################
################################       void       #############################
###############################################################################

# Create a subdictionary for region 1 material data
takeda_materials['void'] = {}

takeda_materials['void']['Absorption XS'] = [0.0000465132, 0.00132890]

takeda_materials['void']['Total XS'] = [0.0128407, 0.0120676]

takeda_materials['void']['Scattering XS'] = [0.01277, 0.0000240997, 0.00, 0.0107387]

takeda_materials['void']['Fission XS'] = [0.0, 0.0]

takeda_materials['void']['Nu Fission XS'] = [0.0, 0.0]

takeda_materials['void']['Chi'] = [1.0, 0.0]
