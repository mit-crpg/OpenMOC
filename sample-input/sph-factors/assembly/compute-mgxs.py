import copy
import numpy as np
import openmc.mgxs
from infermc.materialize import differentiate_mgxs_lib


###############################################################################
#               Build MGXS Library on Distribcell Domains
###############################################################################

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.100.h5')
su = openmc.Summary('summary.h5')
sp.link_with_summary(su)

# Instantiate a 2-group EnergyGroups object
groups = openmc.mgxs.EnergyGroups()
groups.group_edges = np.array([0., 0.625e-6, 20.])

# Initialize an 2-group MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(su.openmc_geometry)
mgxs_lib.energy_groups = groups
mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'distribcell'
mgxs_lib.correction = None
mgxs_lib.build_library()
mgxs_lib.load_from_statepoint(sp)

# Store library and its MGXS objects in a pickled binary file
mgxs_lib.dump_to_file(filename='base', directory='mgxs')


###############################################################################
#                            Store MGXS Libraries
###############################################################################

## Cell-averaged MGXS library
# Get cell-averaged OpenMC MGXS library from distribcell-averaged library
cell_mgxs_lib = mgxs_lib.get_subdomain_avg_library()

# Store library and its MGXS objects in a pickled binary file
cell_mgxs_lib.dump_to_file(filename='cell-avg', directory='mgxs')


## Cell instance (region)-averaged MGXS library
# Differentiate pin cells by cell instance or "region"
opencg_geometry = copy.deepcopy(mgxs_lib.opencg_geometry)
regions_to_regions = np.arange(1, opencg_geometry.num_regions+1)
new_regions = opencg_geometry.differentiate(regions_to_regions)

# Differentiate OpenMC MGXS library for new region-differentiated geometry
region_mgxs_lib = differentiate_mgxs_lib(mgxs_lib, new_regions, opencg_geometry)

# Differentiate OpenMC MGXS library for new region-differentiated geometry
region_mgxs_lib = differentiate_mgxs_lib(mgxs_lib, new_regions, opencg_geometry)

# Store library and its MGXS objects in a pickled binary file
region_mgxs_lib.dump_to_file(filename='region-avg', directory='mgxs')
