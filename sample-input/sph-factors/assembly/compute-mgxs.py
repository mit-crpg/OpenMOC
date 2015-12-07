import copy
import numpy as np
import openmc.mgxs
from openmc.statepoint import StatePoint
from openmc.summary import Summary
from infermc.materialize import differentiate_mgxs_lib

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
mgxs_lib.mgxs_types = ['transport', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'distribcell'
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


## OpenCG Local Neighor Symmetry-averaged MGXS library
# Differentiate fuel pin cells with OpenCG Local Neighbor Symmetry (LNS)
opencg_geometry = copy.deepcopy(mgxs_lib.opencg_geometry)
opencg_geometry.build_neighbors()
opencg_geometry.count_neighbors()

regions_to_lns = np.zeros(opencg_geometry.num_regions, dtype=np.int)
for region in range(opencg_geometry.num_regions):
    coords = opencg_geometry.find_region(region)
    regions_to_lns[region] = opencg_geometry.get_neighbors_hash(region)

# Differentiate OpenMC MGXS library for new region-differentiated geometry
new_regions = opencg_geometry.differentiate(regions_to_lns)
lns_mgxs_lib = \
    differentiate_mgxs_lib(mgxs_lib, new_regions, opencg_geometry)

# Store library and its MGXS objects in a pickled binary file
lns_mgxs_lib.dump_to_file(filename='lns-avg', directory='mgxs')


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
