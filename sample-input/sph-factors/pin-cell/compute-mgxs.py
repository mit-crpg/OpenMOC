import openmc.mgxs
from openmc.statepoint import StatePoint
from openmc.summary import Summary

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.100.h5')
su = openmc.Summary('summary.h5')
sp.link_with_summary(su)

# Instantiate a 16-group EnergyGroups object
groups = openmc.mgxs.EnergyGroups()
groups.group_edges = [0., 0.03e-6, 0.058e-6, 0.14e-6, 0.28e-6, 0.35e-6, 
                      0.625e-6, 0.85e-6, 0.972e-6, 1.02e-6, 1.097e-6, 
                      1.15e-6, 1.3e-6, 4.e-6, 5.53e-3, 821.e-3, 20.]

# Initialize an 2-group MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(su.openmc_geometry)
mgxs_lib.energy_groups = groups
mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'cell'
mgxs_lib.correction = None
mgxs_lib.build_library()
mgxs_lib.load_from_statepoint(sp)

# Store library and its MGXS objects in a pickled binary file
mgxs_lib.dump_to_file(filename='mgxs', directory='.')
