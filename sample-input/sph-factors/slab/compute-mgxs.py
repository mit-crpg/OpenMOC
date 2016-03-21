import openmc
import openmc.mgxs

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.100.h5')
su = openmc.Summary('summary.h5')
sp.link_with_summary(su)

# Instantiate a 1-group EnergyGroups object
groups = openmc.mgxs.EnergyGroups()
groups.group_edges = [0., 20.]

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
