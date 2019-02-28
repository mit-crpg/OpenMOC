import openmc
import openmc.mgxs

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.100.h5')

# Instantiate a 1-group EnergyGroups object
groups = openmc.mgxs.EnergyGroups()
groups.group_edges = [0., 20e6]

# Initialize an 2-group MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(sp.summary.geometry)
mgxs_lib.energy_groups = groups
mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'cell'
mgxs_lib.correction = None
mgxs_lib.build_library()
mgxs_lib.load_from_statepoint(sp)

# Store library and its MGXS objects in a pickled binary file
mgxs_lib.dump_to_file(filename='mgxs', directory='.')
