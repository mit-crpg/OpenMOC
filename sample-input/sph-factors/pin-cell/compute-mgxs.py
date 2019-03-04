import openmc
import openmc.mgxs

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.100.h5')

# Instantiate a 16-group EnergyGroups object
groups = openmc.mgxs.EnergyGroups()
groups.group_edges = [0., 0.03, 0.058, 0.14, 0.28, 0.35,
                      0.625, 0.85, 0.972, 1.02, 1.097,
                      1.15, 1.3, 4., 5.53e3, 821.e3, 20.e6]

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
