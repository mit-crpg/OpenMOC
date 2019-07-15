import numpy as np
import openmc
import openmc.mgxs
import openmc.openmoc_compatible


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate some Materials and register the appropriate Nuclides
fuel = openmc.Material(name='3.1 w/o enriched UO2')
fuel.set_density('sum')
fuel.add_nuclide('U235', 7.18132E-4,)
fuel.add_nuclide('U238', 2.21546E-2)
fuel.add_nuclide('O16', 4.57642E-2)

clad = openmc.Material(name='Cladding')
clad.set_density('sum')
clad.add_element('Zr', 4.98349e-2)

water = openmc.Material(name='Borated Water')
water.set_density('sum')
water.add_nuclide('H1', 4.41459E-2)
water.add_nuclide('O16', 2.20729E-2)
water.add_nuclide('B10', 9.52537E-6)
water.add_nuclide('B11', 3.83408E-5)
water.add_s_alpha_beta(name='HH2O')

# Instantiate a MaterialsFile, add Materials
materials_file = openmc.Materials([fuel, clad, water])
materials_file.make_isotropic_in_lab()
materials_file.default_xs = '71c'

# Export to "materials.xml"
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Create bounding surfaces
min_x = openmc.XPlane(boundary_type='reflective', x0=0.0)
max_x = openmc.XPlane(boundary_type='reflective', x0=5.0)
min_y = openmc.YPlane(boundary_type='reflective', y0=0.0)
max_y = openmc.YPlane(boundary_type='reflective', y0=10.0)
min_z = openmc.ZPlane(boundary_type='reflective', z0=0.0)
max_z = openmc.ZPlane(boundary_type='reflective', z0=10.0)

# Create material interfacial surfaces
left = openmc.XPlane(surface_id=1, x0=2.0)
right = openmc.XPlane(surface_id=2, x0=2.4)

# Create a Universe to encapsulate the 1D slab
slab_universe = openmc.Universe(name='1D slab')

# Create fuel Cells
slab_width = 2. / 40.
for i in range(40):
    fuel_cell = openmc.Cell(name='fuel {}'.format(i))
    fuel_cell.fill = fuel
    fuel_cell.add_surface(
        halfspace=+1, surface=openmc.XPlane(x0=slab_width*i))
    fuel_cell.add_surface(
        halfspace=-1, surface=openmc.XPlane(x0=slab_width * (i+1)))
    slab_universe.add_cell(fuel_cell)

# Create fuel Cells
slab_width = 0.4 / 20.
for i in range(20):
    clad_cell = openmc.Cell(name='clad {}'.format(i))
    clad_cell.fill = clad
    clad_cell.add_surface(
        halfspace=+1, surface=openmc.XPlane(x0=2.+slab_width*i))
    clad_cell.add_surface(
        halfspace=-1, surface=openmc.XPlane(x0=2.+slab_width * (i+1)))
    slab_universe.add_cell(clad_cell)

# Create water Cells
slab_width = 2.6 / 40.
for i in range(40):
    water_cell = openmc.Cell(name='water {}'.format(i))
    water_cell.fill = water
    water_cell.add_surface(
        halfspace=+1, surface=openmc.XPlane(x0=2.4+slab_width * i))
    water_cell.add_surface(
        halfspace=-1, surface=openmc.XPlane(x0=2.4+slab_width * (i+1)))
    slab_universe.add_cell(water_cell)

# Create root Cell
root_cell = openmc.Cell(name='root cell')
root_cell.fill = slab_universe

# Add boundary planes
root_cell.add_surface(halfspace=+1, surface=min_x)
root_cell.add_surface(halfspace=-1, surface=max_x)
root_cell.add_surface(halfspace=+1, surface=min_y)
root_cell.add_surface(halfspace=-1, surface=max_y)
root_cell.add_surface(halfspace=+1, surface=min_z)
root_cell.add_surface(halfspace=-1, surface=max_z)

# Create root Universe
root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# Create Geometry and set root Universe
openmc_geometry = openmc.Geometry()
openmc_geometry.root_universe = root_universe
openmc_geometry.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Construct uniform initial source distribution over fissionable zones
lower_left = [0., 0., 0.]
upper_right = [5., 5., 10.]
source = openmc.source.Source(space=openmc.stats.Box(lower_left, upper_right))
source.only_fissionable = True

# Instantiate a Settings collection
settings_file = openmc.Settings()
settings_file.batches = 100
settings_file.inactive = 10
settings_file.particles = 10000
settings_file.output = {'tallies': False, 'summary': True}
settings_file.source = source
settings_file.sourcepoint_write = False
settings_file.export_to_xml()


###############################################################################
#                     Exporting to OpenMC plots.xml File
###############################################################################

# Instantiate a Plot
plot = openmc.Plot()
plot.origin = [2.5, 2.5, 5.]
plot.width = [5., 5.]
plot.pixels = [250, 250]
plot.color = 'cell'

# Instantiate a Plots collection and export to "plots.xml"
plot_file = openmc.Plots([plot])
plot_file.export_to_xml()


###############################################################################
#                        Create OpenMC MGXS Library
###############################################################################

# Instantiate a 1-group EnergyGroups object
groups = openmc.mgxs.EnergyGroups()
groups.group_edges = [0., 20e6]

# Initialize an MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(openmc_geometry, by_nuclide=False)
mgxs_lib.energy_groups = groups
mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'cell'
mgxs_lib.correction = None
mgxs_lib.domains = openmc_geometry.get_all_material_cells().values()
mgxs_lib.build_library()

# Create a "tallies.xml" file for the MGXS Library
tallies_file = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
tallies_file.export_to_xml()


###############################################################################
#                         Run OpenMC Simulation
###############################################################################

# Run OpenMC
openmc.run()
