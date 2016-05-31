import openmc
import openmc.mgxs


###############################################################################
#                 Exporting to OpenMC materials.xml File
###############################################################################

# Instantiate some Materials and register the appropriate Nuclides
uo2 = openmc.Material(name='UO2 Fuel')
uo2.set_density('g/cm3', 10.29769)
uo2.add_nuclide('U-235', 5.5815e-4)
uo2.add_nuclide('U-238', 2.2408e-2)
uo2.add_nuclide('O-16', 4.5829e-2)

helium = openmc.Material(name='Helium')
helium.set_density('g/cm3', 0.001598)
helium.add_nuclide('He-4', 2.4044e-4)

zircaloy = openmc.Material(name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_nuclide('O-16', 3.0743e-4)
zircaloy.add_nuclide('Fe-56', 1.3610e-4)
zircaloy.add_nuclide('Zr-90', 2.1827e-2)

borated_water = openmc.Material(name='Borated Water')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_nuclide('B-10', 8.0042e-6)
borated_water.add_nuclide('B-11', 3.2218e-5)
borated_water.add_nuclide('H-1', 4.9457e-2)
borated_water.add_nuclide('O-16', 2.4672e-2)
borated_water.add_s_alpha_beta('HH2O', '71t')

# Instantiate a MaterialsFile, register all Materials, and export to XML
materials_file = openmc.Materials([uo2, helium, zircaloy, borated_water])
materials_file.default_xs = '71c'
materials_file.make_isotropic_in_lab()
materials_file.export_to_xml()


###############################################################################
#                 Exporting to OpenMC geometry.xml File
###############################################################################

# Instantiate ZCylinder surfaces
fuel_or = openmc.ZCylinder(x0=0, y0=0, R=0.39218, name='Fuel OR')
clad_ir = openmc.ZCylinder(x0=0, y0=0, R=0.40005, name='Clad IR')
clad_or = openmc.ZCylinder(x0=0, y0=0, R=0.45720, name='Clad OR')
min_x = openmc.XPlane(x0=-0.62992, name='min x')
max_x = openmc.XPlane(x0=+0.62992, name='max x')
min_y = openmc.YPlane(y0=-0.62992, name='min y')
max_y = openmc.YPlane(y0=+0.62992, name='max y')
min_z = openmc.ZPlane(z0=-0.62992, name='min z')
max_z = openmc.ZPlane(z0=+0.62992, name='max z')

min_x.boundary_type = 'reflective'
max_x.boundary_type = 'reflective'
min_y.boundary_type = 'reflective'
max_y.boundary_type = 'reflective'
min_z.boundary_type = 'reflective'
max_z.boundary_type = 'reflective'

# Instantiate Cells
fuel = openmc.Cell(name='fuel')
gap = openmc.Cell(name='gap')
clad = openmc.Cell(name='clad')
water = openmc.Cell(name='water')
root_cell = openmc.Cell(name='root')

# Use surface half-spaces to define regions
fuel.region = -fuel_or
gap.region = +fuel_or & -clad_ir
clad.region = +clad_ir & -clad_or
water.region = +clad_or
root_cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

# Instantiate Universe
pin = openmc.Universe(name='pin cell universe')
root_univ = openmc.Universe(universe_id=0, name='root universe')

# Register fills with Cells
fuel.fill = uo2
gap.fill = helium
clad.fill = zircaloy
water.fill = borated_water
root_cell.fill = pin

# Register Cells with Universe
pin.add_cells([fuel, gap, clad, water])
root_univ.add_cell(root_cell)

# Instantiate a Geometry and register the root Universe
geometry = openmc.Geometry()
geometry.root_universe = root_univ
geometry.export_to_xml()


###############################################################################
#                   Exporting to OpenMC settings.xml File
###############################################################################

# Construct uniform initial source distribution over fissionable zones
lower_left = [-0.62992, -0.62992, -0.62992]
upper_right = [+0.62992, +0.62992, +0.62992]
source = openmc.source.Source(space=openmc.stats.Box(lower_left, upper_right))
source.only_fissionable = True

# Instantiate a Settings collection
settings_file = openmc.Settings()
settings_file.batches = 100
settings_file.inactive = 10
settings_file.particles = 100000
settings_file.output = {'tallies': False, 'summary': True}
settings_file.source = source
settings_file.sourcepoint_write = False
settings_file.export_to_xml()


###############################################################################
#                        Create OpenMC MGXS Library
###############################################################################

# Instantiate a 16-group EnergyGroups object
groups = openmc.mgxs.EnergyGroups()
groups.group_edges = [0., 0.03e-6, 0.058e-6, 0.14e-6, 0.28e-6, 0.35e-6, 
                      0.625e-6, 0.85e-6, 0.972e-6, 1.02e-6, 1.097e-6, 
                      1.15e-6, 1.3e-6, 4.e-6, 5.53e-3, 821.e-3, 20.]

# Initialize an MGXS Library for OpenMOC
mgxs_lib = openmc.mgxs.Library(geometry)
mgxs_lib.energy_groups = groups
mgxs_lib.mgxs_types = ['total', 'nu-fission', 'nu-scatter matrix', 'chi']
mgxs_lib.domain_type = 'cell'
mgxs_lib.correction = None
mgxs_lib.build_library()

# Create a "tallies.xml" file for the MGXS Library
tallies_file = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
tallies_file.export_to_xml()


###############################################################################
#                         Run OpenMC Simulation
###############################################################################

# Run OpenMC
openmc.run(output=True, mpi_procs=4, threads=1)
