import openmc
import openmc.mgxs
import os

# Parameters to use this file to generate libraries for all the MGXS tests
filenames = ['mgxs_isotropic',
            'mgxs_transport_corrected',
            'mgxs_consistent',
            'mgxs_consistent_nuscatter',
            'mgxs_materials',
            'mgxs_by_nuclide',
            'mgxs_angular_legendre']

for filename in filenames:
    ###############################################################################
    #                 Exporting to OpenMC materials.xml File
    ###############################################################################

    # Select scattering model
    scattering = 'isotropic'
    if 'corrected' in filename:
        scattering = 'anisotropic'

    # Instantiate some Materials and register the appropriate Nuclides
    uo2 = openmc.Material(name='UO2 Fuel')
    uo2.set_density('g/cm3', 10.29769)
    uo2.add_nuclide('U235', 5.5815e-4)
    uo2.add_nuclide('U238', 2.2408e-2)
    uo2.add_nuclide('O16', 4.5829e-2)

    helium = openmc.Material(name='Helium')
    helium.set_density('g/cm3', 0.001598)
    helium.add_nuclide('He4', 2.4044e-4)

    zircaloy = openmc.Material(name='Zircaloy 4')
    zircaloy.set_density('g/cm3', 6.55)
    zircaloy.add_nuclide('O16', 3.0743e-4)
    zircaloy.add_nuclide('Fe56', 1.3610e-4)
    zircaloy.add_nuclide('Zr90', 2.1827e-2)

    borated_water = openmc.Material(name='Borated Water')
    borated_water.set_density('g/cm3', 0.740582)
    borated_water.add_nuclide('B10', 8.0042e-6)
    borated_water.add_nuclide('B11', 3.2218e-5)
    borated_water.add_nuclide('H1', 4.9457e-2)
    borated_water.add_nuclide('O16', 2.4672e-2)
    borated_water.add_s_alpha_beta('HH2O')

    # Instantiate a MaterialsFile, register all Materials, and export to XML
    materials_file = openmc.Materials([uo2, helium, zircaloy, borated_water])
    materials_file.default_xs = '71c'
    if scattering == 'isotropic':
        materials_file.make_isotropic_in_lab()
    materials_file.export_to_xml()


    ###############################################################################
    #                 Exporting to OpenMC geometry.xml File
    ###############################################################################

    # Instantiate ZCylinder surfaces
    fuel_or = openmc.ZCylinder(x0=0, y0=0, r=0.39218, name='Fuel OR')
    clad_ir = openmc.ZCylinder(x0=0, y0=0, r=0.40005, name='Clad IR')
    clad_or = openmc.ZCylinder(x0=0, y0=0, r=0.45720, name='Clad OR')
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
    settings_file.batches = 50
    settings_file.inactive = 0
    settings_file.particles = 1000
    settings_file.output = {'tallies': False, 'summary': True}
    settings_file.source = source
    settings_file.sourcepoint_write = False
    settings_file.export_to_xml()


    ###############################################################################
    #                        Create OpenMC MGXS Library
    ###############################################################################

    # Instantiate a 16-group EnergyGroups object
    groups = openmc.mgxs.EnergyGroups()
    groups.group_edges = [0., 0.03, 0.058, 0.14, 0.28, 0.35,
                          0.625, 0.85, 0.972, 1.02, 1.097,
                          1.15, 1.3, 4., 5.53e3, 821.e3, 20.e6]

    # Select scattering matrix type
    scatter_matrix = 'scatter matrix'
    if 'nuscatter' in filename:
        scatter_matrix = 'nu-' + scatter_matrix
    if 'consistent' in filename:
        scatter_matrix = 'consistent ' + scatter_matrix

    # Initialize an MGXS Library for OpenMOC
    mgxs_lib = openmc.mgxs.Library(geometry)
    mgxs_lib.energy_groups = groups
    if scattering == 'isotropic':
        mgxs_lib.mgxs_types = ['total', 'nu-fission', scatter_matrix, 'chi']
        mgxs_lib.correction = None
    elif scattering == 'anisotropic':
        mgxs_lib.mgxs_types = ['nu-transport', 'nu-fission', scatter_matrix, 'chi']
        mgxs_lib.correction = 'P0'
    mgxs_lib.domain_type = 'cell'
    if 'materials' in filename:
        mgxs_lib.domain_type = 'material'
    if 'nuclide' in filename:
        mgxs_lib.by_nuclide = True
    if 'legendre' in filename:
        mgxs_lib.legendre_order = 3
    mgxs_lib.build_library()

    # Create a "tallies.xml" file for the MGXS Library
    tallies_file = openmc.Tallies()
    mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
    tallies_file.export_to_xml()


    ###############################################################################
    #                  Run OpenMC Simulation and generate library
    ###############################################################################

    # Run OpenMC
    openmc.run()

    # Load the library from the last statepoint file
    statepoint = 'statepoint.50.h5'
    sp = openmc.StatePoint(statepoint)
    mgxs_lib.load_from_statepoint(sp)

    # Store library and its MGXS objects in a pickled binary file
    mgxs_lib.dump_to_file(filename=filename, directory='.')

    # Delete statepoint and summary
    os.system("rm summary.h5 "+statepoint)

# Clean directory
os.system("rm *xml")
