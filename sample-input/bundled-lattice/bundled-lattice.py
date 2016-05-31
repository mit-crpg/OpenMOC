import openmoc

###############################################################################
#                          Main Simulation Parameters
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
#                            Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

left = openmoc.XPlane(x=-34.0, name='left')
right = openmoc.XPlane(x=34.0, name='right')
top = openmoc.YPlane(y=-34.0, name='top')
bottom = openmoc.YPlane(y=34.0, name='bottom')
boundaries = [left, right, top, bottom]
for boundary in boundaries: boundary.setBoundaryType(openmoc.REFLECTIVE)

zcylinders = list()
radii = [0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
for r in radii: zcylinders.append(openmoc.ZCylinder(x=0.0, y=0.0, radius=r))


###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

# Create a list of 12 Cell instances
cells = [openmoc.Cell() for i in range(12)]

# Append 3 CellFills for the assemblies and full core
assembly1 = openmoc.Cell(name='assembly 1')
assembly2 = openmoc.Cell(name='assembly 2')
root_cell = openmoc.Cell(name='full core')

# Create fuel/moderator by adding the appropriate Surfaces and Materials
cells[0].addSurface(halfspace=-1, surface=zcylinders[0])
cells[1].addSurface(halfspace=+1, surface=zcylinders[0])
cells[2].addSurface(halfspace=-1, surface=zcylinders[1])
cells[3].addSurface(halfspace=+1, surface=zcylinders[1])
cells[4].addSurface(halfspace=-1, surface=zcylinders[2])
cells[5].addSurface(halfspace=+1, surface=zcylinders[2])
cells[6].addSurface(halfspace=-1, surface=zcylinders[3])
cells[7].addSurface(halfspace=+1, surface=zcylinders[3])
cells[8].addSurface(halfspace=-1, surface=zcylinders[4])
cells[9].addSurface(halfspace=+1, surface=zcylinders[4])
cells[10].addSurface(halfspace=-1, surface=zcylinders[5])
cells[11].addSurface(halfspace=+1, surface=zcylinders[5])

cells[0].setFill(materials['UO2'])
cells[1].setFill(materials['Water'])
cells[2].setFill(materials['UO2'])
cells[3].setFill(materials['Water'])
cells[4].setFill(materials['UO2'])
cells[5].setFill(materials['Water'])
cells[6].setFill(materials['UO2'])
cells[7].setFill(materials['Water'])
cells[8].setFill(materials['UO2'])
cells[9].setFill(materials['Water'])
cells[10].setFill(materials['UO2'])
cells[11].setFill(materials['Water'])

# Add the boundary Planes to the "root" Cell
root_cell.addSurface(halfspace=+1, surface=boundaries[0])
root_cell.addSurface(halfspace=-1, surface=boundaries[1])
root_cell.addSurface(halfspace=+1, surface=boundaries[2])
root_cell.addSurface(halfspace=-1, surface=boundaries[3])


###############################################################################
#                            Creating Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

u1 = openmoc.Universe(name='pin 1')
u2 = openmoc.Universe(name='pin 2')
u3 = openmoc.Universe(name='pin 3')
u4 = openmoc.Universe(name='pin 4')
u5 = openmoc.Universe(name='pin 5')
u6 = openmoc.Universe(name='pin 6')
u7 = openmoc.Universe(name='2x2 lattice')
u8 = openmoc.Universe(name='2x2 lattice')
root_universe = openmoc.Universe(name='root universe')

# Add the appropriate Cells to each Universe
u1.addCell(cells[0])
u1.addCell(cells[1])
u2.addCell(cells[2])
u2.addCell(cells[3])
u3.addCell(cells[4])
u3.addCell(cells[5])
u4.addCell(cells[6])
u4.addCell(cells[7])
u5.addCell(cells[8])
u5.addCell(cells[9])
u6.addCell(cells[10])
u6.addCell(cells[11])
u7.addCell(assembly1)
u8.addCell(assembly2)
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating 4 x 4 core of 17 x 17 assemblies...')

# 1st 17x17 assembly
a1 = openmoc.Lattice(name='assembly 1')
a1.setWidth(width_x=1.0, width_y=1.0)
a1.setUniverses([[
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1]]])

# 2nd 17x17 assembly
a2 = openmoc.Lattice(name='assembly 2')
a2.setWidth(width_x=1.0, width_y=1.0)
a2.setUniverses([[
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4]]])

# 4x4 core
core = openmoc.Lattice(name='full core')
core.setWidth(width_x=17.0, width_y=17.0)
core.setUniverses([[[u7, u8, u7, u8],
                    [u8, u7, u8, u7],
                    [u7, u8, u7, u8],
                    [u8, u7, u8, u7]]])

assembly1.setFill(a1)
assembly2.setFill(a2)
root_cell.setFill(core)


###############################################################################
#                         Creating the Geometry
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(root_universe)


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, opts.num_azim,
                                         opts.azim_spacing)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setNumThreads(opts.num_omp_threads)
solver.setConvergenceThreshold(opts.tolerance)
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry, gridsize=500)
openmoc.plotter.plot_cells(geometry, gridsize=500)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=500)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

openmoc.log.py_printf('TITLE', 'Finished')
