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

boundary = openmoc.RectangularPrism(4., 4.)
boundary.setBoundaryType(openmoc.PERIODIC)

zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.3, name='fuel pin')


###############################################################################
#                            Creating Grid Spacer
###############################################################################

outer_box = openmoc.RectangularPrism(1.00, 1.00)
inner_box = openmoc.RectangularPrism(0.96, 0.96)
inner_complement = openmoc.Complement()
inner_complement.addNode(inner_box)

grid_spacer_region = openmoc.Intersection()
grid_spacer_region.addNode(outer_box)
grid_spacer_region.addNode(inner_complement)


###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

fuel = openmoc.Cell(name='fuel')
fuel.setNumRings(3)
fuel.setNumSectors(8)
fuel.setFill(materials['UO2'])
fuel.addSurface(halfspace=-1, surface=zcylinder)

moderator = openmoc.Cell(name='amoderator')
fuel.setNumSectors(8)
moderator.setFill(materials['Water'])
moderator.setRegion(inner_box)
moderator.addSurface(halfspace=+1, surface=zcylinder)

grid_spacer = openmoc.Cell(name='mox grid spacer')
grid_spacer.setFill(materials['MOX-4.3%'])
grid_spacer.setRegion(grid_spacer_region)

root_cell = openmoc.Cell(name='root cell')
root_cell.setRegion(boundary)


###############################################################################
#                            Creating Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

pin = openmoc.Universe(name='pin cell')
root_universe = openmoc.Universe(name='root universe')

pin.addCell(fuel)
pin.addCell(moderator)
pin.addCell(grid_spacer)
root_universe.addCell(root_cell)


###############################################################################
#                            Creating Lattices
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating 4 x 4 lattice w/ grid spacers...')

lattice = openmoc.Lattice(name='4x4 lattice w/ grid spacers')
lattice.setWidth(width_x=1.0, width_y=1.0)
lattice.setUniverses([[[pin, pin, pin, pin],
                       [pin, pin, pin, pin],
                       [pin, pin, pin, pin],
                       [pin, pin, pin, pin]]])
root_cell.setFill(lattice)


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

openmoc.plotter.plot_segments(track_generator)
openmoc.plotter.plot_materials(geometry, gridsize=500)
openmoc.plotter.plot_cells(geometry, gridsize=500)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=500, centroids=True)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

openmoc.log.py_printf('TITLE', 'Finished')
