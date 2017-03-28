import openmoc

###############################################################################
#                           Main Simulation Parameters
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')


###############################################################################
#                              Creating Materials
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = openmoc.materialize.load_from_hdf5('c5g7-mgxs.h5', '../')


###############################################################################
#                              Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.8, name='pin')

boundary = openmoc.RectangularPrism(4., 4.)
boundary.setBoundaryType(openmoc.REFLECTIVE)


###############################################################################
#                                Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

fuel = openmoc.Cell(name='fuel')
fuel.setFill(materials['UO2'])
fuel.addSurface(halfspace=-1, surface=zcylinder)

moderator = openmoc.Cell(name='moderator')
moderator.setFill(materials['Water'])
moderator.addSurface(halfspace=+1, surface=zcylinder)

root_cell = openmoc.Cell(name='root cell')
root_cell.setRegion(boundary)


###############################################################################
#                            Creating Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

pincell = openmoc.Universe(name='pin cell')
root_universe = openmoc.Universe(name='root universe')

pincell.addCell(fuel)
pincell.addCell(moderator)
root_universe.addCell(root_cell)


###############################################################################
#                             Creating Lattices
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating simple 2 x 2 lattice...')

lattice = openmoc.Lattice(name='2x2 lattice')
lattice.setWidth(width_x=2.0, width_y=2.0)
lattice.setUniverses([[[pincell, pincell], [pincell, pincell]]])

root_cell.setFill(lattice)


###############################################################################
#                            Creating the Geometry
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
#                              Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry, gridsize=50)
openmoc.plotter.plot_cells(geometry, gridsize=50)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=50)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

openmoc.log.py_printf('TITLE', 'Finished')
