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

left = openmoc.XPlane(x=-2.0, name='left')
right = openmoc.XPlane(x=2.0, name='right')
top = openmoc.YPlane(y=2.0, name='top')
bottom = openmoc.YPlane(y=-2.0, name='bottom')

left.setBoundaryType(openmoc.REFLECTIVE)
right.setBoundaryType(openmoc.REFLECTIVE)
top.setBoundaryType(openmoc.REFLECTIVE)
bottom.setBoundaryType(openmoc.REFLECTIVE)


###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

fuel = openmoc.Cell(name='fuel')
fuel.setFill(materials['UO2'])
fuel.addSurface(halfspace=+1, surface=left)
fuel.addSurface(halfspace=-1, surface=right)
fuel.addSurface(halfspace=+1, surface=bottom)
fuel.addSurface(halfspace=-1, surface=top)


###############################################################################
#                            Creating Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

root_universe = openmoc.Universe(name='root universe')
root_universe.addCell(fuel)


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

track_generator = openmoc.TrackGenerator(geometry, 128, .05)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setNumThreads(opts.num_omp_threads)
solver.setConvergenceThreshold(1e-7)
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
