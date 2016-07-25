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


##############################################################################
#                            Running a Simulation
###############################################################################


solver = openmoc.MCSolver()
#solver.setConvergenceThreshold(opts.tolerance)
solver.setGeometry(geometry)
solver.initialize()
solver.setNumThreads(opts.num_omp_threads)
solver.computeEigenvalue(10000000,1,7)
#solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

'''
openmoc.plotter.plot_quadrature(solver)
openmoc.plotter.plot_tracks(track_generator)
openmoc.plotter.plot_segments(track_generator)
'''
openmoc.plotter.plot_flat_source_regions(geometry)
openmoc.plotter.plot_cells(geometry)
openmoc.plotter.plot_materials(geometry)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])
openmoc.plotter.plot_energy_fluxes(solver, fsrs=range(geometry.getNumFSRs()))

openmoc.log.py_printf('TITLE', 'Finished')
