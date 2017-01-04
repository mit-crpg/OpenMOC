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

materials = openmoc.materialize.load_from_hdf5('single-group-materials', './')


###############################################################################
#                            Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

#zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=1.0, name='pin')
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
fuel.setFill(materials['single'])
#fuel.addSurface(halfspace=-1, surface=zcylinder)
#fuel.addSurface(halfspace=+1, surface=zcylinder)
fuel.addSurface(halfspace=+1, surface=left)
fuel.addSurface(halfspace=-1, surface=right)
fuel.addSurface(halfspace=+1, surface=bottom)
fuel.addSurface(halfspace=-1, surface=top)

'''
moderator = openmoc.Cell(name='moderator')
moderator.setFill(materials['Water'])
moderator.addSurface(halfspace=+1, surface=zcylinder)
moderator.addSurface(halfspace=+1, surface=left)
moderator.addSurface(halfspace=-1, surface=right)
moderator.addSurface(halfspace=+1, surface=bottom)
moderator.addSurface(halfspace=-1, surface=top)
'''

###############################################################################
#                            Creating Universes
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

root_universe = openmoc.Universe(name='root universe')
root_universe.addCell(fuel)
#root_universe.addCell(moderator)


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
solver.computeEigenvalue(100000,10000,7)
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
