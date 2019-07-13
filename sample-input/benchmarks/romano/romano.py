import numpy
import openmoc

###############################################################################
#                          Main Simulation Parameters
###############################################################################

opts = openmoc.options.Options()

openmoc.log.set_log_level('NORMAL')

openmoc.log.py_printf('TITLE', 'Simulating HW3 from Fall 2010 22.212...')


###############################################################################
#                            Creating Materials
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating materials...')

fuel = openmoc.Material(name='fuel')
moderator = openmoc.Material(name='moderator')

fuel.setNumEnergyGroups(1)
moderator.setNumEnergyGroups(1)

fuel.setSigmaT(numpy.array([0.452648699]))
fuel.setSigmaF(numpy.array([0.0414198575]))
fuel.setNuSigmaF(numpy.array([0.0994076580]))
fuel.setSigmaS(numpy.array([0.38259177]))
fuel.setChi(numpy.array([1.0]))

moderator.setSigmaT(numpy.array([0.841545641]))
moderator.setSigmaF(numpy.array([0.0]))
moderator.setNuSigmaF(numpy.array([0.0]))
moderator.setSigmaS(numpy.array([0.837794542]))
moderator.setChi(numpy.array([1.0]))


###############################################################################
#                            Creating Surfaces
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating surfaces...')

zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.4)
left = openmoc.XPlane(x=-0.635)
right = openmoc.XPlane(x=0.635)
top = openmoc.YPlane(y=0.635)
bottom = openmoc.YPlane(y=-0.635)

left.setBoundaryType(openmoc.REFLECTIVE)
right.setBoundaryType(openmoc.REFLECTIVE)
top.setBoundaryType(openmoc.REFLECTIVE)
bottom.setBoundaryType(openmoc.REFLECTIVE)


###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating cells...')

fuel_cell = openmoc.Cell(name='fuel')
fuel_cell.setFill(fuel)
fuel_cell.addSurface(halfspace=-1, surface=zcylinder)

moderator_cell = openmoc.Cell(name='moderator')
moderator_cell.setFill(moderator)
moderator_cell.addSurface(halfspace=+1, surface=zcylinder)
moderator_cell.addSurface(halfspace=+1, surface=left)
moderator_cell.addSurface(halfspace=-1, surface=right)
moderator_cell.addSurface(halfspace=+1, surface=bottom)
moderator_cell.addSurface(halfspace=-1, surface=top)


###############################################################################
#                             Creating Cells
###############################################################################

openmoc.log.py_printf('NORMAL', 'Creating universes...')

root_universe = openmoc.Universe(name='root universe')
root_universe.addCell(fuel_cell)
root_universe.addCell(moderator_cell)


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

openmoc.log.py_printf('TITLE', 'Finished')
