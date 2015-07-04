from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = Options()

num_threads = options.getNumThreads()
azim_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()
num_polar = 2
polar_spacing = 1.0

log.set_log_level('NORMAL')


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
#                            Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(x=0.0, y=0.0, radius=1.0, name='pin')
xmin = XPlane(x=-2.0, name='xmin')
ymin = YPlane(y=-2.0, name='ymin')
zmin = ZPlane(z=-2.0, name='zmin')
xmax = XPlane(x=2.0, name='xmax')
ymax = YPlane(y=2.0, name='ymax')
zmax = ZPlane(z=2.0, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(REFLECTIVE)
ymax.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(REFLECTIVE)

###############################################################################
#                             Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = Cell(name='fuel')
fuel.setFill(materials['UO2'])
fuel.addSurface(halfspace=-1, surface=circle)
fuel.addSurface(halfspace=+1, surface=zmin)
fuel.addSurface(halfspace=-1, surface=zmax)

moderator = Cell(name='moderator')
moderator.setFill(materials['Water'])
moderator.addSurface(halfspace=+1, surface=circle)
moderator.addSurface(halfspace=+1, surface=xmin)
moderator.addSurface(halfspace=-1, surface=xmax)
moderator.addSurface(halfspace=+1, surface=ymin)
moderator.addSurface(halfspace=-1, surface=ymax)
moderator.addSurface(halfspace=+1, surface=zmin)
moderator.addSurface(halfspace=-1, surface=zmax)

###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

root_universe = Universe(name='root universe')
root_universe.addCell(fuel)
root_universe.addCell(moderator)


###############################################################################
#                         Creating the Geometry
###############################################################################


log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, num_polar, azim_spacing, \
                                 polar_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()


###############################################################################
#                             Generating Plots
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_tracks(track_generator)
#plotter.plot_segments(track_generator)
#plotter.plot_materials(geometry)
#plotter.plot_cells(geometry)
#plotter.plot_flat_source_regions(geometry)
#plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])
#plotter.plot_energy_fluxes(solver, fsrs=range(geometry.getNumFSRs()))

log.py_printf('TITLE', 'Finished')
