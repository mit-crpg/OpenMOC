from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#                           Main Simulation Parameters 
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

log.set_log_level('NORMAL')


###############################################################################
#                              Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
#                              Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circle = Circle(x=0.0, y=0.0, radius=0.8, name='pin')
left = XPlane(x=-2.0, name='left')
right = XPlane(x=2.0, name='right')
top = YPlane(y=2.0, name='top')
bottom = YPlane(y=-2.0, name='bottom')

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#                                Creating Cells
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = list()
cells.append(CellBasic(name='fuel'))
cells.append(CellBasic(name='moderator'))
cells.append(CellFill(name='root cell'))

cells[0].setMaterial(materials['UO2'])
cells[1].setMaterial(materials['Water'])

cells[0].addSurface(halfspace=-1, surface=circle)
cells[1].addSurface(halfspace=+1, surface=circle)
cells[2].addSurface(halfspace=+1, surface=left)
cells[2].addSurface(halfspace=-1, surface=right)
cells[2].addSurface(halfspace=+1, surface=bottom)
cells[2].addSurface(halfspace=-1, surface=top)


###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

pincell = Universe(name='pin cell')
root = Universe(name='root universe')

pincell.addCell(cells[0])
pincell.addCell(cells[1])
root.addCell(cells[2])


###############################################################################
#                             Creating Lattices
###############################################################################

log.py_printf('NORMAL', 'Creating simple 2 x 2 lattice...')

lattice = Lattice(name='2x2 lattice')
lattice.setWidth(width_x=2.0, width_y=2.0)
lattice.setUniverses([[pincell, pincell], [pincell, pincell]])

cells[2].setFill(lattice)


###############################################################################
#                            Creating the Geometry
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root)
geometry.initializeFlatSourceRegions()


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

solver = CPUSolver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()


###############################################################################
#                              Generating Plots
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_tracks(track_generator)
plotter.plot_segments(track_generator)
plotter.plot_materials(geometry, gridsize=50)
plotter.plot_cells(geometry, gridsize=50)
plotter.plot_flat_source_regions(geometry, gridsize=50)
#plotter.plot_fluxes(geometry, solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
