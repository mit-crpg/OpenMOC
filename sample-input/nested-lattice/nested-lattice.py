from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

log.set_log_level('NORMAL')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circles = list()
planes = list()
planes.append(XPlane(x=-2.0, name='left'))
planes.append(XPlane(x=2.0, name='right'))
planes.append(YPlane(y=-2.0, name='top'))
planes.append(YPlane(y=2.0, name='bottom'))
circles.append(Circle(x=0.0, y=0.0, radius=0.4, name='large pin'))
circles.append(Circle(x=0.0, y=0.0, radius=0.3, name='medium pin'))
circles.append(Circle(x=0.0, y=0.0, radius=0.2, name='small pin'))
for plane in planes: plane.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = list()
cells = list()
cells.append(CellBasic(name='large pin fuel', rings=3, sectors=8))
cells.append(CellBasic(name='large pin moderator', sectors=8))
cells.append(CellBasic(name='medium pin fuel', rings=3, sectors=8))
cells.append(CellBasic(name='medium pin moderator', sectors=8))
cells.append(CellBasic(name='small pin fuel', rings=3, sectors=8))
cells.append(CellBasic(name='small pin moderator', sectors=8))
cells.append(CellFill(name='lattice cell'))
cells.append(CellFill(name='root cell'))

cells[0].setMaterial(materials['UO2'])
cells[1].setMaterial(materials['Water'])
cells[2].setMaterial(materials['UO2'])
cells[3].setMaterial(materials['Water'])
cells[4].setMaterial(materials['UO2'])
cells[5].setMaterial(materials['Water'])

cells[0].addSurface(halfspace=-1, surface=circles[0])
cells[1].addSurface(halfspace=+1, surface=circles[0])
cells[2].addSurface(halfspace=-1, surface=circles[1])
cells[3].addSurface(halfspace=+1, surface=circles[1])
cells[4].addSurface(halfspace=-1, surface=circles[2])
cells[5].addSurface(halfspace=+1, surface=circles[2])

cells[7].addSurface(halfspace=+1, surface=planes[0])
cells[7].addSurface(halfspace=-1, surface=planes[1])
cells[7].addSurface(halfspace=+1, surface=planes[2])
cells[7].addSurface(halfspace=-1, surface=planes[3])


###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

pin1 = Universe(name='large pin cell')
pin2 = Universe(name='medium pin cell')
pin3 = Universe(name='small pin cell')
assembly = Universe(name='2x2 lattice')
root = Universe(name='root universe')

pin1.addCell(cells[0])
pin1.addCell(cells[1])
pin2.addCell(cells[2])
pin2.addCell(cells[3])
pin3.addCell(cells[4])
pin3.addCell(cells[5])
assembly.addCell(cells[6])
root.addCell(cells[7])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating nested 2 x 2 lattices...')

# 2x2 assembly
lattice = Lattice(name='2x2 lattice')
lattice.setWidth(width_x=1.0, width_y=1.0)
lattice.setUniverses([[pin1, pin2], [pin1, pin3]])
cells[6].setFill(lattice)

# 2x2 core
core = Lattice(name='2x2 core')
core.setWidth(width_x=2.0, width_y=2.0)
core.setUniverses([[assembly, assembly], [assembly, assembly]])
cells[7].setFill(core)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root)
geometry.initializeFlatSourceRegions()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_tracks(track_generator)
#plotter.plot_segments(track_generator)
#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
#plotter.plot_fluxes(geometry, solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
