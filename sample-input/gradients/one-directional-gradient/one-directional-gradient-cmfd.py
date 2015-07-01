import numpy
from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
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

log.py_printf('TITLE', 'Simulating a one group homogeneous, one directional'
    ' gradient...')

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

basic_material = Material(name='1-group infinite medium')
basic_material.setNumEnergyGroups(1)
basic_material.setSigmaA(numpy.array([0.069389522]))
basic_material.setSigmaF(numpy.array([0.0414198575]))
basic_material.setNuSigmaF(numpy.array([0.0994076580]))
basic_material.setSigmaS(numpy.array([0.383259177]))
basic_material.setChi(numpy.array([1.0]))
basic_material.setSigmaT(numpy.array([0.452648699]))


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

L = 100.0
left = XPlane(x=-L/2, name='left')
right = XPlane(x=L/2, name='right')
top = YPlane(y=L/2, name='top')
bottom = YPlane(y=-L/2, name='bottom')

left.setBoundaryType(VACUUM)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fill = Cell(name='fill')
fill.setFill(basic_material)

root_cell = Cell(name='root cell')
root_cell.addSurface(halfspace=+1, surface=left)
root_cell.addSurface(halfspace=-1, surface=right)
root_cell.addSurface(halfspace=+1, surface=bottom)
root_cell.addSurface(halfspace=-1, surface=top)


###############################################################################
###########################    Creating Universes   ###########################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

fill_universe = Universe(name='homogeneous fill cell')
fill_universe.addCell(fill)

root_universe = Universe(name='root universe')
root_universe.addCell(root_cell)

###############################################################################
############################    Creating Lattices   ###########################
###############################################################################

log.py_printf('NORMAL', 'Creating 100 x 1 lattice...')

num_cells_x = 100
num_cells_y = 1
lattice = Lattice(name='MxN lattice')
lattice.setWidth(width_x=L/num_cells_x, width_y=L/num_cells_y)
lattice.setUniverses([[fill_universe] * num_cells_x]*num_cells_y)
root_cell.setFill(lattice)

###############################################################################
###########################   Creating CMFD Mesh    ###########################
###############################################################################
log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = Cmfd()
cmfd.setMOCRelaxationFactor(0.6)
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(51,1)
cmfd.setGroupStructure([1])

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()

###############################################################################
############################    Generating Plots   ############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=100)
plotter.plot_cells(geometry, gridsize=100)
plotter.plot_flat_source_regions(geometry, gridsize=100)
plotter.plot_spatial_fluxes(solver, energy_groups=[1])

log.py_printf('TITLE', 'Finished')
