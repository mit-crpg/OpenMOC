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
azim_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()
num_polar = 6
polar_spacing = 0.05
log.set_log_level('NORMAL')
set_line_length(120)

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

xmin = XPlane(x=-0.5, name='xmin')
xmax = XPlane(x= 0.5, name='xmax')
ymin = YPlane(y=-0.5, name='ymin')
ymax = YPlane(y= 0.5, name='ymax')
zmin = ZPlane(z=-0.5, name='zmin')
zmax = ZPlane(z= 0.5, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
ymax.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(REFLECTIVE)

###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = CellBasic(name='fuel')
fuel.setMaterial(materials['UO2'])

moderator = CellBasic(name='moderator')
moderator.setMaterial(materials['Water'])

root_cell = CellFill(name='root cell')
root_cell.addSurface(halfspace=+1, surface=xmin)
root_cell.addSurface(halfspace=-1, surface=xmax)
root_cell.addSurface(halfspace=+1, surface=ymin)
root_cell.addSurface(halfspace=-1, surface=ymax)
root_cell.addSurface(halfspace=+1, surface=zmin)
root_cell.addSurface(halfspace=-1, surface=zmax)


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

fue_univ = Universe(name='homogeneous fue cell')
fue_univ.addCell(fuel)

mod_univ = Universe(name='homogeneous mod cell')
mod_univ.addCell(moderator)

root_universe = Universe(name='root universe')
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple 2 x 2 lattice...')

lattice = Lattice(name='2x2 lattice')
lattice.setWidth(width_x=0.5, width_y=0.5, width_z=0.5)
lattice.setUniverses3D([[[fue_univ, mod_univ],
                         [mod_univ, mod_univ]],
                        [[mod_univ, mod_univ],
                         [mod_univ, mod_univ]]])

root_cell.setFill(lattice)

###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = Cmfd()
cmfd.setMOCRelaxationFactor(0.6)
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(2,2,2)
#cmfd.setGroupStructure([1,4,8])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()

quad = EqualAnglePolarQuad()
quad.setNumPolarAngles(num_polar)


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, num_polar, azim_spacing, polar_spacing)
track_generator.setQuadrature(quad)
track_generator.setNumThreads(num_threads)
#track_generator.setSolve2D()
track_generator.setZLevel(0.1)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.useExponentialIntrinsic()
solver.convergeSource(max_iters)
solver.printTimerReport()


#plotter.plot_tracks(track_generator)
#plotter.plot_tracks_3d(track_generator)
#plotter.plot_segments_3d(track_generator)
#plotter.plot_segments(track_generator)

log.py_printf('TITLE', 'Finished')
