import sys
sys.path.append('..')
from cube import geometry, root_cell, left, right, top, bottom
import openmoc

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = openmoc.options.Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

openmoc.log.set_log_level('NORMAL')
openmoc.log.py_printf('TITLE', \
  'Simulating a one group homogeneous, one directional gradient...')

###############################################################################
###########################   Creating CMFD Mesh    ###########################
###############################################################################
openmoc.log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setMOCRelaxationFactor(0.6)
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(25,1)

###############################################################################
#########################   Load the Cubic Geometry   #########################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Importing cubic geometry...')

left.setBoundaryType(openmoc.VACUUM)
right.setBoundaryType(openmoc.REFLECTIVE)
top.setBoundaryType(openmoc.REFLECTIVE)
bottom.setBoundaryType(openmoc.REFLECTIVE)

geometry.setCmfd(cmfd)

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = openmoc.TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.setConvergenceThreshold(tolerance)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()

###############################################################################
############################    Generating Plots   ############################
###############################################################################

openmoc.log.py_printf('NORMAL', 'Plotting data...')

openmoc.plotter.plot_materials(geometry, gridsize=100)
openmoc.plotter.plot_cells(geometry, gridsize=100)
openmoc.plotter.plot_flat_source_regions(geometry, gridsize=100)
openmoc.plotter.plot_spatial_fluxes(solver, energy_groups=[1])

openmoc.log.py_printf('TITLE', 'Finished')
