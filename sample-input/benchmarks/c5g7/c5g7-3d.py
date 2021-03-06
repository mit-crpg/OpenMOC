import openmoc
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options
from lattices import lattices, universes, cells, surfaces

axial_refines = 2

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.num_omp_threads
azim_spacing = options.azim_spacing
num_azim = options.num_azim
polar_spacing = options.polar_spacing
num_polar = options.num_polar
tolerance = options.tolerance
max_iters = options.max_iters

###############################################################################
##########################   Create Core Lattice  #############################
###############################################################################

cells['Root'].addSurface(+1, surfaces['Root Small z-min'])
cells['Root'].addSurface(-1, surfaces['Root Small z-max'])

uu = universes['UO2 Unrodded Assembly']
ur = universes['UO2 Rodded Assembly']
mu = universes['MOX Unrodded Assembly']
mr = universes['MOX Rodded Assembly']
ru = universes['Reflector Unrodded Assembly']
rr = universes['Reflector Rodded Assembly']
ri = universes['Reflector Right Assembly']
rb = universes['Reflector Bottom Assembly']
rc = universes['Reflector Corner Assembly']

# 3 x 3 x 10 core to represent 3D core
lattices['Root'].setWidth(width_x=21.42, width_y=21.42, width_z=21.42/axial_refines)
lattices['Root'].setUniverses([[[ru, ru, ri],
                                  [ru, ru, ri],
                                  [rb, rb, rc]]] * axial_refines +
                                [[[uu, mu, ri],
                                  [mu, uu, ri],
                                  [rb, rb, rc]]] * 9 * axial_refines)

###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(51,51,10*axial_refines)
cmfd.setGroupStructure([[1,2,3],[4,5,6,7]])
cmfd.setCentroidUpdateOn(False)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(universes['Root'])
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

quad = openmoc.EqualAnglePolarQuad()
quad.setNumPolarAngles(num_polar)

track_generator = openmoc.TrackGenerator3D(geometry, num_azim, num_polar,
                                           azim_spacing, polar_spacing)
track_generator.setQuadrature(quad)
track_generator.setNumThreads(num_threads)
track_generator.setSegmentFormation(openmoc.OTF_STACKS)
track_generator.setSegmentationZones([-32.13, 32.13])
track_generator.generateTracks()

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()

###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=500, plane='xy')
plotter.plot_materials(geometry, gridsize=500, plane='xz', offset=-10.0)
plotter.plot_materials(geometry, gridsize=500, plane='yz')
plotter.plot_cells(geometry, gridsize=500)
plotter.plot_flat_source_regions(geometry, gridsize=500, plane='xy')
plotter.plot_flat_source_regions(geometry, gridsize=500, plane='xz')
plotter.plot_flat_source_regions(geometry, gridsize=500, plane='yz')
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7],
                            gridsize=500, plane='xy', offset=0.)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7],
                            gridsize=500, plane='xz', offset=0.)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7],
                            gridsize=500, plane='yz', offset=0.)

log.py_printf('TITLE', 'Finished')
