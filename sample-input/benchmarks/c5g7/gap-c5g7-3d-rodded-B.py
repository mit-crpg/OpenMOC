import openmoc
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options
from openmoc import process
from lattices import lattices, universes, cells, surfaces
from surfaces import gap

axial_refines = 6

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.num_omp_threads
azim_spacing = options.azim_spacing
num_azim = options.num_azim
polar_spacing = 1.2
num_polar = options.num_polar
tolerance = options.tolerance
max_iters = options.max_iters

###############################################################################
##########################   Create Core Lattice  #############################
###############################################################################

cells['Root'].addSurface(+1, surfaces['Root Small z-min'])
cells['Root'].addSurface(-1, surfaces['Root Small z-max'])
cells['Gap Root'].addSurface(+1, surfaces['Root Small z-min'])
cells['Gap Root'].addSurface(-1, surfaces['Root Small z-max'])

uu = universes['Gap UO2 Unrodded Assembly']
ur = universes['Gap UO2 Rodded Assembly']
mu = universes['Gap MOX Unrodded Assembly']
mr = universes['Gap MOX Rodded Assembly']
rr = universes['Gap Reflector Rodded Assembly']
ri = universes['Gap Reflector Right Assembly']
rb = universes['Gap Reflector Bottom Assembly']
rc = universes['Gap Reflector Corner Assembly']

# 3 x 3 x 9 core to represent 3D core
lattices['Root'].setWidth(width_x=21.42+2*gap, width_y=21.42+2*gap, width_z=7.14/axial_refines)
lattices['Root'].setUniverses(  [[[rr, rr, ri],
                                  [rr, rr, ri],
                                  [rb, rb, rc]]] * 3 * axial_refines +
                                [[[ur, mr, ri],
                                  [mr, uu, ri],
                                  [rb, rb, rc]]] * 2 * axial_refines +
                                [[[ur, mu, ri],
                                  [mu, uu, ri],
                                  [rb, rb, rc]]] * 2 * axial_refines +
                                [[[uu, mu, ri],
                                  [mu, uu, ri],
                                  [rb, rb, rc]]] * 2 * axial_refines)

###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(51,51,9*axial_refines)
cmfd.setCentroidUpdateOn(False)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = openmoc.Geometry()
geometry.setRootUniverse(universes['Gap Root'])
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

quad = openmoc.EqualWeightPolarQuad()
quad.setNumPolarAngles(num_polar)

track_generator = openmoc.TrackGenerator3D(geometry, num_azim, num_polar,
                                           azim_spacing, polar_spacing)
track_generator.setQuadrature(quad)
track_generator.setNumThreads(num_threads)
track_generator.setSegmentFormation(openmoc.OTF_STACKS)
track_generator.setSegmentationZones([-32.13, -10.71, 10.71, 32.13])
track_generator.generateTracks()

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()

mesh = process.Mesh()
mesh.dimension = [34,34,3]
mesh.lower_left = [-32.13, -10.71, -32.13]
mesh.upper_right = [10.71, 32.13, 10.71]
mesh.width = [1.26, 1.26, 14.28]
fission_rates = mesh.tally_fission_rates(solver, volume='integrated')
for k in range(3):
  print('Z = ' + str(k))
  for i in range(34):
    msg = ''
    for j in range(34):
      msg += str(fission_rates[i][j][k])
      msg += ' '
    print(msg)
  print('...')

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
