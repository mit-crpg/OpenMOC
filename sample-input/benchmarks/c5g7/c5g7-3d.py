from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options
from assemblies import *

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.getNumThreads()
azim_spacing = options.getAzimSpacing()
num_azim = options.getNumAzimAngles()
polar_spacing = options.getPolarSpacing()
num_polar = options.getNumPolarAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()


# 3 x 3 x 9 core to represent 3D core
lattices.append(Lattice(name='Full Geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42, width_z=7.14)
lattices[-1].setUniverses3D([[[assembly_rfl_unrod    , assembly_rfl_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod    , assembly_rfl_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_rfl_unrod    , assembly_rfl_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod    , assembly_rfl_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_rfl_unrod    , assembly_rfl_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod    , assembly_rfl_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_uo2_unrod    , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_uo2_unrod    , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_uo2_unrod    , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_uo2_unrod    , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_uo2_unrod    , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]],
                             [[assembly_uo2_unrod    , assembly_mox_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_mox_unrod    , assembly_uo2_unrod    , assembly_rfl_unrod_rgt],
                              [assembly_rfl_unrod_btm, assembly_rfl_unrod_btm, assembly_rfl_unrod_cnr]]])

root_cell.setFill(lattices[-1])


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = Cmfd()
cmfd.setMOCRelaxationFactor(1.0)
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(51,51,9)
cmfd.setGroupStructure([1,4,8])
cmfd.setKNearest(4)


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

quad = EqualAnglePolarQuad()
quad.setNumPolarAngles(num_polar)

track_generator = TrackGenerator(geometry, num_azim, num_polar, azim_spacing,
                                 polar_spacing)
track_generator.setQuadrature(quad)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(track_generator)
solver.setConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)
solver.computeEigenvalue(max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
#plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7],
#                            gridsize=500, plane='xy', offset=0.)
#plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7],
#                            gridsize=500, plane='xz', offset=0.)
#plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7],
#                            gridsize=500, plane='yz', offset=0.)


log.py_printf('TITLE', 'Finished')
