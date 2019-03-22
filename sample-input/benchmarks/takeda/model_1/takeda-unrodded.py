import openmoc
import openmoc.plotter as plotter
import openmoc.process as process
from openmoc.options import Options
from lattices import lattices, universes, cells
import numpy as np

refines = 1

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

opts = Options()

###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

c = universes['Core']
v = universes['Void']
a = universes['Control Rod']
r = universes['Reflector']

lattices['Root'].setWidth(width_x=5.0/refines, width_y=5.0/refines,
                          width_z=5.0/refines)
lattices['Root'].setUniverses(
    [[np.repeat([r, r, r, r, r], refines).tolist()] * 4 * refines +
     [np.repeat([r, r, r, v, r], refines).tolist()] * refines] * 2 * refines +
    [[np.repeat([r, r, r, r, r], refines).tolist()] * 2 * refines +
     [np.repeat([c, c, c, r, r], refines).tolist()] * 2 * refines +
     [np.repeat([c, c, c, v, r], refines).tolist()] * refines] * 3 * refines)

###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

cmfd = openmoc.Cmfd()
cmfd.setSORRelaxationFactor(1.5)
cmfd.setLatticeStructure(5*refines, 5*refines, 5*refines)
cmfd.setCentroidUpdateOn(False)

###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

geometry = openmoc.Geometry()
geometry.setRootUniverse(universes['Root'])
geometry.setCmfd(cmfd)
geometry.initializeFlatSourceRegions()

###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

quad = openmoc.EqualAnglePolarQuad()
quad.setNumPolarAngles(opts.num_polar)

track_generator = openmoc.TrackGenerator3D(geometry, opts.num_azim,
                                           opts.num_polar, opts.azim_spacing,
                                           opts.polar_spacing)
track_generator.setQuadrature(quad)
track_generator.setNumThreads(opts.num_omp_threads)
track_generator.setSegmentFormation(openmoc.OTF_STACKS)
track_generator.generateTracks()

###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = openmoc.CPUSolver(track_generator)
solver.setConvergenceThreshold(opts.tolerance)
solver.setNumThreads(opts.num_omp_threads)
solver.computeEigenvalue(opts.max_iters)
solver.printTimerReport()

###############################################################################
############################   Generating Plots   #############################
###############################################################################

process.compute_fission_rates(solver, use_hdf5=False)

plotter.plot_materials(geometry, gridsize=500, plane='xy')
plotter.plot_materials(geometry, gridsize=500, plane='xz', offset=-10.0)
plotter.plot_materials(geometry, gridsize=500, plane='yz')
plotter.plot_cells(geometry, gridsize=500)
plotter.plot_flat_source_regions(geometry, gridsize=500, plane='xy')
plotter.plot_flat_source_regions(geometry, gridsize=500, plane='xz')
plotter.plot_flat_source_regions(geometry, gridsize=500, plane='yz')
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2],
                            gridsize=500, plane='xy', offset=0.)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2],
                            gridsize=500, plane='xz', offset=0.)
plotter.plot_spatial_fluxes(solver, energy_groups=[1,2],
                            gridsize=500, plane='yz', offset=0.)

