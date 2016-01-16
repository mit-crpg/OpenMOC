import openmoc
import openmoc.plotter as plotter
import openmoc.process as process
from openmoc.options import Options
from universes import universes, cells

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
refines = 1

###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

c = universes['Core']
v = universes['Void']
a = universes['Control Rod']
r = universes['Reflector']

lattice = openmoc.Lattice(name='5x5 lattice')
lattice.setWidth(width_x=5.0/refines, width_y=5.0/refines, width_z=5.0/refines)
template = [[[r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, a, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, r, r],
             [r, r, r, a, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [c, c, c, r, r],
             [c, c, c, r, r],
             [c, c, c, a, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [c, c, c, r, r],
             [c, c, c, r, r],
             [c, c, c, a, r]],
            [[r, r, r, r, r],
             [r, r, r, r, r],
             [c, c, c, r, r],
             [c, c, c, r, r],
             [c, c, c, a, r]]]


# refine lattice
for k in range(5):
  for j in range(5):
    a = template[k][j]
    template[k][j] = sum([[a[i]]*refines for i in range(len(a))], [])
  a = template[k]
  template[k] = sum([[a[i]]*refines for i in range(len(a))], [])
a = template
template = sum([[a[i]]*refines for i in range(len(a))], [])

lattice.setUniverses3D(template)
cells['Root'].setFill(lattice)


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

cmfd = openmoc.Cmfd()
cmfd.setMOCRelaxationFactor(1.0)
cmfd.setSORRelaxationFactor(1.5)
cmfd.setOpticallyThick(True)
cmfd.setLatticeStructure(5*refines, 5*refines, 5*refines)
cmfd.setCentroidUpdateOn(False)
cmfd.setKNearest(1)


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
quad.setNumPolarAngles(num_polar)

track_generator = openmoc.TrackGenerator(geometry, num_azim, num_polar, azim_spacing,
                                         polar_spacing)
track_generator.setQuadrature(quad)
track_generator.setNumThreads(num_threads)
track_generator.setOTF()
track_generator.setSegmentationHeights([0.1])
#track_generator.setGlobalZMesh()
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

process.compute_fission_rates(solver, use_hdf5=False)
process.compute_material_fluxes(solver, use_hdf5=False)

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
