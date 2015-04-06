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

left = XPlane(x=-34.0, name='left')
right = XPlane(x=34.0, name='right')
top = YPlane(y=-34.0, name='top')
bottom = YPlane(y=34.0, name='bottom')
boundaries = [left, right, top, bottom]
for boundary in boundaries: boundary.setBoundaryType(REFLECTIVE)

circles = list()
radii = [0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
for r in radii: circles.append(Circle(x=0.0, y=0.0, radius=r))


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

# Create a list of 10 CellBasic instances
cells = [CellBasic() for i in range(12)]

# Append 3 CellFills for the assemblies and full core
assembly1 = CellFill(name='assembly 1')
assembly2 = CellFill(name='assembly 2')
root_cell = CellFill(name='full core')

# Create fuel/moderator by adding the appropriate Surfaces and Materials
cells[0].addSurface(halfspace=-1, surface=circles[0])
cells[1].addSurface(halfspace=+1, surface=circles[0])
cells[2].addSurface(halfspace=-1, surface=circles[1])
cells[3].addSurface(halfspace=+1, surface=circles[1])
cells[4].addSurface(halfspace=-1, surface=circles[2])
cells[5].addSurface(halfspace=+1, surface=circles[2])
cells[6].addSurface(halfspace=-1, surface=circles[3])
cells[7].addSurface(halfspace=+1, surface=circles[3])
cells[8].addSurface(halfspace=-1, surface=circles[4])
cells[9].addSurface(halfspace=+1, surface=circles[4])
cells[10].addSurface(halfspace=-1, surface=circles[5])
cells[11].addSurface(halfspace=+1, surface=circles[5])

cells[0].setMaterial(materials['UO2'])
cells[1].setMaterial(materials['Water'])
cells[2].setMaterial(materials['UO2'])
cells[3].setMaterial(materials['Water'])
cells[4].setMaterial(materials['UO2'])
cells[5].setMaterial(materials['Water'])
cells[6].setMaterial(materials['UO2'])
cells[7].setMaterial(materials['Water'])
cells[8].setMaterial(materials['UO2'])
cells[9].setMaterial(materials['Water'])
cells[10].setMaterial(materials['UO2'])
cells[11].setMaterial(materials['Water'])

# Add the boundary Planes to the "root" Cell
root_cell.addSurface(halfspace=+1, surface=boundaries[0])
root_cell.addSurface(halfspace=-1, surface=boundaries[1])
root_cell.addSurface(halfspace=+1, surface=boundaries[2])
root_cell.addSurface(halfspace=-1, surface=boundaries[3])


###############################################################################
#                            Creating Universes
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

u1 = Universe(name='pin 1')
u2 = Universe(name='pin 2')
u3 = Universe(name='pin 3')
u4 = Universe(name='pin 4')
u5 = Universe(name='pin 5')
u6 = Universe(name='pin 6')
u7 = Universe(name='2x2 lattice')
u8 = Universe(name='2x2 lattice')
root_universe = Universe(name='root universe')

# Add the appropriate Cells to each Universe
u1.addCell(cells[0])
u1.addCell(cells[1])
u2.addCell(cells[2])
u2.addCell(cells[3])
u3.addCell(cells[4])
u3.addCell(cells[5])
u4.addCell(cells[6])
u4.addCell(cells[7])
u5.addCell(cells[8])
u5.addCell(cells[9])
u6.addCell(cells[10])
u6.addCell(cells[11])
u7.addCell(assembly1)
u8.addCell(assembly2)
root_universe.addCell(root_cell)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating 4 x 4 core of 17 x 17 assemblies...')

# 1st 17x17 assembly
a1 = Lattice(name='assembly 1')
a1.setWidth(width_x=1.0, width_y=1.0)
a1.setUniverses([
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1],
    [u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2, u3, u2],
    [u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1, u2, u1]])

# 2nd 17x17 assembly
a2 = Lattice(name='assembly 2')
a2.setWidth(width_x=1.0, width_y=1.0)
a2.setUniverses([
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4],
    [u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5, u6, u5],
    [u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4, u5, u4]])

# 4x4 core
core = Lattice(name='full core')
core.setWidth(width_x=17.0, width_y=17.0)
core.setUniverses([[u7, u8, u7, u8],
                   [u8, u7, u8, u7],
                   [u7, u8, u7, u8],
                   [u8, u7, u8, u7]])

assembly1.setFill(a1)
assembly2.setFill(a2)
root_cell.setFill(core)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
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
#plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
