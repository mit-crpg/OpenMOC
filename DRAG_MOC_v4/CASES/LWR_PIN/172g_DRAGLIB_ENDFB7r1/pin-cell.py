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
num_azim=64
track_spacing=0.05
log.set_log_level('NORMAL')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('XS_isot_pin.py')

uo2_id = materials['Fuel'].getId()
water_id = materials['Moderator'].getId()
clad_id = materials['Cladding'].getId()


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circle1 = Circle(x=0.0, y=0.0, radius=0.412660)
circle2 = Circle(x=0.0, y=0.0, radius=0.474364)

left = XPlane(x=-0.63236)
right = XPlane(x=0.63236)
top = YPlane(y=0.63236)
bottom = YPlane(y=-0.63236)

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic(universe=1, material=uo2_id, rings=3, sectors=8))
cells.append(CellBasic(universe=1, material=clad_id, sectors=8))
cells.append(CellBasic(universe=1, material=water_id, sectors=8))
cells.append(CellFill(universe=0, universe_fill=3))

cells[0].addSurface(halfspace=-1, surface=circle1)
cells[1].addSurface(halfspace=+1, surface=circle1)
cells[1].addSurface(halfspace=-1, surface=circle2)
cells[2].addSurface(halfspace=+1, surface=circle2)


cells[3].addSurface(halfspace=+1, surface=left)
cells[3].addSurface(halfspace=-1, surface=right)
cells[3].addSurface(halfspace=+1, surface=bottom)
cells[3].addSurface(halfspace=-1, surface=top)


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating simple pin cell lattice...')

lattice = Lattice(id=3, width_x=1.26472, width_y=1.26472)
lattice.setLatticeCells([[1]])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
for material in materials.values(): geometry.addMaterial(material)
for cell in cells: geometry.addCell(cell)
geometry.addLattice(lattice)

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

solver = ThreadPrivateSolver(geometry, track_generator)
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
