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

log.py_printf('TITLE', 'Simulating the LRA Benchmark Problem...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from py...')

materials = materialize.materialize('LRA-materials.py')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

planes = []
planes.append(XPlane(x=-82.5))
planes.append(XPlane(x=82.5))
planes.append(YPlane(y=-82.5))
planes.append(YPlane(y=82.5))
planes[0].setBoundaryType(REFLECTIVE)
planes[1].setBoundaryType(VACUUM)
planes[2].setBoundaryType(REFLECTIVE)
planes[3].setBoundaryType(VACUUM)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = []
cells.append(CellBasic())
cells.append(CellBasic())
cells.append(CellBasic())
cells.append(CellBasic())
cells.append(CellBasic())
cells.append(CellBasic())
cells.append(CellFill())
cells.append(CellFill())
cells.append(CellFill())
cells.append(CellFill())
cells.append(CellFill())
cells.append(CellFill())
cells.append(CellFill())

cells[0].setMaterial(materials['region_1'])
cells[1].setMaterial(materials['region_2'])
cells[2].setMaterial(materials['region_3'])
cells[3].setMaterial(materials['region_4'])
cells[4].setMaterial(materials['region_5'])
cells[5].setMaterial(materials['region_6'])

cells[12].addSurface(halfspace=+1, surface=planes[0])
cells[12].addSurface(halfspace=-1, surface=planes[1])
cells[12].addSurface(halfspace=+1, surface=planes[2])
cells[12].addSurface(halfspace=-1, surface=planes[3])


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

universes = list()
universes.append(Universe(name='region 1'))
universes.append(Universe(name='region 2'))
universes.append(Universe(name='region 3'))
universes.append(Universe(name='region 4'))
universes.append(Universe(name='region 5'))
universes.append(Universe(name='region 6'))

universes.append(Universe(name='assembly 1'))
universes.append(Universe(name='assembly 2'))
universes.append(Universe(name='assembly 3'))
universes.append(Universe(name='assembly 4'))
universes.append(Universe(name='assembly 5'))
universes.append(Universe(name='assembly 6'))
universes.append(Universe(name='core'))
root = Universe(name='root universe')

universes[0].addCell(cells[0])
universes[1].addCell(cells[1])
universes[2].addCell(cells[2])
universes[3].addCell(cells[3])
universes[4].addCell(cells[4])
universes[5].addCell(cells[5])
universes[6].addCell(cells[6])
universes[7].addCell(cells[7])
universes[8].addCell(cells[8])
universes[9].addCell(cells[9])
universes[10].addCell(cells[10])
universes[11].addCell(cells[11])
root.addCell(cells[12])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating LRA lattices...')

assembly1 = Lattice(name='assembly 1')
assembly1.setWidth(width_x=1.5, width_y=1.5)
template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

for i in range(10):
  for j in range(10):
    template[i][j] = universes[template[i][j]-1]
assembly1.setUniverses(template)


assembly2 = Lattice(name='assembly 2')
assembly2.setWidth(width_x=1.5, width_y=1.5)
template = [[2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]]

for i in range(10):
  for j in range(10):
    template[i][j] = universes[template[i][j]-1]
assembly2.setUniverses(template)


assembly3 = Lattice(name='assembly 2')
assembly3.setWidth(width_x=1.5, width_y=1.5)
template = [[3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]]

for i in range(10):
  for j in range(10):
    template[i][j] = universes[template[i][j]-1]
assembly3.setUniverses(template)


assembly4 = Lattice(name='assembly 4')
assembly4.setWidth(width_x=1.5, width_y=1.5)
template = [[4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]]

for i in range(10):
  for j in range(10):
    template[i][j] = universes[template[i][j]-1]
assembly4.setUniverses(template)

assembly5 = Lattice(name='assembly 5')
assembly5.setWidth(width_x=1.5, width_y=1.5)
template = [[5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]]

for i in range(10):
  for j in range(10):
    template[i][j] = universes[template[i][j]-1]
assembly5.setUniverses(template)


assembly6 = Lattice(name='assembly 6')
assembly6.setWidth(width_x=1.5, width_y=1.5)
template = [[6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6]]

for i in range(10):
  for j in range(10):
    template[i][j] = universes[template[i][j]-1]
assembly6.setUniverses(template)


core = Lattice(name='core')
core.setWidth(width_x=15.0, width_y=15.0)
template = [[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
            [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12],
            [9, 9, 9, 9, 9, 9, 9, 12, 12, 12, 12],
            [9, 9, 9, 9, 9, 9, 9, 10, 12, 12, 12],
            [8, 7, 7, 7, 7, 8, 8, 11, 11, 12, 12],
            [8, 7, 7, 7, 7, 8, 8, 11, 11, 12, 12],
            [7, 7, 7, 7, 7, 7, 7, 9, 9, 12, 12],
            [7, 7, 7, 7, 7, 7, 7, 9, 9, 12, 12],
            [7, 7, 7, 7, 7, 7, 7, 9, 9, 12, 12],
            [7, 7, 7, 7, 7, 7, 7, 9, 9, 12, 12],
            [8, 7, 7, 7, 7, 8, 8, 9, 9, 12, 12]]

for i in range(11):
  for j in range(11):
    template[i][j] = universes[template[i][j]-1]
core.setUniverses(template)


cells[6].setFill(assembly1)
cells[7].setFill(assembly2)
cells[8].setFill(assembly3)
cells[9].setFill(assembly4)
cells[10].setFill(assembly5)
cells[11].setFill(assembly6)
cells[12].setFill(core)


###############################################################################
##########################     Creating Cmfd mesh    ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating Cmfd mesh...')

cmfd = Cmfd()
cmfd.setLatticeStructure(110,110)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setCmfd(cmfd)
geometry.setRootUniverse(root)
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
solver.setSourceConvergenceThreshold(tolerance)
solver.setNumThreads(num_threads)
solver.convergeSource(max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_tracks(track_generator)
#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
#plotter.plot_fluxes(geometry, solver, energy_groups=[1,2,3,4,5,6,7])
#plotter.plot_mesh_fluxes(mesh, energy_groups=[1,2,3,4,5,6,7])
#plotter.plot_cmfd_cells(geometry, cmfd, gridsize=500)

log.py_printf('TITLE', 'Finished')

