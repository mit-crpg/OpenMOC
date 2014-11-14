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

log.py_printf('TITLE', 'Simulating the OECD\'s C5G7 Benchmark Problem...')


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../../c5g7-materials.h5')


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

circles = list()
planes = list()
planes.append(XPlane(x=-32.13, name='left'))
planes.append(XPlane(x=32.13, name='right'))
planes.append(YPlane(y=-32.13, name='bottom'))
planes.append(YPlane(y=32.13, name='top'))
circles.append(Circle(x=0., y=0., radius=0.54))
circles.append(Circle(x=0., y=0., radius=0.58))
circles.append(Circle(x=0., y=0., radius=0.62))
planes[0].setBoundaryType(REFLECTIVE)
planes[1].setBoundaryType(VACUUM)
planes[2].setBoundaryType(VACUUM)
planes[3].setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

cells = list()

# UO2 pin cells
cells.append(CellBasic(rings=3, sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells[0].setMaterial(materials['UO2'])
cells[1].setMaterial(materials['Water'])
cells[2].setMaterial(materials['Water'])
cells[3].setMaterial(materials['Water'])
cells[0].addSurface(-1, circles[0])
cells[1].addSurface(+1, circles[0])
cells[1].addSurface(-1, circles[1])
cells[2].addSurface(+1, circles[1])
cells[2].addSurface(-1, circles[2])
cells[3].addSurface(+1, circles[2])

# 4.3% MOX pin cells
cells.append(CellBasic(rings=3, sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells[4].setMaterial(materials['MOX-4.3%'])
cells[5].setMaterial(materials['Water'])
cells[6].setMaterial(materials['Water'])
cells[7].setMaterial(materials['Water'])
cells[4].addSurface(-1, circles[0])
cells[5].addSurface(+1, circles[0])
cells[5].addSurface(-1, circles[1])
cells[6].addSurface(+1, circles[1])
cells[6].addSurface(-1, circles[2])
cells[7].addSurface(+1, circles[2])

# 7% MOX pin cells
cells.append(CellBasic(rings=3, sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells[8].setMaterial(materials['MOX-7%'])
cells[9].setMaterial(materials['Water'])
cells[10].setMaterial(materials['Water'])
cells[11].setMaterial(materials['Water'])
cells[8].addSurface(-1, circles[0])
cells[9].addSurface(+1, circles[0])
cells[9].addSurface(-1, circles[1])
cells[10].addSurface(+1, circles[1])
cells[10].addSurface(-1, circles[2])
cells[11].addSurface(+1, circles[2])

# 8.7% MOX pin cells
cells.append(CellBasic(rings=3, sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells[12].setMaterial(materials['MOX-8.7%'])
cells[13].setMaterial(materials['Water'])
cells[14].setMaterial(materials['Water'])
cells[15].setMaterial(materials['Water'])
cells[12].addSurface(-1, circles[0])
cells[13].addSurface(+1, circles[0])
cells[13].addSurface(-1, circles[1])
cells[14].addSurface(+1, circles[1])
cells[14].addSurface(-1, circles[2])
cells[15].addSurface(+1, circles[2])

# Fission chamber pin cells
cells.append(CellBasic(rings=3, sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells[16].setMaterial(materials['Fission Chamber'])
cells[17].setMaterial(materials['Water'])
cells[18].setMaterial(materials['Water'])
cells[19].setMaterial(materials['Water'])
cells[16].addSurface(-1, circles[0])
cells[17].addSurface(+1, circles[0])
cells[17].addSurface(-1, circles[1])
cells[18].addSurface(+1, circles[1])
cells[18].addSurface(-1, circles[2])
cells[19].addSurface(+1, circles[2])

# Guide tube pin cells
cells.append(CellBasic(rings=3, sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells.append(CellBasic(sectors=8))
cells[20].setMaterial(materials['Guide Tube'])
cells[21].setMaterial(materials['Water'])
cells[22].setMaterial(materials['Water'])
cells[23].setMaterial(materials['Water'])
cells[20].addSurface(-1, circles[0])
cells[21].addSurface(+1, circles[0])
cells[21].addSurface(-1, circles[1])
cells[22].addSurface(+1, circles[1])
cells[22].addSurface(-1, circles[2])
cells[23].addSurface(+1, circles[2])

cells.append(CellBasic(name='moderator'))
cells[24].setMaterial(materials['Water'])
cells.append(CellFill(name='assembly 1'))
cells.append(CellFill(name='assembly 2'))
cells.append(CellFill(name='semi-finely spaced moderator'))
cells.append(CellFill(name='right reflector'))
cells.append(CellFill(name='bottom corner reflector'))
cells.append(CellFill(name='bottom reflector'))

cells.append(CellFill(name='full geometry'))
cells[-1].addSurface(+1, planes[0])
cells[-1].addSurface(-1, planes[1])
cells[-1].addSurface(+1, planes[2])
cells[-1].addSurface(-1, planes[3])


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

universes = list()
universes.append(Universe(name='uo2 pin cell'))
universes.append(Universe(name='mox-4.3% pin cell'))
universes.append(Universe(name='mox-7% pin cell'))
universes.append(Universe(name='mox-8.7% pin cell'))
universes.append(Universe(name='guide tube'))
universes.append(Universe(name='fission chamber'))
universes.append(Universe(name='moderator'))
universes.append(Universe(name='assembly 1'))
universes.append(Universe(name='assembly 2'))
universes.append(Universe(name='semi-finely spaced moderator'))
universes.append(Universe(name='right reflector'))
universes.append(Universe(name='bottom corner reflector'))
universes.append(Universe(name='bottom reflector'))
root = Universe(name='root universe')

universes[0].addCell(cells[0])
universes[0].addCell(cells[1])
universes[0].addCell(cells[2])
universes[0].addCell(cells[3])
universes[1].addCell(cells[4])
universes[1].addCell(cells[5])
universes[1].addCell(cells[6])
universes[1].addCell(cells[7])
universes[2].addCell(cells[8])
universes[2].addCell(cells[9])
universes[2].addCell(cells[10])
universes[2].addCell(cells[11])
universes[3].addCell(cells[12])
universes[3].addCell(cells[13])
universes[3].addCell(cells[14])
universes[3].addCell(cells[15])
universes[4].addCell(cells[16])
universes[4].addCell(cells[17])
universes[4].addCell(cells[18])
universes[4].addCell(cells[19])
universes[5].addCell(cells[20])
universes[5].addCell(cells[21])
universes[5].addCell(cells[22])
universes[5].addCell(cells[23])
universes[6].addCell(cells[24])
universes[7].addCell(cells[25])
universes[8].addCell(cells[26])
universes[9].addCell(cells[27])
universes[10].addCell(cells[28])
universes[11].addCell(cells[29])
universes[12].addCell(cells[30])
root.addCell(cells[31])


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating lattices...')

lattices = list()

# Top left, bottom right 17 x 17 assemblies
lattices.append(Lattice(name='assembly 1'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
            [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 6, 1, 1, 6, 1, 1, 5, 1, 1, 6, 1, 1, 6, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
            [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]-1]
lattices[-1].setUniverses(template)

# Top right, bottom left 17 x 17 assemblies
lattices.append(Lattice(name='assembly 2'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
            [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
            [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
            [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
            [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 6, 4, 4, 6, 4, 4, 5, 4, 4, 6, 4, 4, 6, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
            [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
            [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
            [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
            [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]]
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]-1]
lattices[-1].setUniverses(template)

# Sliced up water cells - semi finely spaced
lattices.append(Lattice(name='semi-finely spaced moderator'))
lattices[-1].setWidth(width_x=0.126, width_y=0.126)
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
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses(template)

# Sliced up water cells - right side of geometry
lattices.append(Lattice(name='right reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6]]
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses(template)

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='bottom corner reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]]
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses(template)

# Sliced up water cells for bottom of geometry
lattices.append(Lattice(name='bottom reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]]
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses(template)

# 4 x 4 core to represent two bundles and water
lattices.append(Lattice(name='full geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42)
lattices[-1].setUniverses([[universes[7], universes[8], universes[10]],
                           [universes[8], universes[7], universes[10]],
                           [universes[12], universes[12], universes[11]]])

cells[25].setFill(lattices[0])
cells[26].setFill(lattices[1])
cells[27].setFill(lattices[2])
cells[28].setFill(lattices[3])
cells[29].setFill(lattices[4])
cells[30].setFill(lattices[5])
cells[31].setFill(lattices[6])


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
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

log.py_printf('TITLE', 'Finished')
