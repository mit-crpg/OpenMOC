import openmoc
from surfaces import surfaces
import openmoc.materialize as materialize

###############################################################################
#                           Creating Materials
###############################################################################

materials = openmoc.materialize.load_from_hdf5('LRA-mgxs.h5', './')

###############################################################################
#                   Create dictionary of all cells
###############################################################################

cells = {}

# Instantiate Cells
cells['Region 1']          = openmoc.Cell(name='Region 1')
cells['Region 2']          = openmoc.Cell(name='Region 2')
cells['Region 3']          = openmoc.Cell(name='Region 3')
cells['Region 4']          = openmoc.Cell(name='Region 4')
cells['Region 5']          = openmoc.Cell(name='Region 5')
cells['Region 1 Assembly'] = openmoc.Cell(name='Region 1 Assembly')
cells['Region 2 Assembly'] = openmoc.Cell(name='Region 2 Assembly')
cells['Region 3 Assembly'] = openmoc.Cell(name='Region 3 Assembly')
cells['Region 4 Assembly'] = openmoc.Cell(name='Region 4 Assembly')
cells['Region 5 Assembly'] = openmoc.Cell(name='Region 5 Assembly')
cells['Root']              = openmoc.Cell(name='Root')

# Register Materials with Cells
cells['Region 1'].setFill(materials['Region 1'])
cells['Region 2'].setFill(materials['Region 2'])
cells['Region 3'].setFill(materials['Region 3'])
cells['Region 4'].setFill(materials['Region 4'])
cells['Region 5'].setFill(materials['Region 5'])

# Add surfaces to root cell
cells['Root'].addSurface(halfspace=+1, surface=surfaces['x-min'])
cells['Root'].addSurface(halfspace=-1, surface=surfaces['x-max'])
cells['Root'].addSurface(halfspace=+1, surface=surfaces['y-min'])
cells['Root'].addSurface(halfspace=-1, surface=surfaces['y-max'])
