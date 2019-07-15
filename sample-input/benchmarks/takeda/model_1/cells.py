import openmoc
import openmoc.materialize as materialize
from surfaces import surfaces

###############################################################################
#                           Creating Materials
###############################################################################

materials = materialize.load_from_hdf5('takeda-mgxs.h5', '')

###############################################################################
#                   Create dictionary of all cells
###############################################################################

cells = {}

# Instantiate Cells
cells['Core']        = openmoc.Cell(name='Core')
cells['Control Rod'] = openmoc.Cell(name='Control Rod')
cells['Void']        = openmoc.Cell(name='Void')
cells['Reflector']   = openmoc.Cell(name='Reflector')
cells['Root']        = openmoc.Cell(name='Root')

# Register Materials with Cells
cells['Core'].setFill(materials['Core'])
cells['Void'].setFill(materials['Void'])
cells['Control Rod'].setFill(materials['Control Rod'])
cells['Reflector'].setFill(materials['Reflector'])

# Add surfaces to root cell
cells['Root'].addSurface(halfspace=+1, surface=surfaces['Root x-min'])
cells['Root'].addSurface(halfspace=-1, surface=surfaces['Root x-max'])
cells['Root'].addSurface(halfspace=+1, surface=surfaces['Root y-min'])
cells['Root'].addSurface(halfspace=-1, surface=surfaces['Root y-max'])
cells['Root'].addSurface(halfspace=+1, surface=surfaces['Root z-min'])
cells['Root'].addSurface(halfspace=-1, surface=surfaces['Root z-max'])
