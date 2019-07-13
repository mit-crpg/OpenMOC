import openmoc
from universes import universes, cells

###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

lattices = {}

# Instantiate Lattices
lattices['Root'] = openmoc.Lattice()

# Fill cells with lattices
cells['Root'].setFill(lattices['Root'])
