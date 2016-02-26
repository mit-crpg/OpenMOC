import openmoc
from universes import universes, cells

###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

# The number of uniform refinements to perform on the assembly lattices
refines = 10

lattices = {}

# Instantiate Lattices
lattices['Region 1 Assembly'] = openmoc.Lattice()
lattices['Region 2 Assembly'] = openmoc.Lattice()
lattices['Region 3 Assembly'] = openmoc.Lattice()
lattices['Region 4 Assembly'] = openmoc.Lattice()
lattices['Region 5 Assembly'] = openmoc.Lattice()
lattices['Root'] = openmoc.Lattice()

u1 = universes['Region 1']
u2 = universes['Region 2']
u3 = universes['Region 3']
u4 = universes['Region 4']
u5 = universes['Region 5']

lattices['Region 1 Assembly'].setWidth(15.0/refines, 15.0/refines)
lattices['Region 1 Assembly'].setUniverses([[[u1] * refines] * refines])

lattices['Region 2 Assembly'].setWidth(15.0/refines, 15.0/refines)
lattices['Region 2 Assembly'].setUniverses([[[u2] * refines] * refines])

lattices['Region 3 Assembly'].setWidth(15.0/refines, 15.0/refines)
lattices['Region 3 Assembly'].setUniverses([[[u3] * refines] * refines])

lattices['Region 4 Assembly'].setWidth(15.0/refines, 15.0/refines)
lattices['Region 4 Assembly'].setUniverses([[[u4] * refines] * refines])

lattices['Region 5 Assembly'].setWidth(15.0/refines, 15.0/refines)
lattices['Region 5 Assembly'].setUniverses([[[u5] * refines] * refines])

u1 = universes['Region 1 Assembly']
u2 = universes['Region 2 Assembly']
u3 = universes['Region 3 Assembly']
u4 = universes['Region 4 Assembly']
u5 = universes['Region 5 Assembly']

lattices['Root'].setWidth(15.0, 15.0)
template = [[u5] * 11] * 2 + \
           [[u3] * 7 + [u5] * 4] + \
           [[u3] * 7 + [u4] + [u5] * 3] + \
           [[u2] + [u1] * 4 + [u2] * 2 + [u3] * 2 + [u5] * 2] * 2 + \
           [[u1] * 7 + [u3] * 2 + [u5] * 2] * 4 + \
           [[u2] + [u1] * 4 + [u2] * 2 + [u3] * 2 + [u5] * 2]
lattices['Root'].setUniverses([template])

# Fill cells with lattices
cells['Region 1 Assembly'].setFill(lattices['Region 1 Assembly'])
cells['Region 2 Assembly'].setFill(lattices['Region 2 Assembly'])
cells['Region 3 Assembly'].setFill(lattices['Region 3 Assembly'])
cells['Region 4 Assembly'].setFill(lattices['Region 4 Assembly'])
cells['Region 5 Assembly'].setFill(lattices['Region 5 Assembly'])
cells['Root'].setFill(lattices['Root'])
