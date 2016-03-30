import openmoc
from cells import cells

###############################################################################
#                            Creating Universes
###############################################################################

universes = {}

universes['Region 1']          = openmoc.Universe(name='Region 1')
universes['Region 2']          = openmoc.Universe(name='Region 2')
universes['Region 3']          = openmoc.Universe(name='Region 3')
universes['Region 4']          = openmoc.Universe(name='Region 4')
universes['Region 5']          = openmoc.Universe(name='Region 5')
universes['Region 1 Assembly'] = openmoc.Universe(name='Region 1 Assembly')
universes['Region 2 Assembly'] = openmoc.Universe(name='Region 2 Assembly')
universes['Region 3 Assembly'] = openmoc.Universe(name='Region 3 Assembly')
universes['Region 4 Assembly'] = openmoc.Universe(name='Region 4 Assembly')
universes['Region 5 Assembly'] = openmoc.Universe(name='Region 5 Assembly')
universes['Root']              = openmoc.Universe(name='Root')

# Add cells to universes
universes['Region 1']         .addCell(cells['Region 1'])
universes['Region 2']         .addCell(cells['Region 2'])
universes['Region 3']         .addCell(cells['Region 3'])
universes['Region 4']         .addCell(cells['Region 4'])
universes['Region 5']         .addCell(cells['Region 5'])
universes['Region 1 Assembly'].addCell(cells['Region 1 Assembly'])
universes['Region 2 Assembly'].addCell(cells['Region 2 Assembly'])
universes['Region 3 Assembly'].addCell(cells['Region 3 Assembly'])
universes['Region 4 Assembly'].addCell(cells['Region 4 Assembly'])
universes['Region 5 Assembly'].addCell(cells['Region 5 Assembly'])
universes['Root']             .addCell(cells['Root'])
