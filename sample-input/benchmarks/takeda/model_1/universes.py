import openmoc
from cells import cells

###############################################################################
#                            Creating Universes
###############################################################################

universes = {}

universes['Core']        = openmoc.Universe(name='Core')
universes['Reflector']   = openmoc.Universe(name='Reflector')
universes['Control Rod'] = openmoc.Universe(name='Control Rod')
universes['Void']        = openmoc.Universe(name='Void')
universes['Root']        = openmoc.Universe(name='Root')

# Add cells to universes
universes['Core']       .addCell(cells['Core'])
universes['Control Rod'].addCell(cells['Control Rod'])
universes['Reflector']  .addCell(cells['Reflector'])
universes['Void']       .addCell(cells['Void'])
universes['Root']       .addCell(cells['Root'])
