import openmoc
import openmoc.materialize as materialize
from surfaces import surfaces
from surfaces import gap

fuel_rings      = 5
moderator_rings = 0
num_sectors     = 4

###############################################################################
###########################   Creating Materials   ############################
###############################################################################

materials = materialize.load_from_hdf5('c5g7-mgxs.h5', '../../')

###############################################################################
######################   Creating Cells and Universes   #######################
###############################################################################

cells = {}

# Instantiate Cells
cells['Root']                        = openmoc.Cell()
cells['Gap Root']                    = openmoc.Cell()
cells['UO2']                         = openmoc.Cell()
cells['MOX 4.3%']                    = openmoc.Cell()
cells['MOX 7.0%']                    = openmoc.Cell()
cells['MOX 8.7%']                    = openmoc.Cell()
cells['Guide Tube']                  = openmoc.Cell()
cells['Fission Chamber']             = openmoc.Cell()
cells['Control Rod']                 = openmoc.Cell()
cells['Moderator']                   = openmoc.Cell()
cells['Reflector']                   = openmoc.Cell()
cells['Moderator in Pin']            = openmoc.Cell()
cells['Refined Reflector Pin']       = openmoc.Cell()
cells['Refined Reflector Mesh']      = openmoc.Cell()
cells['UO2 Unrodded Assembly']       = openmoc.Cell()
cells['UO2 Rodded Assembly']         = openmoc.Cell()
cells['MOX Unrodded Assembly']       = openmoc.Cell()
cells['MOX Rodded Assembly']         = openmoc.Cell()
cells['Reflector Unrodded Assembly'] = openmoc.Cell()
cells['Reflector Rodded Assembly']   = openmoc.Cell()
cells['Reflector Right Assembly']    = openmoc.Cell()
cells['Reflector Bottom Assembly']   = openmoc.Cell()
cells['Reflector Corner Assembly']   = openmoc.Cell()
cells['Reflector Assembly']          = openmoc.Cell()

cells['Gap UO2 Unrodded Assembly']       = openmoc.Cell()
cells['Gap UO2 Rodded Assembly']         = openmoc.Cell()
cells['Gap MOX Unrodded Assembly']       = openmoc.Cell()
cells['Gap MOX Rodded Assembly']         = openmoc.Cell()
cells['Gap Reflector Rodded Assembly']   = openmoc.Cell()
cells['Gap Reflector Right Assembly']    = openmoc.Cell()
cells['Gap Reflector Bottom Assembly']   = openmoc.Cell()
cells['Gap Reflector Corner Assembly']   = openmoc.Cell()

# Set material fills
cells['UO2']             .setFill(materials['UO2'])
cells['MOX 4.3%']        .setFill(materials['MOX-4.3%'])
cells['MOX 7.0%']        .setFill(materials['MOX-7%'])
cells['MOX 8.7%']        .setFill(materials['MOX-8.7%'])
cells['Guide Tube']      .setFill(materials['Guide Tube'])
cells['Fission Chamber'] .setFill(materials['Fission Chamber'])
cells['Control Rod']     .setFill(materials['Control Rod'])
cells['Moderator']       .setFill(materials['Water'])
cells['Reflector']       .setFill(materials['Water'])
cells['Moderator in Pin'].setFill(materials['Water'])

# Set rings and sectors
'''
cells['UO2']             .setNumRings(fuel_rings)
cells['MOX 4.3%']        .setNumRings(fuel_rings)
cells['MOX 7.0%']        .setNumRings(fuel_rings)
cells['MOX 8.7%']        .setNumRings(fuel_rings)
'''
cells['Guide Tube']      .setNumRings(fuel_rings)
cells['Fission Chamber'] .setNumRings(fuel_rings)
cells['Control Rod']     .setNumRings(fuel_rings)
cells['Moderator in Pin'].setNumRings(fuel_rings)
'''
cells['Moderator']       .setNumRings(moderator_rings)
'''
cells['UO2']             .setNumSectors(num_sectors)
cells['MOX 4.3%']        .setNumSectors(num_sectors)
cells['MOX 7.0%']        .setNumSectors(num_sectors)
cells['MOX 8.7%']        .setNumSectors(num_sectors)
cells['Guide Tube']      .setNumSectors(num_sectors)
cells['Fission Chamber'] .setNumSectors(num_sectors)
cells['Control Rod']     .setNumSectors(num_sectors)
cells['Moderator in Pin'].setNumSectors(num_sectors)
cells['Moderator']       .setNumSectors(8)

# Add surfaces
cells['UO2']             .addSurface(-1, surfaces['Fuel Cylinder'])
cells['MOX 4.3%']        .addSurface(-1, surfaces['Fuel Cylinder'])
cells['MOX 7.0%']        .addSurface(-1, surfaces['Fuel Cylinder'])
cells['MOX 8.7%']        .addSurface(-1, surfaces['Fuel Cylinder'])
cells['Guide Tube']      .addSurface(-1, surfaces['Fuel Cylinder'])
cells['Fission Chamber'] .addSurface(-1, surfaces['Fuel Cylinder'])
cells['Control Rod']     .addSurface(-1, surfaces['Fuel Cylinder'])
cells['Moderator in Pin'].addSurface(-1, surfaces['Fuel Cylinder'])
cells['Moderator']       .addSurface(+1, surfaces['Fuel Cylinder'])
cells['Root']            .addSurface(+1, surfaces['Root x-min'])
cells['Root']            .addSurface(-1, surfaces['Root x-max'])
cells['Root']            .addSurface(+1, surfaces['Root y-min'])
cells['Root']            .addSurface(-1, surfaces['Root y-max'])

cells['Gap Root']        .addSurface(+1, surfaces['Gap Root x-min'])
cells['Gap Root']        .addSurface(-1, surfaces['Gap Root x-max'])
cells['Gap Root']        .addSurface(+1, surfaces['Gap Root y-min'])
cells['Gap Root']        .addSurface(-1, surfaces['Gap Root y-max'])