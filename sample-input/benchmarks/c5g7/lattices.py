import openmoc
from universes import universes, cells, surfaces

###############################################################################
#########################   Set Simulation Param   ############################
###############################################################################

reflector_refines = 5

###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

lattices = {}

# Instantiate Lattices
lattices['Refined Reflector Mesh']      = openmoc.Lattice()
lattices['Reflector Unrodded Assembly'] = openmoc.Lattice()
lattices['Reflector Rodded Assembly']   = openmoc.Lattice()
lattices['Reflector Right Assembly']    = openmoc.Lattice()
lattices['Reflector Bottom Assembly']   = openmoc.Lattice()
lattices['Reflector Corner Assembly']   = openmoc.Lattice()
lattices['UO2 Unrodded Assembly']       = openmoc.Lattice()
lattices['UO2 Rodded Assembly']         = openmoc.Lattice()
lattices['MOX Unrodded Assembly']       = openmoc.Lattice()
lattices['MOX Rodded Assembly']         = openmoc.Lattice()
lattices['Root']                        = openmoc.Lattice()

# Abbreviate universes that will fill lattices
u = universes['UO2']
m = universes['MOX 4.3%']
o = universes['MOX 7.0%']
x = universes['MOX 8.7%']
g = universes['Guide Tube']
f = universes['Fission Chamber']
c = universes['Control Rod']
p = universes['Moderator Pin']
r = universes['Reflector']
a = universes['Refined Reflector Mesh']

# Sliced up water cells - semi finely spaced
width_xy = 1.26 / reflector_refines
lattices['Refined Reflector Mesh'].setWidth\
    (width_x=width_xy, width_y=width_xy, width_z=100.)
template = [[[r] * reflector_refines] * reflector_refines]
lattices['Refined Reflector Mesh'].setUniverses(template)


# UO2 unrodded 17 x 17 assemblies
lattices['UO2 Unrodded Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[[u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, g, u, u, g, u, u, g, u, u, u, u, u],
             [u, u, u, g, u, u, u, u, u, u, u, u, u, g, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, g, u, u, g, u, u, g, u, u, g, u, u, g, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, g, u, u, g, u, u, f, u, u, g, u, u, g, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, g, u, u, g, u, u, g, u, u, g, u, u, g, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, g, u, u, u, u, u, u, u, u, u, g, u, u, u],
             [u, u, u, u, u, g, u, u, g, u, u, g, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u]]]

lattices['UO2 Unrodded Assembly'].setUniverses(template)

# UO2 rodded 17 x 17 assemblies
lattices['UO2 Rodded Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[[u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, c, u, u, c, u, u, c, u, u, u, u, u],
             [u, u, u, c, u, u, u, u, u, u, u, u, u, c, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, c, u, u, c, u, u, c, u, u, c, u, u, c, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, c, u, u, c, u, u, f, u, u, c, u, u, c, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, c, u, u, c, u, u, c, u, u, c, u, u, c, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, c, u, u, u, u, u, u, u, u, u, c, u, u, u],
             [u, u, u, u, u, c, u, u, c, u, u, c, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
             [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u]]]

lattices['UO2 Rodded Assembly'].setUniverses(template)

# MOX unrodded 17 x 17 assemblies
lattices['MOX Unrodded Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[[m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m],
             [m, o, o, o, o, o, o, o, o, o, o, o, o, o, o, o, m],
             [m, o, o, o, o, g, o, o, g, o, o, g, o, o, o, o, m],
             [m, o, o, g, o, x, x, x, x, x, x, x, o, g, o, o, m],
             [m, o, o, o, x, x, x, x, x, x, x, x, x, o, o, o, m],
             [m, o, g, x, x, g, x, x, g, x, x, g, x, x, g, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, g, x, x, g, x, x, f, x, x, g, x, x, g, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, g, x, x, g, x, x, g, x, x, g, x, x, g, o, m],
             [m, o, o, o, x, x, x, x, x, x, x, x, x, o, o, o, m],
             [m, o, o, g, o, x, x, x, x, x, x, x, o, g, o, o, m],
             [m, o, o, o, o, g, o, o, g, o, o, g, o, o, o, o, m],
             [m, o, o, o, o, o, o, o, o, o, o, o, o, o, o, o, m],
             [m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m]]]

lattices['MOX Unrodded Assembly'].setUniverses(template)

# MOX rodded 17 x 17 assemblies
lattices['MOX Rodded Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[[m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m],
             [m, o, o, o, o, o, o, o, o, o, o, o, o, o, o, o, m],
             [m, o, o, o, o, c, o, o, c, o, o, c, o, o, o, o, m],
             [m, o, o, c, o, x, x, x, x, x, x, x, o, c, o, o, m],
             [m, o, o, o, x, x, x, x, x, x, x, x, x, o, o, o, m],
             [m, o, c, x, x, c, x, x, c, x, x, c, x, x, c, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, c, x, x, c, x, x, f, x, x, c, x, x, c, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, o, x, x, x, x, x, x, x, x, x, x, x, o, o, m],
             [m, o, c, x, x, c, x, x, c, x, x, c, x, x, c, o, m],
             [m, o, o, o, x, x, x, x, x, x, x, x, x, o, o, o, m],
             [m, o, o, c, o, x, x, x, x, x, x, x, o, c, o, o, m],
             [m, o, o, o, o, c, o, o, c, o, o, c, o, o, o, o, m],
             [m, o, o, o, o, o, o, o, o, o, o, o, o, o, o, o, m],
             [m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m]]]

lattices['MOX Rodded Assembly'].setUniverses(template)

# Reflector unrodded 17 x 17 assemblies
lattices['Reflector Unrodded Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[[p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, g, p, p, g, p, p, g, p, p, p, p, p],
             [p, p, p, g, p, p, p, p, p, p, p, p, p, g, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, g, p, p, g, p, p, g, p, p, g, p, p, g, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, g, p, p, g, p, p, f, p, p, g, p, p, g, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, g, p, p, g, p, p, g, p, p, g, p, p, g, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, g, p, p, p, p, p, p, p, p, p, g, p, p, p],
             [p, p, p, p, p, g, p, p, g, p, p, g, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p]]]

lattices['Reflector Unrodded Assembly'].setUniverses(template)

# Reflector rodded 17 x 17 assemblies
lattices['Reflector Rodded Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[[p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, c, p, p, c, p, p, c, p, p, p, p, p],
             [p, p, p, c, p, p, p, p, p, p, p, p, p, c, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, c, p, p, c, p, p, c, p, p, c, p, p, c, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, c, p, p, c, p, p, f, p, p, c, p, p, c, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, c, p, p, c, p, p, c, p, p, c, p, p, c, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, c, p, p, p, p, p, p, p, p, p, c, p, p, p],
             [p, p, p, p, p, c, p, p, c, p, p, c, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p],
             [p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p, p]]]

lattices['Reflector Rodded Assembly'].setUniverses(template)

# Reflector right 17 x 17 assemblies
lattices['Reflector Right Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[[a] * 11 + [r] * 6] * 17]
lattices['Reflector Right Assembly'].setUniverses(template)

# Reflector bottom 17 x 17 assemblies
lattices['Reflector Bottom Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[a] * 17] * 11
template += [[r] * 17] * 6
template = [template]
lattices['Reflector Bottom Assembly'].setUniverses(template)

# Reflector corner 17 x 17 assemblies
lattices['Reflector Corner Assembly'].setWidth(width_x=1.26, width_y=1.26, width_z=100.)
template = [[a] * 11 + [r] * 6] * 11
template += [[r] * 17] * 6
template = [template]
lattices['Reflector Corner Assembly'].setUniverses(template)


# Fill cells with lattices
cells['Refined Reflector Mesh']     .setFill(lattices['Refined Reflector Mesh'])
cells['UO2 Unrodded Assembly']      .setFill(lattices['UO2 Unrodded Assembly'])
cells['UO2 Rodded Assembly']        .setFill(lattices['UO2 Rodded Assembly'])
cells['MOX Unrodded Assembly']      .setFill(lattices['MOX Unrodded Assembly'])
cells['MOX Rodded Assembly']        .setFill(lattices['MOX Rodded Assembly'])
cells['Reflector Unrodded Assembly'].setFill(lattices['Reflector Unrodded Assembly'])
cells['Reflector Rodded Assembly']  .setFill(lattices['Reflector Rodded Assembly'])
cells['Reflector Right Assembly']   .setFill(lattices['Reflector Right Assembly'])
cells['Reflector Bottom Assembly']  .setFill(lattices['Reflector Bottom Assembly'])
cells['Reflector Corner Assembly']  .setFill(lattices['Reflector Corner Assembly'])
cells['Root']                       .setFill(lattices['Root'])
