import openmoc

###############################################################################
#                 Create dictionary of common surfaces
###############################################################################

surfaces = {}

# Instantiate surfaces
surfaces['x-min'] = openmoc.XPlane(x=-82.5, name='x-min')
surfaces['y-min'] = openmoc.YPlane(y=-82.5, name='y-min')
surfaces['x-max'] = openmoc.XPlane(x= 82.5, name='x-max')
surfaces['y-max'] = openmoc.YPlane(y= 82.5, name='y-max')

surfaces['x-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['y-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['x-max'].setBoundaryType(openmoc.VACUUM)
surfaces['y-max'].setBoundaryType(openmoc.VACUUM)
