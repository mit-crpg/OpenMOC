import openmoc

###############################################################################
#                 Create dictionary of common surfaces
###############################################################################

surfaces = {}

# Instantiate surfaces
surfaces['Root x-min'] = openmoc.XPlane(x=-82.5, name='Root x-min')
surfaces['Root y-min'] = openmoc.YPlane(y=-82.5, name='Root y-min')
surfaces['Root x-max'] = openmoc.XPlane(x= 82.5, name='Root x-max')
surfaces['Root y-max'] = openmoc.YPlane(y= 82.5, name='Root y-max')

surfaces['Root x-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root y-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root x-max'].setBoundaryType(openmoc.VACUUM)
surfaces['Root y-max'].setBoundaryType(openmoc.VACUUM)
