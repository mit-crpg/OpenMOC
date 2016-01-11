import openmoc

###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

surfaces = {}

# Instantiate surfaces
surfaces['Root x-min']    = openmoc.XPlane(x=-32.13, name='Root x-min')
surfaces['Root x-max']    = openmoc.XPlane(x= 32.13, name='Root x-max')
surfaces['Root y-min']    = openmoc.YPlane(y=-32.13, name='Root y-min')
surfaces['Root y-max']    = openmoc.YPlane(y= 32.13, name='Root y-max')
surfaces['Root z-min']    = openmoc.ZPlane(z=-32.13, name='Root z-min')
surfaces['Root z-max']    = openmoc.ZPlane(z= 32.13, name='Root z-max')
surfaces['Fuel Cylinder'] = openmoc.Circle(x=0.0, y=0.0, radius=0.54, name='Fuel Cylinder')

surfaces['Root x-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root x-max'].setBoundaryType(openmoc.VACUUM)
surfaces['Root y-min'].setBoundaryType(openmoc.VACUUM)
surfaces['Root y-max'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root z-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root z-max'].setBoundaryType(openmoc.VACUUM)
