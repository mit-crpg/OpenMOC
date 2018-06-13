import openmoc

###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

surfaces = {}

# Instantiate surfaces
surfaces['Root x-min']       = openmoc.XPlane(x=-32.13, name='Root x-min')
surfaces['Root x-max']       = openmoc.XPlane(x= 32.13, name='Root x-max')
surfaces['Root y-min']       = openmoc.YPlane(y=-32.13, name='Root y-min')
surfaces['Root y-max']       = openmoc.YPlane(y= 32.13, name='Root y-max')
surfaces['Root Small z-min'] = openmoc.ZPlane(z=-32.13, name='Root Small z-min')
surfaces['Root Small z-max'] = openmoc.ZPlane(z= 32.13, name='Root Small z-max')
surfaces['Root Big z-min']   = openmoc.ZPlane(z=-107.1, name='Root Big z-min')
surfaces['Root Big z-max']   = openmoc.ZPlane(z= 107.1, name='Root Big z-max')
surfaces['Fuel Cylinder']    = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.54, name='Fuel Cylinder')

surfaces['Root x-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root x-max'].setBoundaryType(openmoc.VACUUM)
surfaces['Root y-min'].setBoundaryType(openmoc.VACUUM)
surfaces['Root y-max'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root Small z-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root Small z-max'].setBoundaryType(openmoc.VACUUM)
surfaces['Root Big z-min'].setBoundaryType(openmoc.REFLECTIVE)
surfaces['Root Big z-max'].setBoundaryType(openmoc.VACUUM)
