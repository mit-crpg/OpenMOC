from openmoc import *
import openmoc.log as log
import openmoc.materialize as materialize


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

xmin = XPlane(x=-32.13, name='xmin')
xmax = XPlane(x= 32.13, name='xmax')
ymin = YPlane(y=-32.13, name='ymin')
ymax = YPlane(y= 32.13, name='ymax')
zmin = ZPlane(z=-1.0, name='zmin')
zmax = ZPlane(z= 1.0, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(VACUUM)
ymin.setBoundaryType(VACUUM)
ymax.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(REFLECTIVE)

# Create Circles for the fuel as well as to discretize the moderator into rings
fuel_radius = Circle(x=0.0, y=0.0, radius=0.54)
moderator_inner_radius = Circle(x=0.0, y=0.0, radius=0.58)
moderator_outer_radius = Circle(x=0.0, y=0.0, radius=0.62)
