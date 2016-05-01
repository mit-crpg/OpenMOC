from openmoc import *
import openmoc.log as log
import openmoc.materialize as materialize


###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

xmin = XPlane(x=-10.71, name='xmin')
xmax = XPlane(x= 10.71, name='xmax')
ymin = YPlane(y=-10.71, name='ymin')
ymax = YPlane(y= 10.71, name='ymax')
zmin = ZPlane(z=-32.13, name='zmin')
zmax = ZPlane(z= 32.13, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
ymax.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
zmax.setBoundaryType(VACUUM)

# Create ZCylinders for the fuel as well as to discretize the moderator into rings
fuel_radius = ZCylinder(x=0.0, y=0.0, radius=0.54)
moderator_inner_radius = ZCylinder(x=0.0, y=0.0, radius=0.58)
moderator_outer_radius = ZCylinder(x=0.0, y=0.0, radius=0.62)
