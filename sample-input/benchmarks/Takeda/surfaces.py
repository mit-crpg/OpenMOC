from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter



###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

xmin = XPlane(x=-12.5, name='xmin')
xmax = XPlane(x= 12.5, name='xmax')
ymin = YPlane(y=-12.5, name='ymin')
ymax = YPlane(y= 12.5, name='ymax')
zmin = ZPlane(z=-12.5, name='zmin')
zmax = ZPlane(z= 12.5, name='zmax')

xmin.setBoundaryType(REFLECTIVE)
ymin.setBoundaryType(REFLECTIVE)
zmin.setBoundaryType(REFLECTIVE)
xmax.setBoundaryType(VACUUM)
ymax.setBoundaryType(VACUUM)
zmax.setBoundaryType(VACUUM)
