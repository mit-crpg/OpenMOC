from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
from openmoc.options import Options
import h5py as h5
import numpy as np

###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()

log.set_log_level('NORMAL')


###############################################################################
###########################   Creating Isotopes    ############################
###############################################################################

log.py_printf('NORMAL', 'Importing isotope data from HDF5...')

mgxs = h5.File('mgxs-8-groups.h5')
mgxs = mgxs['material']

# O-16 in Fuel
o16f = Isotope(name='O-16 in Fuel')
o16f.setNumEnergyGroups(8)
xs = mgxs['material 10000']['O-16']
o16f.setSigmaT(np.array(xs['transport']['average'][:]))
o16f.setSigmaA(np.array(xs['absorption']['average'][:]))
o16f.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))

# U-235
u235 = Isotope(name='U-235')
u235.setNumEnergyGroups(8)
xs = mgxs['material 10000']['U-235']
u235.setSigmaT(np.array(xs['transport']['average'][:]))
u235.setSigmaA(np.array(xs['absorption']['average'][:]))
u235.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))
u235.setSigmaF(np.array(xs['fission']['average'][:]))
u235.setNuSigmaF(np.array(xs['nu-fission']['average'][:]))
u235.setChi(np.array(xs['chi']['average'][:]))

# U-238
u238 = Isotope(name='U-238')
u238.setNumEnergyGroups(8)
xs = mgxs['material 10000']['U-238']
u238.setSigmaT(np.array(xs['transport']['average'][:]))
u238.setSigmaA(np.array(xs['absorption']['average'][:]))
u238.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))
u238.setSigmaF(np.array(xs['fission']['average'][:]))
u238.setNuSigmaF(np.array(xs['nu-fission']['average'][:]))
u238.setChi(np.array(xs['chi']['average'][:]))

# H-1
h1 = Isotope(name='H-1')
h1.setNumEnergyGroups(8)
xs = mgxs['material 10003']['H-1']
h1.setSigmaT(np.array(xs['transport']['average'][:]))
h1.setSigmaA(np.array(xs['absorption']['average'][:]))
h1.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))


# O-16 in Moderator
o16m = Isotope(name='O-16 in Moderator')
o16m.setNumEnergyGroups(8)
xs = mgxs['material 10003']['O-16']
o16m.setSigmaT(np.array(xs['transport']['average'][:]))
o16m.setSigmaA(np.array(xs['absorption']['average'][:]))
o16m.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))


# B-10
b10 = Isotope(name='B-10')
b10.setNumEnergyGroups(8)
xs = mgxs['material 10003']['B-10']
b10.setSigmaT(np.array(xs['transport']['average'][:]))
b10.setSigmaA(np.array(xs['absorption']['average'][:]))
b10.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))

# B-11
b11 = Isotope(name='B-11')
b11.setNumEnergyGroups(8)
xs = mgxs['material 10003']['B-11']
b11.setSigmaT(np.array(xs['transport']['average'][:]))
b11.setSigmaA(np.array(xs['absorption']['average'][:]))
b11.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))



# Zr-90
zr90 = Isotope(name='Zr-90')
zr90.setNumEnergyGroups(8)
xs = mgxs['material 10005']['Zr-90']
zr90.setSigmaT(np.array(xs['transport']['average'][:]))
zr90.setSigmaA(np.array(xs['absorption']['average'][:]))
zr90.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))


# Zr-91
zr91 = Isotope(name='Zr-91')
zr91.setNumEnergyGroups(8)
xs = mgxs['material 10005']['Zr-91']
zr91.setSigmaT(np.array(xs['transport']['average'][:]))
zr91.setSigmaA(np.array(xs['absorption']['average'][:]))
zr91.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))


# Zr-92
zr92 = Isotope(name='Zr-92')
zr92.setNumEnergyGroups(8)
xs = mgxs['material 10005']['Zr-92']
zr92.setSigmaT(np.array(xs['transport']['average'][:]))
zr92.setSigmaA(np.array(xs['absorption']['average'][:]))
zr92.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))


# Zr-94
zr94 = Isotope(name='Zr-94')
zr94.setNumEnergyGroups(8)
xs = mgxs['material 10005']['Zr-94']
zr94.setSigmaT(np.array(xs['transport']['average'][:]))
zr94.setSigmaA(np.array(xs['absorption']['average'][:]))
zr94.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))


# Zr-96
zr96 = Isotope(name='Zr-96')
zr96.setNumEnergyGroups(8)
xs = mgxs['material 10005']['Zr-96']
zr96.setSigmaT(np.array(xs['transport']['average'][:]))
zr96.setSigmaA(np.array(xs['absorption']['average'][:]))
zr96.setSigmaS(np.array(xs['scatter matrix']['average'][:]).flatten('F'))


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating materials...')

# Fuel
fuel_matl = IsoMaterial(name='Fuel')
fuel_matl.addIsotope(u235,3.7503e-4)
fuel_matl.addIsotope(u238,2.2625e-2)
fuel_matl.addIsotope(o16f,4.5895e-2)

# Moderator
mod_matl = IsoMaterial(name='Moderator')
mod_matl.addIsotope(b10,8.0042e-6)
mod_matl.addIsotope(b11,3.2218e-5)
mod_matl.addIsotope(h1,4.9457e-2)
mod_matl.addIsotope(o16m,2.4672e-2)

# Clad
clad_matl = IsoMaterial(name='Clad')
clad_matl.addIsotope(zr90,2.1827e-2)
clad_matl.addIsotope(zr91,4.7600e-3)
clad_matl.addIsotope(zr92,7.2758e-3)
clad_matl.addIsotope(zr94,7.3734e-3)
clad_matl.addIsotope(zr96,1.1879e-3)



###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

p = 1.25984
pin_circle = Circle(x=0.0, y=0.0, radius=0.39218, name='pin')
clad_circle = Circle(x=0.0, y=0.0, radius=0.45720, name='clad')
left = XPlane(x=-p/2, name='left')
right = XPlane(x=p/2, name='right')
top = YPlane(y=p/2, name='top')
bottom = YPlane(y=-p/2, name='bottom')

left.setBoundaryType(REFLECTIVE)
right.setBoundaryType(REFLECTIVE)
top.setBoundaryType(REFLECTIVE)
bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

fuel = CellBasic(name='fuel')
fuel.setMaterial(fuel_matl)
fuel.addSurface(halfspace=-1, surface=pin_circle)

clad = CellBasic(name='clad')
clad.setMaterial(clad_matl)
clad.addSurface(halfspace=+1, surface=pin_circle)
clad.addSurface(halfspace=-1, surface=clad_circle)

moderator = CellBasic(name='moderator')
moderator.setMaterial(mod_matl)
moderator.addSurface(halfspace=+1, surface=clad_circle)
moderator.addSurface(halfspace=+1, surface=left)
moderator.addSurface(halfspace=-1, surface=right)
moderator.addSurface(halfspace=+1, surface=bottom)
moderator.addSurface(halfspace=-1, surface=top)


###############################################################################
###########################   Creating Universes   ############################
###############################################################################

log.py_printf('NORMAL', 'Creating universes...')

root_universe = Universe(name='root universe')
root_universe.addCell(fuel)
root_universe.addCell(clad)
root_universe.addCell(moderator)


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

solver = CPUSolver(geometry, track_generator)
solver.setNumThreads(num_threads)
solver.setSourceConvergenceThreshold(tolerance)
solver.convergeSource(max_iters)
solver.printTimerReport()


###############################################################################
############################   Generating Plots   #############################
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

#plotter.plot_tracks(track_generator)
#plotter.plot_segments(track_generator)
#plotter.plot_materials(geometry, gridsize=500)
#plotter.plot_cells(geometry, gridsize=500)
#plotter.plot_flat_source_regions(geometry, gridsize=500)
#plotter.plot_fluxes(geometry, solver, energy_groups=[1,2,3,4,5,6,7])

log.py_printf('TITLE', 'Finished')
