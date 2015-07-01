from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#                          Main Simulation Parameters
###############################################################################

options = Options()

num_threads = options.getNumThreads()
track_spacing = options.getTrackSpacing()
num_azim = options.getNumAzimAngles()
tolerance = options.getTolerance()
max_iters = options.getMaxIterations()
num_modes = 6

log.set_log_level('NORMAL')

log.py_printf('TITLE', 'Computing %d eigenmodes', num_modes)


###############################################################################
#                            Creating Materials
###############################################################################

log.py_printf('NORMAL', 'Importing materials data from HDF5...')

materials = materialize.materialize('../c5g7-materials.h5')


###############################################################################
#                            Creating Surfaces
###############################################################################

log.py_printf('NORMAL', 'Creating surfaces...')

left = XPlane(x=-64.26, name='left')
right = XPlane(x=64.26, name='right')
top = YPlane(y=64.26, name='top')
bottom = YPlane(y=-64.26, name='bottom')
left.setBoundaryType(VACUUM)
right.setBoundaryType(VACUUM)
top.setBoundaryType(VACUUM)
bottom.setBoundaryType(VACUUM)
boundaries = [left, right, top, bottom]

# Create Circles for the fuel as well as to discretize the moderator into rings
fuel_radius = Circle(x=0.0, y=0.0, radius=0.54)
moderator_inner_radius = Circle(x=0.0, y=0.0, radius=0.62)
moderator_outer_radius = Circle(x=0.0, y=0.0, radius=0.58)


###############################################################################
#                        Creating Cells and Universes
###############################################################################

log.py_printf('NORMAL', 'Creating cells...')

# Moderator rings
moderator_ring1 = Cell()
moderator_ring2 = Cell()
moderator_ring3 = Cell()
moderator_ring1.setNumSectors(8)
moderator_ring2.setNumSectors(8)
moderator_ring3.setNumSectors(8)
moderator_ring1.setFill(materials['Water'])
moderator_ring2.setFill(materials['Water'])
moderator_ring3.setFill(materials['Water'])
moderator_ring1.addSurface(+1, fuel_radius)
moderator_ring1.addSurface(-1, moderator_inner_radius)
moderator_ring2.addSurface(+1, moderator_inner_radius)
moderator_ring2.addSurface(-1, moderator_outer_radius)
moderator_ring3.addSurface(+1, moderator_outer_radius)

# UO2 pin cell
uo2_cell = Cell()
uo2_cell.setNumRings(3)
uo2_cell.setNumSectors(8)
uo2_cell.setFill(materials['UO2'])
uo2_cell.addSurface(-1, fuel_radius)

uo2 = Universe(name='UO2')
uo2.addCell(uo2_cell)
uo2.addCell(moderator_ring1)
uo2.addCell(moderator_ring2)
uo2.addCell(moderator_ring3)

# 4.3% MOX pin cell
mox43_cell = Cell()
mox43_cell.setNumRings(3)
mox43_cell.setNumSectors(8)
mox43_cell.setFill(materials['MOX-4.3%'])
mox43_cell.addSurface(-1, fuel_radius)

mox43 = Universe(name='MOX-4.3%')
mox43.addCell(mox43_cell)
mox43.addCell(moderator_ring1)
mox43.addCell(moderator_ring2)
mox43.addCell(moderator_ring3)

# 7% MOX pin cell
mox7_cell = Cell()
mox7_cell.setNumRings(3)
mox7_cell.setNumSectors(8)
mox7_cell.setFill(materials['MOX-7%'])
mox7_cell.addSurface(-1, fuel_radius)

mox7 = Universe(name='MOX-7%')
mox7.addCell(mox7_cell)
mox7.addCell(moderator_ring1)
mox7.addCell(moderator_ring2)
mox7.addCell(moderator_ring3)

# 8.7% MOX pin cell
mox87_cell = Cell()
mox87_cell.setNumRings(3)
mox87_cell.setNumSectors(8)
mox87_cell.setFill(materials['MOX-8.7%'])
mox87_cell.addSurface(-1, fuel_radius)

mox87 = Universe(name='MOX-8.7%')
mox87.addCell(mox87_cell)
mox87.addCell(moderator_ring1)
mox87.addCell(moderator_ring2)
mox87.addCell(moderator_ring3)

# Fission chamber pin cell
fission_chamber_cell = Cell()
fission_chamber_cell.setNumRings(3)
fission_chamber_cell.setNumSectors(8)
fission_chamber_cell.setFill(materials['Fission Chamber'])
fission_chamber_cell.addSurface(-1, fuel_radius)

fission_chamber = Universe(name='Fission Chamber')
fission_chamber.addCell(fission_chamber_cell)
fission_chamber.addCell(moderator_ring1)
fission_chamber.addCell(moderator_ring2)
fission_chamber.addCell(moderator_ring3)

# Guide tube pin cell
guide_tube_cell = Cell()
guide_tube_cell.setNumRings(3)
guide_tube_cell.setNumSectors(8)
guide_tube_cell.setFill(materials['Guide Tube'])
guide_tube_cell.addSurface(-1, fuel_radius)

guide_tube = Universe(name='Guide Tube')
guide_tube.addCell(guide_tube_cell)
guide_tube.addCell(moderator_ring1)
guide_tube.addCell(moderator_ring2)
guide_tube.addCell(moderator_ring3)

# Reflector
reflector_cell = Cell(name='moderator')
reflector_cell.setFill(materials['Water'])

reflector = Universe(name='Reflector')
reflector.addCell(reflector_cell)

# Cells
assembly1_cell = Cell(name='Assembly 1')
assembly2_cell = Cell(name='Assembly 2')
refined_reflector_cell = Cell(name='Semi-Finely Spaced Reflector')
right_reflector_cell = Cell(name='Right Reflector')
left_reflector_cell = Cell(name='Left Reflector')
bottom_right_corner_reflector_cell = Cell(name='Bottom Right Corner Reflector')
top_right_corner_reflector_cell = Cell(name='Top Right Corner Reflector')
bottom_left_corner_reflector_cell = Cell(name='Bottom Left Corner Reflector')
top_left_corner_reflector_cell = Cell(name='Top Left Corner Reflector')
bottom_reflector_cell = Cell(name='Bottom Reflector')
top_reflector_cell = Cell(name='Top Reflector')

assembly1 = Universe(name='Assembly 1')
assembly2 = Universe(name='Assembly 2')
refined_reflector = Universe(name='Semi-Finely Spaced Moderator')
right_reflector = Universe(name='Right Reflector')
left_reflector = Universe(name='Left Reflector')
bottom_right_corner_reflector = Universe(name='Bottom Right Corner Reflector')
top_right_corner_reflector = Universe(name='Top Right Corner Reflector')
bottom_left_corner_reflector = Universe(name='Bottom Left Corner Reflector')
top_left_corner_reflector = Universe(name='Top Left Corner Reflector')
bottom_reflector = Universe(name='Bottom Reflector')
top_reflector = Universe(name='Top Reflector')

assembly1.addCell(assembly1_cell)
assembly2.addCell(assembly2_cell)
refined_reflector.addCell(refined_reflector_cell)
right_reflector.addCell(right_reflector_cell)
left_reflector.addCell(left_reflector_cell)
bottom_right_corner_reflector.addCell(bottom_right_corner_reflector_cell)
top_right_corner_reflector.addCell(top_right_corner_reflector_cell)
bottom_left_corner_reflector.addCell(bottom_left_corner_reflector_cell)
top_left_corner_reflector.addCell(top_left_corner_reflector_cell)
bottom_reflector.addCell(bottom_reflector_cell)
top_reflector.addCell(top_reflector_cell)

# Root Cell/Universe
root_cell = Cell(name='Full Geometry')
root_cell.addSurface(+1, boundaries[0])
root_cell.addSurface(-1, boundaries[1])
root_cell.addSurface(-1, boundaries[2])
root_cell.addSurface(+1, boundaries[3])

root_universe = Universe(name='Root Universe')
root_universe.addCell(root_cell)


###############################################################################
#                             Creating Lattices
###############################################################################

log.py_printf('NORMAL', 'Creating lattices...')

lattices = list()

# Top left, bottom right 17 x 17 assemblies
lattices.append(Lattice(name='Assembly 1'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1],
            [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1],
            [1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

universes = {1 : uo2, 2 : guide_tube, 3 : fission_chamber}
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses(template)
assembly1_cell.setFill(lattices[-1])

# Top right, bottom left 17 x 17 assemblies
lattices.append(Lattice(name='Assembly 2'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
            [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
            [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
            [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 5, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
            [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
            [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
            [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
universes = {1 : mox43, 2 : mox7, 3 : mox87,
             4 : guide_tube, 5 : fission_chamber}
for i in range(17):
  for j in range(17):
    template[i][j] = universes[template[i][j]]
lattices[-1].setUniverses(template)
assembly2_cell.setFill(lattices[-1])

# Sliced up water cells - semi finely spaced
lattices.append(Lattice(name='Semi-Finely Spaced Reflector'))
lattices[-1].setWidth(width_x=0.126, width_y=0.126)
template = [[reflector] * 10] * 10
lattices[-1].setUniverses(template)
refined_reflector_cell.setFill(lattices[-1])

# Sliced up water cells - right side of geometry
lattices.append(Lattice(name='Right Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 17
lattices[-1].setUniverses(template)
right_reflector_cell.setFill(lattices[-1])

# Sliced up water cells - right side of geometry
lattices.append(Lattice(name='Left Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 6 +  [refined_reflector] * 11] * 17
lattices[-1].setUniverses(template)
left_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Bottom Right Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 11 + [reflector] * 6] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
bottom_right_corner_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Top Right Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 17] * 6
template += [[refined_reflector] * 11 + [reflector] * 6] * 11
lattices[-1].setUniverses(template)
top_right_corner_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Bottom Left Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 6 + [refined_reflector] * 11] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
bottom_left_corner_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom corner of geometry
lattices.append(Lattice(name='Top Left Corner Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 17] * 6
template += [[reflector] * 6 + [refined_reflector] * 11] * 11
lattices[-1].setUniverses(template)
top_left_corner_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for bottom of geometry
lattices.append(Lattice(name='Bottom Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[refined_reflector] * 17] * 11
template += [[reflector] * 17] * 6
lattices[-1].setUniverses(template)
bottom_reflector_cell.setFill(lattices[-1])

# Sliced up water cells for top of geometry
lattices.append(Lattice(name='Top Reflector'))
lattices[-1].setWidth(width_x=1.26, width_y=1.26)
template = [[reflector] * 17] * 6
template += [[refined_reflector] * 17] * 11
lattices[-1].setUniverses(template)
top_reflector_cell.setFill(lattices[-1])

# 4 x 4 core to represent two bundles and water
lattices.append(Lattice(name='Full Geometry'))
lattices[-1].setWidth(width_x=21.42, width_y=21.42)
lattices[-1].setUniverses([
     [top_left_corner_reflector   , top_reflector   , top_reflector   , top_reflector   , top_reflector   , top_right_corner_reflector    ],
     [left_reflector              , assembly1       , assembly2       , assembly2       , assembly1       , right_reflector               ],
     [left_reflector              , assembly2       , assembly1       , assembly1       , assembly2       , right_reflector               ],
     [left_reflector              , assembly2       , assembly1       , assembly1       , assembly2       , right_reflector               ],
     [left_reflector              , assembly1       , assembly2       , assembly2       , assembly1       , right_reflector               ],
     [bottom_left_corner_reflector, bottom_reflector, bottom_reflector, bottom_reflector, bottom_reflector, bottom_right_corner_reflector]])
root_cell.setFill(lattices[-1])


###############################################################################
#                         Creating the Geometry
###############################################################################

log.py_printf('NORMAL', 'Creating geometry...')

geometry = Geometry()
geometry.setRootUniverse(root_universe)
geometry.initializeFlatSourceRegions()


###############################################################################
#                          Creating the TrackGenerator
###############################################################################

log.py_printf('NORMAL', 'Initializing the track generator...')

track_generator = TrackGenerator(geometry, num_azim, track_spacing)
track_generator.setNumThreads(num_threads)
track_generator.generateTracks()


###############################################################################
#                            Running a Simulation
###############################################################################

from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import eigs
import numpy as np

solver = CPUSolver(track_generator)
solver.setNumThreads(num_threads)
solver.initializeMemory()

# Initialize operators, counters
a_count = 0
m_count = 0

def A_operator(flux):
    global solver
    global a_count
    import numpy
    import copy
    
    # Do not pass imaginary numbers to SWIG
    flux = numpy.real(flux).astype(numpy.float64)
    flux_old = copy.copy(flux)
    
    solver.scatterTransportSweep(flux)
    
    a_count += 1
    log.py_printf('NORMAL', "Performed A operator sweep number %d", a_count)
    
    return flux_old - flux

def M_operator(flux):
    global solver
    global m_count
    import numpy
    # Do not pass imaginary numbers to SWIG
    flux = numpy.real(flux).astype(numpy.float64)
    
    solver.fissionTransportSweep(flux)
    
    m_count += 1
    log.py_printf('NORMAL', "Performed M operator sweep number %d", m_count)
    
    return flux
    
A = LinearOperator( (solver.getOperatorSize(), solver.getOperatorSize() ), matvec=A_operator, dtype='float64')
M = LinearOperator( (solver.getOperatorSize(), solver.getOperatorSize() ), matvec=M_operator, dtype='float64')

def F_operator(flux):
    global A
    global M
    global tolerance
    from scipy.sparse.linalg import gmres
    
    flux = M*flux
    flux, x = gmres(A, flux, tol=tolerance/10) # Note, gmres must converge beyond tolerance to work.
    
    return flux
    
F = LinearOperator( (solver.getOperatorSize(), solver.getOperatorSize() ), matvec=F_operator, dtype='float64')

vals, vecs = eigs(F, k=num_modes, tol=tolerance)

log.py_printf('NORMAL', "The eigenvalues are: %s", str(vals))

def plot_eigenmodes(vecs):
    for i in range(vecs.shape[1]):
        # Convert it into a form that SWIG will be happy with
        vec = np.squeeze(np.ascontiguousarray(vecs[:,i]))
        vec = np.real(vec).astype(np.float64)
        
        if(i == 0):
            # Ensure the primary eigenvalue is positive
            vec = np.abs(vec)
        
        # Insert it into OpenMOC
        solver.putFlux(vec)
        # Switch folders
        plotter.subdirectory = "/plots_eig" + str(i+1) + "/"
        # Plot
        plotter.plot_spatial_fluxes(solver, energy_groups=[1,2,3,4,5,6,7])
        

###############################################################################
#                             Generating Plots
###############################################################################

log.py_printf('NORMAL', 'Plotting data...')

plotter.plot_materials(geometry, gridsize=500)
plotter.plot_cells(geometry, gridsize=500)
plotter.plot_flat_source_regions(geometry, gridsize=500)
plot_eigenmodes(vecs)

log.py_printf('TITLE', 'Finished')
