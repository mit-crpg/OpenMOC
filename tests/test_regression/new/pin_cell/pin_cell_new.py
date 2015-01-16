from openmoc import *
import openmoc.log as log
import openmoc.plotter as plotter
import openmoc.materialize as materialize
from openmoc.options import Options


###############################################################################
#######################   Main Simulation Parameters   ########################
###############################################################################

class simulation_object():

    def set_up_options(self):

        self.options = Options()
        self._num_threads = self.options.getNumThreads()
        self._track_spacing = self.options.getTrackSpacing()
        self._num_azim = self.options.getNumAzimAngles()
        self._tolerance = self.options.getTolerance()
        self._max_iters = self.options.getMaxIterations()

        ## Change this to True if you want to plot when you run this simulation
        ## (will NOT plot if you run the test function)
        self._plot = False


###############################################################################
###########################   Creating Materials   ############################
###############################################################################

    def create_materials(self):
        
        self._materials = materialize.materialize('../c5g7-materials.h5')
        self._uo2_id = self._materials['UO2'].getId()
        self._water_id = self._materials['Water'].getId()

###############################################################################
###########################   Creating Surfaces   #############################
###############################################################################

    def create_surfaces(self):

        self._circle = Circle(x=0.0, y=0.0, radius=1.0)
        self._left = XPlane(x=-2.0)
        self._right = XPlane(x=2.0)
        self._top = YPlane(y=2.0)
        self._bottom = YPlane(y=-2.0)
        
        self._left.setBoundaryType(REFLECTIVE)
        self._right.setBoundaryType(REFLECTIVE)
        self._top.setBoundaryType(REFLECTIVE)
        self._bottom.setBoundaryType(REFLECTIVE)


###############################################################################
#############################   Creating Cells   ##############################
###############################################################################


    def create_cells(self):
        
        cells = []
        cells.append(CellBasic(universe=1, material=self._uo2_id))
        cells.append(CellBasic(universe=1, material=self._water_id))
        cells.append(CellFill(universe=0, universe_fill=2))

        cells[0].addSurface(halfspace=-1, surface=self._circle)
        cells[1].addSurface(halfspace=+1, surface=self._circle)
        cells[2].addSurface(halfspace=+1, surface=self._left)
        cells[2].addSurface(halfspace=-1, surface=self._right)
        cells[2].addSurface(halfspace=+1, surface=self._bottom)

        cells[2].addSurface(halfspace=-1, surface=self._top)

        self._cells = cells


###############################################################################
###########################   Creating Lattices   #############################
###############################################################################

    def create_lattices(self):

        lattice = Lattice(id=2, width_x=4.0, width_y=4.0)

        lattice.setLatticeCells([[1]])

        self._lattice = lattice


###############################################################################
##########################   Creating the Geometry   ##########################
###############################################################################

    def create_geometry(self):

        geometry = Geometry()
        for material in self._materials.values(): geometry.addMaterial(material)
        for cell in self._cells: geometry.addCell(cell)
        geometry.addLattice(self._lattice)

        geometry.initializeFlatSourceRegions()
        self._geometry = geometry


###############################################################################
########################   Creating the TrackGenerator   ######################
###############################################################################

    def create_track_generator(self):
        track_generator = TrackGenerator(self._geometry, self._num_azim, self._track_spacing)
        track_generator.generateTracks()
        self._track_generator = track_generator


###############################################################################
###########################   Running a Simulation   ##########################
###############################################################################

    def run_simulation_full(self):
        solver = ThreadPrivateSolver(self._geometry, self._track_generator)
        solver.setNumThreads(self._num_threads)
        solver.setSourceConvergenceThreshold(self._tolerance)
        solver.convergeSource(self._max_iters)
        solver.printTimerReport()
        self._solver = solver

    def run_simulation_test(self):

        ## NEW ERROR -- ThreadPrivateSolver is not defined ##

        
        solver = ThreadPrivateSolver(self._geometry, self._track_generator)
        solver.setNumThreads(self._num_threads)
        solver.setSourceConvergenceThreshold(self._tolerance)
        solver.convergeSource(self._max_iters)
        Keff = solver.getKeff()
        
        self._solver = solver
        self._Keff = Keff
        return Keff


###############################################################################
############################   Generating Plots   #############################
###############################################################################

    def plot_data(self):
        
        log.py_printf('NORMAL', 'Plotting data...')

        plotter.plot_tracks(self._track_generator)
        plotter.plot_segments(self._track_generator)
        plotter.plot_materials(self._geometry, gridsize=500)
        plotter.plot_cells(self._geometry, gridsize=500)
        plotter.plot_flat_source_regions(self._geometry, gridsize=500)
        plotter.plot_fluxes(self._geometry, self._solver, energy_groups=[1,2,3,4,5,6,7])

###############################################################################
############################   Full Simulations   #############################
###############################################################################

    def run_full_simulation(self):

        ## run entire simulation, w/ print statements, and plotting if desired

        ## set up options
        self.setUpOptions()
        log.set_log_level('NORMAL')
        
        ## create materials
        log.py_printf('NORMAL', 'Importing materials data from HDF5...')
        self.createMaterials()

        ## create surfaces
        log.py_printf('NORMAL', 'Creating surfaces...')
        self.createSurfaces()

        ## create cells
        log.py_printf('NORMAL', 'Creating cells...')
        self.createCells()

        ## create lattices
        log.py_printf('NORMAL', 'Creating simple pin cell lattice...')
        self.createLattices()

        ## create geometry
        log.py_printf('NORMAL', 'Creating geometry...')
        self.createGeometry()

        ## create track generator
        log.py_printf('NORMAL', 'Initializing the track generator...')
        self.createTrackGenerator()

        ## run simulation
        self.runSimulationFull()

        ## plot?
        if self._plot:
            self.plotData()

        log.py_printf('TITLE', 'Finished')

    def run_test_simulation(self):

        ## run entire simulation, w/ print statements

        ## set up options
        self.setUpOptions()
        log.set_log_level('ERROR')
        
        ## create materials
        self.create_materials()

        ## create surfaces
        self.create_surfaces()

        ## create cells
        self.create_cells()

        ## create lattices
        self.create_lattices()

        ## create geometry
        self.create_geometry()

        ## create track generator
        self.create_track_generator()

        ## run simulation
        Keff = self.run_simulation_test()
        return Keff
        

if __name__ == '__main__':

    ## run simulation w/ print statements
    
    simulation = simulation_object()
    simulation.run_full_simulation()
