##def generic_test_setup(sysargs, filename, Keff, materials_filename=False, log_level='ERROR'):
##
##    ###############################################################################
##    #######################   Main Simulation Parameters   ########################
##    ###############################################################################
##
##    sys.argv = sysargs
##
##    options = Options()
##
##    num_threads = options.getNumThreads()
##    track_spacing = options.getTrackSpacing()
##    num_azim = options.getNumAzimAngles()
##    tolerance = options.getTolerance()
##    max_iters = options.getMaxIterations()
##
##    log.set_log_level(log_level)
##
##    ###############################################################################
##    ###########################   Creating Materials   ############################
##    ###############################################################################
##
##    log.py_printf('NORMAL', 'Creating materials...')
##
##    if materials_filename:
##        materials = materialize.materialize(materials_filename)
##
##    for material in materials:
##
##        if material == 'UO2':
##            uo2_id = materials['UO2'].getId()
##        elif material == 'MOX-4.3%':
##            mox43_id = materials['MOX-4.3%'].getId()
##        elif material == 'MOX-7%':
##            mox7_id = materials['MOX-7%'].getId()
##        elif material == 'MOX-8.7%':
##            mox87_id = materials['MOX-8.7%'].getId()
##        elif material == 'Guide Tube':
##            mox43_id = materials['Guide Tube'].getId()
##        elif material == 'Fission Chamber':
##            mox43_id = materials['Fission Chamber'].getId()
##        elif material == 'Water':
##            mox43_id = materials['Water'].getId()
##
##    ###############################################################################
##    ###########################   Creating Surfaces   #############################
##    ###############################################################################
##
##    log.py_printf('NORMAL', 'Creating surfaces...')
##
##    circle = Circle(x=0.0, y=0.0, radius=0.4)
##    left = XPlane(x=-0.635)
##    right = XPlane(x=0.635)
##    top = YPlane(y=0.635)
##    bottom = YPlane(y=-0.635)
##
##    left.setBoundaryType(REFLECTIVE)
##    right.setBoundaryType(REFLECTIVE)
##    top.setBoundaryType(REFLECTIVE)
##    bottom.setBoundaryType(REFLECTIVE)
##
##
##    ###############################################################################
##    #############################   Creating Cells   ##############################
##    ###############################################################################
##
##    log.py_printf('NORMAL', 'Creating cells...')
##
##    cells = []
##    cells.append(CellBasic(universe=1, material=1))
##    cells.append(CellBasic(universe=1, material=2))
##    cells.append(CellFill(universe=0, universe_fill=2))
##
##    cells[0].addSurface(halfspace=-1, surface=circle)
##    cells[1].addSurface(halfspace=+1, surface=circle)
##    cells[2].addSurface(halfspace=+1, surface=left)
##    cells[2].addSurface(halfspace=-1, surface=right)
##    cells[2].addSurface(halfspace=+1, surface=bottom)
##    cells[2].addSurface(halfspace=-1, surface=top)
##
##
##    ###############################################################################
##    ###########################   Creating Lattices   #############################
##    ###############################################################################
##
##    log.py_printf('NORMAL', 'Creating simple pin cell lattice...')
##
##    lattice = Lattice(id=2, width_x=1.27, width_y=1.27)
##    lattice.setLatticeCells([[1]])
##
##
##    ###############################################################################
##    ##########################   Creating the Geometry   ##########################
##    ###############################################################################
##
##    log.py_printf('NORMAL', 'Creating geometry...')
##
##    geometry = Geometry()
##    geometry.addMaterial(fuel)
##    geometry.addMaterial(moderator)
##    geometry.addCell(cells[0])
##    geometry.addCell(cells[1])
##    geometry.addCell(cells[2])
##    geometry.addLattice(lattice)
##
##    geometry.initializeFlatSourceRegions()
##
##
##    ###############################################################################
##    ########################   Creating the TrackGenerator   ######################
##    ###############################################################################
##
##    log.py_printf('NORMAL', 'Initializing the track generator...')
##
##    track_generator = TrackGenerator(geometry, num_azim, track_spacing)
##    track_generator.generateTracks()
##
##
##    ###############################################################################
##    ###########################   Running a Simulation   ##########################
##    ###############################################################################
##
##    solver = ThreadPrivateSolver(geometry, track_generator)
##    solver.setNumThreads(num_threads)
##    solver.setSourceConvergenceThreshold(tolerance)
##    solver.convergeSource(max_iters)
##    solver.printTimerReport()
##
##    log.py_printf('TITLE', 'Finished')
##
##    Keff = solver.getKeff()
##    TotalTime =  solver.getTotalTime()
##    NumPolarAngles = solver.getNumPolarAngles()
##    PolarQuadratureType = solver.getPolarQuadratureType()
##    NumFSRs = g.getNumFSRs()
##    SourceConvergenceThreshold = solver.getSourceConvergenceThreshold()
##    NumGroups = g.getNumEnergyGroups()
##    NumTracks = track_generator.getNumTracks()
##    NumSegments = track_generator.getNumSegments()
##
##    return [Keff, TotalTime, NumPolarAngles, PolarQuadratureType, NumFSRs,
##            SourceConvergenceThreshold, NumGroups, NumTracks, NumSegments]
##    

import getopt, sys
import unittest

class TestBenchmarkFramework(unittest.TestCase):

    ## parse command-line arguments for filename, Keff
    
    opts, args = getopt.getopt(sys.argv[1:], 'f:k:',['filename=', 'Keff='])
    except getopt.GetoptError as err:
        py_printf('WARNING', str(err))
        pass

    for opt, arg in opts:

        if opt in ('-f', '--filename'):
            filename = arg

        elif opt in ('-k', '--Keff'):
            Keff_actual = arg

    ## could add other quantities here - things that we might know ahead of time + try to compare
    ## to results. not sure.

    import filename


    ## TODO:

        ## need to make materials

        ## create surfaces

        ## create cells

        ## create lattices

        ## create geometry

    ## create trackgenerator

    track_generator = TrackGenerator(geometry, num_azim, track_spacing)
    track_generator.generateTracks()

    ## run test
        
    solver = ThreadPrivateSolver(geometry, track_generator)
    solver.setNumThreads(num_threads)
    solver.setSourceConvergenceThreshold(tolerance)
    solver.convergeSource(max_iters)

    ## return results

        Keff = solver.getKeff()
    TotalTime =  solver.getTotalTime()
    NumPolarAngles = solver.getNumPolarAngles()
    PolarQuadratureType = solver.getPolarQuadratureType()
    NumFSRs = g.getNumFSRs()
    SourceConvergenceThreshold = solver.getSourceConvergenceThreshold()
    NumGroups = g.getNumEnergyGroups()
    NumTracks = track_generator.getNumTracks()
    NumSegments = track_generator.getNumSegments()

    self.assertEqual(Keff_actual, Keff)
##    return [Keff, TotalTime, NumPolarAngles, PolarQuadratureType, NumFSRs,
##            SourceConvergenceThreshold, NumGroups, NumTracks, NumSegments]

