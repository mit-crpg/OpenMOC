import getopt, sys
import unittest

class test_benchmark_framework(unittest.TestCase):

    def test_Keff(self):
        ## parse command-line arguments for filename, Keff

        filename = False
        
        opts, args = getopt.getopt(sys.argv[1:], 'f:k:',['filename=', 'Keff='])
    ##    except getopt.GetoptError as err:
    ##        py_printf('WARNING', str(err))
    ##        pass

        for opt, arg in opts:

            if opt in ('-f', '--filename'):
                filename = arg

            elif opt in ('-k', '--Keff'):
                Keff_actual = float(arg)

        if filename:
            file_name = __import__(filename)
            ## CHECK THE SYNTAX/STRUCTURE OF THIS
            ## Does this import as desired?

        else:
            log.py_printf('ERROR', 'No filename given')

        sys.argv = [filename]
        ## add way for other options to be added here
        
        test_sim = file_name.simulation_object()
        Keff = test_sim.run_test_simulation()

    ##    NumPolarAngles = solver.getNumPolarAngles()
    ##    PolarQuadratureType = solver.getPolarQuadratureType()
    ##    NumFSRs = g.getNumFSRs()
    ##    SourceConvergenceThreshold = solver.getSourceConvergenceThreshold()
    ##    NumGroups = g.getNumEnergyGroups()
    ##    NumTracks = track_generator.getNumTracks()
    ##    NumSegments = track_generator.getNumSegments()

        self.assertTrue(abs(Keff - Keff_actual) < 0.000001)
        ## maybe change threshold here

        
    ##    return [Keff, TotalTime, NumPolarAngles, PolarQuadratureType, NumFSRs,
    ##            SourceConvergenceThreshold, NumGroups, NumTracks, NumSegments]

suite = unittest.TestLoader().loadTestsFromTestCase(test_benchmark_framework)
unittest.TextTestRunner(verbosity=2).run(suite)
