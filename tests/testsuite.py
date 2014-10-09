import testlog
import testmaterialize
import testmaterializeh5
import testoptions
import unittest


## options for if we want to run all
fullsuite = False

def buildQuickTestSuite():

    quickTestSuite = unittest.TestLoader().loadTestsFromTestCase(testlog.TestLogLevel)
    quickTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testoptions.TestDefaultInit))
    quickTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testoptions.TestCustomInit))
    quickTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testmaterialize.TestPyFiles))
##    quickTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testmaterialize.TestMatPySuite))
##    quickTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testmaterializh5e.TestMatH5Suite))

    ## above 2 lines still throwing the Class error -- fix
    
    return quickTestSuite

def buildFullTestSuite():

    fullTestSuite = buildQuickTestSuite()
    fullTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase("TESTNAMEGOESHERE"))

    ## TO ADD: REGRESSION TESTS (slower)
    ## First - rework regression tests

if fullsuite:
    TestSuite = buildFullTestSuite()

else:
    TestSuite = buildQuickTestSuite()

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(TestSuite)

