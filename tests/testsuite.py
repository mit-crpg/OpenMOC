## import statements don't work bc the tests are in individual folders now -- fix

import testLog.testLog as testlog
import testMaterialize.testMaterialize as testmaterialize
import testOptions.testOptions as testoptions
import testRegression
import unittest

## TODO: fix

## options for if we want to run all
fullsuite = False

def buildQuickTestSuite():

    quickTestSuite = unittest.TestLoader().loadTestsFromTestCase(testlog.TestLogLevel)
    quickTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testoptions.TestDefaultInit))
    quickTestSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testoptions.TestCustomInit))
    quickTestSuite.addTest(testmaterialize.MaterializeSuite)
    
    return quickTestSuite

def buildFullTestSuite():

    fullTestSuite = buildQuickTestSuite()
    fullTestSuite.addTest()

    ## TO ADD: REGRESSION TESTS (slower)
    ## First - rework regression tests

if fullsuite:
    TestSuite = buildFullTestSuite()

else:
    TestSuite = buildQuickTestSuite()

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(TestSuite)

