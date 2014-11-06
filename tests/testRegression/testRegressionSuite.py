import openmoc
import unittest

print 'starting import'
import benchmarks.homogeneousOneGroup.testH1G as H1G
import benchmarks.homogeneousTwoGroup.testH2G as H2G
import benchmarks.LRA.testLRA as LRA
import benchmarks.romano.testRomano as romano
import benchmarks.testc5g7.testC5g7 as C5G7

import new.pinCell.pincellnew as pinCellNew
print 'imports done'


def buildRegressionSuite():

    regressionSuite = unittest.TestSuite()
#    regressionSuite.addTest(H1G.testH1G)
#    regressionSuite.addTest(H2G.testH2G)
    regressionSuite.addTest(LRA.testLRA)
    regressionSuite.addTest(romano.testRomano)
    regressionSuite.addTest(C5G7.testC5G7)


    return regressionSuite

RegressionSuite = buildRegressionSuite()


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(RegressionSuite)
