import openmoc
import unittest

import benchmarks.homogeneousOneGroup.testH1G as H1G
import benchmarks.homogeneousTwoGroup as H2G
import benchmarks.LRA as LRA
import benchmarks.romano as romano
import benchmarks.testC5g7 as testC5g7

import new.pinCell as pinCellNew


def buildRegressionSuite():

    regressionSuite = unittest.TestSuite()
    regressionSuite.addTest(H1G.testH1G)

    return regressionSuite

RegressionSuite = buildRegressionSuite()


if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(RegressionSuite)
