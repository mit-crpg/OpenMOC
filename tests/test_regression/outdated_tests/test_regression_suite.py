import openmoc
import unittest
import time

## output = open("regression_failed_results.txt", "a+")
## create output file that the test runner will append all failed/error results to
## but first -- create test runner

## consider using glob module -- automatically add tests from benchmarks directory

## currently: importing regression test suites as modules

import benchmarks.homogeneous_one_group.test_homogeneous_one_group as H1G
import benchmarks.homogeneous_two_group.test_homogeneous_two_groups as H2G
import benchmarks.LRA.test_LRA as LRA
import benchmarks.romano.test_romano as romano
import benchmarks.c5g7_cmfd.test_c5g7_cmfd as c5g7_cmfd
print 'managed to import all'

## import new.pin_cell.pin_cell_new as pin_cell_new


        
def build_regression_suite():

    regression_suite = unittest.TestSuite()
#    regression_suite.addTest(H1G.test_h1g)
#    regression_suite.addTest(H2G.test_h2g)
#    regression_suite.addTest(LRA.test_LRA)
#    regression_suite.addTest(romano.test_romano)
    regression_suite.addTest(c5g7_cmfd.test_c5g7_cmfd)

    print 'added all tests to suite'

    return regression_suite


regression_suite = build_regression_suite()

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(regression_suite)
