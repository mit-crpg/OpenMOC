## THIS WORKS AS OF 1/15/15

import unittest


import test_logger.test_log as test_log
import test_materialize.test_materialize as test_materialize
import test_options.test_options as test_options

from test_regression.regression_test_runner import *
from test_regression.test_regression_suite import *


## options for if we want to run all regression tests
fullsuite = True

def build_quick_test_suite():

    quick_test_suite = unittest.TestLoader().loadTestsFromTestCase(test_log.test_log_level)
    quick_test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(test_options.test_default_init))
    quick_test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(test_options.test_custom_init))

## check if adding test_options.test_options_suite gets same result as adding 2 tests!

    quick_test_suite.addTest(test_materialize.materialize_suite)
    
    return quick_test_suite

test_suite = build_quick_test_suite()

if __name__ == '__main__' and fullsuite:
    unittest.TextTestRunner(verbosity=2).run(test_suite)
    regression_suite.run_tests()

elif __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(test_suite)

