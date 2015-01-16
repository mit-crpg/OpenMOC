## import statements don't work bc the tests are in individual folders now -- fix

import test_logger.test_log as test_log
import test_materialize.test_materialize as test_materialize
import test_options.test_options as test_options
import test_regression.test_regression_suite as test_regression
import unittest


## options for if we want to run all
fullsuite = True

def build_quick_test_suite():

    quick_test_suite = unittest.TestLoader().loadTestsFromTestCase(test_log.test_log_level)
    quick_test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(test_options.test_default_init))
    quick_test_suite.addTest(unittest.TestLoader().loadTestsFromTestCase(test_options.test_custom_init))

## check if adding test_options.test_options_suite gets same result as adding 2 tests!

    quick_test_suite.addTest(test_materialize.materialize_suite)
    
    return quick_test_suite

def build_full_test_suite():

    full_test_suite = build_quick_test_suite()
    full_test_suite.addTest(test_regression.regression_suite)

    return full_test_suite

if fullsuite:
    test_suite = build_full_test_suite()

else:
    test_suite = build_quick_test_suite()

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(test_suite)

