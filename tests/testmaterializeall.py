## Run all materialize tests.

import testmaterialize as pytest
import testmaterializeh5 as h5test
import unittest

## collect tests from testmaterialize.py

print '************************************************************'
print '************************************************************'
print '************************************************************'
print type(pytest.PyTestSuite)
print '************************************************************'
print '************************************************************'
print '************************************************************'
MaterializeSuite = (unittest.TestLoader().loadTestsFromTestCase(h5test.MaterializeH5Suite))

MaterializeSuite.addTests(TestLoader().loadTestsFromTestCase(pytestsuite))

## collect tests from testmaterializeh5.py

if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(MaterializeSuite)


## NOT WORKING: says it needs a class as input to loadTestsFromTestCase() 
