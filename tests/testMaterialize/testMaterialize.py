## Run all materialize tests.

import testMaterializePython as pytest
import testMaterializeH5 as h5test
import unittest

## collect tests from testmaterialize.py

MaterializeSuite = (h5test.TestMatH5Suite)

MaterializeSuite.addTests(pytest.TestMatPySuite)

## collect tests from testmaterializeh5.py

if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(MaterializeSuite)
