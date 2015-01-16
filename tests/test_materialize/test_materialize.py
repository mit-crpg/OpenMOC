## Run all materialize tests.

import test_materialize_python as pytest
import test_materialize_H5 as h5test
import unittest

## collect tests from testmaterialize.py

materialize_suite = (h5test.test_materialize_H5_suite)

materialize_suite.addTests(pytest.test_materials_py_suite)

## collect tests from testmaterializeh5.py

if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(materialize_suite)
