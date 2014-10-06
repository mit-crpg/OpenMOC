## Run all materialize tests.

import testmaterialize
import testmaterializeh5
import unittest

## collect tests from testmaterialize.py

MaterializeSuite = unittest.TestLoader().loadTestsFromTestCase(MaterializePySuite)

## collect tests from testmaterializeh5.py

MaterializeSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(MaterializeH5Suite))

if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(MaterializeSuite)



## Tests run on their own right now
## TODO: aggregate tests into one suite
