import unittest
import openmoc
import log

class TestLogLevel(unittest.TestCase):

    ## Tests the setters / getters

    ## List of levels and their associated value:
    ## 0: DEBUG
    ## 1: INFO
    ## 2: NORMAL
    ## 3: SEPARATOR (currently sets to HEADER via log.cpp)
    ## 4: HEADER
    ## 5: TITLE
    ## 6: WARNING
    ## 7: CRITICAL
    ## 8: RESULT
    ## 9: UNITTEST
    ## 10: ERROR

    def testDebugLevel(self):
        openmoc.set_log_level('DEBUG')
        self.assertEqual(0, openmoc.get_log_level())

    def testInfoLevel(self):
        openmoc.set_log_level('INFO')
        self.assertEqual(1, openmoc.get_log_level())

    def testNormalLevel(self):
        openmoc.set_log_level('NORMAL')
        self.assertEqual(2, openmoc.get_log_level())

    def testSeparatorLevel(self):
        openmoc.set_log_level('SEPARATOR')
        self.assertEqual(3, openmoc.get_log_level())

    def testHeaderLevel(self):
        # CURRENTLY FAILS
        openmoc.set_log_level('HEADER')
        self.assertEqual(4, openmoc.get_log_level())

    def testTitleLevel(self):
        openmoc.set_log_level('TITLE')
        self.assertEqual(5, openmoc.get_log_level())

    def testWarningLevel(self):
        openmoc.set_log_level('WARNING')
        self.assertEqual(6, openmoc.get_log_level())

    def testCriticalLevel(self):
        openmoc.set_log_level('CRITICAL')
        self.assertEqual(7, openmoc.get_log_level())

    def testResultLevel(self):
        openmoc.set_log_level('RESULT')
        self.assertEqual(8, openmoc.get_log_level())

    def testUnittestLevel(self):
        openmoc.set_log_level('UNITTEST')
        self.assertEqual(9, openmoc.get_log_level())

    def testErrorLevel(self):
        openmoc.set_log_level('ERROR')
        self.assertEqual(10, openmoc.get_log_level())




suite = unittest.TestLoader().loadTestsFromTestCase(TestLogLevel)

unittest.TextTestRunner(verbosity=2).run(suite)
