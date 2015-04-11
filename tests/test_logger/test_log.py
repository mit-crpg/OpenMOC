import unittest
import openmoc
import openmoc.log as log

class test_log_level(unittest.TestCase):

    ## Tests the setters / getters

    ## List of levels and their associated value:
    ## 0: DEBUG
    ## 1: INFO
    ## 2: NORMAL
    ## 3: SEPARATOR
    ## 4: HEADER
    ## 5: TITLE
    ## 6: WARNING
    ## 7: CRITICAL
    ## 8: RESULT
    ## 9: UNITTEST
    ## 10: ERROR

    def test_debug_level(self):
        openmoc.set_log_level('DEBUG')
        self.assertEqual(0, openmoc.get_log_level())

    def test_info_level(self):
        openmoc.set_log_level('INFO')
        self.assertEqual(1, openmoc.get_log_level())

    def test_normal_level(self):
        openmoc.set_log_level('NORMAL')
        self.assertEqual(2, openmoc.get_log_level())

    def test_separator_level(self):
        openmoc.set_log_level('SEPARATOR')
        self.assertEqual(3, openmoc.get_log_level())

    def test_header_level(self):
        openmoc.set_log_level('HEADER')
        self.assertEqual(4, openmoc.get_log_level())

    def test_title_level(self):
        openmoc.set_log_level('TITLE')
        self.assertEqual(5, openmoc.get_log_level())

    def test_warning_level(self):
        openmoc.set_log_level('WARNING')
        self.assertEqual(6, openmoc.get_log_level())

    def test_critical_level(self):
        openmoc.set_log_level('CRITICAL')
        self.assertEqual(7, openmoc.get_log_level())

    def test_result_level(self):
        openmoc.set_log_level('RESULT')
        self.assertEqual(8, openmoc.get_log_level())

    def test_unittest_level(self):
        openmoc.set_log_level('UNITTEST')
        self.assertEqual(9, openmoc.get_log_level())

    def test_error_level(self):
        openmoc.set_log_level('ERROR')
        self.assertEqual(10, openmoc.get_log_level())




test_log_suite = unittest.TestLoader().loadTestsFromTestCase(test_log_level)

if __name__ == '__main__':

    unittest.TextTestRunner(verbosity=2).run(test_log_suite)
