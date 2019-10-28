import unittest
import numpy

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc


class TestLogParameters(unittest.TestCase):

    def test_log_level(self):

        # Test setting log levels with the openmoc logs
        openmoc.set_log_level(openmoc.DEBUG)
        self.assertEqual(openmoc.get_log_level(), 0)

        openmoc.set_log_level(openmoc.INFO)
        self.assertEqual(openmoc.get_log_level(), 1)

        openmoc.set_log_level(openmoc.NORMAL)
        self.assertEqual(openmoc.get_log_level(), 3)

        openmoc.set_log_level(openmoc.SEPARATOR)
        self.assertEqual(openmoc.get_log_level(), 5)

        openmoc.set_log_level(openmoc.HEADER)
        self.assertEqual(openmoc.get_log_level(), 6)

        openmoc.set_log_level(openmoc.TITLE)
        self.assertEqual(openmoc.get_log_level(), 7)

        openmoc.set_log_level(openmoc.WARNING)
        self.assertEqual(openmoc.get_log_level(), 8)

        openmoc.set_log_level(openmoc.CRITICAL)
        self.assertEqual(openmoc.get_log_level(), 10)

        openmoc.set_log_level(openmoc.RESULT)
        self.assertEqual(openmoc.get_log_level(), 11)

        openmoc.set_log_level(openmoc.UNITTEST)
        self.assertEqual(openmoc.get_log_level(), 12)

        openmoc.set_log_level(openmoc.ERROR)
        self.assertEqual(openmoc.get_log_level(), 13)

        openmoc.set_log_level("DEBUG")
        self.assertEqual(openmoc.get_log_level(), 0)

        openmoc.set_log_level("INFO")
        self.assertEqual(openmoc.get_log_level(), 1)

        openmoc.set_log_level("NORMAL")
        self.assertEqual(openmoc.get_log_level(), 3)

        openmoc.set_log_level("SEPARATOR")
        self.assertEqual(openmoc.get_log_level(), 5)

        openmoc.set_log_level("HEADER")
        self.assertEqual(openmoc.get_log_level(), 6)

        openmoc.set_log_level("TITLE")
        self.assertEqual(openmoc.get_log_level(), 7)

        openmoc.set_log_level("WARNING")
        self.assertEqual(openmoc.get_log_level(), 8)

        openmoc.set_log_level("CRITICAL")
        self.assertEqual(openmoc.get_log_level(), 10)

        openmoc.set_log_level("RESULT")
        self.assertEqual(openmoc.get_log_level(), 11)

        openmoc.set_log_level("UNITTEST")
        self.assertEqual(openmoc.get_log_level(), 12)

        openmoc.set_log_level("ERROR")
        self.assertEqual(openmoc.get_log_level(), 13)

        # Setting log level with an integer
        openmoc.set_log_level(3)
        self.assertEqual(openmoc.get_log_level(), 3)

    def test_log_characters(self):

        openmoc.set_separator_character('z')
        self.assertEqual(openmoc.get_separator_character(), 'z')

        openmoc.set_header_character('w')
        self.assertEqual(openmoc.get_header_character(), 'w')

        openmoc.set_title_character('u')
        self.assertEqual(openmoc.get_title_character(), 'u')

if __name__ == '__main__':
    unittest.main()
