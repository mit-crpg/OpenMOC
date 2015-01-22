import time
import glob
import string
import importlib
import imp
import sys

from regression_test_runner import *

start = time.clock()
# open file for failed test results and information
output = open("regression_suite_failed_results.txt", "a+")

# add all relevant folders to the path for importing
input_folders = glob.glob('benchmarks/*')[2:]
for folder in input_folders:
    sys.path.insert(0, folder)

# make list of all test files -- then another list with just the file names
# without benchmarks/*/ and .py
test_files_with_path = glob.glob('benchmarks/*/test_*.py')
test_files = []
for test_file in test_files_with_path:

    # removes benchmarks/ and .py from string
    test_file_module_with_directory = test_file[11:-3]
    
    # file name only (removes its directory)
    test_file_module_name = test_file_module_with_directory[string.find(test_file_module_with_directory,'/')+1:]

    # import the module
    test_files.append(test_file_module_name)
    importlib.import_module(test_file_module_name)

# dictionary mapping module names (strings) to modules
modules_dict = sys.modules
        
def build_regression_suite():

    regression_suite = regression_test_suite([], output)
    
    for test_file in test_files:

        test_module = modules_dict[test_file]
        regression_suite.add_test(test_module.test_list)
    
    print 'added all tests to suite'
    return regression_suite

regression_suite = build_regression_suite()

if __name__ == '__main__':
    regression_suite.run_tests()
    print 'Ran all tests in', time.clock() - start, 'seconds'
