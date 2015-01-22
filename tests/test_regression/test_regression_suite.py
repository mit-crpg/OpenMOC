import time
import glob
import string
import importlib
import imp
import sys
import os

if __name__ == '__main__':
    run_from_regression_suite = True
else:
    run_from_regression_suite = False

from regression_test_runner import *

output = open("regression_suite_failed_results.txt", "a") # create output file in case of failures

start = time.clock()

# open file for failed test results and information
output = open("regression_suite_failed_results.txt", "a+")

##current_directory = os.path.dirname(os.path.realpath(__file__))
##sys.path.insert(0,current_directory)

# add all relevant folders to the path for importing
if run_from_regression_suite:
    input_folders = glob.glob('benchmarks/*')[2:]
    
if not run_from_regression_suite:
    input_folders = glob.glob('test_regression/benchmarks/*')[2:]

for folder in input_folders:
    sys.path.insert(0, folder)

# make list of all test files (and directories)
test_files = []

if run_from_regression_suite:
    test_files_with_path = glob.glob('benchmarks/*/test_*.py')

if not run_from_regression_suite:
    test_files_with_path = glob.glob('test_regression/benchmarks/*/test_*.py')

# string processing to remove directories, .py
for test_file in test_files_with_path:
    if run_from_regression_suite:
        # removes benchmarks/ and .py from string
        test_file_module_with_directory = test_file[11:-3]
    if not run_from_regression_suite:
        test_file_module_with_directory = test_file[27:-3]
    
    # file name only (removes its directory)
    test_file_module_name = test_file_module_with_directory[string.find(test_file_module_with_directory,'/')+1:]

    # import the module
    test_files.append(test_file_module_name)
    importlib.import_module(test_file_module_name)

# dictionary mapping module names (strings) to modules
modules_dict = sys.modules


def build_regression_suite():

    reg_suite = regression_test_suite([], output)
    
    for test_file in test_files:

        test_module = modules_dict[test_file]
        reg_suite.add_test(test_module.test_list)
    
    return reg_suite

regression_suite = build_regression_suite()

if __name__ == '__main__':
    regression_suite.run_tests()
