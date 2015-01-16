import openmoc
import unittest
import time

from regression_test_runner import *

output = open("regression_failed_results_NEW.txt", "a+")

## consider using glob module -- automatically add tests from benchmarks directory

## IMPORTS ARE SLOW BECAUSE they actually do the calculations when it's imported now.
## So each import takes 2x as long as it takes to just run the code. Very slow.
## But the tests themselves are instant.

print 'starting imports'
start = time.clock()
####import benchmarks.homogeneous_one_group.test_homogeneous_one_group as H1G
####import benchmarks.homogeneous_two_group.test_homogeneous_two_groups as H2G
import benchmarks.LRA.test_LRA_new as LRA
elapsed_LRA = time.clock() - start
print 'Time so far: ', elapsed_LRA
####import benchmarks.romano.test_romano as romano
import benchmarks.c5g7_cmfd.test_c5g7_cmfd_new as c5g7_cmfd
elapsed_tot = time.clock() - start
print 'Time so far: ',elapsed_tot

## import new.pin_cell.pin_cell_new as pin_cell_new
        
def build_regression_suite():

    regression_suite = regression_test_suite([], output)
    print 'made initial empty suite'
    regression_suite.add_test(c5g7_cmfd.test_c5g7_cmfd)
    regression_suite.add_test(LRA.test_LRA)

    print 'added all tests to suite'

    return regression_suite

print 'about to build suite'
regression_suite = build_regression_suite()

if __name__ == '__main__':
    regression_suite.run_tests()
