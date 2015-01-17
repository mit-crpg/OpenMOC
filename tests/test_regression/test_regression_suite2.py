import openmoc
import unittest
import time

from regression_test_runner import *

output = open("regression_suite_failed_results.txt", "a+")

## consider using glob module -- automatically add tests from benchmarks directory

## IMPORTS ARE SLOW BECAUSE they actually do the calculations when it's imported now.
## So each import takes 2x as long as it takes to just run the code. Very slow.
## But the tests themselves are instant.

print 'starting imports'
start = time.clock()

import benchmarks.homogeneous_one_group.test_homogeneous_one_group_new as H1G
import benchmarks.homogeneous_two_group.test_homogeneous_two_groups_new as H2G
import benchmarks.LRA.test_LRA_new as LRA
import benchmarks.romano.test_romano_new as romano
import benchmarks.c5g7_cmfd.test_c5g7_cmfd_new as c5g7_cmfd
import benchmarks.c5g7.test_c5g7_new as c5g7

elapsed_tot = time.clock() - start
print 'Time so far: ',elapsed_tot
        
def build_regression_suite():

    regression_suite = regression_test_suite([], output)
    regression_suite.add_test(c5g7_cmfd.test_c5g7_cmfd)
    regression_suite.add_test(c5g7.test_c5g7)    
    regression_suite.add_test(LRA.test_LRA)
    regression_suite.add_test(H1G.test_H1G)
    regression_suite.add_test(H2G.test_H2G)
    regression_suite.add_test(romano.test_romano)
    print 'added all tests to suite'

    return regression_suite

regression_suite = build_regression_suite()

if __name__ == '__main__':
    regression_suite.run_tests()
    print 'total time:', time.clock() - start
