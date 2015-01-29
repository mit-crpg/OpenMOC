import time
import multiprocessing as mp
import platform
import sys

print_times = True
# variable for if you want it to print time taken for each individual test
# this could be a command-line thing

class regression_test_case():

    """Single regression test case to be handled as part of a suite.

    Attributes:
      test_type = string; 'Keff', 'pin powers'
      benchmark = string; 'LRA', 'c5g7-cmfd'
      benchmark_value = float; known value
      error_margin = float; acceptable margin for correct answers
      filename = string; filename of the benchmark file to be acted on, i.e. 'LRA.py'
      setup_func = variable indicating the setup function (part of the test file),
          i.e. setup_LRA, setup_romano, etc.
      num_threads = string or int; the number of threads the test is for.
          Should be either 1 or 'DEFAULT'.
    """

    def __init__(self, test_type, benchmark, benchmark_value, error_margin, filename, num_threads):

        self.test_type = test_type
        self.benchmark = benchmark
        self.benchmark_value = benchmark_value
        self.error_margin = error_margin
        self.filename = filename        
 #       self.setup_func = setup_func
        self.num_threads = num_threads # this should be either 'DEFAULT' or 1

    def set_up_test(self, test_mod):

        """Calculates the Keff and Keff 1 thread of the test. This is set up
        so that the tests are set up and then run individually instead of waiting 30 min
        for all the simulations to run in the import statements.
        """

        filename = self.filename
        num_threads = self.num_threads

        module_dict = sys.modules
        test_module = module_dict[test_mod]

        if num_threads == 'DEFAULT':
            self.calculated_value = test_module.setup([filename])
        elif type(num_threads) == int:
            self.calculated_value = test_module.setup([filename, '--num-omp-threads', str(num_threads)])
        else:
            print 'ERROR: unexpected type for num_threads (expected int or str, got', type(num_threads)


class regression_test_suite():

    """A suite including a list of tests to run and a desired output file.

    Attributes:
      tests = list of tuples: test[0] is regression_test_case instance,
          test[1] is the name of the module whose setup routine will be used 
      output = variable referring to output filename
      none_failed = boolean referring to whether or not a test has failed thus
          far (defaults to True - meaning no test has failed thus far)

    """

    def __init__(self, tests, output, none_failed=True):
        self.tests = tests
        self.none_failed = none_failed
        self.output = output
        self.num_failed = 0

    def add_test(self, test):
        # handles both single tests and lists of test objects
        self.tests.extend(test)

    def get_tests(self):
        return self.tests

    def test_failed(self):
        self.none_failed = False
        self.num_failed += 1

    def get_num_tests(self):
        return len(self.tests)

    def get_num_failed(self):
        return self.num_failed

    def end_file(self):

        # closes output file visually
        
        output = self.output
        output.write('------------------TESTS COMPLETE------------------'+"\n")
        output.write('\n')

    def close_file(self):

        # closes output file
        
        output = self.output
        output.close()

    def run_tests(self):

        start_time = time.clock()
        print '------------------BEGINNING REGRESSION TESTS------------------'

        # runs all tests and adds to output file; closes file when done
        for test_tuple in self.get_tests():

            test = test_tuple[0]
            test_mod = test_tuple[1]
            test_start_time = time.clock()
            test.set_up_test(test_mod)
            if not run_regression_test(test, self.output, test_start_time, self.none_failed):
                self.test_failed()
            

        print '------------------TESTS COMPLETE------------------'
        duration = time.clock() - start_time
        print 'Completed '+str(self.get_num_tests())+' tests in '+str(duration)+' seconds.'
        print 'Tests passed:', self.get_num_tests() - self.get_num_failed()
        print 'Tests failed:', self.get_num_failed()
        
        if not self.none_failed:
            self.end_file()
            
        self.close_file()

def run_regression_test(test, output,test_start_time, none_failed=True):

    # test = a test object
    # output = a variable referring to the file that is open and being written to
    # none_failed: a boolean referring to whether or not any tests have failed yet.
        # True: no tests have failed so far
        # False: 1+ tests have failed so far

    test_type = test.test_type              # string indicating what type of test, i.e. Keff, pin powers
    benchmark = test.benchmark              # string indicating which benchmark is being tested, i.e. LRA, romano
    benchmark_value = test.benchmark_value  # float indicating the known benchmark value
    error_margin = test.error_margin        # float indicating the margin of error to accept for passing
    calculated_value = test.calculated_value # float indicating value found w/ default num threads
    num_threads = test.num_threads

    if test_type != 'Keff':
        print 'ERROR: Unrecognized test type'

    ## fail
    if abs(calculated_value - benchmark_value) >= error_margin:
        write_failed_results(output, test, none_failed)
        this_test_passed = False
        
    ## pass
    elif abs(calculated_value - benchmark_value) < error_margin:
        write_passed_results(test_type, benchmark,num_threads)
        this_test_passed = True

    if print_times:
        print 'Test completed in', time.clock() - test_start_time, 'seconds.'
    
    return this_test_passed

def write_failed_results(output, test, none_failed=True):

    test_type = test.test_type
    benchmark = test.benchmark
    benchmark_value = test.benchmark_value
    error_margin = test.error_margin
    calculated_value = test.calculated_value
    num_threads = test.num_threads

    if num_threads == 'DEFAULT':
        num_threads = mp.cpu_count()

    if none_failed:
        output.write("--------------------------FAILURE--------------------------"+"\n")
        output.write("Date: "+str(time.strftime("%d/%m/%Y"))+", Time: "+str(time.strftime("%H:%M:%S"))+"\n")
        output.write("Platform: "+platform.platform()+"\n")
        output.write("Python version: "+platform.python_version()+"\n")
        output.write("System: "+platform.system()+"\n")
        output.write("Default num_threads: "+str(mp.cpu_count())+"\n")

    output.write("\n")

    output.write("Test Failed: "+benchmark+" "+test_type+', number of threads: '+str(num_threads)+"\n")
    output.write("Benchmark "+test_type+" is "+str(round(benchmark_value,5))+", Keff found was "+str(round(calculated_value,5))+"\n")

    print "------------------FAILURE------------------"
    print "Test Failed: "+benchmark+" "+test_type+', number of threads: '+str(num_threads)    
    print "Benchmark "+test_type+" is "+str(round(benchmark_value,5))+", Keff found was "+str(round(calculated_value,5))
    print "-------------------------------------------"

def write_passed_results(test_type, benchmark, num_threads):
    
    if num_threads == 'DEFAULT':
        num_threads = mp.cpu_count()
        
    print benchmark + " "+test_type+' with '+str(num_threads)+' threads ... ok'
    

