def __init__():
    pass

class regression_test_case():

    def __init__(self, test_type, benchmark, benchmark_value, calculated_value, calculated_value_1t, error_margin):

        self.test_type = test_type
        self.benchmark = benchmark
        self.benchmark_value = benchmark_value
        self.calculated_value = calculated_value
        self.calculated_value_1t = calculated_value_1t
        self.error_margin = error_margin

class regression_test_suite():

    def __init__(self, tests, output, none_failed=True):

    ## tests = list of (custom) test objects
    ## output = variable referring to output filename
    ## none_failed = boolean referring to whether or not a test has failed thus far in the suite

        self.tests = tests
        self.none_failed = none_failed
        self.output = output

    def add_test(self, test):
        self.tests.append(test)

    def get_tests(self):
        return self.tests

    def test_failed(self):
        self.none_failed = False

    def close_file(self):
        output = self.output
        output.close()

    def run_tests(self):
        
        for test in self.get_tests():

            if not run_regression_test(test, self.output, self.none_failed):
                self.test_failed()

        self.close_file()
            

def run_regression_test(test, output, none_failed=True):

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
    calculated_1t_value = test.calculated_value_1t # float indicating value found w/ 1 thread

    if test_type != 'Keff':
        print 'ERROR: Unrecognized test type'

    ## Test 1: Keff with default num threads

    ## fail
    if abs(calculated_value - benchmark_value) >= error_margin:
        write_failed_results(output, test_type, benchmark, benchmark_value, calculated_value, none_failed)
        none_failed = False
        
    ## pass
    elif abs(calculated_value - benchmark_value) < error_margin:
        write_passed_results(test_type, benchmark)

    ## Test 2: same test_type test, but with 1 thread

    ## fail
    if abs(calculated_1t_value - benchmark_value) >= error_margin:
        write_failed_results(output, test_type+" with 1 thread", benchmark, benchmark_value, calculated_1t_value, none_failed)
        none_failed = False

    ## pass
    elif abs(calculated_1t_value - benchmark_value) < error_margin:
        write_passed_results(test_type+" with 1 thread", benchmark)

    return none_failed

def write_failed_results(output, test_type, benchmark, benchmark_value, calculated_value,none_failed=True):

    if none_failed:
        import time
        output.write("------------------FAILURE------------------"+"\n")
        output.write("Date: "+str(time.strftime("%d/%m/%Y"))+", Time: "+str(time.strftime("%H:%M:%S"))+"\n")

    output.write("Test Failed: "+benchmark+" "+test_type+ "\n")
    output.write("Benchmark "+test_type+" is "+str(benchmark_value)+", Keff found was "+str(calculated_value)+"\n")
    output.write("\n")

    print "------------------FAILURE------------------"
    print "Test Failed: "+benchmark+" "+test_type
    print "Benchmark "+test_type+" is "+str(benchmark_value)+", Keff found was "+str(calculated_value)

def write_passed_results(test_type, benchmark):

    print benchmark + " "+test_type+' ... ok'

    
##        elif test_type == 'pin powers':
##            pin_powers_benchmark = test._pin_powers_benchmark


##        if test_type == 'Keff':

    




