def run_regression_test(self, test):

# tests should have markers (i.e. what benchmark, what is it testing)
# issue: tests currently written for unittest (assertAlmostEqual) and i don't know how to access info from within
    # individual tests - may need to restructure so that each test case is only one comparison
    # and that comparison can be generalized in this runner


    test_type = test._test_type # string indicating what type of test, i.e. Keff, pin powers
    benchmark = test._benchmark # string indicating which benchmark is being tested, i.e. LRA, romano
    
    if test_type == 'Keff':
        
        benchmark_value = test._Keff_benchmark
        calculated_value = test._Keff
        error_margin = test._error_margin
        
    elif test_type == 'pin powers':

        benchmark_value = test._pin_powers_benchmark
        calculated_value = test._pin_powers
        error_margin = test._error_margin # or 0?

    else:
        print 'ERROR: Unrecognized test type'
        break

    if abs(Keff - Keff_benchmark) >= error_margin:
        if not test_failed:
            import time
            output.write("------------------------------------------------------"+"\n")
            output.write("Date: "+str(time.strftime("%d/%m/%Y"))+", Time: "+str(time.strftime("%H:%M:%S"))+"\n")
            test_failed = True
        write_failed_results(test_type, benchmark, benchmark_value, calculated_value)

    else:
        pass

        # TODO:
        # ADD RESULTS WHEN TEST PASSES
        

def write_failed_results(test_type, benchmark, benchmark_value, calculated_value):
        output.write("Test Failed: "+benchmark+" "+test_type+ "\n")
        output.write("Benchmark "+test_type+" is "+str(benchmark_value)+", Keff found was "+str(self._calculated_value)+"\n")
        output.write("\n")

    
##        elif test_type == 'pin powers':
##            pin_powers_benchmark = test._pin_powers_benchmark


##        if test_type == 'Keff':

    




