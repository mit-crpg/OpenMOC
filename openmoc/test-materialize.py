## Functions to test:
# materialize(filename)

    # test the error for when filename is NOT a string (should print an error)
    # test with H5/HDF5 data files:
        # check if filename has 'Energy Groups' attribute (an int/fl)
        # material_names = list(f) :: list version of the file
        # for name, should print INFO-level: Importing material (strname)
        # for each, generates a material instance and assigns the num Energy Groups to it
        # checks name of each and sets known values

        # forms dict: materials[name] = new_material

        # test checkSigmaT() ?

    # test .py files
        # same basic tests (loaded differently)
        # data['Energy Groups'] = num groups

        # data = name of file --> dataset LOOK THAT UP
        # later, data = data['Materials']

        # data.keys() = names of materials
        # for each key, there is a dictionary as the stored value (created in the openmoc.Material

    # check output: list of material instances with correct info




    # gonna need 3 files to test on (plus a non-str to test error)
    # one .py, one .hdf5, and one diff type

import materialize
import unittest
import imp
import openmoc

class TestPyFiles(unittest.TestCase):

    ## Test the materialize.py module for when the input is a
    ## .py file.

    @classmethod
    def setUpClass(cls):
        
        # Stores the input filename, the imported file,
        # and the output dictionary in the test class.
        
        cls._input = 'materials-test.py'
        cls._import = imp.load_source(cls._input, cls._input).dataset
        cls._output = materialize.materialize(cls._input)
    
    def testEnergyGroup(self):

        # Asserts the number of energy groups is 7 (from file)
        
        data = self._import
        num_groups = data['Energy Groups']
        self.assertEqual(num_groups, 7)

    def testNumOfMaterials(self):

        # Asserts the number of materials is the same in the data
        # and in the processed list
        
        data = self._import
        data = data['Materials']
        output = self._output
        self.assertEqual(len(data), len(output))

    def testNamesOfMaterials(self):

        # Makes sure every material in the imported file
        # is present in the output file.
        
        list_of_results = []
        list_of_true = []
        data = self._import
        data = data['Materials']
        output = self._output
        for material in data.keys():
            list_of_results.append(material in output)
            list_of_true.append(True)
        self.assertEqual(list_of_results, list_of_true)

class TestUO2(TestPyFiles):

    ## Test the assigned values to UO2 match those in the input file.

    @classmethod
    def setUpClass(cls):
        cls._input = 'materials-test.py'
        cls._import = imp.load_source(cls._input, cls._input).dataset
        cls._output = materialize.materialize(cls._input)
        data = cls._import
        cls._data = data['Materials']
        cls._UO2_input = cls._data['UO2']
        cls._UO2_output = cls._output['UO2']

    def testUO2TotalXS(self):

        ## this is just wrong. 
        
##        UO2_input = self._UO2_input
##        UO2_output = self._UO2_output
##        print UO2_output.getSigmaT() # this is a swig object
##        print UO2_output._swig_getattr(float) # this is not right. not sure how to use wrappers
##        self.assertEqual(UO2_input['Total XS'], UO2_output.getSigmaT())


        UO2_input = self._UO2_input
        UO2_output_orig = self._UO2_output
        UO2_output = self._UO2_output
        UO2_TotalXS_input = UO2_input['Total XS']
        list_of_results = []
        list_of_true = []

        for energy_group in xrange(len(UO2_TotalXS_input)):
            # this would set the values one by one (by energy group)
            UO2_output.setSigmaTByGroup(UO2_TotalXS_input[energy_group], energy_group)

        # could write individual tests for each energy level, or an aggregate test:
            list_of_true.append(True)

            ##############################
            ##############################
            ##############################
            ## FIX THIS WHEN getSigmaTByGroup() has been implemented!!!
            list_of_results.append(UO2_output.getSigmaTByGroup(energy_group) = UO2_TotalXS_input[energy_group])
        
        self.assertEqual(list_of_results, list_of_true)

            
suite = unittest.TestLoader().loadTestsFromTestCase(TestPyFiles)
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUO2))

unittest.TextTestRunner(verbosity=2).run(suite)

