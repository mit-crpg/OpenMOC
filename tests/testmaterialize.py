## Currently imports 'materials-test.py' which is a copy of 'c5g7-materials.py'
## and tests the materialize module by comparing the inputted values for all
## parameters to the values stored when the Material instance is created.

## Test cases from c5g7-materials:
##  UO2
##  MOX-4.3%
##  MOX-7%
##  MOX-8.7%
##  Fission Chamber

## Can also add tests for Guide Tube, Water, and Control Rod if that would be useful.
## Add water (don't need all fuels)

## make sure we get an error if we try energy group = 0
## test isFissionable

## This DOES NOT currently test the following modules because they aren't included in the test input file:
##  getDifCoefByGroup, getBucklingByGroup, getDifHat, getDifTilde

## This DOES NOT test the case where the filename is not a string

## This DOES NOT test the case where the filetype is unsupported.


    # test with H5/HDF5 data files:
        # check if filename has 'Energy Groups' attribute (an int/fl)
        # material_names = list(f) :: list version of the file
        # for name, should print INFO-level: Importing material (strname)
        # for each, generates a material instance and assigns the num Energy Groups to it
        # checks name of each and sets known values

import openmoc.materialize as materialize
import unittest
import imp
import sys
import optparse

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



##################################################################################
################################   Test Case: UO2   ##############################
##################################################################################

class TestUO2(unittest.TestCase):

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

        UO2_input = self._UO2_input
        UO2_output = self._UO2_output
        UO2_TotalXS_input = UO2_input['Total XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(UO2_TotalXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaT = UO2_TotalXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getSigmaTByGroup(energy_group) - SigmaT) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testUO2ScatteringXS(self):

        # The input is a list of num_groups**2 values that should be mapped as follows:
        # list[0] = getSigmaSByGroup(1,1)
        # list[1] = getSigmaSByGroup(2,1)
        # list[2] = getSigmaSByGroup(3,1), etc.


        UO2_input = self._UO2_input
        UO2_output = self._UO2_output
        UO2_Scattering_input = UO2_input['Scattering XS']
        list_of_results = []
        list_of_true = []
        num_energy_groups = int((len(UO2_Scattering_input))**0.5)

        #print len(UO2_Scattering_input)
        #print UO2_Scattering_input
        
        input_list_index = 0
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(1, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = UO2_Scattering_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(UO2_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

                # debugging code
                if abs(UO2_output.getSigmaSByGroup(group1, group2) - SigmaS) >= 0.05:
                    print 'PROBLEM AREA!!!!', group1, group2
                    print 'We expected:', SigmaS, 'We got:', UO2_output.getSigmaSByGroup(group1, group2)

        len(list_of_results)
        self.assertEqual(list_of_results, list_of_true)

    def testUO2FissionXS(self):

        UO2_input = self._UO2_input
        UO2_output = self._UO2_output
        UO2_FissionXS_input = UO2_input['Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(UO2_FissionXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaF = UO2_FissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getSigmaFByGroup(energy_group) - SigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testUO2NuFissionXS(self):

        UO2_input = self._UO2_input
        UO2_output = self._UO2_output
        UO2_NuFissionXS_input = UO2_input['Nu Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(UO2_NuFissionXS_input)+1):

            # identify what was passed in for SigmaT here
            NuSigmaF = UO2_NuFissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getNuSigmaFByGroup(energy_group) - NuSigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testUO2Chi(self):

        UO2_input = self._UO2_input
        UO2_output = self._UO2_output
        UO2_Chi_input = UO2_input['Chi']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(UO2_Chi_input)+1):

            # identify what was passed in for SigmaT here
            Chi = UO2_Chi_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getChiByGroup(energy_group) - Chi) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testUO2AbsorptionXS(self):

        UO2_input = self._UO2_input
        UO2_output = self._UO2_output
        UO2_AbsorptionXS_input = UO2_input['Absorption XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaA to the assigned value.
        for energy_group in xrange(1, len(UO2_AbsorptionXS_input)+1):

            # identify what was passed in for SigmaA here
            SigmaA = UO2_AbsorptionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getSigmaAByGroup(energy_group) - SigmaA) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testUO2IsFissionable(self):
        UO2_output = self._UO2_output
        self.assertTrue(UO2_output.isFissionable())


#############################################################################################
################################   Test Case: Fission Chamber  ##############################
#############################################################################################


class TestFissionChamber(unittest.TestCase):

    ## Test the assigned values to UO2 match those in the input file.

    @classmethod
    def setUpClass(cls):
        cls._input = 'materials-test.py'
        cls._import = imp.load_source(cls._input, cls._input).dataset
        cls._output = materialize.materialize(cls._input)
        data = cls._import
        cls._data = data['Materials']
        cls._FissionChamber_input = cls._data['Fission Chamber']
        cls._FissionChamber_output = cls._output['Fission Chamber']

    def testFCTotalXS(self):

        FC_input = self._FissionChamber_input
        FC_output = self._FissionChamber_output
        FC_TotalXS_input = FC_input['Total XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(FC_TotalXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaT = FC_TotalXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(FC_output.getSigmaTByGroup(energy_group) - SigmaT) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testFCScatteringXS(self):

        # The input is a list of num_groups**2 values that should be mapped as follows:
        # list[0] = getSigmaSByGroup(1,1)
        # list[1] = getSigmaSByGroup(2,1)
        # list[2] = getSigmaSByGroup(3,1), etc.

        FC_input = self._FissionChamber_input
        FC_output = self._FissionChamber_output
        FC_ScatteringXS_input = FC_input['Scattering XS']
        list_of_results = []
        list_of_true = []
        
        num_energy_groups = int((len(FC_ScatteringXS_input))**0.5)
        
        input_list_index = 0
        
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(1, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = FC_ScatteringXS_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(FC_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

                # debugging code
##                if abs(MOX_output.getSigmaSByGroup(group1, group2) - SigmaS) >= 0.05:
##                    print 'PROBLEM AREA!!!!', group1, group2

        len(list_of_results)
        self.assertEqual(list_of_results, list_of_true)

    def testFCFissionXS(self):

        FC_input = self._FissionChamber_input
        FC_output = self._FissionChamber_output
        FC_FissionXS_input = FC_input['Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(FC_FissionXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaF = FC_FissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(FC_output.getSigmaFByGroup(energy_group) - SigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testFCNuFissionXS(self):

        FC_input = self._FissionChamber_input
        FC_output = self._FissionChamber_output
        FC_NuFissionXS_input = FC_input['Nu Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(FC_NuFissionXS_input)+1):

            # identify what was passed in for SigmaT here
            NuSigmaF = FC_NuFissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(FC_output.getNuSigmaFByGroup(energy_group) - NuSigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testFCChi(self):

        FC_input = self._FissionChamber_input
        FC_output = self._FissionChamber_output
        FC_Chi_input = FC_input['Chi']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(FC_Chi_input)+1):

            # identify what was passed in for SigmaT here
            Chi = FC_Chi_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(FC_output.getChiByGroup(energy_group) - Chi) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testFCAbsorptionXS(self):

        FC_input = self._FissionChamber_input
        FC_output = self._FissionChamber_output
        FC_AbsorptionXS_input = FC_input['Absorption XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaA to the assigned value.
        for energy_group in xrange(1, len(FC_AbsorptionXS_input)+1):

            # identify what was passed in for SigmaA here
            SigmaA = FC_AbsorptionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(FC_output.getSigmaAByGroup(energy_group) - SigmaA) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testFCIsFissionable(self):
        FC_output = self._FissionChamber_output
        self.assertTrue(FC_output.isFissionable())


################################################################################
#############################   Test Case: Water  ##############################
################################################################################


class TestWater(unittest.TestCase):

    ## Test the assigned values to UO2 match those in the input file.

    @classmethod
    def setUpClass(cls):
        cls._input = 'materials-test.py'
        cls._import = imp.load_source(cls._input, cls._input).dataset
        cls._output = materialize.materialize(cls._input)
        data = cls._import['Materials']
        cls._Water_input = data['Water']
        cls._Water_output = cls._output['Water']

    def testWaterTotalXS(self):

        Water_input = self._Water_input
        Water_output = self._Water_output
        Water_TotalXS_input = Water_input['Total XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(Water_TotalXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaT = Water_TotalXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(Water_output.getSigmaTByGroup(energy_group) - SigmaT) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testWaterScatteringXS(self):

        # The input is a list of num_groups**2 values that should be mapped as follows:
        # list[0] = getSigmaSByGroup(1,1)
        # list[1] = getSigmaSByGroup(2,1)
        # list[2] = getSigmaSByGroup(3,1), etc.

        Water_input = self._Water_input
        Water_output = self._Water_output
        Water_ScatteringXS_input = Water_input['Scattering XS']
        list_of_results = []
        list_of_true = []
        
        num_energy_groups = int((len(Water_ScatteringXS_input))**0.5)
        
        input_list_index = 0
        
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(1, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = Water_ScatteringXS_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(Water_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

        self.assertEqual(list_of_results, list_of_true)

    def testWaterFissionXS(self):

        Water_input = self._Water_input
        Water_output = self._Water_output
        Water_FissionXS_input = Water_input['Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(Water_FissionXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaF = Water_FissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(Water_output.getSigmaFByGroup(energy_group) - SigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testWaterNuFissionXS(self):

        Water_input = self._Water_input
        Water_output = self._Water_output
        Water_NuFissionXS_input = Water_input['Nu Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(Water_NuFissionXS_input)+1):

            # identify what was passed in for SigmaT here
            NuSigmaF = Water_NuFissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(Water_output.getNuSigmaFByGroup(energy_group) - NuSigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testWaterChi(self):

        Water_input = self._Water_input
        Water_output = self._Water_output
        Water_Chi_input = Water_input['Chi']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(Water_Chi_input)+1):

            # identify what was passed in for SigmaT here
            Chi = Water_Chi_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(Water_output.getChiByGroup(energy_group) - Chi) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testWaterAbsorptionXS(self):

        Water_input = self._Water_input
        Water_output = self._Water_output
        Water_AbsorptionXS_input = Water_input['Absorption XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaA to the assigned value.
        for energy_group in xrange(1, len(Water_AbsorptionXS_input)+1):

            # identify what was passed in for SigmaA here
            SigmaA = Water_AbsorptionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(Water_output.getSigmaAByGroup(energy_group) - SigmaA) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testWaterIsFissionable(self):
        Water_output = self._Water_output
        self.assertFalse(Water_output.isFissionable())

class TestMaterialFile(unittest.TestCase):

    ## Runs basic tests (num energy groups, num materials, names of materials) on an input file
    ## and the Material() object created from that file.

    ## TODO: Create a framework that would also test the values for each material
    ## in the input file.

    def __init__(filename):
        self._input = filename
        
    @classmethod
    def setUpClass(cls):

        cls._import = imp.load_source(cls._input, cls._input).dataset
        cls._output = materialize.materialize(cls._input)

    def testEnergyGroup(self):

        # Asserts the number of energy groups is the same as the lengths of the cross-section lists
        
        data = self._import
        num_groups = data['Energy Groups']
        list_of_results = []
        list_of_true = []

        for material in data['Materials']:
            for attribute in material:
                if attribute != 'Scattering XS':
                    list_of_results.append(len(attribute) == num_groups)
                list_of_true.append(True)
        
        self.assertEqual(list_of_results, list_of_true)

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


##MaterializePySuite = unittest.TestLoader().loadTestsFromTestCase(TestPyFiles)
##MaterializePySuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUO2))
##MaterializePySuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFissionChamber))
##MaterializePySuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestWater))

TestMatPySuite = unittest.TestSuite()
TestMatPySuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestPyFiles))
TestMatPySuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUO2))
TestMatPySuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFissionChamber))
TestMatPySuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestWater))
#PyTestSuite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestMaterialFile))


## Currently all the _getValueByGroup methods are broken -- update

## this is just playing around with command line stuff. please ignore for now
##
##p = optparse.OptionParser()
##p.add_option('--custom', '-c')
##opts, arguments = p.parse_args()
##if '--custom' in opts or '-c' in opts:
##    print 'success'
##    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestMaterialFile))
##    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestCustomMaterials))


## if sys.argv says to use custom one:
##  for attr in input_attributes:
##      add test case

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(TestMatPySuite)

