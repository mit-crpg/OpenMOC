## Currently imports 'materials-test.h5' (which is a copy of 'c5g7-materials.h5')
## and tests the materialize module by comparing the inputted values for all
## parameters to the values stored when the Material instance is created.

## Test cases from c5g7-materials:
##  UO2
##  MOX-4.3%
##  MOX-7%
##  MOX-8.7%
##  Fission Chamber

## Can also add tests for Guide Tube, Water, and Control Rod if that would be useful.

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

import materialize
import unittest
import imp
import h5py

class TestH5Files(unittest.TestCase):

    ## Test the materialize.py module for when the input is a
    ## .py file.

    @classmethod
    def setUpClass(cls):
        
        # Stores the input filename, the imported file,
        # and the output dictionary in the test class.
        
        cls._input = 'materials-test.h5'
        cls._import = h5py.File(cls._input)
        cls._output = materialize.materialize(cls._input)
    
    def testEnergyGroup(self):

        # Asserts the number of energy groups is 7 (from file)
        
        data = self._import
        num_groups = data.attrs['Energy Groups']
        self.assertEqual(num_groups, 7)

    def testNumOfMaterials(self):

        # Asserts the number of materials is the same in the data
        # and in the processed list
        
        data = self._import
        data = list(data)
        output = self._output
        self.assertEqual(len(data), len(output))

    def testNamesOfMaterials(self):

        # Makes sure every material in the imported file
        # is present in the output file.
        
        list_of_results = []
        list_of_true = []
        data = self._import
        data = list(data)
        output = self._output
        for material in data:
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
        cls._input = 'materials-test.h5'
        cls._import = h5py.File(cls._input)
        cls._output = materialize.materialize(cls._input)
        data = cls._import
        cls._data = data
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
        low = 1
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(low, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = UO2_Scattering_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(UO2_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

                # if we just got to num_energy_groups, bump up low
                if group1 == num_energy_groups:
                    low += 1

                # debugging code
##                if abs(UO2_output.getSigmaSByGroup(group1, group2) - SigmaS) >= 0.05:
##                    print 'PROBLEM AREA!!!!', group1, group2

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


#######################################################################################
################################   Test Case: MOX-4.3%   ##############################
#######################################################################################

class TestMOX43(unittest.TestCase):

    ## Test the assigned values to UO2 match those in the input file.

    @classmethod
    def setUpClass(cls):

        cls._input = 'materials-test.h5'
        cls._import = h5py.File(cls._input)
        cls._output = materialize.materialize(cls._input)
        data = cls._import
        cls._data = data
        cls._MOX_input = cls._data['MOX-4.3%']
        cls._MOX_output = cls._output['MOX-4.3%']


    def testMOX43TotalXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_TotalXS_input = MOX_input['Total XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_TotalXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaT = MOX_TotalXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaTByGroup(energy_group) - SigmaT) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX43ScatteringXS(self):

        # The input is a list of num_groups**2 values that should be mapped as follows:
        # list[0] = getSigmaSByGroup(1,1)
        # list[1] = getSigmaSByGroup(2,1)
        # list[2] = getSigmaSByGroup(3,1), etc.

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_ScatteringXS_input = MOX_input['Scattering XS']
        list_of_results = []
        list_of_true = []
        
        num_energy_groups = int((len(MOX_ScatteringXS_input))**0.5)

        #print len(UO2_Scattering_input)
        #print UO2_Scattering_input
        
        input_list_index = 0
        low = 1
        
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(low, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = MOX_ScatteringXS_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(MOX_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

                # if we just got to num_energy_groups, bump up low
                if group1 == num_energy_groups:
                    low += 1

                # debugging code
##                if abs(MOX_output.getSigmaSByGroup(group1, group2) - SigmaS) >= 0.05:
##                    print 'PROBLEM AREA!!!!', group1, group2

        len(list_of_results)
        self.assertEqual(list_of_results, list_of_true)

    def testMOX43FissionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_FissionXS_input = MOX_input['Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_FissionXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaF = MOX_FissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaFByGroup(energy_group) - SigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX43NuFissionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_NuFissionXS_input = MOX_input['Nu Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_NuFissionXS_input)+1):

            # identify what was passed in for SigmaT here
            NuSigmaF = MOX_NuFissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getNuSigmaFByGroup(energy_group) - NuSigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX43Chi(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_Chi_input = MOX_input['Chi']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_Chi_input)+1):

            # identify what was passed in for SigmaT here
            Chi = MOX_Chi_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getChiByGroup(energy_group) - Chi) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX43AbsorptionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_AbsorptionXS_input = MOX_input['Absorption XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaA to the assigned value.
        for energy_group in xrange(1, len(MOX_AbsorptionXS_input)+1):

            # identify what was passed in for SigmaA here
            SigmaA = MOX_AbsorptionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaAByGroup(energy_group) - SigmaA) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

#####################################################################################
################################   Test Case: MOX-7%   ##############################
#####################################################################################


class TestMOX7(unittest.TestCase):

    ## Test the assigned values to UO2 match those in the input file.

    @classmethod
    def setUpClass(cls):
        cls._input = 'materials-test.h5'
        cls._import = h5py.File(cls._input)
        cls._output = materialize.materialize(cls._input)
        data = cls._import
        cls._MOX_input = cls._data['MOX-7%']
        cls._MOX_output = cls._output['MOX-7%']

    def testMOX7TotalXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_TotalXS_input = MOX_input['Total XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_TotalXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaT = MOX_TotalXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaTByGroup(energy_group) - SigmaT) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX7ScatteringXS(self):

        # The input is a list of num_groups**2 values that should be mapped as follows:
        # list[0] = getSigmaSByGroup(1,1)
        # list[1] = getSigmaSByGroup(2,1)
        # list[2] = getSigmaSByGroup(3,1), etc.

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_ScatteringXS_input = MOX_input['Scattering XS']
        list_of_results = []
        list_of_true = []
        
        num_energy_groups = int((len(MOX_ScatteringXS_input))**0.5)

        #print len(UO2_Scattering_input)
        #print UO2_Scattering_input
        
        input_list_index = 0
        low = 1
        
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(low, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = MOX_ScatteringXS_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(MOX_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

                # if we just got to num_energy_groups, bump up low
                if group1 == num_energy_groups:
                    low += 1

                # debugging code
##                if abs(MOX_output.getSigmaSByGroup(group1, group2) - SigmaS) >= 0.05:
##                    print 'PROBLEM AREA!!!!', group1, group2

        len(list_of_results)
        self.assertEqual(list_of_results, list_of_true)

    def testMOX7FissionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_FissionXS_input = MOX_input['Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_FissionXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaF = MOX_FissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaFByGroup(energy_group) - SigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX7NuFissionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_NuFissionXS_input = MOX_input['Nu Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_NuFissionXS_input)+1):

            # identify what was passed in for SigmaT here
            NuSigmaF = MOX_NuFissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getNuSigmaFByGroup(energy_group) - NuSigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX7Chi(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_Chi_input = MOX_input['Chi']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_Chi_input)+1):

            # identify what was passed in for SigmaT here
            Chi = MOX_Chi_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getChiByGroup(energy_group) - Chi) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX7AbsorptionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_AbsorptionXS_input = MOX_input['Absorption XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaA to the assigned value.
        for energy_group in xrange(1, len(MOX_AbsorptionXS_input)+1):

            # identify what was passed in for SigmaA here
            SigmaA = MOX_AbsorptionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaAByGroup(energy_group) - SigmaA) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)


#######################################################################################
################################   Test Case: MOX-8.7%   ##############################
#######################################################################################


class TestMOX87(unittest.TestCase):

    ## Test the assigned values to UO2 match those in the input file.

    @classmethod
    def setUpClass(cls):
        cls._input = 'materials-test.h5'
        cls._import = h5py.File(cls._input)
        cls._output = materialize.materialize(cls._input)
        data = cls._import
        cls._MOX_input = cls._data['MOX-8.7%']
        cls._MOX_output = cls._output['MOX-8.7%']

    def testMOX87TotalXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_TotalXS_input = MOX_input['Total XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_TotalXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaT = MOX_TotalXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaTByGroup(energy_group) - SigmaT) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX87ScatteringXS(self):

        # The input is a list of num_groups**2 values that should be mapped as follows:
        # list[0] = getSigmaSByGroup(1,1)
        # list[1] = getSigmaSByGroup(2,1)
        # list[2] = getSigmaSByGroup(3,1), etc.

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_ScatteringXS_input = MOX_input['Scattering XS']
        list_of_results = []
        list_of_true = []
        
        num_energy_groups = int((len(MOX_ScatteringXS_input))**0.5)

        #print len(UO2_Scattering_input)
        #print UO2_Scattering_input
        
        input_list_index = 0
        low = 1
        
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(low, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = MOX_ScatteringXS_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(MOX_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

                # if we just got to num_energy_groups, bump up low
                if group1 == num_energy_groups:
                    low += 1

                # debugging code
##                if abs(MOX_output.getSigmaSByGroup(group1, group2) - SigmaS) >= 0.05:
##                    print 'PROBLEM AREA!!!!', group1, group2

        len(list_of_results)
        self.assertEqual(list_of_results, list_of_true)

    def testMOX87FissionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_FissionXS_input = MOX_input['Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_FissionXS_input)+1):

            # identify what was passed in for SigmaT here
            SigmaF = MOX_FissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaFByGroup(energy_group) - SigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX87NuFissionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_NuFissionXS_input = MOX_input['Nu Fission XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_NuFissionXS_input)+1):

            # identify what was passed in for SigmaT here
            NuSigmaF = MOX_NuFissionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(UO2_output.getNuSigmaFByGroup(energy_group) - NuSigmaF) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX87Chi(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_Chi_input = MOX_input['Chi']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for energy_group in xrange(1, len(MOX_Chi_input)+1):

            # identify what was passed in for SigmaT here
            Chi = MOX_Chi_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getChiByGroup(energy_group) - Chi) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

    def testMOX87AbsorptionXS(self):

        MOX_input = self._MOX_input
        MOX_output = self._MOX_output
        MOX_AbsorptionXS_input = MOX_input['Absorption XS']
        list_of_results = []
        list_of_true = []

        # for each energy group, compare the value that SHOULD be SigmaA to the assigned value.
        for energy_group in xrange(1, len(MOX_AbsorptionXS_input)+1):

            # identify what was passed in for SigmaA here
            SigmaA = MOX_AbsorptionXS_input[energy_group - 1]
            
            list_of_true.append(True)
            list_of_results.append(abs(MOX_output.getSigmaAByGroup(energy_group) - SigmaA) < 0.05)
        
        self.assertEqual(list_of_results, list_of_true)

#############################################################################################
################################   Test Case: Fission Chamber  ##############################
#############################################################################################


class TestFissionChamber(unittest.TestCase):

    ## Test the assigned values to UO2 match those in the input file.

    @classmethod
    def setUpClass(cls):
        cls._input = 'materials-test.h5'
        cls._import = h5py.File(cls._input)
        cls._output = materialize.materialize(cls._input)
        data = cls._import
        cls._MOX_input = cls._data['Fission Chamber']
        cls._MOX_output = cls._output['Fission Chamber']

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
        low = 1
        
        # for each energy group, compare the value that SHOULD be SigmaT to the assigned value.
        for group2 in xrange(1, num_energy_groups + 1):

            for group1 in xrange(low, num_energy_groups + 1):
            
                # identify what was passed in for SigmaS at this point
                SigmaS = FC_ScatteringXS_input[input_list_index]
                input_list_index += 1           
                list_of_true.append(True)
                list_of_results.append(abs(FC_output.getSigmaSByGroup(group1, group2) - SigmaS) < 0.05)

                # if we just got to num_energy_groups, bump up low
                if group1 == num_energy_groups:
                    low += 1

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
            list_of_results.append(abs(MOX_output.getChiByGroup(energy_group) - Chi) < 0.05)
        
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



        
suite = unittest.TestLoader().loadTestsFromTestCase(TestH5Files)
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestUO2))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestMOX43))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestMOX7))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestMOX87))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestFissionChamber))

unittest.TextTestRunner(verbosity=2).run(suite)
