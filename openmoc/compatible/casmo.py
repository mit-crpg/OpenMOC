##
# @file casmo.py
# @package openmoc.compatible.casmo
# @brief The parsing module provides utility functions to parse in data 
#        necessary to construct assembly geometries in OpenMOC
# @author Davis Tran (dvtran@mit.edu)
# @date April 24, 2014

import numpy
import h5py
import os
import openmoc.log as log

##
# @class casmo.py 'openmoc/compatible/casmo.py'
# @brief Contains data parsed from casmo output file
class Casmo(object):

  ##
  # @brief Casmo object class constructor
  def __init__(self):
    self._assembly_name = None
    self._filename = None
    self._directory = None
    self._is_symmetric = True
    self._energy_groups = None
    self._num_micro_regions = None
    self._fuel_pin_rad = None
    self._lattice_pitch = None
    self._siga = None
    self._sigd = None
    self._sigt = None
    self._sigf = None
    self._signf = None
    self._sigs = None
    self._chi = None
    self._width = None
    self._min_microregions = None
    self._max_microregions = None
    self._kinf = None
    self._pin_powers = None
    self._cell_types = {}
    self._cell_type_array = None 
    self._string_cell_type_array = None
    self._average_cross_sections = None


  ##
  # @brief Returns assembly type as string
  # @return assembly type (string)
  def getAssemblyName(self): 
    return self._assembly_name

  ##
  # @brief Sets assembly type
  # @param assembly_name a string that indicates assembly type
  def setAssemblyName(self, assembly_name): 
    self._assembly_name = assembly_name

  ##
  # @brief Returns name of casmo output file to be parsed
  # @return name of casmo output file to be parsed
  def getFilename(self): 
    return self._filename

  ##
  # @brief Sets file name of casmo output file to be parsed
  # @param filename the name of the casmo output file to be parsed (string)
  def setFilename(self, filename): 
    self._filename = filename

  ##
  # @brief Returns directory of casmo output file being parsed
  # @return directory of casmo output file being parsed
  def getDirectory(self):
    return self._directory

  ##
  # @brief Sets directory of casmo output file to be parsed
  # @param directory directory of the casmo output file to be parsed (string)
  def setDirectory(self, directory): 
    self._directory = directory
    
  ##
  # @brief Returns whether the assembly for the casmo output file is symmetric 
  # @return True if symmetric, else False
  def isSymmetric(self):
    return self._is_symmetric

  ##
  # @brief Sets whether the assembly for the casmo output file is symmetric
  # @param is_symmetric boolean indicating whether the geometry is symmetric
  def setSymmetric(self, is_symmetric):
    self._is_symmetric = is_symmetric

  ##
  # @brief Checks to see if assembly for casmo output file is symmetric
  # @param f casmo output file
  def checkSymmetry(self, f):
    sym_counter = 0
    for sym_line in f: 
      if 'LPI' in sym_line: 
        sym_counter += 1      
        continue
      if sym_counter ==1: 
        sym_tokens = sym_line.split()
        if len(sym_tokens) > 2:
          self._is_symmetric = False
          break
        else: 
          self._is_symmetric = True
          break
  

  ##
  # @brief This method parses the casmo output file for the number of
  #        energy groups
  # @return number of energy groups directly from casmo output file
  def parseEnergyGroups(self):
    f = open(self._directory + self._filename,'r')
    for line in f:
      if '[Usage Note]' in line:
        tokens = line.split()
        energy_groups = int(tokens[5])
        break
    f.close()
    return energy_groups
  
  ##
  # @brief Returns number of energy groups
  # @return number of energy groups
  def getEnergyGroups(self): 
    return self._energy_groups

  ##
  # @brief Sets number of energy groups
  # @param energy_groups number of energy groups (int)
  def setEnergyGroups(self, energy_groups):
    self._energy_groups = energy_groups

  ##
  # @brief parses and sets number of energy groups from casmo output file
  def importEnergyGroups(self): 
    self.setEnergyGroups(self.parseEnergyGroups())

  ##
  # @brief This method parses the casmo output file for the number of 
  #        microregions in the assembly
  # @return number of microregions directly from casmo output file
  def parseNumRegions(self):
    f = open(self._directory + self._filename, 'r')
    
    #check for symmetry
    self.checkSymmetry(f)
          
    counter = 0
    newcounter = 0
    num_micro_regions = 0
    if self._is_symmetric:
      for line in f:
        if 'Micro-region number ' in line:
          counter += 1
          continue
        if counter == 1:
          tokens = line.split()
          num_micro_regions = int(tokens[1])
          break
    else:
      for newline in f:
        if '--- ---- ---------------  ------------    ' in newline:
          newcounter += 1
          continue
        if newcounter == 1:
          newtokens = newline.split()
          num_micro_regions = int(newtokens[0])
          break
      
    
    f.close()
    return num_micro_regions

  ##
  # @brief Returns number of microregions in assembly
  # @return number of microregions
  def getNumRegions(self): 
    return self._num_micro_regions

  ##
  # @brief Sets the number of microregions
  # @param num_micro_regions the number of microregions in the assembly
  def setNumRegions(self, num_micro_regions):
    self._num_micro_regions = num_micro_regions

  ##
  # @brief parses and sets number of microregions from casmo output file
  def importNumRegions(self):
    self.setNumRegions(self.parseNumRegions())

  ##
  # @brief This method parses the casmo output file for the thermally
  #        expanded fuel pin radii
  # @return fuel pin radii (float)
  def parseFuelPinRadii(self):
    f = open(self._directory + self._filename, 'r')
    for line in f:
      if 'Average fuel pellet diam.' in line:
        tokens = line.split()
        diameter = tokens[5]
        break
    f.close()
    E = diameter.index('E')
    radii = (0.5 * float(diameter[0:E]) * 10 ** int(diameter[E+1:]))
    return radii

  ##
  # @brief Returns fuel pin radii of the assembly
  # @return fuel pin radii (float)
  def getFuelPinRadii(self): 
    return self._fuel_pin_rad

  ##
  # @brief Sets fuel pin radii of the assembly
  # @param fuel_pin_rad fuel pin radii to be set for assembly (float)
  def setFuelPinRadii(self, fuel_pin_rad):
    self._fuel_pin_rad = fuel_pin_rad

  ##
  # @brief parses and sets fuel pin radii of the assembly
  def importFuelPinRadii(self): 
    self.setFuelPinRadii(self.parseFuelPinRadii())

  ##
  # @brief This method parses the casmo output file for the thermally
  #        expanded lattice pitch
  # @return lattice pitch (float)
  def parseLatticePitch(self):
    f = open(self._directory + self._filename, 'r')
    for line in f:
      if 'Bundle pitch' in line:
        tokens = line.split()
        raw_pitch = tokens[3]
        break
    f.close()
    E = raw_pitch.index('E')
    pitch = (float(raw_pitch[0:E]) * 10 ** int(raw_pitch[E+1:]))
    return pitch

  ##
  # @brief Returns lattice pitch of the assembly
  # @return lattice pitch (float)
  def getLatticePitch(self): 
    return self._lattice_pitch

  ##
  # @brief Sets lattice pitch of the assembly
  # @param lattice_pitch lattice pitch to be set for assembly (float)
  def setLatticePitch(self, lattice_pitch):
    self._lattice_pitch = lattice_pitch

  ##
  # @brief parses and sets lattice pitch of the assembly
  def importLatticePitch(self): 
    self.setLatticePitch(self.parseLatticePitch())    

  ##
  # @brief This method parses the casmo output file for the materials 
  #        cross sections for every microregion in the assembly
  # @param xs_name the name of cross section type (string in all CAPS)
  # @return numpy array of cross sections
  def parseXS(self, xs_name):

    # Parses for cross sections that are not the scattering matrix
    if xs_name != 'SIGS' and xs_name!='CHI':
      xs_array = numpy.zeros((self._num_micro_regions, self._energy_groups))
      f = open(self._directory + self._filename, 'r')
      counter = 0
      for line in f:
        if xs_name in line:
          tokens = line.split()
          xs_array[counter, :] = [float(xs) for xs in tokens[2:2+self._energy_groups]]
          counter += 1
        if counter == self._num_micro_regions:
          break
      f.close()

    # Parses for scattering matrix cross sections
    if xs_name == 'SIGS':
      xs_array = numpy.zeros((self._num_micro_regions, self._energy_groups, self._energy_groups))
      f = open(self._directory + self._filename, 'r')
      cur_region = 0
      cur_group = 0
      for line in f:
        if xs_name in line:
          words = line.split()
          xs_array[cur_region, cur_group, :] = [float(xs) for xs in words[2:2+self._energy_groups]]
          cur_group += 1
        if cur_group == self._energy_groups:
          cur_region += 1
          cur_group = 0
        if cur_region == self._num_micro_regions:
          break
    f.close()
    return xs_array

  ##
  # @brief Returns a specific cross section numpy array
  # @param xs_name the name of a type of cross section (string)
  # @return a cross section numpy array
  def getXS(self, xs_name):
    '''Retrieves cross-section attribute.'''

    if xs_name == 'SIGA':
      return self._siga
    if xs_name == 'SIGD':
      return self._sigd
    if xs_name == 'SIGT':
      return self._sigt
    if xs_name == 'SIGF':
      return self._sigf
    if xs_name == 'SIGNF':
      return self._signf
    if xs_name == 'SIGS':
      return self._sigs
    if xs_name == 'CHI':
      return self._chi

  ##
  # @brief Sets a specific cross section
  # @param xs_name the name of a type of cross section (string)
  # @param xs_array a numpy array of cross section values
  def setXS(self, xs_name, xs_array):

    if xs_name == 'SIGA':
      self._siga = xs_array
    if xs_name == 'SIGD':
      self._sigd = xs_array
    if xs_name == 'SIGT':
      self._sigt = xs_array
    if xs_name == 'SIGF':
      self._sigf = xs_array
    if xs_name == 'SIGNF':
      self._signf = xs_array
    if xs_name == 'SIGS':
      self._sigs = xs_array
    if xs_name == 'CHI':
      self._chi = xs_array

  ##
  # @brief parses and sets a specific cross section type from casmo ouput file
  # @param xs_name the name of a type of cross section (string)
  def importXS(self, xs_name):
    self.setXS(xs_name, self.parseXS(xs_name))

  ##
  # @brief calls importXS for all types of cross sections needed by OpenMOC
  def importAllXS(self):
    xs_list = ['SIGA', 'SIGD', 'SIGT', 'SIGF', 'SIGNF', 'SIGS']
    for xs_name in xs_list:
      self.importXS(xs_name)

  ##
  # @brief This method parses the casmo output file for the dimensions of
  #        the assembly. The width equals the number of fuel pins in a row
  #        or column of an assembly.
  # @return width of the assembly
  def parseWidth(self):
    half_width = -1
    f = open(self._directory + self._filename, 'r')
    
    #check for symmetry
    self.checkSymmetry(f)
    
    
    for line in f:
      if 'Layout' in line:
        half_width += 1
        continue
      if half_width>=0 and line == '\n':
        break
      if half_width>=0:
        half_width += 1
    f.close()
    if self._is_symmetric:
      return half_width*2-1
    else:
      return half_width
      
      
  ##
  # @brief Returns width of the assembly
  # @return width of the assembly (int)
  def getWidth(self):
    return self._width

  ##
  # @brief Sets width of the assembly
  # @param width the width to be set for the assembly
  def setWidth(self, width):
    self._width = width

  ##
  # @brief parses and sets a width of assembly from casmo ouput file
  def importWidth(self):
    self.setWidth(self.parseWidth())

  ##
  # @brief This method parses the casmo output file for microregion ranges
  #        and returns a tuple of two numpy arrays, one is the minimum values
  #        and the other is the maximum values of those microregion ranges, each
  #        each located within its specific macroregion
  # @return numpy array tuple (min microregion values, max microregion values)
  def parseMicroregions(self):

    half_width = (self._width+1)/2
    min_array = numpy.zeros((self._width,self._width), dtype=numpy.int32)
    max_array = numpy.zeros((self._width,self._width), dtype=numpy.int32)
    min_quadrant4 = numpy.zeros((half_width,half_width), dtype=numpy.int32)
    max_quadrant4 = numpy.zeros((half_width,half_width), dtype=numpy.int32)
    min_values = []
    max_values = []
    
   
    

    f = open(self._directory + self._filename, 'r')
    counter = 0
    
    #check for symmetry
    self.checkSymmetry(f)
    
    if self._is_symmetric:
      for line in f:
        if counter >= 1 and '1_________' in line:
          break
        if 'Micro-region' in line:
          counter += 1
          continue
        if counter >= 1:
          tokens = line.split()
          for index, token in enumerate(tokens):
            token = token.strip('*')
            token = token.strip('-')
            if index%2 ==0:
              min_quadrant4[counter-1, index/2] = float(token)
              min_quadrant4[index/2, counter-1] = float(token)
            else:
              max_quadrant4[counter-1, (index-1)/2] = float(token)
              max_quadrant4[(index-1)/2, counter-1] = float(token)
          counter += 1
      f.close()

      min_array[(half_width-1):,(half_width-1):] = min_quadrant4
      min_array[(half_width-1):, 0:(half_width)] = numpy.fliplr(min_quadrant4)
      min_array[0:(half_width), (half_width-1):] = numpy.flipud(min_quadrant4)
      min_array[0:(half_width), 0:(half_width)] = numpy.flipud(numpy.fliplr(min_quadrant4))

      max_array[(half_width-1):,(half_width-1):] = max_quadrant4
      max_array[(half_width-1):, 0:(half_width)] = numpy.fliplr(max_quadrant4)
      max_array[0:(half_width), (half_width-1):] = numpy.flipud(max_quadrant4)
      max_array[0:(half_width), 0:(half_width)] = numpy.flipud(numpy.fliplr(max_quadrant4))
    else:
      counter = 0
      for line in f:
        if 'Micro :' in line:
          newline = line.lstrip('Micro :')
          xline = newline.translate(None, '-')
          tokens = xline.split()
          for index, token in enumerate(tokens):
            if index%2 ==0:
              min_values.append(token)
            else:
              max_values.append(token)
      for index, value in enumerate(min_values):
        min_array[int(counter)/int(self._width), index%self._width] = float(value)
        counter += 1
        continue
      
      counter = 0
      for index, value in enumerate(max_values):
        max_array[int(counter)/int(self._width), index%self._width] = float(value)
        counter += 1
        continue
      f.close()
    return min_array, max_array

  ##
  # @brief Returns numpy array of minimum values of microregion range within
  #        each macroregion
  # @return numpy array of minimum values of microregion ranges
  def getMinMicroregions(self):
    return self._min_microregions

  ##
  # @brief Sets minimum values of microregion ranges within each macroregion
  # @param min_array numpy array of minimum values of microregion ranges
  def setMinMicroregions(self, min_array):
    self._min_microregions = min_array

  ##
  # @brief Returns numpy array of maximum values of microregion ranges within
  #        each macroregion
  # @return numpy array of maximum values of microregion ranges
  def getMaxMicroregions(self):
    return self._max_microregions

  ##
  # @brief Sets minimum values of microregion ranges within each macroregion
  # @param max_array numpy array of minimum values of microregion ranges
  def setMaxMicroregions(self, max_array):
    self._max_microregions = max_array

  ##
  # @brief parses and sets microregion value numpy arrays
  def importMicroregions(self):
      self.setMinMicroregions(self.parseMicroregions()[0])
      self.setMaxMicroregions(self.parseMicroregions()[1])

  ##
  # @brief This method parses the casmo output file for reference eigenvalue
  # @return reference eigenvalue of assembly (float)
  def parseKinf(self):
    f = open(self._directory + self._filename, 'r')

    for line in f:
      if 'k-infinity' in line:
        tokens = line.split()
        kinf = float(tokens[2])
        break
    f.close()
    return kinf

  ##
  # @brief Returns reference eigenvalue of assembly from casmo output file
  # @return reference eigenvalue of assembly (float)
  def getKinf(self):
    return self._kinf

  ##
  # @brief Sets reference eigenvalue of assembly
  # @param kinf the reference eigenvalue to be set for the assembly
  def setKinf(self, kinf):
    self._kinf = kinf

  ##
  # @brief parses and sets eigenvalue of assembly
  def importKinf(self):
    self.setKinf(self.parseKinf())

  ##
  # @brief This method parses the casmo output file for reference pin powers
  # @return numpy array of float-valued reference pin powers of assembly
  def parsePinPowers(self):

    f = open(self._directory + self._filename, 'r')

    half_width = (self._width+1)/2
    pin_power_array = numpy.zeros((self._width,self._width), dtype=numpy.float32)
    quadrant4 = numpy.zeros((half_width,half_width), dtype=numpy.float32)

    counter = 0
    
    #check for symmetry
    self.checkSymmetry(f)
    
    for line in f:
    
      if counter >= 1 and line == '\n':
        break
        
      if 'Power Distribution' in line:
        counter += 1
        continue
        
      if self._is_symmetric:
        if counter >= 1:
          powers = line.split()
          for index, power in enumerate(powers):
            power = power.strip('*')
            quadrant4[counter-1, index] = float(power)
            quadrant4[index, counter-1] = float(power)
          counter += 1
        # Arranges section of pin powers into larger array by symmetry
        pin_power_array[(half_width-1):,(half_width-1):] = quadrant4
        pin_power_array[(half_width-1):, 0:(half_width)] = numpy.fliplr(quadrant4)
        pin_power_array[0:(half_width), (half_width-1):] = numpy.flipud(quadrant4)
        pin_power_array[0:(half_width), 0:(half_width)] = numpy.flipud(numpy.fliplr(quadrant4))
        
      else:
        if counter >= 1:
          powers = line.split()
          
          for index, power in enumerate(powers):
            power = power.strip('*')
            pin_power_array[counter-1, index] = float(power)
          counter+=1
          
    f.close()
    return pin_power_array

  ##
  # @brief Returns reference pin powers of assembly from casmo output file
  # @return numpy array of float valued reference pin powers of assembly
  def getPinPowers(self):
    return self._pin_powers

  ##
  # @brief Sets reference pin powers of assembly
  # @param pin_power_array numpy array of float-valued reference pin powers
  def setPinPowers(self, pin_power_array):
    self._pin_powers = pin_power_array

  ##
  # @brief parses and sets pin powers of assembly
  def importPinPowers(self):
    self.setPinPowers(self.parsePinPowers())

  ##
  # @brief Returns dictionary of cell type associated with each id number
  # @return dictionary cell types by id number, int-->string
  def getCellTypes(self):
    return self._cell_types

  ##
  # @brief Sets a cell type and cell type id key-value pair
  # @param cell_types_id id number for a certain cell type (int)
  # @param name name of a specific cell type associated an id number (string)
  def setCellType(self, cell_types_id, name):
    self._cell_types[cell_types_id] = name

  ##
  # @brief This method parses the casmo output file for the type of material in
  #        each cell
  # @return numpy array of int-valued cell types
  def parseCellTypeArray(self):

    half_width = (self._width+1)/2
    full_width = self._width
    cell_type_array = numpy.zeros((full_width,full_width), dtype=numpy.int32)
    quadrant4 = numpy.zeros((half_width,half_width), dtype=numpy.int32)

    counter = 0
    f = open(self._directory + self._filename, 'r')
    
    #check for symmetry
    self.checkSymmetry(f)
    
    for line in f:
      if counter >=1 and line == '\n':
        break
      if 'Layout' in line:
        counter += 1
        continue
      if counter >= 1:
        cell_types = line.split()
        for index, cell_type in enumerate(cell_types):
          cell_type = cell_type.strip('*')
          if self._is_symmetric:
            quadrant4[counter-1, index] = int(cell_type)
          else:
            cell_type_array[counter-1, index] = int(cell_type)
        counter += 1
    f.close()
    if self._is_symmetric:
      # Arranges section of cell types into larger array by symmetry
      cell_type_array[(half_width-1):,(half_width-1):] = quadrant4
      cell_type_array[(half_width-1):, 0:(half_width)] = numpy.fliplr(quadrant4)
      cell_type_array[0:(half_width), (half_width-1):] = numpy.flipud(quadrant4)
      cell_type_array[0:(half_width), 0:(half_width)] = numpy.flipud(numpy.fliplr(quadrant4))

    cell_type_array[half_width-1,half_width-1] = 2
    return cell_type_array

  ##
  # @brief Returns array of cell type ids for assembly
  # @return array of cell types for every cell in assembly
  def getCellTypeArray(self):
    return self._cell_type_array

  ##
  # @brief Sets array of cell type ids for assembly
  # @param cell_type_array numpy array of int-valued cell type ids
  def setCellTypeArray(self, cell_type_array):
    self._cell_type_array = cell_type_array

  ##
  # @brief parses and sets cell type ids for assembly
  def importCellTypeArray(self):
    self.setCellTypeArray(self.parseCellTypeArray())

  ##
  # @brief This method converts the numerical cell type array to strings that
  #        indicate the cell type in clearer language
  # @return numpy array of cell types as strings
  def stringCellTypeArray(self):

    #id of 1 corresponds to fuel (string of fuel)
    #id of 2 corresponds to guide tube (string of gt)
    #id of 3 corresponds to burnable poison (string of bp)
    string_cell_type_array = numpy.zeros((self._width,self._width), dtype=numpy.str)


    for i, row in enumerate(self._cell_type_array):
      for j, cell in enumerate(row):
        if self._cell_type_array[i,j] in self._cell_types.keys():
          string_cell_type_array[i,j] = self._cell_types[self._cell_type_array[i,j]]
        else:
          log.py_printf('WARNING', 'Cell type id %d does not exist. Call'
          ' setCellTypes to set cell name for id.', self._cell_type_array[i,j])

    return string_cell_type_array

  ##
  # @brief Returns array of cell types as strings for assembly
  # @return array of cell types as strings for assembly
  def getStringCellTypeArray(self):
    return self._string_cell_type_array

  ##
  # @brief Sets array of cell types as strings for assembly
  # @param string_cell_type_array array of cell types as strings
  def setStringCellTypeArray(self, string_cell_type_array):
    self._string_cell_type_array = string_cell_type_array

  ##
  # @brief This method calls the Casmo import methods necessary to construct
  #        the geometry of an assembly in OpenMOC
  # @param filename filename of casmo output file to be parsed
  # @param directory directory of casmo output file to be parsed
  def importFromCasmo(self, filename, directory):
    self._filename = filename
    self._directory = directory
    self.importEnergyGroups()
    self.importNumRegions()
    self.importAllXS()
    self.importWidth()
    self.importMicroregions()
    self.importKinf()
    self.importPinPowers()
    self.importCellTypeArray()
    self.importFuelPinRadii()
    self.importLatticePitch()

  ##
  # @brief This method exports all data contained within member variables
  #        of the Casmo object to an hdf5 data file, data sets expect arrays
  # @param filename filename of hdf5 data file
  # @param directory directory where hdf5 data file will be stored
  def export(self, directory = 'casmo-data/', filename = 'casmo-data.h5'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    f = h5py.File(directory + filename, 'w')
    f.attrs['Energy Groups'] = self._energy_groups
    f.attrs['Assembly Width'] = self._width
    f.attrs['Num Microregions'] = self._num_micro_regions
    f.attrs['Fuel Pin Radii'] = self._fuel_pin_rad
    f.attrs['Lattice Pitch'] = self._lattice_pitch
    big_data = f.create_group('Casmo Data')
    big_data.create_dataset('K-Infinity', data=self._kinf)
    big_data.create_dataset('Total XS', data=self._sigt)
    big_data.create_dataset('Absorption XS', data=self._siga)
    big_data.create_dataset('Fission XS', data=self._sigf)
    big_data.create_dataset('Nu Fission XS', data=self._signf)
    big_data.create_dataset('Scattering XS', data=self._sigs)
    big_data.create_dataset('Dif Coefficient', data=self._sigd)
    big_data.create_dataset('Chi', data=self._chi)
    big_data.create_dataset('Pin Powers', data=self._pin_powers)
    big_data.create_dataset('Cell Types', data=self._cell_type_array)
    big_data.create_dataset('String Cell Types', data=self._string_cell_type_array)
    big_data.create_dataset('Min Microregions', data=self._min_microregions)
    big_data.create_dataset('Max Microregions', data=self._max_microregions)
    f.close()

  ##
  # @brief This method imports data from an hdf5 data file and assigns it
  #        to the corresponding member variables
  # @param filename filename of hdf5 data file
  # @param directory directory where hdf5 data file is stored
  def importFromHDF5(self, directory = 'casmo-data/', filename = 'casmo-data.h5'):

    f = h5py.File(directory + filename, 'r')
    self._directory = directory    
    self._filename = filename
    self._energy_groups = f.attrs['Energy Groups']
    self._kinf = f.attrs['K-Infinity']
    self._width = f.attrs['Assembly Width']
    self._num_micro_regions = f.attrs['Num Microregions']
    self._fuel_pin_rad = f.attrs['Fuel Pin Radii']
    self._lattice_pitch = f.attrs['Lattice Pitch']
    self._sigt = f['Casmo Data']['Total XS'][...]
    self._siga = f['Casmo Data']['Absorption XS'][...]
    self._sigf = f['Casmo Data']['Fission XS'][...]
    self._signf = f['Casmo Data']['Nu Fission XS'][...]
    self._sigs = f['Casmo Data']['Scattering XS'][...]
    self._sigd = f['Casmo Data']['Dif Coefficient'][...]
    self._chi = f['Casmo Data']['Chi'][...]
    self._pin_powers = f['Casmo Data']['Pin Powers'][...]
    self._cell_type_array = f['Casmo Data']['Cell Types'][...]
    self._min_microregions = f['Casmo Data']['Min Microregions'][...]
    self._max_microregions = f['Casmo Data']['Max Microregions'][...]
    f.close()

  ##
  # @brief This method exports only cross sectional arrays contained within 
  #        member variables of the Casmo object to an hdf5 data file
  # @param assembly_name name of assembly for materials being exported
  # @param directory directory where hdf5 data file will be stored
  def exportAllXSToHDF5(self, assembly_name, directory = 'casmo-data'):
    if not os.path.exists(directory):
      os.makedirs(directory)
    f = h5py.File(directory + '/' + assembly_name + '-all-materials.hdf5','w')
    f.attrs['Energy Groups'] = self._energy_groups
    for region in range(self._num_micro_regions):
      material = f.create_group('microregion-' + str((region + 1)))
      material.create_dataset('Total XS', data=self._sigt[region, :])
      material.create_dataset('Absorption XS', data=self._siga[region, :])
      material.create_dataset('Fission XS', data=self._sigf[region, :])
      material.create_dataset('Nu Fission XS', data=self._signf[region, :])
      material.create_dataset('Scattering XS', data=numpy.ravel(self._sigs[region, :, :]))
      material.create_dataset('Dif Coefficient', data=self._sigd[region, :])
      material.create_dataset('Chi', data=self._chi[region, :])
    f.close()


  ##
  # @brief This method exports average cross sectional arrays contained within 
  #        member variables of the Casmo object to an hdf5 data file
  # @param assembly_name name of assembly for materials being exported
  # @param directory directory where hdf5 data file will be stored
  def exportAvgXSToHDF5(self, assembly_name, directory = 'casmo-data'):
  
    #check if cross sections have been computed
    if len(self._average_cross_sections) == 0: 
      log.py_printf('WARNING', 'Average Cross Sections do not exist. Call'
      ' averageXSGenerator to compute them.')
      
    else:
    
      #create/set directory in which to store hdf5 file
      if not os.path.exists(directory):
        os.makedirs(directory)
      f = h5py.File(directory + '/' + assembly_name + '-avg-materials.hdf5','w')
      f.attrs['Energy Groups'] = self._energy_groups
      
      #create an hdf5 dataset to store each average cross section
      for material in self._average_cross_sections.keys():
        material_group = f.create_group(material)
        for xs_type in self._average_cross_sections[material].keys():
          material_group.create_dataset(xs_type,data=self._average_cross_sections[material][xs_type])
      f.close()
   

  ##
  # @brief This method determines the average materials based on average cross 
  #        parsed from the output file
  def averageXSGenerator(self):

    materials = ['fuel','water','cladding','helium']
    
    #check for burnable poisons
    if 'b' in self._string_cell_type_array:
      materials.extend(['bp','ss304'])
      
    #create dictionary of variables
    variable_dict = {'Absorption XS':self._siga,'Dif Coefficient':self._sigd,
      'Total XS':self._sigt,'Fission XS':self._sigf,'Nu Fission XS':self._signf,
      'Scattering XS':self._sigs,'Chi':self._chi}
      
    #create dictionary of values  
    val_dict = {}
    
    #compute average cross section for each material 
    for material in materials:
      val_dict[material] = {}
      for xs_type in variable_dict.keys():
        val_dict[material][xs_type] = []
        
    for i in range(len(self._string_cell_type_array)):
      for j in range(len(self._string_cell_type_array[i])):
        for xs_type in variable_dict.keys():
        
          #if pin cell is guide tube
          if self._string_cell_type_array[i][j]=='g':
            val_dict['water'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]-1])
            val_dict['cladding'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]])
            for k in range(self._min_microregions[i][j]+1,self._max_microregions[i][j]):
              val_dict['water'][xs_type].append(variable_dict[xs_type][k])
              
          #if pin cell is fuel
          elif self._string_cell_type_array[i][j]=='f':
            val_dict['fuel'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]-1])
            val_dict['helium'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]])
            val_dict['cladding'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]+1])
            for k in range(self._min_microregions[i][j]+2,self._max_microregions[i][j]):
              val_dict['water'][xs_type].append(variable_dict[xs_type][k])
              
          #if pin cell is burnable poison
          elif self._string_cell_type_array[i][j]=='b':
            val_dict['helium'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]-1])
            val_dict['ss304'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]])
            val_dict['helium'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]+1])
            val_dict['bp'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]+2])
            val_dict['helium'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]+3])
            val_dict['ss304'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]+4])
            val_dict['water'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]+5])
            val_dict['cladding'][xs_type].append(variable_dict[xs_type][self._min_microregions[i][j]+6])
            for k in range(self._min_microregions[i][j]+7,self._max_microregions[i][j]):
              val_dict['water'][xs_type].append(variable_dict[xs_type][k])
      avg_dict = {}         
      
    #add avg cross sections to dictionary   
    for material in materials:
      avg_dict[material] = {}
      for xs_type in variable_dict.keys():
        avg_dict[material][xs_type] = []
        for group in range(self._energy_groups):
          numerator = sum([e[group] for e in val_dict[material][xs_type]])
          denominator = float(len(val_dict[material][xs_type]))
          if xs_type == 'Scattering XS':
            avg_dict[material][xs_type].extend(numerator/denominator)
          else:
            avg_dict[material][xs_type].append(numerator/denominator)
    self._average_cross_sections = avg_dict
