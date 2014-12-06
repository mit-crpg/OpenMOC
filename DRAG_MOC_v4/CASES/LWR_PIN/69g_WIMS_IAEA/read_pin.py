""" This is an example on how to use DRAG_MOC.py
This file contains the neccesary info for computing OpenMOC material 
files both in python and HDF5 format                                 """

#----THIS PARTICULAR EXAMPLE CORRESPONDS TO THE UAM MOX PWR CASE-------#

#!/usr/local/bin/python
from DRAG_MOC import *
from os import uname
read_ss = DRAG_MOC()


#------------PARAMETERS TO BE DEFINED BY THE USER----------------------#

#----------------------------No. 1-------------------------------------#
read_ss.groups=69
#----------------------------No. 2-------------------------------------#
read_ss.libraryFile="library_ss"
#----------------------------No. 3-------------------------------------#
""" Names defined by the user for the different mixtures         """

read_ss.mixtures=['Moderator','Cladding','Fuel']
read_ss.mixtures_h5=['water','clad','uo2']

#----------------------------No. 4-------------------------------------#
""" Strings of the mixture number identification tag contained in DRAGON
E.g. If in the DRAGON input deck the first mixture is defined as MIX 1, 
here it will be 0001                                                 """

read_ss.mix_drag=['0001', '0002', '0003']

#----------------------------No. 5-------------------------------------#
""" Isotopes as defined exactly in the DRAGON input deck 
(EACH STRING SHOULD CONTAIN 6 TERMS!!!)                              """

read_ss.isotopes=['H1H2O ','O16_1 ',           
                  'Zr91  ','O16_2 ','Fe56  ','Cr52  ',
                  'O16_3 ','U235  ','U238  ']
		  

#----------------------------No. 6-------------------------------------#
""" Atomic densities as exactly in the DRAGON input deck 
(the order of the different elements of this array should match
 the above array)                                                    """
#--------------Atomic densities for the pin case-----------------------#
read_ss.atom_dens=['4.76690E-2','2.38345E-2',
                   '4.18621E-2','3.06711E-4','1.47624E-4','7.54987E-5',
                   '4.49355E-2','7.39237E-4','2.17285E-2']		   

                  
#---------------------------FILES NAMES--------------------------------#
read_ss.isotFileName_py_to_h5="XS_isot_pin_hdf5.py"
read_ss.isotFileName_h5="XS_isot_pin_hdf5.h5"
read_ss.trancFileName_py_to_h5="XS_tranc_pin_hdf5.py"
read_ss.trancFileName_h5="XS_tranc_pin_hdf5.h5"
read_ss.isotFileName_py="XS_isot_pin.py"
read_ss.trancFileName_py="XS_tranc_pin.py"
read_ss.matlabFileName="XS_matlab_pin_file.m"

read_ss.DRAGtoMOC()

