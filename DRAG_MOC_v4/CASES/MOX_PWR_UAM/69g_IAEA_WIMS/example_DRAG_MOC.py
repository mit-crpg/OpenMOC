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
read_ss.libraryFile="/home/augusto/library_ss"
#----------------------------No. 3-------------------------------------#
""" Variables defined by the user for the different mixtures         """

read_ss.mixtures=['Moderator','Clad_cold','Clad_hot','MOX98','MOX65','MOX37']
read_ss.mixtures_h5=['water','clad_c','clad_h','mox98','mox65','mox37']

#----------------------------No. 4-------------------------------------#
""" Strings of the mixture number identification tag contained in DRAGON
E.g. If in the DRAGON input deck the first mixture is defined as MIX 1, 
here it will be 0001                                                 """

read_ss.mix_drag=['0001', '0002', '0003', '0004', '0005', '0006']

#----------------------------No. 5-------------------------------------#
""" Isotopes as defined exactly in the DRAGON input deck 
(EACH STRING SHOULD CONTAIN 6 TERMS!!!)                              """

#read_ss.isotopes=['H1H2O ','O161','B10   ','B11   ',           
#                  'Zr91  ','O16   ','Fe56  ','Cr52  ',
#		  'Zr91  ','O16   ','Fe56  ','Cr52  ',
#	          'O16   ','U235  ','U238  ','Pu238 ','Pu239 ','Pu240 ','Pu241 ','Pu242 ','AM241 ',
#		  'O16   ','U235  ','U238  ','Pu238 ','Pu239 ','Pu240 ','Pu241 ','Pu242 ','AM241 ',
#		  'O16   ','U235  ','U238  ','Pu238 ','Pu239 ','Pu240 ','Pu241 ','Pu242 ','AM241 ']
read_ss.isotopes=['H1H2O ','O161  ','B10   ','B11   ',           
                  'Zr912 ','O162  ','Fe562 ','Cr522 ',
		  'Zr913 ','O163  ','Fe563 ','Cr523 ',
	          'O164  ','U2354 ','U2384 ','Pu2384','Pu2394','Pu2404','Pu2414','Pu2424','AM2414',
		  'O165  ','U2355 ','U2385 ','Pu2385','Pu2395','Pu2405','Pu2415','Pu2425','AM2415',
		  'O166  ','U2356 ','U2386 ','Pu2386','Pu2396','Pu2406','Pu2416','Pu2426','AM2416']
		  

#----------------------------No. 6-------------------------------------#
""" Atomic densities as exactly in the DRAGON input deck 
(the order of the different elements of this array should match
 the above array)                                                    """
#--------------Atomic densities for the assembly case------------------#
read_ss.atom_dens=['4.70000E-2','2.34000E-2','1.01600E-5','4.09000E-5',
                   '4.18621E-2','3.06711E-4','1.47624E-4','6.24587E-5',
		   '4.18621E-2','3.06711E-4','1.47624E-4','6.24587E-5',
	           '4.62000E-2','5.34900E-5','2.08000E-2','5.64600E-5','1.22600E-3','5.64300E-4','1.91800E-4','1.75500E-4','2.89900E-5',
		   '4.62100E-2','5.34900E-5','2.15600E-2','3.75000E-5','8.14100E-4','3.74900E-4','1.27400E-4','1.16500E-4','1.92600E-5',
		   '4.62200E-2','5.34900E-5','2.22200E-2','2.10900E-5','4.57900E-4','2.10800E-4','7.16500E-5','6.55400E-5','1.08300E-5']
		   

                  
#---------------------------FILES NAMES--------------------------------#
read_ss.isotFileName_py_to_h5="XS_isot_assm_hdf5.py"
read_ss.isotFileName_h5="XS_isot_assm_hdf5.h5"
read_ss.trancFileName_py_to_h5="XS_tranc_assm_hdf5.py"
read_ss.trancFileName_h5="XS_tranc_assm_hdf5.h5"
read_ss.isotFileName_py="XS_isot_assm.py"
read_ss.trancFileName_py="XS_tranc_assm.py"
read_ss.matlabFileName="XS_matlab_assm_file.m"

read_ss.DRAGtoMOC()

