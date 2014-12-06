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
""" Variables defined by the user for the different mixtures         """

read_ss.mixtures=['Moderator','Tube1','Tube2','Tube3','Tube4','Tube5',
                  'Fuel','BA']
read_ss.mixtures_h5=['water','clad_1','clad_2','clad_3','clad_4','clad_5',
                     'uo2','uo2Gd']

#----------------------------No. 4-------------------------------------#
""" Strings of the mixture number identification tag contained in DRAGON
E.g. If in the DRAGON input deck the first mixture is defined as MIX 1, 
here it will be 0001                                                 """

read_ss.mix_drag=['0001', '0002', '0003', '0004', '0005', '0006', 
                  '0007', '0008']

#----------------------------No. 5-------------------------------------#
""" Isotopes as defined exactly in the DRAGON input deck 
(EACH STRING SHOULD CONTAIN 6 TERMS!!!)                              """

read_ss.isotopes=['H1H2O1','O161  ','BNat1 ',           
                  'Zr912 ','O162  ','Fe562 ','Cr522 ',
		  'Zr913 ','O163  ','Fe563 ','Cr523 ','H1H2O3','Ni583 ','MoNat3','Al273 ','Mn553 ','BNat3 ',
                  'Zr914 ','O164  ','Fe564 ','Cr524 ',
                  'Zr915 ','O165  ','Fe565 ','Cr525 ','H1H2O5','Ni585 ','MoNat5','Al275 ','Mn555 ','BNat5 ',
                  'Zr916 ','O166  ','Fe566 ','Cr526 ','H1H2O6','Ni586 ','MoNat6','Al276 ','Mn556 ','BNat6 ',
                  'O167  ','U2357 ','U2387 ',
		  'O168  ','U2358 ','U2388 ','Gd1568','Gd1578']
		  

#----------------------------No. 6-------------------------------------#
""" Atomic densities as exactly in the DRAGON input deck 
(the order of the different elements of this array should match
 the above array)                                                    """
#--------------Atomic densities for the assembly case------------------#
read_ss.atom_dens=['4.76690E-2','2.38345E-2','2.38103E-5',
                   '4.18621E-2','3.06711E-4','1.47624E-4','7.54987E-5',
      '8.92427E-4','2.32646E-2','4.45845E-5','4.79927E-5','4.65292E-2','1.13521E-4','4.03755E-6','2.35231E-6','4.15901E-7','2.32761E-5',
	           '3.92175E-2','2.87335E-4','1.38298E-4','7.07291E-5',
      '4.18372E-4','2.35673E-2','2.09013E-5','2.24991E-5','4.71346E-2','5.32188E-5','1.89281E-6','1.10277E-6','1.94976E-7','2.35598E-5',
      '3.92583E-4','2.35838E-2','1.96130E-5','2.11122E-5','4.71676E-2','4.99383E-5','1.77614E-6','1.03479E-6','1.82957E-7','2.35753E-5',
                   '4.49355E-2','7.39237E-4','2.17285E-2',
                   '4.49355E-2','7.39237E-4','2.17285E-2','5.62500E-4','4.17900E-4']		   

                  
#---------------------------FILES NAMES--------------------------------#
read_ss.isotFileName_py_to_h5="XS_isot_assm_hdf5.py"
read_ss.isotFileName_h5="XS_isot_assm_hdf5.h5"
read_ss.trancFileName_py_to_h5="XS_tranc_assm_hdf5.py"
read_ss.trancFileName_h5="XS_tranc_assm_hdf5.h5"
read_ss.isotFileName_py="XS_isot_assm.py"
read_ss.trancFileName_py="XS_tranc_assm.py"
read_ss.matlabFileName="XS_matlab_assm_file.m"

read_ss.DRAGtoMOC()

