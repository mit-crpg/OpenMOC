#-------------------------------------------------------------------------#
""" 
This program reads self-shielded (SS) microscopic cross-sections (XS) 
from ASCII outputs of the module SHI of the DRAGON code, and creates 
material files for OpenMOC.

     AUTHOR: Augusto Hernandez-Solis, October 2014 @ KTH, Sweden
     
In general, from the version 4 of the DRAGON code, two different library 
formats where microscopic XS are tabulated as a function of temperatures 
and background XS can be handled: the WIMS and the DRAGLIB formats. 
More info about those type of libraries can be found respectively in the
following links:

https://www-nds.iaea.org/wimsd/

www.polymtl.ca/merlin/version5.htm
     
                                                                        """
#-------------------------------------------------------------------------#
                                                                     

from operator import sub
from operator import add
import os,time

class DRAG_MOC_Error(Exception):
  """Exception indicating an error in DRAG_MOC"""

class DRAG_MOC:
  def DRAGtoMOC(self): 	
      
   myCwd = os.getcwd()
   if not os.path.isfile(os.path.expandvars(self.libraryFile)):
    raise DRAG_MOC_Error("SS library file "+self.libraryFile+" not found")
   
   """ This array matches the different reaction names (in DRAGON format)  
   that are contained in the SS XS library produced by DRAGON. If other 
   reaction is wanted to be added, it should match the syntax of DRAGON """

   reactions=[
   'NTOT0                                                                           ',
   'SIGS00                                                                          ',
   'SIGS01                                                                          ',
   'SCAT00                                                                          ',
   'NJJS00                                                                          ',
   'IJJS00                                                                          ',
   'SCAT01                                                                          ',
   'NJJS01                                                                          ',
   'IJJS01                                                                          ',
   'CHI                                                                             ',
   'NUSIGF                                                                          ',
   'NFTOT                                                                           ',
   'NG                                                                              ' ]
	   
   reac_py=['NTOT0','SIG00','SIG01','SCAT0','NJS00','IJS00','SCAT1','NJS01'
   ,'IJS01','CHI_X','NUSIF','NFTOT','N_GAM']

   #----------------------------------------------------------------------#
   #-------------LOOP OVER THE NUMBER OF ISOTOPES-------------------------#
   global GROUPS
   GROUPS=self.groups
   reac=0
   mix_aux2 =[]
   isot_aux2=[]
   print("------------Reading the following data---------------------")
   while reac < len(reactions):
    mix=0
    while mix < len(self.mix_drag):
     exec("f=open('%s','r')"%self.libraryFile) 
     while True:
      a=f.readline()
      if len(a)==0:
       break
      if a[8:12]==self.mix_drag[mix]:
       isot_aux=a[0:6]
       a=f.readline()
       a=f.readline()
       if a[0:3]=='AWR':
        flag=0
        while a[0:10]!='->      -2':
         a=f.readline()
	 if a[0:80]==reactions[reac]:
          flag=1
	  if reac==1:
	   mix_aux2.append(self.mix_drag[mix])
	   isot_aux2.append(isot_aux)
          a=f.readline()
          print("-------Mixture %s, Isotope %s, Reaction %s -------\n" 
          %(self.mix_drag[mix],isot_aux,reac_py[reac]))
          exec("vec_%s_%s_%s=[]"%(reac_py[reac],self.mix_drag[mix],isot_aux))
          while a[0]!='-':
        
           b=a.split()
           c= [float(d) for d in b] 
        
           for i in range(len(c)):
            exec("vec_%s_%s_%s.append(c[i])"%(reac_py[reac],self.mix_drag[mix],isot_aux))
	   a=f.readline()
          #
	 #
	#
        if flag==0:
         exec("vec_%s_%s_%s=[0]*GROUPS"%(reac_py[reac],self.mix_drag[mix],isot_aux))
        #
       #	
      #  
     f.close()
     mix=mix+1
    reac=reac+1 
   #-----------------------END OF READING-----------------------------#

   for ind in range(len(mix_aux2)):
    #print("%s,%s"%(mix_aux2[ind], isot_aux2[ind]))
    for ind2 in range(len(self.isotopes)):
     isotope_ind=self.isotopes[ind2]
     if isotope_ind==isot_aux2[ind]:
      atom_dens_ind=self.atom_dens[ind2]
     
    #---------------SCATTERING MATRICES COMPUTATION-------------------#
    """ 
    Calculation of the scattering matrices columns (S0, S1 and SCAT_TRANC)
    This calculation is based on the algorithm that DRAGON has implemented
    to only print non-zero values in case there is no up-scattering. This 
    is why it is necessary to read NJJS00, IJJS00, NJJS01 and IJJS01 from 
    the SS library. For more info about the algorithm, check the DRAGON
    user guide.
    Remember that the WIMS library only defines P1 matrices for moderator
    materials.                                                         """
  
    diff_aux0=[]
    diff_aux1=[]
    for g in range(GROUPS):
   
     exec("col_%s_%s_SCAT0_MTX_%s=[0]*GROUPS"%(g,mix_aux2[ind],isot_aux2[ind]))
     exec("col_%s_%s_SCAT1_MTX_%s=[0]*GROUPS"%(g,mix_aux2[ind],isot_aux2[ind]))
     exec("diff_aux0=map(sub,vec_IJS00_%s_%s,vec_NJS00_%s_%s)"
     %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
     exec("diff_aux1=map(sub,vec_IJS01_%s_%s,vec_NJS01_%s_%s)"
     %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))

     if g==0:
      #---------------------GROUP No. 1---------------------------------#
      from_NJS00=0
      exec("to_NJS00=vec_NJS00_%s_%s[g]"%(mix_aux2[ind],isot_aux2[ind]))
     
      from_NJS01=0
      exec("to_NJS01=vec_NJS01_%s_%s[g]"%(mix_aux2[ind],isot_aux2[ind]))

     elif g==1:
      #----------------------GROUP No. 2--------------------------------#
      exec("from_NJS00=vec_NJS00_%s_%s[g-1]"%(mix_aux2[ind],isot_aux2[ind]))
      exec("to_NJS00=vec_NJS00_%s_%s[g]+vec_NJS00_%s_%s[g-1]"
      %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
      exec("sum_aux0=vec_NJS00_%s_%s[g-1]"%(mix_aux2[ind],isot_aux2[ind])) 
      
      exec("from_NJS01=vec_NJS01_%s_%s[g-1]"%(mix_aux2[ind],isot_aux2[ind]))
      exec("to_NJS01=vec_NJS01_%s_%s[g]+vec_NJS01_%s_%s[g-1]"
      %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
      exec("sum_aux1=vec_NJS01_%s_%s[g-1]"%(mix_aux2[ind],isot_aux2[ind]))
     
     else:
      #-------------------GROUP GT. 2-----------------------------------#
      exec("sum_aux0=sum_aux0+vec_NJS00_%s_%s[g-1]"%(mix_aux2[ind],isot_aux2[ind]))
      from_NJS00=sum_aux0
      exec("to_NJS00=from_NJS00+vec_NJS00_%s_%s[g]"%(mix_aux2[ind],isot_aux2[ind]))

      exec("sum_aux1=sum_aux1+vec_NJS01_%s_%s[g-1]"%(mix_aux2[ind],isot_aux2[ind]))
      from_NJS01=sum_aux1
      exec("to_NJS01=from_NJS01+vec_NJS01_%s_%s[g]"%(mix_aux2[ind],isot_aux2[ind]))
     
     col_aux01=[]
     col_aux11=[]
     exec("col_aux01=vec_SCAT0_%s_%s[int(from_NJS00):int(to_NJS00)]"
     %(mix_aux2[ind],isot_aux2[ind]))
     exec("col_aux11=vec_SCAT1_%s_%s[int(from_NJS01):int(to_NJS01)]"
     %(mix_aux2[ind],isot_aux2[ind]))
     
     col_aux02=[0]*len(col_aux01)
     count=0
     end=len(col_aux01)
     while count < len(col_aux01):
      col_aux02[count] = col_aux01[end-1]
      count=count+1
      end=end-1
          
     col_aux12=[0]*len(col_aux11)
     count=0
     end=len(col_aux11)
     while count < len(col_aux11):
      col_aux12[count] = col_aux11[end-1]
      count=count+1
      end=end-1
      
     #---------P0 AND P1 SCATTERING MATRICES COULMNS COMPUTATION-------#
     count=0
     while count < len(col_aux02):
      exec("col_%s_%s_SCAT0_MTX_%s[int(diff_aux0[g])+count]= \
      col_aux02[count]"%(g,mix_aux2[ind],isot_aux2[ind]))
      count=count+1

     count=0
     while count < len(col_aux12):
      exec("col_%s_%s_SCAT1_MTX_%s[int(diff_aux1[g])+count]= \
      col_aux12[count]"%(g,mix_aux2[ind],isot_aux2[ind]))
      count=count+1
      
     #TRANSPORT-CORRECTION IN THE COLUMNS OF THE MICROSCOPIC SCATTERING
     #MATRIX (i.e. SCAT_MTX = SCAT0 - SCAT1)--------------------------#
     exec("col_%s_%s_mic_SCAT_TRANC_MTX_%s=[]"%(g,mix_aux2[ind],isot_aux2[ind]))
     
     exec("col_%s_%s_mic_SCAT_TRANC_MTX_%s=map(sub,col_%s_%s_SCAT0_MTX_%s, \
     col_%s_%s_SCAT1_MTX_%s)"%(g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind], 
     isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind]))
     
     #MACROSCOPIC SCATTERING MATRIX (TRANSPORT-CORRECTED) AND P1 MATRIX 
     #COLUMNS COMPUTATION---------------------------------------------#
     exec("col_%s_%s_mac_SCAT_TRANC_MTX_%s=[]"%(g,mix_aux2[ind],isot_aux2[ind]))
     
     exec("col_%s_%s_mac_SCAT_TRANC_MTX_%s= \
     [float(atom_dens_ind)*col_%s_%s_mic_SCAT_TRANC_MTX_%s[i] \
     for i in range(len(col_%s_%s_mic_SCAT_TRANC_MTX_%s))]"
     %(g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind]))
     
     exec("col_%s_%s_mac_SCAT_P1_MTX_%s=[]"%(g,mix_aux2[ind],isot_aux2[ind]))
    
     exec("col_%s_%s_mac_SCAT_P1_MTX_%s= \
     [float(atom_dens_ind)*col_%s_%s_SCAT1_MTX_%s[i] \
     for i in range(len(col_%s_%s_SCAT1_MTX_%s))]"
     %(g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind]))


     #MACROSCOPIC SCATTERING MATRIX (ISOTROPIC) COLUMNS COMPUTATION---#
     exec("col_%s_%s_mic_SCAT_ISOT_MTX_%s=[]"%(g,mix_aux2[ind],isot_aux2[ind]))
     
     exec("col_%s_%s_mic_SCAT_ISOT_MTX_%s=col_%s_%s_SCAT0_MTX_%s" \
     %(g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind]))
     
     exec("col_%s_%s_mac_SCAT_ISOT_MTX_%s=[]"%(g,mix_aux2[ind],isot_aux2[ind]))

     exec("col_%s_%s_mac_SCAT_ISOT_MTX_%s= \
     [float(atom_dens_ind)*col_%s_%s_mic_SCAT_ISOT_MTX_%s[i] \
     for i in range(len(col_%s_%s_mic_SCAT_ISOT_MTX_%s))]"
     %(g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind],g,mix_aux2[ind],isot_aux2[ind]))
     #exec("print col_%s_%s_mac_SCAT_ISOT_MTX_%s"%(g,mix_aux2[ind],isot_aux2[ind]))

    #-----CALCULATION OF GROUP MACROSCOPIC XS (EXCEPT CHI)--------------#
   
    exec("n_abs_mic_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_scat_tranc_mic_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_scat_isot_mic_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_tot_isot_mic_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_tot_tranc_mic_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))


    exec("n_abs_mic_%s_%s= \
    map(sub,vec_NTOT0_%s_%s,vec_SIG00_%s_%s)"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind])) 
    exec("n_scat_tranc_mic_%s_%s= \
    map(sub,vec_SIG00_%s_%s,vec_SIG01_%s_%s)"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))  
    exec("n_scat_isot_mic_%s_%s=vec_SIG00_%s_%s" \
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
    exec("n_tot_isot_mic_%s_%s=vec_NTOT0_%s_%s" \
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))  
    exec("n_tot_tranc_mic_%s_%s= \
    map(sub,n_tot_isot_mic_%s_%s,vec_SIG01_%s_%s)"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))

    exec("n_abs_mac_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_scat_tranc_mac_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_scat_isot_mac_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_tot_isot_mac_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_tot_tranc_mac_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_nufis_mac_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("n_fis_mac_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))
    exec("chi_mac_aux_%s_%s=[]"%(mix_aux2[ind],isot_aux2[ind]))

    exec("n_abs_mac_%s_%s= \
    [float(atom_dens_ind)*n_abs_mic_%s_%s[i] \
    for i in range(len(n_abs_mic_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
   
    exec("n_scat_tranc_mac_%s_%s= \
    [float(atom_dens_ind)*n_scat_tranc_mic_%s_%s[i] \
    for i in range(len(n_scat_tranc_mic_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
   
    exec("n_scat_isot_mac_%s_%s= \
    [float(atom_dens_ind)*n_scat_isot_mic_%s_%s[i] \
    for i in range(len(n_scat_isot_mic_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
    exec("n_tot_isot_mac_%s_%s= \
    [float(atom_dens_ind)*n_tot_isot_mic_%s_%s[i] \
    for i in range(len(n_tot_isot_mic_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
    exec("n_tot_tranc_mac_%s_%s= \
    [float(atom_dens_ind)*n_tot_tranc_mic_%s_%s[i] \
    for i in range(len(n_tot_tranc_mic_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
    exec("n_nufis_mac_%s_%s= \
    [float(atom_dens_ind)*vec_NUSIF_%s_%s[i] \
    for i in range(len(vec_NUSIF_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
    exec("n_fis_mac_%s_%s= \
    [float(atom_dens_ind)*vec_NFTOT_%s_%s[i] \
    for i in range(len(vec_NFTOT_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
    
    exec("sum_nufis_mac_%s_%s=0"%(mix_aux2[ind],isot_aux2[ind]))   
   
    for g_aux in range(GROUPS):
     exec("sum_nufis_mac_%s_%s= \
     sum_nufis_mac_%s_%s+n_nufis_mac_%s_%s[g_aux]"
     %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))
    
    exec("chi_mac_aux_%s_%s= \
    [(vec_CHI_X_%s_%s[i])*sum_nufis_mac_%s_%s \
    for i in range(len(vec_CHI_X_%s_%s))]"
    %(mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],isot_aux2[ind],mix_aux2[ind],
    isot_aux2[ind],mix_aux2[ind],isot_aux2[ind]))    
   
   #----------MACROSCOPIC XS COMPUTATION FOR THE DIFFERENT MIXTURES-------#
 
   print "Now computing macroscopic XS for the different mixtures"
   for mix in range(len(self.mixtures)):

    exec("n_tot_tranc_mac_%s=[0]*GROUPS"%(self.mixtures[mix]))
    exec("n_tot_isot_mac_%s=[0]*GROUPS"%(self.mixtures[mix]))
    exec("n_abs_mac_%s=[0]*GROUPS"%(self.mixtures[mix]))
    exec("n_scat_isot_mac_%s=[0]*GROUPS"%(self.mixtures[mix]))
    exec("n_scat_tranc_mac_%s=[0]*GROUPS"%(self.mixtures[mix]))
    exec("n_fis_mac_%s=[0]*GROUPS"%(self.mixtures[mix]))
    exec("n_nufis_mac_%s=[0]*GROUPS"%(self.mixtures[mix]))

    for g in range(GROUPS):
 
     exec("col_%s_mac_SCAT_ISOT_MTX_%s=[0]*GROUPS"%(g,self.mixtures[mix]))
     exec("col_%s_mac_SCAT_TRANC_MTX_%s=[0]*GROUPS"%(g,self.mixtures[mix]))
     exec("col_%s_mac_SCAT_P1_MTX_%s=[0]*GROUPS"%(g,self.mixtures[mix]))


     for isot in range(len(self.isotopes)):

      try:
       exec("n_tot_tranc_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("n_tot_tranc_mac_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))

      try:
       exec("n_tot_isot_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("n_tot_isot_mac_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))

      try:
       exec("n_abs_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("n_abs_mac_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))

      try:
       exec("n_scat_isot_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("n_scat_isot_mac_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))

      try:
       exec("n_scat_tranc_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("n_scat_tranc_mac_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))

      try:
       exec("n_fis_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("n_fis_mac_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))

      try:
       exec("n_nufis_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("n_nufis_mac_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))
       
      try:
       exec("chi_mac_aux_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("chi_mac_aux_%s_%s=[0]*GROUPS"%(self.mix_drag[mix],self.isotopes[isot]))    
       
      try:
       exec("sum_nufis_mac_%s_%s"%(self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("sum_nufis_mac_%s_%s=0.0"%(self.mix_drag[mix],self.isotopes[isot]))    
   
      try:
       exec("col_%s_%s_mac_SCAT_TRANC_MTX_%s"%(g,self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("col_%s_%s_mac_SCAT_TRANC_MTX_%s=[0]*GROUPS"
       %(g,self.mix_drag[mix],self.isotopes[isot]))

      try:
       exec("col_%s_%s_mac_SCAT_ISOT_MTX_%s"%(g,self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("col_%s_%s_mac_SCAT_ISOT_MTX_%s=[0]*GROUPS"
       %(g,self.mix_drag[mix],self.isotopes[isot]))
    
      try:
       exec("col_%s_%s_mac_SCAT_P1_MTX_%s"%(g,self.mix_drag[mix],self.isotopes[isot]))
      except NameError:
       exec("col_%s_%s_mac_SCAT_P1_MTX_%s=[0]*GROUPS"
       %(g,self.mix_drag[mix],self.isotopes[isot]))
   
      for g_aux in range(GROUPS):

       exec("col_%s_mac_SCAT_TRANC_MTX_%s[g_aux]= \
       col_%s_mac_SCAT_TRANC_MTX_%s[g_aux]+col_%s_%s_mac_SCAT_TRANC_MTX_%s[g_aux]"
       %(g,self.mixtures[mix],g,self.mixtures[mix],g,self.mix_drag[mix],self.isotopes[isot]))

       exec("col_%s_mac_SCAT_ISOT_MTX_%s[g_aux]= \
       col_%s_mac_SCAT_ISOT_MTX_%s[g_aux]+col_%s_%s_mac_SCAT_ISOT_MTX_%s[g_aux]"
       %(g,self.mixtures[mix],g,self.mixtures[mix],g,self.mix_drag[mix],self.isotopes[isot]))

       exec("col_%s_mac_SCAT_P1_MTX_%s[g_aux]= \
       col_%s_mac_SCAT_P1_MTX_%s[g_aux]+col_%s_%s_mac_SCAT_P1_MTX_%s[g_aux]"
       %(g,self.mixtures[mix],g,self.mixtures[mix],g,self.mix_drag[mix],self.isotopes[isot]))

      exec("n_tot_tranc_mac_%s[g]= \
      n_tot_tranc_mac_%s[g]+n_tot_tranc_mac_%s_%s[g]"
      %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
   
      exec("n_tot_isot_mac_%s[g]= \
      n_tot_isot_mac_%s[g]+n_tot_isot_mac_%s_%s[g]"
      %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
   
      exec("n_abs_mac_%s[g]= \
      n_abs_mac_%s[g]+n_abs_mac_%s_%s[g]"
      %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
   
      exec("n_scat_isot_mac_%s[g]= \
      n_scat_isot_mac_%s[g]+n_scat_isot_mac_%s_%s[g]"
      %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
   
      exec("n_scat_tranc_mac_%s[g]= \
      n_scat_tranc_mac_%s[g]+n_scat_tranc_mac_%s_%s[g]"
      %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
   
      exec("n_fis_mac_%s[g]= \
      n_fis_mac_%s[g]+n_fis_mac_%s_%s[g]"
      %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
   
      exec("n_nufis_mac_%s[g]= \
      n_nufis_mac_%s[g]+n_nufis_mac_%s_%s[g]"
      %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
     #
    #
   #----------------------------------------------------------------------#

   """ Evaluating the fission spectra (CHI). 
   The spectra is weighted among the different fissionable isotopes that 
   are contained in the fuel part.                                      """

   #----------------------------------------------------------------------#
   for mix in range(len(self.mixtures)):
    if eval("n_nufis_mac_%s[0]"%(self.mixtures[mix]))==0:
     exec("chi_final_%s=[0]*GROUPS"%self.mixtures[mix])
   
    else:
     exec("chi_final_%s=[0]*GROUPS"%self.mixtures[mix])
       
     for g in range(GROUPS):
      
      chi_aux_den=0
      
      for isot in range(len(self.isotopes)):
       exec("chi_final_%s[g]=chi_final_%s[g]+chi_mac_aux_%s_%s[g]"
       %(self.mixtures[mix],self.mixtures[mix],self.mix_drag[mix],self.isotopes[isot]))
       exec("chi_aux_den=chi_aux_den+sum_nufis_mac_%s_%s"
       %(self.mix_drag[mix],self.isotopes[isot]))
      exec("chi_final_%s[g]=chi_final_%s[g]/chi_aux_den"
      %(self.mixtures[mix],self.mixtures[mix]))
     #
     
   #--------------CREATING THE ISOTROPIC FILE (HD5 FORMAT)----------------#

   exec("m=open('%s','w')"%self.isotFileName_py_to_h5)

   exec("print 'Now printing the file for OpenMOC materials %s \
   (isotropic)'%self.isotFileName_py_to_h5")

   m.write('import h5py\n')
   m.write('import numpy\n\n')
   m.write('# Create the file to store multi-groups cross-sections\n')
   m.write("f = h5py.File('%s')\n"%self.isotFileName_h5)
   m.write("f.attrs['Energy Groups'] = %s\n"%GROUPS)

   for mix in range(len(self.mixtures)):
 
     m.write('##########################################################\n')
     m.write('######################      %s      #######################\n' \
     %self.mixtures_h5[mix])
     m.write('###########################################################\n\n')

     m.write('# Create a subgroup for %s materials data\n'%self.mixtures[mix])
     m.write("%s=f.create_group('%s')\n\n"%(self.mixtures_h5[mix],self.mixtures[mix]))
     m.write("sigma_t=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_tot_isot_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_tot_isot_mac_%s[g])"%self.mixtures[mix])))

     m.write("sigma_a=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))  

     m.write("sigma_f=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))

     m.write("nu_sigma_f=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))   

     m.write("chi=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(chi_final_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(chi_final_%s[g])"%self.mixtures[mix]))) 

     m.write("sigma_s=numpy.array([")
 
     for g in range(GROUPS):
      for g_aux in range(GROUPS):
       if g==GROUPS-1 and g_aux==GROUPS-1:
        m.write("%.8f])\n "%(eval("(col_%s_mac_SCAT_ISOT_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g))))
       else:
        m.write("%.8f, "%(eval("(col_%s_mac_SCAT_ISOT_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g)))) 
        
     m.write('\n# Create datasets for each cross-section type\n')
     m.write("%s.create_dataset('Total XS', data=sigma_t)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Absorption XS', data=sigma_a)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Scattering XS', data=sigma_s)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Fission XS', data=sigma_f)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Nu Fission XS', data=nu_sigma_f)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Chi', data=chi)\n\n"%self.mixtures_h5[mix])
     m.write('# Close the hdf5 data file\n')
   m.write('f.close()\n')
   m.close()
   #----------------------------------------------------------------------#
   
   #----------CREATING THE TRANS. CORRECTED FILE (HDF5 FORMAT)------------#  

   exec("m=open('%s','w')"%self.trancFileName_py_to_h5)

   exec("print 'Now printing the file for OpenMOC materials %s \
   (transport corrected)'%self.trancFileName_py_to_h5")

   m.write('import h5py\n')
   m.write('import numpy\n\n')
   m.write('# Create the file to store multi-groups cross-sections\n')
   m.write("f = h5py.File('%s')\n"%self.trancFileName_h5)
   m.write("f.attrs['Energy Groups'] = %s\n"%GROUPS)

   for mix in range(len(self.mixtures)):
 
     m.write('##########################################################\n')
     m.write('######################      %s      #######################\n' \
     %self.mixtures_h5[mix])
     m.write('###########################################################\n\n')

     m.write('# Create a subgroup for %s materials data\n'%self.mixtures[mix])
     m.write("%s=f.create_group('%s')\n\n"%(self.mixtures_h5[mix],self.mixtures[mix]))
     m.write("sigma_t=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_tot_tranc_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_tot_tranc_mac_%s[g])"%self.mixtures[mix])))

     m.write("sigma_a=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))  

     m.write("sigma_f=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))

     m.write("nu_sigma_f=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))   

     m.write("chi=numpy.array([")
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f])\n\n"%(eval("(chi_final_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(chi_final_%s[g])"%self.mixtures[mix]))) 

     m.write("sigma_s=numpy.array([")
 
     for g in range(GROUPS):
      for g_aux in range(GROUPS):
       if g==GROUPS-1 and g_aux==GROUPS-1:
        m.write("%.8f])\n "%(eval("(col_%s_mac_SCAT_TRANC_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g))))
       else:
        m.write("%.8f, "%(eval("(col_%s_mac_SCAT_TRANC_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g)))) 
        
     m.write('\n# Create datasets for each cross-section type\n')
     m.write('\n# Create datasets for each cross-section type\n')
     m.write("%s.create_dataset('Total XS', data=sigma_t)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Absorption XS', data=sigma_a)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Scattering XS', data=sigma_s)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Fission XS', data=sigma_f)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Nu Fission XS', data=nu_sigma_f)\n"%self.mixtures_h5[mix])
     m.write("%s.create_dataset('Chi', data=chi)\n\n"%self.mixtures_h5[mix])
     m.write('# Close the hdf5 data file\n')
   m.write('f.close()\n')
   m.close()
   #-----------------------------------------------------------------------#  
  
   #------------CREATING THE ISOTROPIC FILE (Python FORMAT)----------------#

   exec("m=open('%s','w')"%self.isotFileName_py)

   exec("print 'Now printing the file for OpenMOC materials %s \
   (isotropic)'%self.isotFileName_py")

   m.write('dataset = {}\n')
   m.write("dataset['Energy Groups'] = %s\n"%GROUPS)
   m.write("dataset['Materials'] = {}\n")
   m.write("test_dataset = dataset['Materials']\n\n")

   for mix in range(len(self.mixtures)):
 
     m.write('##########################################################\n')
     m.write('######################      %s      #######################\n' \
     %self.mixtures[mix])
     m.write('###########################################################\n\n')

     m.write('# Create a subdictionary for %s materials data\n'%self.mixtures[mix])
     m.write("test_dataset['%s'] = {}\n\n"%self.mixtures[mix])
     m.write("test_dataset['%s']['Total XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_tot_isot_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_tot_isot_mac_%s[g])"%self.mixtures[mix])))

     m.write("test_dataset['%s']['Absorption XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))  

     m.write("test_dataset['%s']['Fission XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))

     m.write("test_dataset['%s']['Nu Fission XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))   

     m.write("test_dataset['%s']['Chi'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(chi_final_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(chi_final_%s[g])"%self.mixtures[mix]))) 

     m.write("test_dataset['%s']['Scattering XS'] =["%self.mixtures[mix])
 
     for g in range(GROUPS):
      for g_aux in range(GROUPS):
       if g==GROUPS-1 and g_aux==GROUPS-1:
        m.write("%.8f]\n "%(eval("(col_%s_mac_SCAT_ISOT_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g))))
       else:
        m.write("%.8f, "%(eval("(col_%s_mac_SCAT_ISOT_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g)))) 
   m.close()
   #----------------------------------------------------------------------#
   
   #-------CREATING THE TRANS. CORRECTED FILE (Python FORMAT)-------------#  

   exec("m=open('%s','w')"%self.trancFileName_py)

   exec("print 'Now printing the file for OpenMOC materials %s \
   (transport corrected)'%self.trancFileName_py")

   m.write('dataset = {}\n')
   m.write("dataset['Energy Groups'] = %s\n"%GROUPS)
   m.write("dataset['Materials'] = {}\n")
   m.write("test_dataset = dataset['Materials']\n\n")

   for mix in range(len(self.mixtures)):
 
     m.write('##########################################################\n')
     m.write('######################      %s      #######################\n' \
     %self.mixtures[mix])
     m.write('###########################################################\n\n')

     m.write('# Create a subdictionary for %s materials data\n'%self.mixtures[mix])
     m.write("test_dataset['%s'] = {}\n\n"%self.mixtures[mix])
     m.write("test_dataset['%s']['Total XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_tot_tranc_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_tot_tranc_mac_%s[g])"%self.mixtures[mix])))

     m.write("test_dataset['%s']['Absorption XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))  

     m.write("test_dataset['%s']['Fission XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))

     m.write("test_dataset['%s']['Nu Fission XS'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))   

     m.write("test_dataset['%s']['Chi'] =["%self.mixtures[mix])
     for g in range(GROUPS):
      if g==GROUPS-1:
       m.write("%.8f]\n\n"%(eval("(chi_final_%s[g])"%self.mixtures[mix])))
      else:
       m.write("%.8f, "%(eval("(chi_final_%s[g])"%self.mixtures[mix]))) 

     m.write("test_dataset['%s']['Scattering XS'] =["%self.mixtures[mix])
 
     for g in range(GROUPS):
      for g_aux in range(GROUPS):
       if g==GROUPS-1 and g_aux==GROUPS-1:
        m.write("%.8f]\n "%(eval("(col_%s_mac_SCAT_TRANC_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g))))
       else:
        m.write("%.8f, "%(eval("(col_%s_mac_SCAT_TRANC_MTX_%s[%s])"
        %(g_aux,self.mixtures[mix],g)))) 
   m.close()
   #-----------------------------------------------------------------------#
  
   #---------------------CREATING THE MATLAB FILE--------------------------#
   
   exec("m=open('%s','w')"%self.matlabFileName)
   exec("print 'Now printing the file for Matlab/Octave processing %s  ' \
   %self.matlabFileName")
   m.write("%This file contains the different macroscopic XS per mixture that have been processed from DRAGON\n")
   m.write("%As previously warned in the python script, the transport correction is as follows: NTOT_TRANC = NTOT0 - SIGS01 and SCAT_MTX_TRANC = SCAT00-SCAT01\n")
   m.write("%This means that each term of the P0 matrix has been substracted to its correspondent P1 matrix term \n")
   m.write("%This type of transport correction might not be correct and may be disregarded\n")
   m.write("%ONLY THE ISOTROPIC MACROSCOPIC CROSS-SECTIONS PER MIXTURE HAVE BEEN VALIDATED\n")
   m.write("%Nevertheless, since this file contains P0 and P1 matrices of the different mixtures, you can use them to build your own transport corrected XS\n")
   m.write("%You can modify the end of this file in order to add more figures and/or process the XS\n\n")
   m.write("%Augusto Hernandez-Solis, October 2014 @ KTH, Sweden\n\n")
   m.write('clear all\n\n')
   m.write('clc \n\n')
   m.write('close all \n\n') 
   m.write("disp('You can modify the end of the file %s in order to add more figures and/or process the XS')\n"%self.matlabFileName)
   m.write("disp('Type help %s for more info')\n\n"%self.matlabFileName)

   for mix in range(len(self.mixtures)):

    m.write("NTOT0_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(n_tot_isot_mac_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(n_tot_isot_mac_%s[g])"%self.mixtures[mix])))

    m.write("NTOT_TRANC_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(n_tot_tranc_mac_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(n_tot_tranc_mac_%s[g])"%self.mixtures[mix])))

    m.write("ABS_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(n_abs_mac_%s[g])"%self.mixtures[mix])))  

    m.write("NFISTOT_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(n_fis_mac_%s[g])"%self.mixtures[mix])))

    m.write("NUFIS_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(n_nufis_mac_%s[g])"%self.mixtures[mix])))   

    m.write("CHI_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(chi_final_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(chi_final_%s[g])"%self.mixtures[mix]))) 

    m.write("SCAT00_MTX_vec_%s =["%self.mixtures[mix])
 
    for g in range(GROUPS):
     for g_aux in range(GROUPS):
      if g==GROUPS-1 and g_aux==GROUPS-1:
       m.write("%.8f];\n "%(eval("(col_%s_mac_SCAT_ISOT_MTX_%s[%s])"
       %(g_aux,self.mixtures[mix],g))))
      else:
       m.write("%.8f, "%(eval("(col_%s_mac_SCAT_ISOT_MTX_%s[%s])"
       %(g_aux,self.mixtures[mix],g))))

    m.write("SCAT01_MTX_vec_%s =["%self.mixtures[mix])
 
    for g in range(GROUPS):
     for g_aux in range(GROUPS):
      if g==GROUPS-1 and g_aux==GROUPS-1:
       m.write("%.8f];\n "%(eval("(col_%s_mac_SCAT_P1_MTX_%s[%s])"
       %(g_aux,self.mixtures[mix],g))))
      else:
       m.write("%.8f, "%(eval("(col_%s_mac_SCAT_P1_MTX_%s[%s])"
       %(g_aux,self.mixtures[mix],g))))

    m.write("SIGS00_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(n_scat_isot_mac_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(n_scat_isot_mac_%s[g])"%self.mixtures[mix])))

    m.write("SIGS_TRANC_%s =["%self.mixtures[mix])
    for g in range(GROUPS):
     if g==GROUPS-1:
      m.write("%.8f];\n\n"%(eval("(n_scat_tranc_mac_%s[g])"%self.mixtures[mix])))
     else:
      m.write("%.8f, "%(eval("(n_scat_tranc_mac_%s[g])"%self.mixtures[mix])))

    m.write("SIGS01_%s=-1*(SIGS_TRANC_%s-SIGS00_%s);\n\n"
    %(self.mixtures[mix],self.mixtures[mix],self.mixtures[mix]))

    m.write("for i=1:%s\n\n"%GROUPS)

    m.write("    SCAT00_MTX_%s(i,:)=SCAT00_MTX_vec_%s(1+%s*(i-1):%s*(i));\n"
    %(self.mixtures[mix],self.mixtures[mix],GROUPS,GROUPS))
    m.write("    SCAT01_MTX_%s(i,:)=SCAT01_MTX_vec_%s(1+%s*(i-1):%s*(i));\n"
    %(self.mixtures[mix],self.mixtures[mix],GROUPS,GROUPS))
 
    m.write("end\n\n")
   # 
   m.write("grp_172=[1.96403E+07 1.73325E+07 1.49182E+07 1.38403E+07 1.16183E+07 1.00000E+07 8.18731E+06 6.70320E+06 6.06531E+06 5.48812E+06 ...\n")
   m.write("4.49329E+06 3.67879E+06 3.01194E+06 2.46597E+06 2.23130E+06 2.01897E+06 1.65299E+06 1.35335E+06 1.22456E+06 1.10803E+06 ...\n")
   m.write("1.00259E+06 9.07180E+05 8.20850E+05 6.08101E+05 5.50232E+05 4.97871E+05 4.50492E+05 4.07622E+05 3.01974E+05 2.73237E+05 ...\n")
   m.write("2.47235E+05 1.83156E+05 1.22773E+05 1.11090E+05 8.22975E+04 6.73795E+04 5.51656E+04 4.08677E+04 3.69786E+04 2.92830E+04 ...\n")
   m.write("2.73944E+04 2.47875E+04 1.66156E+04 1.50344E+04 1.11378E+04 9.11882E+03 7.46586E+03 5.53085E+03 5.00450E+03 3.52662E+03 ...\n")
   m.write("3.35463E+03 2.24867E+03 2.03468E+03 1.50733E+03 1.43382E+03 1.23410E+03 1.01039E+03 9.14242E+02 7.48518E+02 6.77287E+02 ...\n")
   m.write("4.53999E+02 3.71703E+02 3.04325E+02 2.03995E+02 1.48625E+02 1.36742E+02 9.16609E+01 7.56736E+01 6.79040E+01 5.55951E+01 ...\n")
   m.write("5.15780E+01 4.82516E+01 4.55174E+01 4.01690E+01 3.72665E+01 3.37201E+01 3.05113E+01 2.76077E+01 2.49805E+01 2.26033E+01 ...\n")
   m.write("1.94548E+01 1.59283E+01 1.37096E+01 1.12245E+01 9.90555E+00 9.18981E+00 8.31529E+00 7.52398E+00 6.16012E+00 5.34643E+00 ...\n")
   m.write("5.04348E+00 4.12925E+00 4.00000E+00 3.38075E+00 3.30000E+00 2.76792E+00 2.72000E+00 2.60000E+00 2.55000E+00 2.36000E+00 ...\n")
   m.write("2.13000E+00 2.10000E+00 2.02000E+00 1.93000E+00 1.84000E+00 1.75500E+00 1.67000E+00 1.59000E+00 1.50000E+00 1.47500E+00 ...\n")
   m.write("1.44498E+00 1.37000E+00 1.33750E+00 1.30000E+00 1.23500E+00 1.17000E+00 1.15000E+00 1.12535E+00 1.11000E+00 1.09700E+00 ...\n")
   m.write("1.07100E+00 1.04500E+00 1.03500E+00 1.02000E+00 9.96000E-01 9.86000E-01 9.72000E-01 9.50000E-01 9.30000E-01 9.10000E-01 ...\n")
   m.write("8.60000E-01 8.50000E-01 7.90000E-01 7.80000E-01 7.05000E-01 6.25000E-01 5.40000E-01 5.00000E-01 4.85000E-01 4.33000E-01 ...\n")
   m.write("4.00000E-01 3.91000E-01 3.50000E-01 3.20000E-01 3.14500E-01 3.00000E-01 2.80000E-01 2.48000E-01 2.20000E-01 1.89000E-01 ...\n")
   m.write("1.80000E-01 1.60000E-01 1.40000E-01 1.34000E-01 1.15000E-01 1.00000E-01 9.50000E-02 8.00000E-02 7.70000E-02 6.70000E-02 ...\n")
   m.write("5.80000E-02 5.00000E-02 4.20000E-02 3.50000E-02 3.00000E-02 2.50000E-02 2.00000E-02 1.50000E-02 1.00000E-02 6.90000E-03 ...\n")
   m.write("5.00000E-03 3.00000E-03 1.00000E-05];\n\n")
 
   m.write("grp_69=[1.00000E+07 6.06550E+06 3.67900E+06 2.23100E+06 1.35300E+06 8.21000E+05 5.00000E+05 3.02500E+05 1.83000E+05 ...\n")
   m.write("1.11000E+05 6.73400E+04 4.08500E+04 2.47800E+04 1.50300E+04 9.11800E+03 5.53000E+03 3.51910E+03 2.23945E+03 ...\n")
   m.write("1.42510E+03 9.06899E+02 3.67263E+02 1.48729E+02 7.55014E+01 4.80520E+01 2.77000E+01 1.59680E+01 9.87700E+00 ...\n")
   m.write("4.00000E+00 3.30000E+00 2.60000E+00 2.10000E+00 1.50000E+00 1.30000E+00 1.15000E+00 1.12300E+00 1.09700E+00 ...\n")
   m.write("1.07100E+00 1.04500E+00 1.02000E+00 9.96000E-01 9.72000E-01 9.50000E-01 9.10000E-01 8.50000E-01 7.80000E-01 ...\n")
   m.write("6.25000E-01 5.00000E-01 4.00000E-01 3.50000E-01 3.20000E-01 3.00000E-01 2.80000E-01 2.50000E-01 2.20000E-01 ...\n")
   m.write("1.80000E-01 1.40000E-01 1.00000E-01 8.00000E-02 6.70000E-02 5.80000E-02 5.00000E-02 4.20000E-02 3.50000E-02 ...\n")
   m.write("3.00000E-02 2.50000E-02 2.00000E-02 1.50000E-02 1.00000E-02 5.00000E-03 1.00000E-05];\n\n")

   m.write("figure(1)\n")
   m.write("loglog(grp_%s(1:%s),NFISTOT_%s,'LineWidth',3)\n"%(GROUPS,GROUPS,self.mixtures[mix]))
   m.write("title('NFISTOT Fuel XS vs. Energy')\n")
   m.write("xlabel('Energy (eV)')\n")
   m.write("ylabel('NF (1/cm)')\n")
 
   m.close()
   #----------------------------------------------------------------------#
  #
 #  
#


     
        
     
   
