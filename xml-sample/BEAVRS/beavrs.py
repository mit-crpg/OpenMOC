#!/usr/bin/env python

"""
PWR OpenMC Model Generator

Allows for tweaking of specifications for the PWR OpenMC model, producing:
  geometry.xml
  materials.xml
  settings.xml
  plot.xml
  tallies.xml

"""

from __future__ import division

import copy
import math

from templates import *

############## Material paramters ##############

h2oDens = 0.73986
nominalBoronPPM = 975

############## Geometry paramters ##############

## 2-D or 3-D core
core_D = "2-D"
twoDlower = 220.0
twoDhigher = 230.0

## control rod insertion
controlStep = 0                       # bite position is D bank at 213 steps withdraw (228-213=16)

## pincell parameters
pelletOR        = 0.39218  # 
cladIR          = 0.40005  # 
cladOR          = 0.45720  # 
rodGridSide_tb  = 1.24416  #
rodGridSide_i   = 1.21962  #
guideTubeIR     = 0.56134  # 
guideTubeOR     = 0.60198  # 
guideTubeDashIR = 0.50419  # 
guideTubeDashOR = 0.54610  #
controlPoisonOR = 0.43310  #
controlRodIR    = 0.43688  # this is actually the CR clad only
controlRodOR    = 0.48387  # this is actually the CR clad only
burnabs1        = 0.21400  # 
burnabs2        = 0.23051  # 
burnabs3        = 0.24130  # 
burnabs4        = 0.42672  # 
burnabs5        = 0.43688  # 
burnabs6        = 0.48387  # 
burnabs7        = 0.56134  # 
burnabs8        = 0.60198  # 
instrTubeIR     = 0.43688  # no source                           
instrTubeOR     = 0.48387  # no source                           
plenumSpringOR  = 0.06459  # frapcon source - see beavrs_spring.ods

## lattice parameters
pinPitch        = 1.25984  # 
latticePitch    = 21.50364 #
gridstrapSide   = 21.47270 # 

## <!--- axials copied straight from axials_beavrs.ods

## axial paramters    
lowestExtent         =      0.000  # 
bottomSupportPlate   =     20.000  # arbitrary amount of water below core
topSupportPlate      =     25.000  # guessed
bottomLowerNozzle    =     25.000  # same as topSupportPlate
topLowerNozzle       =     35.160  # approx from ******* NDR of 4.088in for lower nozzle height
bottomFuelRod        =     35.160  # same as topLowerNozzle
topLowerThimble      =     36.007  # approx as 1/3 of inch, this is exact ******** NDR value for bottom thimble
bottomFuelStack      =     36.007  # same as topLowerThimble
activeCoreHeight     =    365.760  # provided by D***
topActiveCore        =    401.767  # bottomFuelStack + activeCoreHeight
botBurnAbs           =     41.087  # approx from ***** NDR of 1.987in for space between bot of BAs and bot of active fuel
   
# from Watts Bar Unit 2 safety analysis report   
topFuelRod           =    423.272  # fig 4.2-3 from Watts Bar Unit 2 safety analysis report
topPlenum            =    421.223  # fig 4.2-3 from Watts Bar Unit 2 safety analysis report
bottomUpperNozzle    =    426.617  # fig 4.2-3 from Watts Bar Unit 2 safety analysis report
topUpperNozzle       =    435.444  # fig 4.2-3 from Watts Bar Unit 2 safety analysis report
                       
highestExtent        =    455.444  # arbitrary amount of water above core
   
# grid z planes (centers provided by D***)  (heights 1.65 top/bottom, 2.25 intermediate)   
grid1Center          =     39.974  # bottomFuelStack + 1.562in
grid1bot             =     37.879  # grid1Center - 1.65/2
grid1top             =     42.070  # grid1Center + 1.65/2
grid2Center          =    102.021  # bottomFuelStack + 25.990in
grid2bot             =     99.164  # grid2Center - 2.25/2
grid2top             =    104.879  # grid2Center + 2.25/2
grid3Center          =    154.218  # bottomFuelStack + 46.540in
grid3bot             =    151.361  # grid3Center - 2.25/2
grid3top             =    157.076  # grid3Center + 2.25/2
grid4Center          =    206.415  # bottomFuelStack + 67.090in
grid4bot             =    203.558  # grid4Center - 2.25/2
grid4top             =    209.273  # grid4Center + 2.25/2
grid5Center          =    258.612  # bottomFuelStack + 87.640in
grid5bot             =    255.755  # grid5Center - 2.25/2
grid5top             =    261.470  # grid5Center + 2.25/2
grid6Center          =    310.809  # bottomFuelStack + 108.190in
grid6bot             =    307.952  # grid6Center - 2.25/2
grid6top             =    313.667  # grid6Center + 2.25/2
grid7Center          =    363.006  # bottomFuelStack + 128.740in
grid7bot             =    360.149  # grid7Center - 2.25/2
grid7top             =    365.864  # grid7Center + 2.25/2
grid8Center          =    414.624  # bottomFuelStack + 149.062in
grid8bot             =    412.529  # grid8Center - 1.65/2
grid8top             =    416.720  # grid8Center + 1.65/2
   
# control rod step heights   
step0H               =     45.079  # chosen to match the step size calculated for intervals between other grids
step36H              =    102.021  # grid2Center
step69H              =    154.218  # grid3Center
step102H             =    206.415  # grid4Center
step135H             =    258.612  # grid5Center
step168H             =    310.809  # grid6Center
step201H             =    363.006  # grid7Center
step228H             =    405.713  # set using calculated step width (27*stepWidth + step201H)
stepWidth            =    1.58173  # calculated from grid center planes

## -- axials copied straight from axials_beavrs.ods -->

## radial paramters
coreBarrelIR       = 187.9600 # Fig 2.1 from NDR 2-8              
coreBarrelOR       = 193.6750 # Fig 2.1 from NDR 2-8              
baffleWidth        =   2.2225 # Fig 2.1 from NDR 2-8              
rpvIR              = 230.0
rpvOR              = 251.9

neutronShieldOR    = 199.3900 # coreBarrelOR+2 find real number    TODO
neutronShield_NWbot_SEtop = "1 {0} 0 0".format(math.tan(math.pi/3))
neutronShield_NWtop_SEbot = "1 {0} 0 0".format(math.tan(math.pi/6))
neutronShield_NEbot_SWtop = "1 {0} 0 0".format(math.tan(-math.pi/3))
neutronShield_NEtop_SWbot = "1 {0} 0 0".format(math.tan(-math.pi/6))

def init_data():
  """All model parameters set here

  Materials, surfaces, cells, and lattices are defined here and automatically written
  to the proper files later.  Dictionary keys need to match those in the templates.

  The order each item is written is the order of appearance as written below

  Notes about core construction:
    The entire axial extent for each pincell is constructed in universes, and then added to fuel assembly lattices
    The fuel assembly lattices are then added to one master core lattice


  """

  matIDs  = set()
  surfIDs = set()
  cellIDs = set()
  univIDs = set()
  plotIDs = set()

  ################## materials ##################

  mo = {'n': 0}
  mats = {}

  from boron_ppm import Water

  h2o = Water()

  borwatmats = h2o.get_water_mats(h2oDens,nominalBoronPPM)
# h2o.print_latex()

  mats['water-nominal'] =   { 'order':    inc_order(mo),
                              'comment':  "Coolant (nominal)",
                              'id':       new_id(matIDs),
                              'density':  '{den:8.6f}'.format(den=h2o.rhoBW), 'unit': 'g/cc',
                              'nuclides': [ {'n':'B-10',   'xs':'71c','woao':'ao','a':h2o.NB10},
                                            {'n':'B-11',   'xs':'71c','woao':'ao','a':h2o.NB11},
                                            {'n':'H-1',    'xs':'71c','woao':'ao','a':h2o.NH1},
                                            {'n':'H-2',    'xs':'71c','woao':'ao','a':h2o.NH2},
                                            {'n':'O-16',   'xs':'71c','woao':'ao','a':h2o.NO16},
                                            {'n':'O-17',   'xs':'71c','woao':'ao','a':h2o.NO17+h2o.NO18},],
                              'sab':      [ {'name': 'lwtr', 'xs': '15t'}]}

  mats['helium'] =          { 'order':    inc_order(mo),
                              'comment':  "Helium for gap",
                              'id':       new_id(matIDs),
                              'density':  0.001598, 'unit': 'g/cc',
                              'nuclides': [ {'n':'He-4',  'xs':'71c','woao':'ao','a':2.4044e-04},]}

  mats['air'] =             { 'order':    inc_order(mo),
                              'comment':  "Air for instrument tubes",
                              'id':       new_id(matIDs),
                              'density':  0.006160, 'unit': 'g/cc',
                              'nuclides': [ {'n':'C-Nat',  'xs':'71c','woao':'ao','a':6.8296E-09},
                                            {'n':'O-16',  'xs':'71c','woao':'ao','a':5.2864e-06},
                                            {'n':'O-17',  'xs':'71c','woao':'ao','a':1.2877E-08},
                                            {'n':'N-14',  'xs':'71c','woao':'ao','a':1.9681e-05},
                                            {'n':'N-15',  'xs':'71c','woao':'ao','a':7.1900e-08},
                                            {'n':'Ar-36', 'xs':'71c','woao':'ao','a':7.9414e-10},
                                            {'n':'Ar-38', 'xs':'71c','woao':'ao','a':1.4915e-10},
                                            {'n':'Ar-40', 'xs':'71c','woao':'ao','a':2.3506e-07},]}

  mats['inconel'] =         { 'order':    inc_order(mo),
                              'comment':  "Inconel 718",
                              'id':       new_id(matIDs),
                              'density':  8.2, 'unit': 'g/cc',
                              'nuclides': [ {'n':'Si-28','xs':'71c','woao':'ao','a':5.6753E-04},
                                            {'n':'Si-29','xs':'71c','woao':'ao','a':2.8831E-05},
                                            {'n':'Si-30','xs':'71c','woao':'ao','a':1.9028E-05},
                                            {'n':'Cr-50','xs':'71c','woao':'ao','a':7.8239E-04},
                                            {'n':'Cr-52','xs':'71c','woao':'ao','a':1.5088E-02},
                                            {'n':'Cr-53','xs':'71c','woao':'ao','a':1.7108E-03},
                                            {'n':'Cr-54','xs':'71c','woao':'ao','a':4.2586E-04},
                                            {'n':'Mn-55','xs':'71c','woao':'ao','a':7.8201E-04},
                                            {'n':'Fe-54','xs':'71c','woao':'ao','a':1.4797E-03},
                                            {'n':'Fe-56','xs':'71c','woao':'ao','a':2.3229E-02},
                                            {'n':'Fe-57','xs':'71c','woao':'ao','a':5.3645E-04},
                                            {'n':'Fe-58','xs':'71c','woao':'ao','a':7.1392E-05},
                                            {'n':'Ni-58','xs':'71c','woao':'ao','a':2.9320E-02},
                                            {'n':'Ni-60','xs':'71c','woao':'ao','a':1.1294E-02},
                                            {'n':'Ni-61','xs':'71c','woao':'ao','a':4.9094E-04},
                                            {'n':'Ni-62','xs':'71c','woao':'ao','a':1.5653E-03},
                                            {'n':'Ni-64','xs':'71c','woao':'ao','a':3.9864E-04},]}

  mats['SS304'] =           { 'order':    inc_order(mo),
                              'comment':  "Stainless Steel 304",
                              'id':       new_id(matIDs),
                              'density':  8.03, 'unit': 'g/cc',
                              'nuclides': [ {'n':'Si-28','xs':'71c','woao':'ao','a':9.5274E-04},
                                            {'n':'Si-29','xs':'71c','woao':'ao','a':4.8400E-05},
                                            {'n':'Si-30','xs':'71c','woao':'ao','a':3.1943E-05},
                                            {'n':'Cr-50','xs':'71c','woao':'ao','a':7.6778E-04},
                                            {'n':'Cr-52','xs':'71c','woao':'ao','a':1.4806E-02},
                                            {'n':'Cr-53','xs':'71c','woao':'ao','a':1.6789E-03},
                                            {'n':'Cr-54','xs':'71c','woao':'ao','a':4.1791E-04},
                                            {'n':'Mn-55','xs':'71c','woao':'ao','a':1.7604E-03},
                                            {'n':'Fe-54','xs':'71c','woao':'ao','a':3.4620E-03},
                                            {'n':'Fe-56','xs':'71c','woao':'ao','a':5.4345E-02},
                                            {'n':'Fe-57','xs':'71c','woao':'ao','a':1.2551E-03},
                                            {'n':'Fe-58','xs':'71c','woao':'ao','a':1.6703E-04},
                                            {'n':'Ni-58','xs':'71c','woao':'ao','a':5.6089E-03},
                                            {'n':'Ni-60','xs':'71c','woao':'ao','a':2.1605E-03},
                                            {'n':'Ni-61','xs':'71c','woao':'ao','a':9.3917E-05},
                                            {'n':'Ni-62','xs':'71c','woao':'ao','a':2.9945E-04},
                                            {'n':'Ni-64','xs':'71c','woao':'ao','a':7.6261E-05},]}
                                            
  mats['carbon steel'] =    { 'order':    inc_order(mo),
                              'comment':  "Carbon Steel ASTM A533 Grade B",
                              'id':       new_id(matIDs),
                              'density':  7.8, 'unit': 'g/cc',
                              'nuclides': [ {'n':'C-Nat','xs':'71c','woao':'ao','a':9.7772E-04},
                                            {'n':'Si-28','xs':'71c','woao':'ao','a':4.2417E-04},
                                            {'n':'Si-29','xs':'71c','woao':'ao','a':2.1548E-05},
                                            {'n':'Si-30','xs':'71c','woao':'ao','a':1.4221E-05},
                                            {'n':'Mn-55','xs':'71c','woao':'ao','a':1.1329E-03},
                                            {'n':'P-31','xs':'71c','woao':'ao','a':3.7913E-05},
                                            {'n':'Mo-92','xs':'71c','woao':'ao','a':3.7965E-05},
                                            {'n':'Mo-94','xs':'71c','woao':'ao','a':2.3725E-05},
                                            {'n':'Mo-96','xs':'71c','woao':'ao','a':4.2875E-05},
                                            {'n':'Mo-97','xs':'71c','woao':'ao','a':2.4573E-05},
                                            {'n':'Mo-98','xs':'71c','woao':'ao','a':6.2179E-05},
                                            {'n':'Mo-100','xs':'71c','woao':'ao','a':2.4856E-05},
                                            {'n':'Fe-54','xs':'71c','woao':'ao','a':4.7714E-03},
                                            {'n':'Fe-56','xs':'71c','woao':'ao','a':7.4900E-02},
                                            {'n':'Fe-57','xs':'71c','woao':'ao','a':1.7298E-03},
                                            {'n':'Fe-58','xs':'71c','woao':'ao','a':2.3020E-04},
                                            {'n':'Ni-58','xs':'71c','woao':'ao','a':2.9965E-04},
                                            {'n':'Ni-60','xs':'71c','woao':'ao','a':1.1543E-04},
                                            {'n':'Ni-61','xs':'71c','woao':'ao','a':5.0175E-06},
                                            {'n':'Ni-62','xs':'71c','woao':'ao','a':1.5998E-05},
                                            {'n':'Ni-64','xs':'71c','woao':'ao','a':4.0742E-06},]}

  mats['zirc'] =            { 'order':    inc_order(mo),
                              'comment':  "Zircaloy-4",
                              'id':       new_id(matIDs),
                              'density':  6.55, 'unit': 'g/cc',
                              'nuclides': [ {'n':'O-16','xs':'71c','woao':'ao','a':3.0743E-04},
                                            {'n':'O-17','xs':'71c','woao':'ao','a':7.4887E-07},
                                            {'n':'Cr-50','xs':'71c','woao':'ao','a':3.2962E-06},
                                            {'n':'Cr-52','xs':'71c','woao':'ao','a':6.3564E-05},
                                            {'n':'Cr-53','xs':'71c','woao':'ao','a':7.2076E-06},
                                            {'n':'Cr-54','xs':'71c','woao':'ao','a':1.7941E-06},
                                            {'n':'Fe-54','xs':'71c','woao':'ao','a':8.6699E-06},
                                            {'n':'Fe-56','xs':'71c','woao':'ao','a':1.3610E-04},
                                            {'n':'Fe-57','xs':'71c','woao':'ao','a':3.1431E-06},
                                            {'n':'Fe-58','xs':'71c','woao':'ao','a':4.1829E-07},
                                            {'n':'Zr-90','xs':'71c','woao':'ao','a':2.1827E-02},
                                            {'n':'Zr-91','xs':'71c','woao':'ao','a':4.7600E-03},
                                            {'n':'Zr-92','xs':'71c','woao':'ao','a':7.2758E-03},
                                            {'n':'Zr-94','xs':'71c','woao':'ao','a':7.3734E-03},
                                            {'n':'Zr-96','xs':'71c','woao':'ao','a':1.1879E-03},
                                            {'n':'Sn-112','xs':'71c','woao':'ao','a':4.6735E-06},
                                            {'n':'Sn-114','xs':'71c','woao':'ao','a':3.1799E-06},
                                            {'n':'Sn-115','xs':'71c','woao':'ao','a':1.6381E-06},
                                            {'n':'Sn-116','xs':'71c','woao':'ao','a':7.0055E-05},
                                            {'n':'Sn-117','xs':'71c','woao':'ao','a':3.7003E-05},
                                            {'n':'Sn-118','xs':'71c','woao':'ao','a':1.1669E-04},
                                            {'n':'Sn-119','xs':'71c','woao':'ao','a':4.1387E-05},
                                            {'n':'Sn-120','xs':'71c','woao':'ao','a':1.5697E-04},
                                            {'n':'Sn-122','xs':'71c','woao':'ao','a':2.2308E-05},
                                            {'n':'Sn-124','xs':'71c','woao':'ao','a':2.7897E-05},]}

  mats['UO2 1.6'] =         { 'order':    inc_order(mo),
                              'comment':  "UO2 Fuel 1.6 w/o",
                              'id':       new_id(matIDs),
                              'density':  10.31362, 'unit': 'g/cc',
                              'nuclides': [ {'n':'U-234','xs':'71c','woao':'ao','a':3.0131E-06},
                                            {'n':'U-235','xs':'71c','woao':'ao','a':3.7503E-04},
                                            {'n':'U-238','xs':'71c','woao':'ao','a':2.2626E-02},
                                            {'n':'O-16','xs':'71c','woao':'ao','a':4.5896E-02},
                                            {'n':'O-17','xs':'71c','woao':'ao','a':1.1180E-04},]}

  mats['UO2 2.4'] =         { 'order':    inc_order(mo),
                              'comment':  "UO2 Fuel 2.4 w/o",
                              'id':       new_id(matIDs),
                              'density':  10.29769, 'unit': 'g/cc',
                              'nuclides': [ {'n':'U-234','xs':'71c','woao':'ao','a':4.4843e-06},
                                            {'n':'U-235','xs':'71c','woao':'ao','a':5.5815e-04},
                                            {'n':'U-238','xs':'71c','woao':'ao','a':2.2408e-02},
                                            {'n':'O-16','xs':'71c','woao':'ao','a':4.5829e-02},
                                            {'n':'O-17','xs':'71c','woao':'ao','a':1.1164E-04},]}

  mats['UO2 3.1'] =         { 'order':    inc_order(mo),
                              'comment':  "UO2 Fuel 3.1 w/o",
                              'id':       new_id(matIDs),
                              'density':  10.30187, 'unit': 'g/cc',
                              'nuclides': [ {'n':'U-234','xs':'71c','woao':'ao','a':5.7988e-06},
                                            {'n':'U-235','xs':'71c','woao':'ao','a':7.2176e-04},
                                            {'n':'U-238','xs':'71c','woao':'ao','a':2.2254e-02},
                                            {'n':'O-16','xs':'71c','woao':'ao','a':4.5851e-02},
                                            {'n':'O-17','xs':'71c','woao':'ao','a':1.1169E-04},]}

  mats['control rod'] =     { 'order':    inc_order(mo),
                              'comment':  "Ag-In-Cd Control Rod",
                              'id':       new_id(matIDs),
                              'density':  10.16, 'unit': 'g/cc',
                              'nuclides': [ {'n':'Ag-107','xs':'71c','woao':'ao','a':2.3523E-02},
                                            {'n':'Ag-109','xs':'71c','woao':'ao','a':2.1854E-02},
                                            {'n':'In-113','xs':'71c','woao':'ao','a':3.4291E-04},
                                            {'n':'In-115','xs':'71c','woao':'ao','a':7.6504E-03},
                                            {'n':'Cd-106','xs':'71c','woao':'ao','a':3.4019E-05},
                                            {'n':'Cd-108','xs':'71c','woao':'ao','a':2.4221E-05},
                                            {'n':'Cd-110','xs':'71c','woao':'ao','a':3.3991E-04},
                                            {'n':'Cd-111','xs':'71c','woao':'ao','a':3.4835E-04},
                                            {'n':'Cd-112','xs':'71c','woao':'ao','a':6.5669E-04},
                                            {'n':'Cd-113','xs':'71c','woao':'ao','a':3.3257E-04},
                                            {'n':'Cd-114','xs':'71c','woao':'ao','a':7.8188E-04},
                                            {'n':'Cd-116','xs':'71c','woao':'ao','a':2.0384E-04},]}

  mats['borosilicate'] =    { 'order':    inc_order(mo),
                              'comment':  "Borosilicate Glass in BA rod",
                              'id':       new_id(matIDs),
                              'density':  2.26, 'unit': 'g/cc',
                              'nuclides': [ {'n':'B-10','xs':'71c','woao':'ao','a':9.6506E-04},
                                            {'n':'B-11','xs':'71c','woao':'ao','a':3.9189E-03},
                                            {'n':'O-16','xs':'71c','woao':'ao','a':4.6511E-02},
                                            {'n':'O-17','xs':'71c','woao':'ao','a':1.1330E-04},
                                            {'n':'Al-27','xs':'71c','woao':'ao','a':1.7352E-03},
                                            {'n':'Si-28','xs':'71c','woao':'ao','a':1.6924E-02},
                                            {'n':'Si-29','xs':'71c','woao':'ao','a':8.5977E-04},
                                            {'n':'Si-30','xs':'71c','woao':'ao','a':5.6743E-04},]}
  
  ################## surfaces ##################

  so = {'n': 0}
  surfs = {}
  surfs['dummy outer'] =        { 'order':   inc_order(so),
                                  'section': comm_t.format("Pincell surfaces"),
                                  'comm':    comm_t.format("dummy outer boundary"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 4000',}
  surfs['pellet OR'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("pellet OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(pelletOR)}
  surfs['plenum spring OR'] =   { 'order':   inc_order(so),
                                  'comm':    comm_t.format("fuel rod plenum spring OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(plenumSpringOR)}
  surfs['clad IR'] =            { 'order':   inc_order(so),
                                  'comm':    comm_t.format("clad IR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(cladIR)}
  surfs['clad OR'] =            { 'order':   inc_order(so),
                                  'comm':    comm_t.format("clad OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(cladOR)}
  surfs['guide tube IR'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("guide tube IR above dashpot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(guideTubeIR)}
  surfs['guide tube OR'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("guide tube OR above dashpot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(guideTubeOR)}
  surfs['GT dashpot IR'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("guide tube IR at dashpot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(guideTubeDashIR)}
  surfs['GT dashpot OR'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("guide tube OR at dashpot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(guideTubeDashOR)}
  surfs['control poison OR'] =  { 'order':   inc_order(so),
                                  'comm':    comm_t.format("control poison OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(controlPoisonOR)}
  surfs['control rod IR'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("control rod clad IR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(controlRodIR)}
  surfs['control rod OR'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("control rod clad OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(controlRodOR)}
  surfs['burnabs rad 1'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 1"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs1)}
  surfs['burnabs rad 2'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 2"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs2)}
  surfs['burnabs rad 3'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 3"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs3)}
  surfs['burnabs rad 4'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 4"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs4)}
  surfs['burnabs rad 5'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 5"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs5)}
  surfs['burnabs rad 6'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 6"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs6)}
  surfs['burnabs rad 7'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 7"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs7)}
  surfs['burnabs rad 8'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("burnable absorber rod inner radius 8"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(burnabs8)}
  surfs['instr tube IR'] =         copy.copy(surfs['burnabs rad 5'])
  surfs['instr tube IR']['dupe'] = True
  surfs['instr tube OR'] =         copy.copy(surfs['burnabs rad 6'])
  surfs['instr tube OR']['dupe'] = True 


  surfs['rod grid box xtop tb'] =            { 'order':   inc_order(so),
                                  'comm':    comm_t.format("X max for grid outside fuel rods in top/bottom spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(rodGridSide_tb/2)}
  surfs['rod grid box xbot tb'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("X min for grid outside fuel rods in top/bottom spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(rodGridSide_tb/2)}
  surfs['rod grid box ytop tb'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("Y max for grid outside fuel rods in top/bottom spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(rodGridSide_tb/2)}
  surfs['rod grid box ybot tb'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("Y min for grid outside fuel rods in top/bottom spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(rodGridSide_tb/2)}

  surfs['rod grid box xtop i'] =            { 'order':   inc_order(so),
                                  'comm':    comm_t.format("X max for grid outside fuel rods in intermediate spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(rodGridSide_i/2)}
  surfs['rod grid box xbot i'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("X min for grid outside fuel rods in intermediate spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(rodGridSide_i/2)}
  surfs['rod grid box ytop i'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("Y max for grid outside fuel rods in intermediate spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(rodGridSide_i/2)}
  surfs['rod grid box ybot i'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("Y min for grid outside fuel rods in intermediate spacers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(rodGridSide_i/2)}

  surfs['lat grid box xtop'] =  { 'order':   inc_order(so),
                                  'comm':    comm_t.format("X max for grid outside fuel assembly"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(gridstrapSide/2)}
  surfs['lat grid box xbot'] =  { 'order':   inc_order(so),
                                  'comm':    comm_t.format("X min for grid outside fuel assembly"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(gridstrapSide/2)}
  surfs['lat grid box ytop'] =  { 'order':   inc_order(so),
                                  'comm':    comm_t.format("Y max for grid outside fuel assembly"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(gridstrapSide/2)}
  surfs['lat grid box ybot'] =  { 'order':   inc_order(so),
                                  'comm':    comm_t.format("Y min for grid outside fuel assembly"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(gridstrapSide/2)}

   

  # lattice surfaces
  surfs['lat box xtop'] =            { 'order':   inc_order(so),
                                  'section': comm_t.format("Lattice surfaces"),
                                  'comm':    comm_t.format("lattice X max"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(17*pinPitch/2)}
  surfs['lat box xbot'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("lattice X min"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(17*pinPitch/2)}
  surfs['lat box ytop'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("lattice Y max"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(17*pinPitch/2)}
  surfs['lat box ybot'] =      { 'order':   inc_order(so),
                                  'comm':    comm_t.format("lattice Y min"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '-{0:0<8.6}'.format(17*pinPitch/2)}
  # axial surfaces
  surfs['lowest extent'] =      { 'order':   inc_order(so),
                                  'section': comm_t.format("Axial surfaces"),
                                  'comm':    comm_t.format("lowest extent"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(lowestExtent)}
  surfs['bot support plate'] =  { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bot support plate"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(bottomSupportPlate)}
  surfs['top support plate'] =  { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top support plate"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(topSupportPlate)}
  surfs['bot lower nozzle'] =   copy.copy(surfs['top support plate'])
  surfs['bot lower nozzle']['dupe'] = True
  surfs['bot fuel rod'] =       { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bot fuel rod"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(bottomFuelRod)}
  surfs['top lower nozzle'] =   copy.copy(surfs['bot fuel rod'])
  surfs['top lower nozzle']['dupe'] = True
  surfs['bot active core'] =    { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bot active core"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(bottomFuelStack)}
  surfs['top lower thimble'] =   copy.copy(surfs['bot active core'])
  surfs['top lower thimble']['dupe'] = True
  surfs['grid1bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 1"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid1bot)}
  surfs['burn abs bot'] =       { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom of burnable absorbers"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(botBurnAbs)}
  surfs['grid1top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 1"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid1top)}
  surfs['dashpot top'] =        { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top dashpot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(step0H)}
  surfs['grid2bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 2"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid2bot)}
  surfs['grid2top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 2"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid2top)}
  surfs['grid3bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 3"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid3bot)}
  surfs['grid3top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 3"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid3top)}
  surfs['grid4bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 4"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid4bot)}
  surfs['grid4top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 4"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid4top)}
  surfs['grid5bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 5"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid5bot)}
  surfs['grid5top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 5"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid5top)}
  surfs['grid6bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 6"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid6bot)}
  surfs['grid6top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 6"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid6top)}
  surfs['grid7bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 7"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid7bot)}
  surfs['grid7top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 7"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid7top)}
  surfs['top active core'] =    { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top active core"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(topActiveCore)}
  surfs['grid8bot'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom grid 8"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid8bot)}
  surfs['grid8top'] =           { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top grid 8"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(grid8top)}
  gridSurfaces = [surfs['grid1bot']['id'],
                  surfs['grid1top']['id'],
                  surfs['grid2bot']['id'],
                  surfs['grid2top']['id'],
                  surfs['grid3bot']['id'],
                  surfs['grid3top']['id'],
                  surfs['grid4bot']['id'],
                  surfs['grid4top']['id'],
                  surfs['grid5bot']['id'],
                  surfs['grid5top']['id'],
                  surfs['grid6bot']['id'],
                  surfs['grid6top']['id'],
                  surfs['grid7bot']['id'],
                  surfs['grid7top']['id'],
                  surfs['grid8bot']['id'],
                  surfs['grid8top']['id'],]
  surfs['top pin plenum'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top fuel rod"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(topPlenum)}
  surfs['top fuel rod'] =       { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top fuel rod"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(topFuelRod)}
  surfs['bot upper nozzle'] =   { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bottom upper nozzle"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(bottomUpperNozzle)}
  surfs['top upper nozzle'] =   { 'order':   inc_order(so),
                                  'comm':    comm_t.format("top upper nozzle"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(topUpperNozzle)}
  surfs['highest extent'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("highest extent"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(highestExtent)}

  # control rod bank surfaces
  controlAxials = control_bank_axials(controlStep)

  surfs['bankA top'] =          { 'order':   inc_order(so),
                                  'section': comm_t.format("Control rod bank surfaces"),
                                  'comm':    comm_t.format("bankA top"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankAtop'])}
  surfs['bankA bot'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bankA bot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankAbot'])}
  surfs['bankB top'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bankB top"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankBtop'])}
  surfs['bankB bot'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bankB bot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankBbot'])}
  surfs['bankC top'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bankC top"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankCtop'])}
  surfs['bankC bot'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bankC bot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankCbot'])}
  surfs['bankD top'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bankD top"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankDtop'])}
  surfs['bankD bot'] =          { 'order':   inc_order(so),
                                  'comm':    comm_t.format("bankD bot"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(controlAxials['bankDbot'])}
                                  
  for b in ['A','B','C','D','E',]:
    surfs['bankS{0} top'.format(b)] =       { 'order':   inc_order(so),
                                    'comm':    comm_t.format("bankS{0} top".format(b)),
                                    'id':      new_id(surfIDs),
                                    'type':    '"z-plane"',
                                    'coeffs':  '{0:0<8.6}'.format(step228H+stepWidth*228)}
    surfs['bankS{0} bot'.format(b)] =       { 'order':   inc_order(so),
                                    'comm':    comm_t.format("bankS{0} bot".format(b)),
                                    'id':      new_id(surfIDs),
                                    'type':    '"z-plane"',
                                    'coeffs':  '{0:0<8.6}'.format(step228H)}

  # outer radial surfaces
  surfs['core barrel IR'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("core barrel IR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(coreBarrelIR)}
  surfs['core barrel OR'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("core barrel OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(coreBarrelOR)}
  surfs['neut shield OR'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("neutron shield OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(neutronShieldOR)}
  surfs['neut shield NWbot SEtop'] =     { 'order':   inc_order(so),
                                  'comm':    comm_t.format("neutron shield planes"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"plane"',
                                  'coeffs':  neutronShield_NWbot_SEtop}
  surfs['neut shield NWtop SEbot'] =     { 'order':   inc_order(so),
                                  'comm':    "",
                                  'id':      new_id(surfIDs),
                                  'type':    '"plane"',
                                  'coeffs':  neutronShield_NWtop_SEbot}
  surfs['neut shield NEbot SWtop'] =     { 'order':   inc_order(so),
                                  'comm':    "",
                                  'id':      new_id(surfIDs),
                                  'type':    '"plane"',
                                  'coeffs':  neutronShield_NEbot_SWtop}
  surfs['neut shield NEtop SWbot'] =     { 'order':   inc_order(so),
                                  'comm':    "",
                                  'id':      new_id(surfIDs),
                                  'type':    '"plane"',
                                  'coeffs':  neutronShield_NEtop_SWbot}
  surfs['RPV IR'] =             { 'order':   inc_order(so),
                                  'comm':    comm_t.format("RPV IR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(rpvIR)}
  surfs['RPV OR'] =             { 'order':   inc_order(so),
                                  'comm':    comm_t.format("RPV OR"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"circle"',
                                  'coeffs':  '0.0 0.0 {0:0<8.6}'.format(rpvOR),
                                  'bc':      'vacuum'}
  surfs['upper bound'] =        { 'order':   inc_order(so),
                                  'comm':    comm_t.format("upper problem boundary"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(highestExtent),
                                  'bc':      'vacuum'}
  surfs['lower bound'] =        { 'order':   inc_order(so),
                                  'comm':    comm_t.format("lower problem boundary"),
                                  'id':      new_id(surfIDs),
                                  'type':    '"z-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(lowestExtent),
                                  'bc':      'vacuum'}
  if core_D == '2-D':
    surfs['upper bound'].update( {  'coeffs':  '{0:0<8.6}'.format(twoDhigher),
                                    'bc':      'reflective'} )
    surfs['lower bound'].update( {  'coeffs':  '{0:0<8.6}'.format(twoDlower),
                                    'bc':      'reflective'} )


  ################## cells ##################

  # (if not a fill cell, set 'fill': None.  if not None, 'mat' is ignored)

  co = {'n': 0}
  cells = {}
  cells['water pin'] =          { 'order':   inc_order(co),
                                  'section': comm_t.format("Empty water pincell universe"),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}
  cells['water pin 2'] =        { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['water pin']['univ'],
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['dummy outer']['id'])}

  # GUIDE TUBE PIN CELLS
  
  
  make_pin('GT empty',"Guide tube pincell universes",co,cells,new_id(univIDs),cellIDs,
            [surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id']],
            [(mats['water-nominal']['id'],"empty guide tube"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('GT empty grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['water-nominal']['id'],"empty guide tube with top/bottom grid"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('GT empty grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['water-nominal']['id'],"empty guide tube with intermediate grid"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')
  make_pin('GT empty nozzle',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['water-nominal']['id'],"empty guide tube SS304 penetration"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('GTd empty',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['GT dashpot IR']['id'],
             surfs['GT dashpot OR']['id'],],
            [(mats['water-nominal']['id'],"empty guide tube at dashpot"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('GTd empty grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['GT dashpot IR']['id'],
             surfs['GT dashpot OR']['id'],],
            [(mats['water-nominal']['id'],"empty guide tube at dashpot with top/bottom grid"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('GTd empty grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['GT dashpot IR']['id'],
             surfs['GT dashpot OR']['id'],],
            [(mats['water-nominal']['id'],"empty guide tube at dashpot with intermediate grid"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')
  make_pin('GTd empty nozzle',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['GT dashpot IR']['id'],
             surfs['GT dashpot OR']['id'],],
            [(mats['water-nominal']['id'],"empty GTd SS304 penetration"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])


  # final combination of all axial pieces for empty guide tube
  stackSurfs = [surfs['bot support plate']['id'],
                surfs['top support plate']['id'],
                surfs['top lower nozzle']['id'],
                surfs['top lower thimble']['id'],
                surfs['grid1bot']['id'],
                surfs['grid1top']['id'],
                surfs['dashpot top']['id'],
                surfs['grid2bot']['id'],
                surfs['grid2top']['id'],
                surfs['grid3bot']['id'],
                surfs['grid3top']['id'],
                surfs['grid4bot']['id'],
                surfs['grid4top']['id'],
                surfs['grid5bot']['id'],
                surfs['grid5top']['id'],
                surfs['grid6bot']['id'],
                surfs['grid6top']['id'],
                surfs['grid7bot']['id'],
                surfs['grid7top']['id'],
                surfs['top active core']['id'],
                surfs['grid8bot']['id'],
                surfs['grid8top']['id'],
                surfs['top pin plenum']['id'],
                surfs['top fuel rod']['id'],
                surfs['bot upper nozzle']['id'],
                surfs['top upper nozzle']['id']]
  make_stack('GT empty stack',co,cells,new_id(univIDs),cellIDs,
            stackSurfs,
            [cells['water pin']['univ'],          # lower plenum
#            cells['GTd empty nozzle']['univ'],   # support plate
#            cells['GTd empty nozzle']['univ'],   # lower nozzle 
             cells['water pin']['univ'],          # support plate
             cells['water pin']['univ'],          # lower nozzle
             cells['GTd empty']['univ'],          # lower thimble
             cells['GTd empty']['univ'],          # dashpot
             cells['GTd empty grid_tb']['univ'],  # dashpot grid 1
             cells['GTd empty']['univ'],          # dashpot
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],    # reg grid 2
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],    # reg grid 3
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],    # reg grid 4
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],    # reg grid 5
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],    # reg grid 6
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],    # reg grid 7
             cells['GT empty']['univ'],           # reg
             cells['GT empty']['univ'],           # pin plenum
             cells['GT empty grid_tb']['univ'],   # pin plenum grid 8
             cells['GT empty']['univ'],           # pin plenum
             cells['GT empty']['univ'],           # upper fuel rod end plug
             cells['GT empty']['univ'],           # space between fuel rod and and upper nozzle
#            cells['GT empty nozzle']['univ'],    # upper nozzle
             cells['water pin']['univ'],          # upper nozzle
             cells['water pin']['univ']])         # upper plenum
  make_stack('GT empty stack instr',co,cells,new_id(univIDs),cellIDs,
            stackSurfs,
            [cells['water pin']['univ'],          # lower plenum
#            cells['GT empty nozzle']['univ'],   # support plate
#            cells['GT empty nozzle']['univ'],   # lower nozzle
             cells['water pin']['univ'],          # support plate
             cells['water pin']['univ'],          # lower nozzle
             cells['GT empty']['univ'],          # lower thimble
             cells['GT empty']['univ'],          # dashpot
             cells['GT empty grid_tb']['univ'],     # dashpot grid 1
             cells['GT empty']['univ'],          # dashpot
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],      # reg grid 2
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],      # reg grid 3
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],      # reg grid 4
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],      # reg grid 5
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],      # reg grid 6
             cells['GT empty']['univ'],           # reg
             cells['GT empty grid_i']['univ'],      # reg grid 7
             cells['GT empty']['univ'],           # reg
             cells['GT empty']['univ'],           # pin plenum
             cells['GT empty grid_tb']['univ'],      # pin plenum grid 8
             cells['GT empty']['univ'],           # pin plenum
             cells['GT empty']['univ'],           # upper fuel rod end plug
             cells['GT empty']['univ'],           # space between fuel rod and and upper nozzle
             cells['water pin']['univ'],          # upper nozzle
#            cells['GT empty nozzle']['univ'],    # upper nozzle
             cells['water pin']['univ']])         # upper plenum

  # INSTRUMENT TUBE PIN CELL

  make_pin('GT instr',"Instrument tube pincell universes",co,cells,new_id(univIDs),cellIDs,
            [surfs['instr tube IR']['id'],
             surfs['instr tube OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id']],
            [(mats['air']['id'],"instr guide tube above dashpot"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('GT instr grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['instr tube IR']['id'],
             surfs['instr tube OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['air']['id'],"instr guide tube at dashpot with top/bottom grid"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('GT instr grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['instr tube IR']['id'],
             surfs['instr tube OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['air']['id'],"instr guide tube at dashpot with intermediate grid"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')
  make_pin('GT instr nozzle',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['instr tube IR']['id'],
             surfs['instr tube OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['air']['id'],"instr guide tube SS304 penetration"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['SS304']['id'],""),])
  make_pin('bare instr',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['instr tube IR']['id'],
             surfs['instr tube OR']['id'],],
            [(mats['air']['id'],"instr guide tube at dashpot"),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])


  # final combination of all axial pieces for instrument tube
  make_stack('GT instr stack',co,cells,new_id(univIDs),cellIDs,
            stackSurfs,
            [cells['bare instr']['univ'],        # lower plenum
#            cells['GT instr nozzle']['univ'],   # support plate
#            cells['GT instr nozzle']['univ'],   # lower nozzle
             cells['bare instr']['univ'],        # support plate
             cells['bare instr']['univ'],        # lower nozzle
             cells['GT instr']['univ'],          # lower thimble
             cells['GT instr']['univ'],          # dashpot
             cells['GT instr grid_tb']['univ'],  # dashpot grid 1
             cells['GT instr']['univ'],          # dashpot
             cells['GT instr']['univ'],          # reg
             cells['GT instr grid_i']['univ'],   # reg grid 2
             cells['GT instr']['univ'],          # reg
             cells['GT instr grid_i']['univ'],   # reg grid 3
             cells['GT instr']['univ'],          # reg
             cells['GT instr grid_i']['univ'],   # reg grid 4
             cells['GT instr']['univ'],          # reg
             cells['GT instr grid_i']['univ'],   # reg grid 5
             cells['GT instr']['univ'],          # reg
             cells['GT instr grid_i']['univ'],   # reg grid 6
             cells['GT instr']['univ'],          # reg
             cells['GT instr grid_i']['univ'],   # reg grid 7
             cells['GT instr']['univ'],          # reg
             cells['GT instr']['univ'],          # pin plenum
             cells['GT instr grid_tb']['univ'],  # pin plenum grid 8
             cells['GT instr']['univ'],          # pin plenum
             cells['GT instr']['univ'],          # upper fuel rod end plug
             cells['GT instr']['univ'],          # space between fuel rod and and upper nozzle
#            cells['GT instr nozzle']['univ'],   # upper nozzle
             cells['bare instr']['univ'],        # upper nozzle
             cells['water pin']['univ']])         # upper plenum


  # CONTROL ROD PIN CELLS

  make_pin('control rod',"Control rod bank pincell universes",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['control rod']['id'],"guide tube above dashpot with control rod"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('control rod grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['control rod']['id'],"guide tube above dashpot with control rod, with top/bottom grid"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('control rod grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['control rod']['id'],"guide tube above dashpot with control rod, with intermediate grid"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')
  make_pin('control rod nozzle',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],],
            [(mats['control rod']['id'],"guide tube above dashpot with control rod, with nozzle"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('control rod blank',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['SS304']['id'],"blank control rod above active region, above nozzle"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('control rod blank grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['SS304']['id'],"guide tube above dashpot with blank control rod, with top/bottom grid"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('control rod blank grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],
             surfs['guide tube IR']['id'],
             surfs['guide tube OR']['id'],],
            [(mats['SS304']['id'],"guide tube above dashpot with blank control rod, with intermediate grid"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')
  make_pin('control rod blank nozzle',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],],
            [(mats['SS304']['id'],"blank control rod above active region, within nozzle"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('control rod blank bare',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],],
            [(mats['SS304']['id'],"blank control rod above active region, within nozzle"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('control rod bare',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['control poison OR']['id'],
             surfs['control rod IR']['id'],
             surfs['control rod OR']['id'],],
            [(mats['control rod']['id'],"guide tube above dashpot with control rod"),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),])


  # for the control rods, we make the axial universe three times, once with the
  # grid everywhere, once without, and once with the nozzle everywhere. Then
  # we can alternate with them as appropriate in the final axial stack

  banks = ['A','B','C','D','SA','SB','SC','SD','SE']
  
  for b in banks:

    # no grid, no nozzle
    make_stack('dummy GT CR bank {0}'.format(b),co,cells,new_id(univIDs),cellIDs,
                [surfs['bot fuel rod']['id'],
                 surfs['dashpot top']['id'],
                 surfs['bank{0} bot'.format(b)]['id'],
                 surfs['bank{0} top'.format(b)]['id'],],
                [cells['water pin']['univ'],
                 cells['GTd empty']['univ'],
                 cells['GT empty']['univ'],
                 cells['control rod']['univ'],
                 cells['control rod blank']['univ'],])
                 
    # top/bottom grid
    make_stack('dummy GT CR bank {0} grid_tb'.format(b),co,cells,new_id(univIDs),cellIDs,
                [surfs['bot fuel rod']['id'],
                 surfs['dashpot top']['id'],
                 surfs['bank{0} bot'.format(b)]['id'],
                 surfs['bank{0} top'.format(b)]['id'],],
                [cells['water pin']['univ'],
                 cells['GTd empty grid_tb']['univ'],
                 cells['GT empty grid_tb']['univ'],
                 cells['control rod grid_tb']['univ'],
                 cells['control rod blank grid_tb']['univ'],])

    # intermediate grid
    make_stack('dummy GT CR bank {0} grid_i'.format(b),co,cells,new_id(univIDs),cellIDs,
                [surfs['bot fuel rod']['id'],
                 surfs['dashpot top']['id'],
                 surfs['bank{0} bot'.format(b)]['id'],
                 surfs['bank{0} top'.format(b)]['id'],],
                [cells['water pin']['univ'],
                 cells['GTd empty grid_i']['univ'],
                 cells['GT empty grid_i']['univ'],
                 cells['control rod grid_i']['univ'],
                 cells['control rod blank grid_i']['univ'],])
    
    # nozzle
    make_stack('dummy GT CR bank {0} nozzle'.format(b),co,cells,new_id(univIDs),cellIDs,
                [surfs['bot fuel rod']['id'],
                 surfs['dashpot top']['id'],
                 surfs['bank{0} bot'.format(b)]['id'],
                 surfs['bank{0} top'.format(b)]['id'],],
                [cells['water pin']['univ'],
                 cells['GTd empty nozzle']['univ'],
                 cells['GT empty nozzle']['univ'],
                 cells['control rod nozzle']['univ'],
                 cells['control rod blank nozzle']['univ'],])

    # bare 
    make_stack('dummy GT CR bank {0} bare'.format(b),co,cells,new_id(univIDs),cellIDs,
                [surfs['bot fuel rod']['id'],
                 surfs['dashpot top']['id'],
                 surfs['bank{0} bot'.format(b)]['id'],
                 surfs['bank{0} top'.format(b)]['id'],],
                [cells['water pin']['univ'],
                 cells['GTd empty nozzle']['univ'],
                 cells['GT empty nozzle']['univ'],
                 cells['control rod bare']['univ'],
                 cells['control rod blank bare']['univ'],])

    # final combination of all axial pieces for control rod bank b
    make_stack('GT CR bank {0}'.format(b),co,cells,new_id(univIDs),cellIDs,
              stackSurfs,
              [cells['water pin']['univ'],                          # lower plenum
               cells['dummy GT CR bank {0} nozzle'.format(b)]['univ'],    # support plate
               cells['dummy GT CR bank {0} nozzle'.format(b)]['univ'],    # lower nozzle
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # lower thimble
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # dashpot
               cells['dummy GT CR bank {0} grid_tb'.format(b)]['univ'],      # dashpot grid 1
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # dashpot
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # reg
               cells['dummy GT CR bank {0} grid_i'.format(b)]['univ'],      # reg grid 2
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # reg
               cells['dummy GT CR bank {0} grid_i'.format(b)]['univ'],      # reg grid 3
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # reg
               cells['dummy GT CR bank {0} grid_i'.format(b)]['univ'],      # reg grid 4
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # reg
               cells['dummy GT CR bank {0} grid_i'.format(b)]['univ'],      # reg grid 5
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # reg
               cells['dummy GT CR bank {0} grid_i'.format(b)]['univ'],      # reg grid 6
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # reg
               cells['dummy GT CR bank {0} grid_i'.format(b)]['univ'],      # reg grid 7
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # reg
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # pin plenum
               cells['dummy GT CR bank {0} grid_tb'.format(b)]['univ'],      # pin plenum grid 8
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # pin plenum
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # upper fuel rod end plug
               cells['dummy GT CR bank {0}'.format(b)]['univ'],           # space between fuel rod and and upper nozzle
#               cells['dummy GT CR bank {0} nozzle'.format(b)]['univ'],    # upper nozzle
               cells['dummy GT CR bank {0} bare'.format(b)]['univ'],    # upper nozzle
               cells['dummy GT CR bank {0} bare'.format(b)]['univ']])                       # upper plenum
  

  # BURNABLE ABSORBER PIN CELLS

  # These suckers don't go all the way down to the bottom of the fuel rods, but
  # they do extend into the dashpot, ending in the middle of grid 1.

  make_pin('burn abs',"Burnable absorber pincell universes",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 1']['id'],
             surfs['burnabs rad 2']['id'],
             surfs['burnabs rad 3']['id'],
             surfs['burnabs rad 4']['id'],
             surfs['burnabs rad 5']['id'],
             surfs['burnabs rad 6']['id'],
             surfs['burnabs rad 7']['id'],
             surfs['burnabs rad 8']['id'],],
            [(mats['air']['id'],"burnable absorber pin"),
             (mats['SS304']['id'],""),
             (mats['air']['id'],""),
             (mats['borosilicate']['id'],""),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('burn abs grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 1']['id'],
             surfs['burnabs rad 2']['id'],
             surfs['burnabs rad 3']['id'],
             surfs['burnabs rad 4']['id'],
             surfs['burnabs rad 5']['id'],
             surfs['burnabs rad 6']['id'],
             surfs['burnabs rad 7']['id'],
             surfs['burnabs rad 8']['id'],],
            [(mats['air']['id'],"burnable absorber pin top/bottom grid"),
             (mats['SS304']['id'],""),
             (mats['air']['id'],""),
             (mats['borosilicate']['id'],""),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('burn abs grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 1']['id'],
             surfs['burnabs rad 2']['id'],
             surfs['burnabs rad 3']['id'],
             surfs['burnabs rad 4']['id'],
             surfs['burnabs rad 5']['id'],
             surfs['burnabs rad 6']['id'],
             surfs['burnabs rad 7']['id'],
             surfs['burnabs rad 8']['id'],],
            [(mats['air']['id'],"burnable absorber pin intermediate grid"),
             (mats['SS304']['id'],""),
             (mats['air']['id'],""),
             (mats['borosilicate']['id'],""),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')
  make_pin('burn abs dashpot',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 1']['id'],
             surfs['burnabs rad 2']['id'],
             surfs['burnabs rad 3']['id'],
             surfs['burnabs rad 4']['id'],
             surfs['burnabs rad 5']['id'],
             surfs['burnabs rad 6']['id'],
             surfs['GT dashpot IR']['id'],
             surfs['GT dashpot OR']['id'],],
            [(mats['air']['id'],"burnable absorber pin dashpot"),
             (mats['SS304']['id'],""),
             (mats['air']['id'],""),
             (mats['borosilicate']['id'],""),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('burn abs dashpot grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 1']['id'],
             surfs['burnabs rad 2']['id'],
             surfs['burnabs rad 3']['id'],
             surfs['burnabs rad 4']['id'],
             surfs['burnabs rad 5']['id'],
             surfs['burnabs rad 6']['id'],
             surfs['GT dashpot IR']['id'],
             surfs['GT dashpot OR']['id'],],
            [(mats['air']['id'],"burnable absorber pin dashpot top/bottomgrid"),
             (mats['SS304']['id'],""),
             (mats['air']['id'],""),
             (mats['borosilicate']['id'],""),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('burn abs dashpot grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 1']['id'],
             surfs['burnabs rad 2']['id'],
             surfs['burnabs rad 3']['id'],
             surfs['burnabs rad 4']['id'],
             surfs['burnabs rad 5']['id'],
             surfs['burnabs rad 6']['id'],
             surfs['GT dashpot IR']['id'],
             surfs['GT dashpot OR']['id'],],
            [(mats['air']['id'],"burnable absorber pin dashpot intermediate grid"),
             (mats['SS304']['id'],""),
             (mats['air']['id'],""),
             (mats['borosilicate']['id'],""),
             (mats['air']['id'],""),
             (mats['SS304']['id'],""),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')
  make_pin('blank burn abs ss',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 6']['id'],
             surfs['burnabs rad 7']['id'],
             surfs['burnabs rad 8']['id']],
            [(mats['SS304']['id'],"blank burnable absorber pin above active poison, inside guide tube"),
             (mats['water-nominal']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('blank burn abs ss bare',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['burnabs rad 6']['id']],
            [(mats['SS304']['id'],"blank burnable absorber pin above active poison, inside nozzle"),
             (mats['water-nominal']['id'],""),])

  # final combination of all axial pieces burnable absorber rod
  stackSurfsBA = [surfs['bot support plate']['id'],
                  surfs['top support plate']['id'],
                  surfs['top lower nozzle']['id'],
                  surfs['top lower thimble']['id'],
                  surfs['grid1bot']['id'],
                  surfs['burn abs bot']['id'],
                  surfs['grid1top']['id'],
                  surfs['dashpot top']['id'],
                  surfs['grid2bot']['id'],
                  surfs['grid2top']['id'],
                  surfs['grid3bot']['id'],
                  surfs['grid3top']['id'],
                  surfs['grid4bot']['id'],
                  surfs['grid4top']['id'],
                  surfs['grid5bot']['id'],
                  surfs['grid5top']['id'],
                  surfs['grid6bot']['id'],
                  surfs['grid6top']['id'],
                  surfs['grid7bot']['id'],
                  surfs['grid7top']['id'],
                  surfs['top active core']['id'],
                  surfs['grid8bot']['id'],
                  surfs['grid8top']['id'],
                  surfs['top pin plenum']['id'],
                  surfs['top fuel rod']['id'],
                  surfs['bot upper nozzle']['id'],
                  surfs['top upper nozzle']['id']]
                
  make_stack('burn abs stack',co,cells,new_id(univIDs),cellIDs,
            stackSurfsBA,
            [cells['water pin']['univ'],          # lower plenum
             cells['water pin']['univ'],   # support plate
             cells['water pin']['univ'],   # lower nozzle
             cells['GTd empty']['univ'],          # lower thimble
             cells['GTd empty']['univ'],          # dashpot
             cells['GTd empty grid_tb']['univ'],     # dashpot grid 1
             cells['burn abs dashpot grid_tb']['univ'],     # dashpot grid 1
             cells['burn abs dashpot']['univ'],          # dashpot
             cells['burn abs']['univ'],           # reg
             cells['burn abs grid_i']['univ'],      # reg grid 2
             cells['burn abs']['univ'],           # reg
             cells['burn abs grid_i']['univ'],      # reg grid 3
             cells['burn abs']['univ'],           # reg
             cells['burn abs grid_i']['univ'],      # reg grid 4
             cells['burn abs']['univ'],           # reg
             cells['burn abs grid_i']['univ'],      # reg grid 5
             cells['burn abs']['univ'],           # reg
             cells['burn abs grid_i']['univ'],      # reg grid 6
             cells['burn abs']['univ'],           # reg
             cells['burn abs grid_i']['univ'],      # reg grid 7
             cells['burn abs']['univ'],           # reg
             cells['blank burn abs ss']['univ'],           # pin plenum
             cells['blank burn abs ss']['univ'],      # pin plenum grid 8
             cells['blank burn abs ss']['univ'],           # pin plenum
             cells['blank burn abs ss']['univ'],           # upper fuel rod end plug
             cells['blank burn abs ss']['univ'],         # space between fuel rod and and upper nozzle
             cells['blank burn abs ss bare']['univ'],         # upper nozzle
             cells['water pin']['univ']])        # upper plenum


  # FUEL PIN CELLS
  
  make_pin('SS pin',"Fuel pincells",co,cells,new_id(univIDs),cellIDs,
            [surfs['clad OR']['id'],],
            [(mats['SS304']['id'],"SS pin for lower thimble and SS penetrations"),
             (mats['water-nominal']['id'],""),])
  make_pin('end plug',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['clad OR']['id'],],
            [(mats['zirc']['id'],"end plug"),
             (mats['water-nominal']['id'],""),])
  make_pin('pin plenum',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['plenum spring OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['inconel']['id'],"pin plenum"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('pin plenum grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['plenum spring OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['inconel']['id'],"pin plenum"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')


  ## 1.6 w/o
  
  make_pin('Fuel 1.6 w/o',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 1.6']['id'],"UO2 Fuel 1.6 w/o"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('Fuel 1.6 w/o grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 1.6']['id'],"UO2 Fuel 1.6 w/o with top/bottom grid"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('Fuel 1.6 w/o grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 1.6']['id'],"UO2 Fuel 1.6 w/o with intermediate grid"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')

  # final combination of all axial pieces for Fuel 1.6 w/o
  make_stack('Fuel 1.6 w/o stack',co,cells,new_id(univIDs),cellIDs,
            stackSurfs,
            [cells['water pin']['univ'],              # lower plenum
             cells['SS pin']['univ'],                 # support plate
             cells['SS pin']['univ'],                 # lower nozzle
             cells['end plug']['univ'],               # lower thimble
             cells['Fuel 1.6 w/o']['univ'],           # dashpot
             cells['Fuel 1.6 w/o grid_tb']['univ'],      # dashpot grid 1
             cells['Fuel 1.6 w/o']['univ'],           # dashpot
             cells['Fuel 1.6 w/o']['univ'],           # reg
             cells['Fuel 1.6 w/o grid_i']['univ'],      # reg grid 2
             cells['Fuel 1.6 w/o']['univ'],           # reg
             cells['Fuel 1.6 w/o grid_i']['univ'],      # reg grid 3
             cells['Fuel 1.6 w/o']['univ'],           # reg
             cells['Fuel 1.6 w/o grid_i']['univ'],      # reg grid 4
             cells['Fuel 1.6 w/o']['univ'],           # reg
             cells['Fuel 1.6 w/o grid_i']['univ'],      # reg grid 5
             cells['Fuel 1.6 w/o']['univ'],           # reg
             cells['Fuel 1.6 w/o grid_i']['univ'],      # reg grid 6
             cells['Fuel 1.6 w/o']['univ'],           # reg
             cells['Fuel 1.6 w/o grid_i']['univ'],      # reg grid 7
             cells['Fuel 1.6 w/o']['univ'],           # reg
             cells['pin plenum']['univ'],             # pin plenum
             cells['pin plenum grid_tb']['univ'],        # pin plenum grid 8
             cells['pin plenum']['univ'],             # pin plenum
             cells['end plug']['univ'],               # upper fuel rod end plug
             cells['water pin']['univ'],              # space between fuel rod and and upper nozzle
             cells['SS pin']['univ'],                 # upper nozzle
             cells['water pin']['univ']])             # upper plenum

  ## 2.4 w/o

  make_pin('Fuel 2.4 w/o',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 2.4']['id'],"UO2 Fuel 2.4 w/o"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('Fuel 2.4 w/o grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 2.4']['id'],"UO2 Fuel 2.4 w/o with top/bottom grid"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('Fuel 2.4 w/o grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 2.4']['id'],"UO2 Fuel 2.4 w/o with intermediate grid"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')

  # final combination of all axial pieces for Fuel 2.4 w/o
  make_stack('Fuel 2.4 w/o stack',co,cells,new_id(univIDs),cellIDs,
            stackSurfs,
            [cells['water pin']['univ'],              # lower plenum
             cells['SS pin']['univ'],                 # support plate
             cells['SS pin']['univ'],                 # lower nozzle
             cells['end plug']['univ'],               # lower thimble
             cells['Fuel 2.4 w/o']['univ'],           # dashpot
             cells['Fuel 2.4 w/o grid_tb']['univ'],      # dashpot grid 1
             cells['Fuel 2.4 w/o']['univ'],           # dashpot
             cells['Fuel 2.4 w/o']['univ'],           # reg
             cells['Fuel 2.4 w/o grid_i']['univ'],      # reg grid 2
             cells['Fuel 2.4 w/o']['univ'],           # reg
             cells['Fuel 2.4 w/o grid_i']['univ'],      # reg grid 3
             cells['Fuel 2.4 w/o']['univ'],           # reg
             cells['Fuel 2.4 w/o grid_i']['univ'],      # reg grid 4
             cells['Fuel 2.4 w/o']['univ'],           # reg
             cells['Fuel 2.4 w/o grid_i']['univ'],      # reg grid 5
             cells['Fuel 2.4 w/o']['univ'],           # reg
             cells['Fuel 2.4 w/o grid_i']['univ'],      # reg grid 6
             cells['Fuel 2.4 w/o']['univ'],           # reg
             cells['Fuel 2.4 w/o grid_i']['univ'],      # reg grid 7
             cells['Fuel 2.4 w/o']['univ'],           # reg
             cells['pin plenum']['univ'],             # pin plenum
             cells['pin plenum grid_tb']['univ'],        # pin plenum grid 8
             cells['pin plenum']['univ'],             # pin plenum
             cells['end plug']['univ'],               # upper fuel rod end plug
             cells['water pin']['univ'],              # space between fuel rod and and upper nozzle
             cells['SS pin']['univ'],                 # upper nozzle
             cells['water pin']['univ']])             # upper plenum

  ## 3.1 w/o

  make_pin('Fuel 3.1 w/o',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 3.1']['id'],"UO2 Fuel 3.1 w/o"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),])
  make_pin('Fuel 3.1 w/o grid_tb',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 3.1']['id'],"UO2 Fuel 3.1 w/o with top/bottom grid"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['inconel']['id'], gridType='tb')
  make_pin('Fuel 3.1 w/o grid_i',"",co,cells,new_id(univIDs),cellIDs,
            [surfs['pellet OR']['id'],
             surfs['clad IR']['id'],
             surfs['clad OR']['id'],],
            [(mats['UO2 3.1']['id'],"UO2 Fuel 3.1 w/o with intermediate grid"),
             (mats['helium']['id'],""),
             (mats['zirc']['id'],""),
             (mats['water-nominal']['id'],""),],
            grid=True, surfs=surfs, gridMat=mats['zirc']['id'], gridType='i')

  # final combination of all axial pieces for Fuel 3.1 w/o
  make_stack('Fuel 3.1 w/o stack',co,cells,new_id(univIDs),cellIDs,
            stackSurfs,
            [cells['water pin']['univ'],              # lower plenum
             cells['SS pin']['univ'],                 # support plate
             cells['SS pin']['univ'],                 # lower nozzle
             cells['end plug']['univ'],               # lower thimble
             cells['Fuel 3.1 w/o']['univ'],           # dashpot
             cells['Fuel 3.1 w/o grid_tb']['univ'],      # dashpot grid 1
             cells['Fuel 3.1 w/o']['univ'],           # dashpot
             cells['Fuel 3.1 w/o']['univ'],           # reg
             cells['Fuel 3.1 w/o grid_i']['univ'],      # reg grid 2
             cells['Fuel 3.1 w/o']['univ'],           # reg
             cells['Fuel 3.1 w/o grid_i']['univ'],      # reg grid 3
             cells['Fuel 3.1 w/o']['univ'],           # reg
             cells['Fuel 3.1 w/o grid_i']['univ'],      # reg grid 4
             cells['Fuel 3.1 w/o']['univ'],           # reg
             cells['Fuel 3.1 w/o grid_i']['univ'],      # reg grid 5
             cells['Fuel 3.1 w/o']['univ'],           # reg
             cells['Fuel 3.1 w/o grid_i']['univ'],      # reg grid 6
             cells['Fuel 3.1 w/o']['univ'],           # reg
             cells['Fuel 3.1 w/o grid_i']['univ'],      # reg grid 7
             cells['Fuel 3.1 w/o']['univ'],           # reg
             cells['pin plenum']['univ'],             # pin plenum
             cells['pin plenum grid_tb']['univ'],        # pin plenum grid 8
             cells['pin plenum']['univ'],             # pin plenum
             cells['end plug']['univ'],               # upper fuel rod end plug
             cells['water pin']['univ'],              # space between fuel rod and and upper nozzle
             cells['SS pin']['univ'],                 # upper nozzle
             cells['water pin']['univ']])             # upper plenum


################## Baffle construction ##################


  lo = {'n': 0}
  latts = {}

  # baffle north

  surfs['baffle surf north'] =  { 'order':   inc_order(so),
                                  'section': comm_t.format("Baffle surfaces"),
                                  'comm':    comm_t.format("chosen for 2x2 baffle lattice, so w={0}".format(baffleWidth)),
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '{0:0<8.5}'.format(latticePitch/4 - baffleWidth)}
  cells['baffle dummy north'] = { 'order':   inc_order(co),
                                  'section': comm_t.format("Baffle cells"),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf north']['id'])}
  cells['baf dummy north 2'] =  { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baffle dummy north']['univ'],
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf north']['id'])}
  latts['baffle north'] =       { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle north"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{0:>3} {0:>3}
{1:>3} {1:>3}
""".format(cells['baffle dummy north']['univ'],cells['water pin']['univ'])}
  cells['baffle north'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("north baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle north']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle south

  surfs['baffle surf south'] =  { 'order':   inc_order(so),
                                  'comm':    "",
                                  'id':      new_id(surfIDs),
                                  'type':    '"y-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(baffleWidth - latticePitch/4)}
  cells['baffle dummy south'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf south']['id'])}
  cells['baf dummy south 2'] =  { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baffle dummy south']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf south']['id'])}
  latts['baffle south'] =       { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle south"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{1:>3} {1:>3}
{0:>3} {0:>3}
""".format(cells['baffle dummy south']['univ'],cells['water pin']['univ'])}
  cells['baffle south'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("south baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle south']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # bafffle east

  surfs['baffle surf east'] =   { 'order':   inc_order(so),
                                  'comm':    "",
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '{0:0<8.5}'.format(latticePitch/4 - baffleWidth)}
  cells['baffle dummy east'] =  { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf east']['id'])}
  cells['baf dummy east 2'] =   { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baffle dummy east']['univ'],
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf east']['id'])}
  latts['baffle east'] =       { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle east"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{1:>3} {0:>3}
{1:>3} {0:>3}
""".format(cells['baffle dummy east']['univ'],cells['water pin']['univ'])}
  cells['baffle east'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("east baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle east']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle west

  surfs['baffle surf west'] =   { 'order':   inc_order(so),
                                  'comm':    "",
                                  'id':      new_id(surfIDs),
                                  'type':    '"x-plane"',
                                  'coeffs':  '{0:0<8.6}'.format(baffleWidth - latticePitch/4)}
  cells['baffle dummy west'] =  { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf west']['id'])}
  cells['baf dummy west 2'] =   { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baffle dummy west']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf west']['id'])}
  latts['baffle west'] =       { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle west"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{0:>3} {1:>3}
{0:>3} {1:>3}
""".format(cells['baffle dummy west']['univ'],cells['water pin']['univ'])}
  cells['baffle west'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("west baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle west']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}


  # baffle NW edges

  cells['baf dummy edges NW'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} -{1}'.format(surfs['baffle surf west']['id'],surfs['baffle surf north']['id'])}
  cells['baf dmy edges NW 2'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges NW']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} {1}'.format(surfs['baffle surf west']['id'],surfs['baffle surf north']['id'])}
  cells['baf dmy edges NW 3'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges NW']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf west']['id'])}
  latts['baffle edges NW'] =    { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle NW edges"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{0:>3} {1:>3}
{2:>3} {3:>3}
""".format(cells['baf dummy edges NW']['univ'],cells['baffle dummy north']['univ'],cells['baffle dummy west']['univ'],cells['water pin']['univ'])}
  cells['baffle edges NW'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("NW edges baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle edges NW']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}


  # baffle NE edges

  cells['baf dummy edges NE'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0} -{1}'.format(surfs['baffle surf north']['id'],surfs['baffle surf east']['id'])}
  cells['baf dmy edges NE 2'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges NE']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} -{1}'.format(surfs['baffle surf north']['id'],surfs['baffle surf east']['id'])}
  cells['baf dmy edges NE 3'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges NE']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf east']['id'])}
  latts['baffle edges NE'] =    { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle NE edges"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{0:>3} {1:>3}
{2:>3} {3:>3}
""".format(cells['baffle dummy north']['univ'],cells['baf dummy edges NE']['univ'],cells['water pin']['univ'],cells['baffle dummy east']['univ'])}
  cells['baffle edges NE'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("NE edges baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle edges NE']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle SW edges

  cells['baf dummy edges SW'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} {1}'.format(surfs['baffle surf south']['id'],surfs['baffle surf west']['id'])}
  cells['baf dmy edges SW 2'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges SW']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0} {1}'.format(surfs['baffle surf south']['id'],surfs['baffle surf west']['id'])}
  cells['baf dmy edges SW 3'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges SW']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf west']['id'])}
  latts['baffle edges SW'] =    { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle SW edges"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{0:>3} {1:>3}
{2:>3} {3:>3}
""".format(cells['baffle dummy west']['univ'],cells['water pin']['univ'],cells['baf dummy edges SW']['univ'],cells['baffle dummy south']['univ'])}
  cells['baffle edges SW'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("NE edges baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle edges SW']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle SE edges

  cells['baf dummy edges SE'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} -{1}'.format(surfs['baffle surf south']['id'],surfs['baffle surf east']['id'])}
  cells['baf dmy edges SE 2'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges SE']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0} -{1}'.format(surfs['baffle surf south']['id'],surfs['baffle surf east']['id'])}
  cells['baf dmy edges SE 3'] = { 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy edges SE']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf east']['id'])}
  latts['baffle edges SE'] =    { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle SE edges"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{0:>3} {1:>3}
{2:>3} {3:>3}
""".format(cells['water pin']['univ'],cells['baffle dummy east']['univ'],cells['baffle dummy south']['univ'],cells['baf dummy edges SE']['univ'])}
  cells['baffle edges SE'] =       { 'order':   inc_order(co),
                                  'comm':    comm_t.format("NE edges baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle edges SE']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle NW corner

  cells['baf dummy corner NW'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0} -{1}'.format(surfs['baffle surf west']['id'],surfs['baffle surf north']['id'])}
  cells['baf dmy corner NW 2'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner NW']['univ'],
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf west']['id'])}
  cells['baf dmy corner NW 3'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner NW']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} -{1}'.format(surfs['baffle surf north']['id'],surfs['baffle surf west']['id'])}
  latts['baffle corner NW'] =   { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle NW corner"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{0:>3} {1:>3}
{1:>3} {1:>3}
""".format(cells['baf dummy corner NW']['univ'],cells['water pin']['univ'])}
  cells['baffle corner NW'] =   { 'order':   inc_order(co),
                                  'comm':    comm_t.format("NW corner baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle corner NW']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle NE corner

  cells['baf dummy corner NE'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} -{1}'.format(surfs['baffle surf east']['id'],surfs['baffle surf north']['id'])}
  cells['baf dmy corner NE 2'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner NE']['univ'],
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf east']['id'])}
  cells['baf dmy corner NE 3'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner NE']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} {1}'.format(surfs['baffle surf north']['id'],surfs['baffle surf east']['id'])}
  latts['baffle corner NE'] =   { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle NE corner"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{1:>3} {0:>3}
{1:>3} {1:>3}
""".format(cells['baf dummy corner NE']['univ'],cells['water pin']['univ'])}
  cells['baffle corner NE'] =   { 'order':   inc_order(co),
                                  'comm':    comm_t.format("NE corner baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle corner NE']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle SE corner

  cells['baf dummy corner SE'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0} {1}'.format(surfs['baffle surf east']['id'],surfs['baffle surf south']['id'])}
  cells['baf dmy corner SE 2'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner SE']['univ'],
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['baffle surf east']['id'])}
  cells['baf dmy corner SE 3'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner SE']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0} {1}'.format(surfs['baffle surf south']['id'],surfs['baffle surf east']['id'])}
  latts['baffle corner SE'] =   { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle SE corner"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{1:>3} {1:>3}
{1:>3} {0:>3}
""".format(cells['baf dummy corner SE']['univ'],cells['water pin']['univ'])}
  cells['baffle corner SE'] =   { 'order':   inc_order(co),
                                  'comm':    comm_t.format("SE corner baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle corner SE']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  # baffle SW corner

  cells['baf dummy corner SW'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0} {1}'.format(surfs['baffle surf west']['id'],surfs['baffle surf south']['id'])}
  cells['baf dmy corner SW 2'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner SW']['univ'],
                                  'mat':     mats['water-nominal']['id'],
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['baffle surf west']['id'])}
  cells['baf dmy corner SW 3'] ={ 'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    cells['baf dummy corner SW']['univ'],
                                  'mat':     mats['SS304']['id'],
                                  'fill':    None,
                                  'surfs':  '-{0} -{1}'.format(surfs['baffle surf south']['id'],surfs['baffle surf west']['id'])}
  latts['baffle corner SW'] =   { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Baffle SW corner"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     2,
                                  'lleft':   -latticePitch/2,
                                  'width':   latticePitch/2,
                                  'univs':   """
{1:>3} {1:>3}
{0:>3} {1:>3}
""".format(cells['baf dummy corner SW']['univ'],cells['water pin']['univ'])}
  cells['baffle corner SW'] =   { 'order':   inc_order(co),
                                  'comm':    comm_t.format("SW corner baffle universe"),
                                  'id':      new_id(cellIDs),
                                  'univ':    new_id(univIDs),
                                  'fill':    latts['baffle corner SW']['id'],
                                  'surfs':  '-{0}'.format(surfs['dummy outer']['id'])}

  baffle = {}
  baffle['bafn_'] = cells['baffle south']['univ']
  baffle['bafs_'] = cells['baffle north']['univ']
  baffle['bafe_'] = cells['baffle west']['univ']
  baffle['bafw_'] = cells['baffle east']['univ']
  baffle['bafnw'] = cells['baffle corner SE']['univ']
  baffle['bafne'] = cells['baffle corner SW']['univ']
  baffle['bafsw'] = cells['baffle corner NE']['univ']
  baffle['bafse'] = cells['baffle corner NW']['univ']
  baffle['bfcnw'] = cells['baffle edges SE']['univ']
  baffle['bfcne'] = cells['baffle edges SW']['univ']
  baffle['bfcsw'] = cells['baffle edges NE']['univ']
  baffle['bfcse'] = cells['baffle edges NW']['univ']




  ################## lattices ##################

  assemblyCells = {} # convenience dictionary holds which cells point to which assembly type

  # commonly needed universes
  gtu = cells['GT empty stack']['univ']
  gti = cells['GT empty stack instr']['univ']
  bas = cells['burn abs stack']['univ']
  ins = cells['GT instr stack']['univ']
  crA = cells['GT CR bank A']['univ']
  crB = cells['GT CR bank B']['univ']
  crC = cells['GT CR bank C']['univ']
  crD = cells['GT CR bank D']['univ']
  crSA = cells['GT CR bank SA']['univ']
  crSB = cells['GT CR bank SB']['univ']
  crSC = cells['GT CR bank SC']['univ']
  crSD = cells['GT CR bank SD']['univ']
  crSE = cells['GT CR bank SE']['univ']

  latDims = { 'dim':17, 'lleft':-17*pinPitch/2, 'width':pinPitch}

  ## 1.6 w/o assemblies

  for cent,comment1 in [(gti,""),(ins," + instr")]:

    if cent == gti:
      sss = "Core Lattice universes"
    else:
      sss = None

    # No BAs
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 1.6 w/o'+comment1,comm="Assembly 1.6 w/o no BAs"+comment1, sect=sss,
                                 univs=pinLattice_t.format(cells['Fuel 1.6 w/o stack']['univ'],
                                                    a=gtu,  b=gtu,  c=gtu,
                                            d=gtu,                          e=gtu,
                                            f=gtu,  g=gtu,  h=gtu,  i=gtu,  j=gtu,
                                            k=gtu,  l=gtu,  m=cent, n=gtu,  o=gtu,
                                            p=gtu,  q=gtu,  r=gtu,  s=gtu,  t=gtu,
                                            u=gtu,                          v=gtu,
                                                    w=gtu,  x=gtu,  y=gtu,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 1.6 w/o'+comment1] = newCells

    for bank,comment2 in [(crA,' + CRA'),(crB,' + CRB'),(crC,' + CRC'),(crD,' + CRD'),
                          (crSB,' + shutB'),(crSC,' + shutC'),(crSD,' + shutD'),(crSE,' + shutE')]:

      newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                    'Fuel 1.6 w/o'+comment1+comment2,comm="Assembly 1.6 w/o"+comment1+comment2,
                                   univs=pinLattice_t.format(cells['Fuel 1.6 w/o stack']['univ'],
                                                      a=bank,  b=bank,  c=bank,
                                              d=bank,                             e=bank,
                                              f=bank,  g=bank,  h=bank,  i=bank,  j=bank,
                                              k=bank,  l=bank,  m=cent,  n=bank,  o=bank,
                                              p=bank,  q=bank,  r=bank,  s=bank,  t=bank,
                                              u=bank,                             v=bank,
                                                      w=bank,  x=bank,  y=bank,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
      assemblyCells['Fuel 1.6 w/o'+comment1+comment2] = newCells



  ## 2.4 w/o assemblies

  for cen,comment1 in [(gti,""),(ins," + instr")]:

    # no BAs
    bank = gtu
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 2.4 w/o'+comment1,comm="Assembly 2.4 w/o no BAs"+comment1,
                                 univs=pinLattice_t.format(cells['Fuel 2.4 w/o stack']['univ'],
                                                    a=bank,  b=bank,  c=bank,
                                            d=bank,                             e=bank,
                                            f=bank,  g=bank,  h=bank,  i=bank,  j=bank,
                                            k=bank,  l=bank,  m=cen,   n=bank,  o=bank,
                                            p=bank,  q=bank,  r=bank,  s=bank,  t=bank,
                                            u=bank,                             v=bank,
                                                    w=bank,  x=bank,  y=bank,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 2.4 w/o'+comment1] = newCells

    # CRD
    bank = crD
    comment2 = ' + CRD'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 2.4 w/o'+comment1+comment2,comm="Assembly 2.4 w/o "+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 2.4 w/o stack']['univ'],
                                                    a=bank,  b=bank,  c=bank,
                                            d=bank,                             e=bank,
                                            f=bank,  g=bank,  h=bank,  i=bank,  j=bank,
                                            k=bank,  l=bank,  m=cen,   n=bank,  o=bank,
                                            p=bank,  q=bank,  r=bank,  s=bank,  t=bank,
                                            u=bank,                             v=bank,
                                                    w=bank,  x=bank,  y=bank,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 2.4 w/o'+comment1+comment2] = newCells

    # 12 BAs
    comment2 = ' + 12BA'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 2.4 w/o'+comment1+comment2,comm="Assembly 2.4 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 2.4 w/o stack']['univ'],
                                                    a=bas,  b=gtu,  c=bas,
                                            d=bas,                          e=bas,
                                            f=bas,  g=gtu,  h=gtu,  i=gtu,  j=bas,
                                            k=gtu,  l=gtu,  m=cen,  n=gtu,  o=gtu,
                                            p=bas,  q=gtu,  r=gtu,  s=gtu,  t=bas,
                                            u=bas,                          v=bas,
                                                    w=bas,  x=gtu,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 2.4 w/o'+comment1+comment2] = newCells

    # 16 BAs
    comment2 = ' + 16BA'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 2.4 w/o'+comment1+comment2,comm="Assembly 2.4 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 2.4 w/o stack']['univ'],
                                                    a=bas,  b=bas,  c=bas,
                                            d=bas,                          e=bas,
                                            f=bas,  g=gtu,  h=gtu,  i=gtu,  j=bas,
                                            k=bas,  l=gtu,  m=cen,  n=gtu,  o=bas,
                                            p=bas,  q=gtu,  r=gtu,  s=gtu,  t=bas,
                                            u=bas,                          v=bas,
                                                    w=bas,  x=bas,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 2.4 w/o'+comment1+comment2] = newCells
    
  ## 3.1 w/o assemblies

  for cen,comment1 in [(gti,""),(ins," + instr")]:

    # no BAs
    bank = gtu
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1,comm="Assembly 3.1 w/o no BAs"+comment1,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bank,  b=bank,  c=bank,
                                            d=bank,                             e=bank,
                                            f=bank,  g=bank,  h=bank,  i=bank,  j=bank,
                                            k=bank,  l=bank,  m=cen,   n=bank,  o=bank,
                                            p=bank,  q=bank,  r=bank,  s=bank,  t=bank,
                                            u=bank,                             v=bank,
                                                    w=bank,  x=bank,  y=bank,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1] = newCells
    
    # shut
    bank = crSA
    comment2 = ' + shutA'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bank,  b=bank,  c=bank,
                                            d=bank,                             e=bank,
                                            f=bank,  g=bank,  h=bank,  i=bank,  j=bank,
                                            k=bank,  l=bank,  m=cen,   n=bank,  o=bank,
                                            p=bank,  q=bank,  r=bank,  s=bank,  t=bank,
                                            u=bank,                             v=bank,
                                                    w=bank,  x=bank,  y=bank,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells
    
    # 20 BAs
    comment2 = ' + 20BA'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bas,  b=bas,  c=bas,
                                            d=bas,                          e=bas,
                                            f=bas,  g=bas,  h=gtu,  i=bas,  j=bas,
                                            k=bas,  l=gtu,  m=cen,  n=gtu,  o=bas,
                                            p=bas,  q=bas,  r=gtu,  s=bas,  t=bas,
                                            u=bas,                          v=bas,
                                                    w=bas,  x=bas,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells

    # 16 BAs
    comment2 = ' + 16BA'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bas,  b=bas,  c=bas,
                                            d=bas,                          e=bas,
                                            f=bas,  g=gtu,  h=gtu,  i=gtu,  j=bas,
                                            k=bas,  l=gtu,  m=cen,  n=gtu,  o=bas,
                                            p=bas,  q=gtu,  r=gtu,  s=gtu,  t=bas,
                                            u=bas,                          v=bas,
                                                    w=bas,  x=bas,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells
    
    # 15 BAs NW
    comment2 = ' + 15BANW'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=gtu,  b=gtu,  c=gtu,
                                            d=gtu,                          e=gtu,
                                            f=gtu,  g=bas,  h=bas,  i=bas,  j=bas,
                                            k=gtu,  l=bas,  m=cen,  n=bas,  o=bas,
                                            p=gtu,  q=bas,  r=bas,  s=bas,  t=bas,
                                            u=gtu,                          v=bas,
                                                    w=bas,  x=bas,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells

    # 15 BAs NE
    comment2 = ' + 15BANE'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=gtu,  b=gtu,  c=gtu,
                                            d=gtu,                          e=gtu,
                                            f=bas,  g=bas,  h=bas,  i=bas,  j=gtu,
                                            k=bas,  l=bas,  m=cen,  n=bas,  o=gtu,
                                            p=bas,  q=bas,  r=bas,  s=bas,  t=gtu,
                                            u=bas,                          v=gtu,
                                                    w=bas,  x=bas,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells
    
    # 15 BAs SW
    comment2 = ' + 15BASW'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bas,  b=bas,  c=bas,
                                            d=gtu,                          e=bas,
                                            f=gtu,  g=bas,  h=bas,  i=bas,  j=bas,
                                            k=gtu,  l=bas,  m=cen,  n=bas,  o=bas,
                                            p=gtu,  q=bas,  r=bas,  s=bas,  t=bas,
                                            u=gtu,                          v=gtu,
                                                    w=gtu,  x=gtu,  y=gtu,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells
    
    # 15 BAs SE
    comment2 = ' + 15BASE'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bas,  b=bas,  c=bas,
                                            d=bas,                          e=gtu,
                                            f=bas,  g=bas,  h=bas,  i=bas,  j=gtu,
                                            k=bas,  l=bas,  m=cen,  n=bas,  o=gtu,
                                            p=bas,  q=bas,  r=bas,  s=bas,  t=gtu,
                                            u=gtu,                          v=gtu,
                                                    w=gtu,  x=gtu,  y=gtu,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells
    
    # 6 BAs N
    comment2 = ' + 6BAN'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=gtu,  b=gtu,  c=gtu,
                                            d=gtu,                          e=gtu,
                                            f=gtu,  g=gtu,  h=gtu,  i=gtu,  j=gtu,
                                            k=gtu,  l=gtu,  m=cen,  n=gtu,  o=gtu,
                                            p=bas,  q=gtu,  r=gtu,  s=gtu,  t=bas,
                                            u=bas,                          v=bas,
                                                    w=bas,  x=gtu,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells
    
    # 6 BAs S
    comment2 = ' + 6BAS'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bas,  b=gtu,  c=bas,
                                            d=bas,                          e=bas,
                                            f=bas,  g=gtu,  h=gtu,  i=gtu,  j=bas,
                                            k=gtu,  l=gtu,  m=cen,  n=gtu,  o=gtu,
                                            p=gtu,  q=gtu,  r=gtu,  s=gtu,  t=gtu,
                                            u=gtu,                          v=gtu,
                                                    w=gtu,  x=gtu,  y=gtu,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells

    # 6 BAs W
    comment2 = ' + 6BAW'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=gtu,  b=gtu,  c=bas,
                                            d=gtu,                          e=bas,
                                            f=gtu,  g=gtu,  h=gtu,  i=gtu,  j=bas,
                                            k=gtu,  l=gtu,  m=cen,  n=gtu,  o=gtu,
                                            p=gtu,  q=gtu,  r=gtu,  s=gtu,  t=bas,
                                            u=gtu,                          v=bas,
                                                    w=gtu,  x=gtu,  y=bas,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells
    
    # 6 BAs E
    comment2 = ' + 6BAE'
    newCells = make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, mats['water-nominal']['id'],
                  'Fuel 3.1 w/o'+comment1+comment2,comm="Assembly 3.1 w/o"+comment1+comment2,
                                 univs=pinLattice_t.format(cells['Fuel 3.1 w/o stack']['univ'],
                                                    a=bas,  b=gtu,  c=gtu,
                                            d=bas,                          e=gtu,
                                            f=bas,  g=gtu,  h=gtu,  i=gtu,  j=gtu,
                                            k=gtu,  l=gtu,  m=cen,  n=gtu,  o=gtu,
                                            p=bas,  q=gtu,  r=gtu,  s=gtu,  t=gtu,
                                            u=bas,                          v=gtu,
                                                    w=bas,  x=gtu,  y=gtu,),
                                 gridSurfs=gridSurfaces, sleeveMats=[mats['SS304']['id'],mats['zirc']['id']],
                                 **latDims)
    assemblyCells['Fuel 3.1 w/o'+comment1+comment2] = newCells

    

  ################## Main Core Lattices ##################

  latts['Main Core'] =          { 'order':   inc_order(lo),
                                  'comm':    comm_t.format("Main Core Lattice"),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular',
                                  'dim':     19,
                                  'lleft':   -19*latticePitch/2,
                                  'width':   latticePitch,
                                  'univs':   coreLattice_t.format(
dummy = cells['water pin']['univ'],
L___1 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
K___1 = cells['Fuel 3.1 w/o + 6BAN lattice']['univ'],
J___1 = cells['Fuel 3.1 w/o lattice']['univ'], # Sp spare location
H___1 = cells['Fuel 3.1 w/o + instr + 6BAN lattice']['univ'],
G___1 = cells['Fuel 3.1 w/o lattice']['univ'],
F___1 = cells['Fuel 3.1 w/o + instr + 6BAN lattice']['univ'],
E___1 = cells['Fuel 3.1 w/o lattice']['univ'],
N___2 = cells['Fuel 3.1 w/o lattice']['univ'],
M___2 = cells['Fuel 3.1 w/o + shutA lattice']['univ'],
L___2 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
K___2 = cells['Fuel 1.6 w/o + CRB lattice']['univ'],
J___2 = cells['Fuel 3.1 w/o + instr + 20BA lattice']['univ'],
H___2 = cells['Fuel 1.6 w/o + CRC lattice']['univ'],
G___2 = cells['Fuel 3.1 w/o + 20BA lattice']['univ'], #really should be 23BA1P
F___2 = cells['Fuel 1.6 w/o + CRB lattice']['univ'],
E___2 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
D___2 = cells['Fuel 3.1 w/o + instr + shutA lattice']['univ'],
C___2 = cells['Fuel 3.1 w/o lattice']['univ'],
P___3 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
N___3 = cells['Fuel 3.1 w/o + instr + 15BANW lattice']['univ'],
M___3 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
L___3 = cells['Fuel 1.6 w/o + shutD lattice']['univ'],
K___3 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
J___3 = cells['Fuel 1.6 w/o + shutB lattice']['univ'],
H___3 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
G___3 = cells['Fuel 1.6 w/o + shutB lattice']['univ'],
F___3 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
E___3 = cells['Fuel 1.6 w/o + shutD lattice']['univ'],
D___3 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
C___3 = cells['Fuel 3.1 w/o + 15BANE lattice']['univ'],
B___3 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
R___4 = cells['Fuel 3.1 w/o lattice']['univ'],
P___4 = cells['Fuel 3.1 w/o + shutA lattice']['univ'],
N___4 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
M___4 = cells['Fuel 2.4 w/o + CRD lattice']['univ'],
L___4 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
K___4 = cells['Fuel 1.6 w/o lattice']['univ'],      #4 secondary source rods here
J___4 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
H___4 = cells['Fuel 1.6 w/o + shutE lattice']['univ'],
G___4 = cells['Fuel 2.4 w/o + instr + 12BA lattice']['univ'],
F___4 = cells['Fuel 1.6 w/o lattice']['univ'],
E___4 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
D___4 = cells['Fuel 2.4 w/o + CRD lattice']['univ'],
C___4 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
B___4 = cells['Fuel 3.1 w/o + shutA lattice']['univ'],
A___4 = cells['Fuel 3.1 w/o lattice']['univ'],
R___5 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
P___5 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
N___5 = cells['Fuel 1.6 w/o + instr + shutC lattice']['univ'],
M___5 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
L___5 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
K___5 = cells['Fuel 2.4 w/o + instr + 12BA lattice']['univ'],
J___5 = cells['Fuel 1.6 w/o lattice']['univ'],
H___5 = cells['Fuel 2.4 w/o + instr + 12BA lattice']['univ'],
G___5 = cells['Fuel 1.6 w/o lattice']['univ'],
F___5 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
E___5 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
D___5 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
C___5 = cells['Fuel 1.6 w/o + shutD lattice']['univ'],
B___5 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
A___5 = cells['Fuel 3.1 w/o lattice']['univ'],
R___6 = cells['Fuel 3.1 w/o + 6BAW lattice']['univ'],
P___6 = cells['Fuel 1.6 w/o + CRB lattice']['univ'],
N___6 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
M___6 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
L___6 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
K___6 = cells['Fuel 1.6 w/o + CRC lattice']['univ'],
J___6 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
H___6 = cells['Fuel 1.6 w/o + CRA lattice']['univ'],
G___6 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
F___6 = cells['Fuel 1.6 w/o + instr + CRC lattice']['univ'],
E___6 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
D___6 = cells['Fuel 1.6 w/o lattice']['univ'],
C___6 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
B___6 = cells['Fuel 1.6 w/o + instr + CRB lattice']['univ'],
A___6 = cells['Fuel 3.1 w/o + 6BAE lattice']['univ'],
R___7 = cells['Fuel 3.1 w/o lattice']['univ'],
P___7 = cells['Fuel 3.1 w/o + instr + 20BA lattice']['univ'],
N___7 = cells['Fuel 1.6 w/o + shutB lattice']['univ'],
M___7 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
L___7 = cells['Fuel 1.6 w/o lattice']['univ'],
K___7 = cells['Fuel 2.4 w/o + instr + 12BA lattice']['univ'],
J___7 = cells['Fuel 1.6 w/o lattice']['univ'],
H___7 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
G___7 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
F___7 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
E___7 = cells['Fuel 1.6 w/o lattice']['univ'],
D___7 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
C___7 = cells['Fuel 1.6 w/o + shutB lattice']['univ'],
B___7 = cells['Fuel 3.1 w/o + 20BA lattice']['univ'],
A___7 = cells['Fuel 3.1 w/o + instr lattice']['univ'], # Sp spare location
R___8 = cells['Fuel 3.1 w/o + instr + 6BAW lattice']['univ'],
P___8 = cells['Fuel 1.6 w/o + CRC lattice']['univ'],
N___8 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
M___8 = cells['Fuel 1.6 w/o + shutE lattice']['univ'],
L___8 = cells['Fuel 2.4 w/o + instr + 12BA lattice']['univ'],
K___8 = cells['Fuel 1.6 w/o + CRA lattice']['univ'],
J___8 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
H___8 = cells['Fuel 1.6 w/o + CRD lattice']['univ'],
G___8 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
F___8 = cells['Fuel 1.6 w/o + instr + CRA lattice']['univ'],
E___8 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
D___8 = cells['Fuel 1.6 w/o + instr + shutE lattice']['univ'],
C___8 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
B___8 = cells['Fuel 1.6 w/o + instr + CRC lattice']['univ'],
A___8 = cells['Fuel 3.1 w/o + 6BAE lattice']['univ'],
R___9 = cells['Fuel 3.1 w/o lattice']['univ'], # Sp spare location
P___9 = cells['Fuel 3.1 w/o + 20BA lattice']['univ'],
N___9 = cells['Fuel 1.6 w/o + shutB lattice']['univ'],
M___9 = cells['Fuel 2.4 w/o + instr + 12BA lattice']['univ'],
L___9 = cells['Fuel 1.6 w/o lattice']['univ'],
K___9 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
J___9 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
H___9 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
G___9 = cells['Fuel 1.6 w/o lattice']['univ'],
F___9 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
E___9 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
D___9 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
C___9 = cells['Fuel 1.6 w/o + shutB lattice']['univ'],
B___9 = cells['Fuel 3.1 w/o + 20BA lattice']['univ'],
A___9 = cells['Fuel 3.1 w/o lattice']['univ'],
R__10 = cells['Fuel 3.1 w/o + 6BAW lattice']['univ'],
P__10 = cells['Fuel 1.6 w/o + instr + CRB lattice']['univ'],
N__10 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
M__10 = cells['Fuel 1.6 w/o lattice']['univ'],
L__10 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
K__10 = cells['Fuel 1.6 w/o + CRC lattice']['univ'],
J__10 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
H__10 = cells['Fuel 1.6 w/o + instr + CRA lattice']['univ'],
G__10 = cells['Fuel 2.4 w/o + instr + 12BA lattice']['univ'],
F__10 = cells['Fuel 1.6 w/o + CRC lattice']['univ'],
E__10 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
D__10 = cells['Fuel 1.6 w/o lattice']['univ'],
C__10 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
B__10 = cells['Fuel 1.6 w/o + CRB lattice']['univ'],
A__10 = cells['Fuel 3.1 w/o + instr + 6BAE lattice']['univ'],
R__11 = cells['Fuel 3.1 w/o lattice']['univ'],
P__11 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
N__11 = cells['Fuel 1.6 w/o + shutD lattice']['univ'],
M__11 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
L__11 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
K__11 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
J__11 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
H__11 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
G__11 = cells['Fuel 1.6 w/o lattice']['univ'],
F__11 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
E__11 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
D__11 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
C__11 = cells['Fuel 1.6 w/o + shutC lattice']['univ'],
B__11 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
A__11 = cells['Fuel 3.1 w/o lattice']['univ'],
P__12 = cells['Fuel 3.1 w/o + instr + shutA lattice']['univ'],
N__12 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
M__12 = cells['Fuel 2.4 w/o + instr + CRD lattice']['univ'],
L__12 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
K__12 = cells['Fuel 1.6 w/o + instr lattice']['univ'],
J__12 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
H__12 = cells['Fuel 1.6 w/o + instr + shutE lattice']['univ'],
G__12 = cells['Fuel 2.4 w/o + 12BA lattice']['univ'],
F__12 = cells['Fuel 1.6 w/o lattice']['univ'],      #4 secondary source rods here
E__12 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
D__12 = cells['Fuel 2.4 w/o + CRD lattice']['univ'],
C__12 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
B__12 = cells['Fuel 3.1 w/o + shutA lattice']['univ'],
P__13 = cells['Fuel 3.1 w/o lattice']['univ'],
N__13 = cells['Fuel 3.1 w/o + 15BASW lattice']['univ'],
M__13 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
L__13 = cells['Fuel 1.6 w/o + shutC lattice']['univ'],
K__13 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
J__13 = cells['Fuel 1.6 w/o + shutB lattice']['univ'],
H__13 = cells['Fuel 2.4 w/o + instr + 16BA lattice']['univ'],
G__13 = cells['Fuel 1.6 w/o + instr + shutB lattice']['univ'],
F__13 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
E__13 = cells['Fuel 1.6 w/o + instr + shutD lattice']['univ'],
D__13 = cells['Fuel 2.4 w/o + 16BA lattice']['univ'],
C__13 = cells['Fuel 3.1 w/o + 15BASE lattice']['univ'],
B__13 = cells['Fuel 3.1 w/o lattice']['univ'],
N__14 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
M__14 = cells['Fuel 3.1 w/o + shutA lattice']['univ'],
L__14 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
K__14 = cells['Fuel 1.6 w/o + CRB lattice']['univ'],
J__14 = cells['Fuel 3.1 w/o + 20BA lattice']['univ'],
H__14 = cells['Fuel 1.6 w/o + instr + CRC lattice']['univ'],
G__14 = cells['Fuel 3.1 w/o + 20BA lattice']['univ'], #really should be 23BA1P
F__14 = cells['Fuel 1.6 w/o + instr + CRB lattice']['univ'],
E__14 = cells['Fuel 3.1 w/o + 16BA lattice']['univ'],
D__14 = cells['Fuel 3.1 w/o + shutA lattice']['univ'],
C__14 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
L__15 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
K__15 = cells['Fuel 3.1 w/o + 6BAS lattice']['univ'],
J__15 = cells['Fuel 3.1 w/o + instr lattice']['univ'],
H__15 = cells['Fuel 3.1 w/o + 6BAS lattice']['univ'],
G__15 = cells['Fuel 3.1 w/o lattice']['univ'], # Sp spare location
F__15 = cells['Fuel 3.1 w/o + 6BAS lattice']['univ'],
E__15 = cells['Fuel 3.1 w/o lattice']['univ'],
**baffle)}



  ################## universe 0 cells ##################

  # the axial pincell universes contained in the lattices include the nozzles and bot support plate
  cells['inside core barrel'] ={ 'order':  inc_order(co),
                                'section': comm_t.format("Main universe cells"),
                                'comm':    comm_t.format("inside core barrel"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'fill':    latts['Main Core']['id'],
                                'surfs':  '-{0} {1} -{2}'.format(surfs['core barrel IR']['id'],
                                                                 surfs['lower bound']['id'],
                                                                 surfs['upper bound']['id'])}
  cells['core barrel'] =      { 'order':   inc_order(co),
                                'comm':    comm_t.format("core barrel"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['SS304']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} {2} -{3}'.format(surfs['core barrel IR']['id'],surfs['core barrel OR']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
                                                                     
  cells['shield panel NW'] =   { 'order':   inc_order(co),
                                'comm':    comm_t.format("neutron shield panel NW"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['SS304']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} {2} -{3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NWbot SEtop']['id'],surfs['neut shield NWtop SEbot']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['shield panel N'] =   { 'order':   inc_order(co),
                                'comm':    "",
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['water-nominal']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} {2} -{3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NWtop SEbot']['id'],surfs['neut shield NEtop SWbot']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['shield panel SE'] = { 'order':   inc_order(co),
                                'comm':    comm_t.format("neutron shield panel SE"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['SS304']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} -{2} {3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NWbot SEtop']['id'],surfs['neut shield NWtop SEbot']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['shield panel E'] =   { 'order':   inc_order(co),
                                'comm':    "",
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['water-nominal']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} {2} {3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NWbot SEtop']['id'],surfs['neut shield NEbot SWtop']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['shield panel NE'] =  { 'order':   inc_order(co),
                                'comm':    comm_t.format("neutron shield panel NE"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['SS304']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} {2} -{3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NEbot SWtop']['id'],surfs['neut shield NEtop SWbot']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['shield panel S'] =   { 'order':   inc_order(co),
                                'comm':    "",
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['water-nominal']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} -{2} {3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NWtop SEbot']['id'],surfs['neut shield NEtop SWbot']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['shield panel SW'] =  { 'order':   inc_order(co),
                                'comm':    comm_t.format("neutron shield panel SW"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['SS304']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} -{2} {3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NEbot SWtop']['id'],surfs['neut shield NEtop SWbot']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['shield panel W'] =   { 'order':   inc_order(co),
                                'comm':    "",
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['water-nominal']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} -{2} -{3} {4} -{5}'.format(surfs['core barrel OR']['id'],surfs['neut shield OR']['id'],
                                                                     surfs['neut shield NWbot SEtop']['id'],surfs['neut shield NEbot SWtop']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}                                                                     
                                                                   
                                                                     
  cells['downcomer'] =        { 'order':   inc_order(co),
                                'comm':    comm_t.format("downcomer"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['water-nominal']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} {2} -{3}'.format(surfs['neut shield OR']['id'],surfs['RPV IR']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}
  cells['rpv'] =              { 'order':   inc_order(co),
                                'comm':    comm_t.format("pressure vessel"),
                                'id':      new_id(cellIDs),
                                'univ':    0,
                                'mat':     mats['carbon steel']['id'],
                                'fill':    None,
                                'surfs':  '{0} -{1} {2} -{3}'.format(surfs['RPV IR']['id'],surfs['RPV OR']['id'],
                                                                     surfs['lower bound']['id'],
                                                                     surfs['upper bound']['id'])}


  return mats,surfs,cells,latts


def write_geometry(surfs,cells,lats,outFile):

  surfItems= surfs.items()
  surfItems.sort(key=lambda d: d[1]['order'])
  cellItems= cells.items()
  cellItems.sort(key=lambda d: d[1]['order'])
  lattItems= lats.items()
  lattItems.sort(key=lambda d: d[1]['order'])

  outStr = """<?xml version="1.0" encoding="UTF-8"?>\n"""
  outStr += "<geometry>\n"
  outStr += "\n"
  outStr += comm_t.format("This file auto-generated by beavrs.py")
  outStr += "\n"
  for surf,surfDict in surfItems:
    if 'dupe' in surfDict: continue
    if 'section' in surfDict:
      outStr += "\n"
      outStr += surfDict['section'] + '\n'
    if 'bc' in surfDict:
      outStr += surf_t_bc.format(**surfDict)
    else:
      outStr += surf_t.format(**surfDict)
  outStr += "\n"
  for cell,cellDict in cellItems:
    if 'dupe' in surfDict: continue
    if 'section' in cellDict:
      outStr += "\n"
      outStr += cellDict['section'] + '\n'
    if cellDict['fill']:
      outStr += cell_t_fill.format(**cellDict)
    else:
      outStr += cell_t_mat.format(**cellDict)
  outStr += "\n"
  for latt,latDict in lattItems:
    if 'dupe' in surfDict: continue
    outStr += "\n"
    outStr += latt_t.format(**latDict)
  outStr += "\n"
  outStr += "</geometry>"

  # write file
  with open(outFile,'w') as fh:
    fh.write(outStr)

def make_pin(name,section,co,cells,univ,cellIDs,radSurfList,matList,
             grid=False, surfs=None, gridMat=None, gridType=None):
  """Populates the cells dictionary with a radial pincell universe
  
    name       - string, name of the cell as well as  outKeys = [name] the comment to be written
    section    - section comment - not written if left as ""
    co         - cells order dictionary
    cells      - cells dictionary
    univ       - new universe id number for this pincell
    cellIDs    - set of already used cell IDs
    radSurfList- list of surface IDs, from inside to outside
    matList    - list of tuples of material IDs and comments, from inside to outside.
                 length should be one more than radSurfList.
                 first material will be the inside of the first surface,
                 last last will be outside the last surface.
                 if comments are "", no comment will be written

    Since box surfaces aren't fully implemented yet, square grids around the
    outside of pins need to be added manually here.

    grid       - flag for whether or not to make the grid
    surfs      - full surfaces dictionary.  This routine uses the 'rod grid box'
                 surfaces
    gridMat    - material for the grid
    gridType   - either 'i' or 'tb' for intermediate or top/bottom
    
  """
  
  cells[name] = { 'order':   inc_order(co),
                  'comm':    "",
                  'id':      new_id(cellIDs),
                  'univ':    univ,
                  'mat':     matList[0][0],
                  'fill':    None,
                  'surfs':   '-{0}'.format(radSurfList[0])}
  if section != "":
    cells[name]['section'] = comm_t.format(section)
  if matList[0][1] != "":
    cells[name]['comm'] = comm_t.format(matList[0][1])

  for i,(matIDcomment,outSurf) in enumerate(zip(matList[1:-1],radSurfList[:-1])):
    cells[name + str(i)] = {  'order':   inc_order(co),
                              'comm':    "",
                              'id':      new_id(cellIDs),
                              'univ':    univ,
                              'mat':     matIDcomment[0],
                              'fill':    None,
                              'surfs':  '{0} -{1}'.format(outSurf,radSurfList[i+1])}
    if matIDcomment[1] != "":
      cells[name + str(i)]['comm'] = comm_t.format(matIDcomment[1])
      
  cells[name + 'last'] = {  'order':   inc_order(co),
                            'comm':    "",
                            'id':      new_id(cellIDs),
                            'univ':    univ,
                            'mat':     matList[-1][0],
                            'fill':    None,
                            'surfs':  '{0}'.format(radSurfList[-1])}


  if grid:
  
    # if box-type surfaces were implemented in openmc, this entire thing could be
    # skipped, and instead the grid material and box surface could be put in
    # radSurfList and matList
    
    if not surfs or not gridMat: raise Exception('need surfs and gridMat for grid')


    cells[name + 'last']['surfs'] = '{0} {1} -{2} {3} -{4}'.format(radSurfList[-1],
                                                                   surfs['rod grid box ybot {0}'.format(gridType)]['id'],
                                                                   surfs['rod grid box ytop {0}'.format(gridType)]['id'],
                                                                   surfs['rod grid box xbot {0}'.format(gridType)]['id'],
                                                                   surfs['rod grid box xtop {0}'.format(gridType)]['id'])
    
    cells[name + ' grid 1'] = {   'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    univ,
                                  'mat':     gridMat,
                                  'fill':    None,
                                  'surfs':  '-{0}'.format(surfs['rod grid box ybot {0}'.format(gridType)]['id'])}
    cells[name + ' grid 2'] = {   'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    univ,
                                  'mat':     gridMat,
                                  'fill':    None,
                                  'surfs':  '{0}'.format(surfs['rod grid box ytop {0}'.format(gridType)]['id'])}
    cells[name + ' grid 3'] = {   'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    univ,
                                  'mat':     gridMat,
                                  'fill':    None,
                                  'surfs':  '{0} -{1} {2}'.format(surfs['rod grid box xtop {0}'.format(gridType)]['id'],
                                                                  surfs['rod grid box ytop {0}'.format(gridType)]['id'],surfs['rod grid box ybot {0}'.format(gridType)]['id'])}
    cells[name + ' grid 4'] = {   'order':   inc_order(co),
                                  'comm':    "",
                                  'id':      new_id(cellIDs),
                                  'univ':    univ,
                                  'mat':     gridMat,
                                  'fill':    None,
                                  'surfs':  '-{0} -{1} {2}'.format(surfs['rod grid box xbot {0}'.format(gridType)]['id'],
                                                                   surfs['rod grid box ytop {0}'.format(gridType)]['id'],surfs['rod grid box ybot {0}'.format(gridType)]['id'])}

def make_stack(name,co,cells,univ,cellIDs,axSurfList,fillList):
  """Populates the cells dictionary with an axial stack universe
  
    name       - string, name of the cell as well as the comment to be written
    co         - cells order dictionary
    cells      - cells dictionary
    univ       - new universe id number for this stack
    cellIDs    - set of already used cell IDs
    axSurfList - list of surface IDs, from botom to top
    fillList   - list of fill universe IDs, from bottom to top.
                 length should be one more than axSurfList.
                 first fill universe will be below the first surface,
                 last fill universe will be above the last surface
  """
  
  cells[name] = { 'order':   inc_order(co),
                  'comm':    comm_t.format(name),
                  'id':      new_id(cellIDs),
                  'univ':    univ,
                  'fill':    fillList[0],
                  'surfs':  '-{0}'.format(axSurfList[0])}

  for i,(fillU,botSurf) in enumerate(zip(fillList[1:-1],axSurfList[:-1])):
    cells[name + str(i)] = {  'order':   inc_order(co),
                              'comm':    "",
                              'id':      new_id(cellIDs),
                              'univ':    univ,
                              'fill':    fillU,
                              'surfs':  '{0} -{1}'.format(botSurf,axSurfList[i+1])}
  cells[name + 'last'] = {  'order':   inc_order(co),
                            'comm':    "",
                            'id':      new_id(cellIDs),
                            'univ':    univ,
                            'fill':    fillList[-1],
                            'surfs':  '{0}'.format(axSurfList[-1])}


def make_assembly(latts, cells, surfs, lo, co, univIDs, cellIDs, water, name,
                  comm=None,sect=None, dim=None,lleft=None,width=None,univs=None,
                  gridSurfs=None,sleeveMats=None):
  """Populates the cells and latts dictionary with an assembly

    The cell universe handle to use will be:  cells[name+' lattice']['univ']
  
    cells      - cells dictionary
    surfs      - surfs dictionary
    lo         - latts order dictionary
    co         - cells order dictionary
    univIDs    - set of already used universe IDs
    cellIDs    - set of already used cell IDs
    name       - string, name of the latt/cell family
    water      - material to put outside the lattice
    comm       - optional comment for the lattice and cell family
    sect       - optional section comment
    dim        - required lattice dimension
    lleft      - required lattice lower_left
    width      - required lattice width
    univs      - required lattice universe string.  Should be made with pinLattice_t
    gridSurfs  - required list of grid surface ids, from bottom to top
    sleeveMats - required materials for grid sleeves as [topbot,intermediate]

    returns list of all the created cell keys

  """

  name = str(name)

  # first make the lattice

  latts[name] =                 { 'order':   inc_order(lo),
                                  'id':      new_id(univIDs),
                                  'type':    'rectangular'}
  if comm:
    latts[name]['comm'] = comm_t.format(comm)
  else:
    latts[name]['comm'] = ""
  if sect:
    latts[name]['section'] = comm_t.format(sect)

  for key in ['dim','lleft','width','univs']:
    if not locals()[key]:
      raise Exception('make_assembly requires {0}'.format(key))
    else:
      latts[name][key] = locals()[key]
  
  # add lattice to bounding cell
  cells[name+' lattice'] =    { 'order':   inc_order(co),
                                'id':      new_id(cellIDs),
                                'univ':    new_id(univIDs),
                                'fill':    latts[name]['id'],
                                'surfs':  '-{0} {1} -{2} {3}'.format(surfs['lat box xtop']['id'],surfs['lat box xbot']['id'],surfs['lat box ytop']['id'],surfs['lat box ybot']['id'])}
  if comm:
    cells[name+' lattice']['comm'] = comm_t.format(comm)
  else:
    cells[name+' lattice']['comm'] = ""
  if sect:
    cells[name+' lattice']['section'] = comm_t.format(sect)
  
  
  # make axial all cells for outside of assembly

  # !! all of this would be greatly simplified if box-type surfaces were implemented in openmc !!
  
  outKeys = []  # we only keep track of this to facilitate certain plots/tallies
  
  # first bottom part
  cells[name+' lattice 2'] =  { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '-{0} -{1}'.format(surfs['lat box ybot']['id'],gridSurfs[0])}
  outKeys.append(cells[name+' lattice 2']['id'])
  cells[name+' lattice 3'] =  { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '{0} -{1}'.format(surfs['lat box ytop']['id'],gridSurfs[0])}
  outKeys.append(cells[name+' lattice 3']['id'])
  cells[name+' lattice 4'] =  { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '{0} -{1} {2} -{3}'.format(surfs['lat box xtop']['id'],
                                                                surfs['lat box ytop']['id'],surfs['lat box ybot']['id'],
                                                                gridSurfs[0])}
  outKeys.append(cells[name+' lattice 4']['id'])
  cells[name+' lattice 5'] =  { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '-{0} -{1} {2} -{3}'.format(surfs['lat box xbot']['id'],
                                                                 surfs['lat box ytop']['id'],surfs['lat box ybot']['id'],
                                                                 gridSurfs[0])}
  outKeys.append(cells[name+' lattice 5']['id'])

  gridnum = 0
  
  # all middle cells
  for i,botGridSurf in enumerate(gridSurfs[:-1]):
  
  
    if i%2 == 0:
      # make gridsleeve cells

      gridnum += 1
      
      if gridnum == 1 or gridnum == 8:
        smat = sleeveMats[0]
      else:
        smat = sleeveMats[1]
      
      outSurfsKey = 'lat grid box '
      inSurfsKey = 'lat box '
           
      cellKey = name+' lattice grid {0}'.format(6+i*4)
      cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                          'id':      new_id(cellIDs),
                          'univ':    cells[name+' lattice']['univ'],
                          'fill':    None,
                          'mat':     smat,
                          'surfs':  '-{0} {1} {2} -{3} {4} -{5}'.format(surfs[inSurfsKey+'ybot']['id'],surfs[outSurfsKey+'ybot']['id'],
                                                                        surfs[outSurfsKey+'xbot']['id'],surfs[outSurfsKey+'xtop']['id'],
                                                                        botGridSurf,gridSurfs[i+1])}
      outKeys.append(cells[cellKey]['id'])
      cellKey = name+' lattice grid {0}'.format(7+i*4)
      cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                          'id':      new_id(cellIDs),
                          'univ':    cells[name+' lattice']['univ'],
                          'fill':    None,
                          'mat':     smat,
                          'surfs':  '-{0} {1} {2} -{3} {4} -{5}'.format(surfs[outSurfsKey+'ytop']['id'],surfs[inSurfsKey+'ytop']['id'],
                                                                        surfs[outSurfsKey+'xbot']['id'],surfs[outSurfsKey+'xtop']['id'],
                                                                        botGridSurf,gridSurfs[i+1])}
      outKeys.append(cells[cellKey]['id'])
      cellKey = name+' lattice grid {0}'.format(8+i*4)
      cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                          'id':      new_id(cellIDs),
                          'univ':    cells[name+' lattice']['univ'],
                          'fill':    None,
                          'mat':     smat,
                          'surfs':  '-{0} {1} {2} -{3} {4} -{5}'.format(surfs[inSurfsKey+'xbot']['id'],surfs[outSurfsKey+'xbot']['id'],
                                                                        surfs[inSurfsKey+'ybot']['id'],surfs[inSurfsKey+'ytop']['id'],
                                                                        botGridSurf,gridSurfs[i+1])}
      outKeys.append(cells[cellKey]['id'])
      cellKey = name+' lattice grid {0}'.format(9+i*4)
      cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                          'id':      new_id(cellIDs),
                          'univ':    cells[name+' lattice']['univ'],
                          'fill':    None,
                          'mat':     smat,
                          'surfs':  '-{0} {1} {2} -{3} {4} -{5}'.format(surfs[outSurfsKey+'xtop']['id'],surfs[inSurfsKey+'xtop']['id'],
                                                                        surfs[inSurfsKey+'ybot']['id'],surfs[inSurfsKey+'ytop']['id'],
                                                                        botGridSurf,gridSurfs[i+1])}
      outKeys.append(cells[cellKey]['id'])
      
    else:
      outSurfsKey = 'lat box '
  
  
    cellKey = name+' lattice {0}'.format(6+i*4)
    cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                        'id':      new_id(cellIDs),
                        'univ':    cells[name+' lattice']['univ'],
                        'fill':    None,
                        'mat':     water,
                        'surfs':  '-{0} {1} -{2}'.format(surfs[outSurfsKey+'ybot']['id'],botGridSurf,gridSurfs[i+1])}
    outKeys.append(cells[cellKey]['id'])
    cellKey = name+' lattice {0}'.format(7+i*4)
    cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                        'id':      new_id(cellIDs),
                        'univ':    cells[name+' lattice']['univ'],
                        'fill':    None,
                        'mat':     water,
                        'surfs':  '{0} {1} -{2}'.format(surfs[outSurfsKey+'ytop']['id'],botGridSurf,gridSurfs[i+1])}
    outKeys.append(cells[cellKey]['id'])
    cellKey = name+' lattice {0}'.format(8+i*4)
    cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                        'id':      new_id(cellIDs),
                        'univ':    cells[name+' lattice']['univ'],
                        'fill':    None,
                        'mat':     water,
                        'surfs':  '{0} -{1} {2} {3} -{4}'.format(surfs[outSurfsKey+'xtop']['id'],
                                                        surfs[outSurfsKey+'ytop']['id'],surfs[outSurfsKey+'ybot']['id'],
                                                        botGridSurf,gridSurfs[i+1])}
    outKeys.append(cells[cellKey]['id'])
    cellKey = name+' lattice {0}'.format(9+i*4)
    cells[cellKey] =  { 'order':   inc_order(co), 'comm':"",
                        'id':      new_id(cellIDs),
                        'univ':    cells[name+' lattice']['univ'],
                        'fill':    None,
                        'mat':     water,
                        'surfs':  '-{0} -{1} {2} {3} -{4}'.format(surfs[outSurfsKey+'xbot']['id'],
                                                         surfs[outSurfsKey+'ytop']['id'],surfs[outSurfsKey+'ybot']['id'],
                                                         botGridSurf,gridSurfs[i+1])}
    outKeys.append(cells[cellKey]['id'])
    
  # top part
  cells[name+' lat last 1'] = { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '-{0} {1}'.format(surfs['lat box ybot']['id'],gridSurfs[-1])}
  outKeys.append(cells[name+' lat last 1']['id'])
  cells[name+' lat last 2'] = { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '{0} {1}'.format(surfs['lat box ytop']['id'],gridSurfs[-1])}
  outKeys.append(cells[name+' lat last 2']['id'])
  cells[name+' lat last 3'] = { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '{0} -{1} {2} {3}'.format(surfs['lat box xtop']['id'],
                                                                surfs['lat box ytop']['id'],surfs['lat box ybot']['id'],
                                                                gridSurfs[-1])}
  outKeys.append(cells[name+' lat last 3']['id'])
  cells[name+' lat last 4'] = { 'order':   inc_order(co), 'comm':"",
                                'id':      new_id(cellIDs),
                                'univ':    cells[name+' lattice']['univ'],
                                'fill':    None,
                                'mat':     water,
                                'surfs':  '-{0} -{1} {2} {3}'.format(surfs['lat box xbot']['id'],
                                                                 surfs['lat box ytop']['id'],surfs['lat box ybot']['id'],
                                                                 gridSurfs[-1])}
  outKeys.append(cells[name+' lat last 4']['id'])
  
  
  
  return outKeys # we only keep track of this to facilitate certain plots/tallies


def control_bank_axials(step):
  """Given the total number of steps withdrawn, returns a dictionary of control rod axial surfaces

  Starting from all out, we follow this sequence:
  First D moves in alone, until it gets to 113 steps withdrawn
  Now D and C move together until C gets to 113 steps withdrawn (D is all the way in by C at 115)
  Now C and B move together until B gets to 113 steps withdrawn (C is all the way in by B at 115)
  Now B and A move together until A gets to 0 steps withdrawn (B is all the way in by A at 115)

  Assuming only movement of each control rod bank by one step, in total this sequence yields 574 unique positions,
  which we specify with step.  If step=1, all control rods are out of the core, and if step=574, all are fully inserted.

  For example for BOL normal operation if you want the D bank at the bite position of 184 steps withdrawn,
  set step=228-184+1=45

  """

  if step < 0 or step > 574:
    raise Exception("Invalid control bank step specification {0}".format(step))

  bankDstep = max(0.0,228-step)
  bankCstep = (bankDstep<113)*max(0.0,228-step+113  +3) + (bankDstep>=113)*228
  bankBstep = (bankCstep<113)*max(0.0,228-step+113*2+5) + (bankCstep>=113)*228
  bankAstep = (bankBstep<113)*max(0.0,228-step+113*3+7) + (bankBstep>=113)*228

  bankAbot = step0H + stepWidth*bankAstep
  bankBbot = step0H + stepWidth*bankBstep
  bankCbot = step0H + stepWidth*bankCstep
  bankDbot = step0H + stepWidth*bankDstep

  bankAtop = bankAbot + stepWidth*228
  bankBtop = bankBbot + stepWidth*228
  bankCtop = bankCbot + stepWidth*228
  bankDtop = bankDbot + stepWidth*228

  return locals()


def new_id(set_,range_=None):
  """Returns a new id unique to set_, within desired range_, and adds it to that set

      range_ should be a tuple or list where range_[0] = lower bound, inclusive
                                             range_[1] = upper bound, inclusive
  """

  if not range_:
    range_ = [1,1000000]

  id_ = range_[0]

  while id_ in set_:
    id_ += 1
    if id_ > range_[1]:
      raise Exception("no more unique ids in range {0}".format(range_))

  set_.add(id_)

  return id_


def inc_order(o):
  o['n'] += 1
  return o['n']


def main():

  mats,surfs,cells,latts = init_data()

  write_geometry(surfs,cells,latts,"geometry.xml")


if __name__ == "__main__":
  main()
