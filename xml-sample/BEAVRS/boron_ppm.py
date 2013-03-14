from __future__ import division

## Molecular Masses

# isotopic molar masses
MH1    =   1.0078250
MH2    =   2.0141018
MHe4   =   4.0026032542
MB10   =  10.0129370
MB11   =  11.0093054
MC12   =  12.0000000
MC13   =  13.003354838
MO16   =  15.9949146196
MO17   =  16.9991317
MO18   =  17.999161
MN14   =  14.003074005
MN15   =  15.000108898
MSi28  =  27.976926532
MSi29  =  28.97649470
MSi30  =  29.97377017
MP31   =  30.9737616
MAl27  =  26.9815386
MAr36  =  35.96754511
MAr38  =  37.9627324
MAr40  =  39.962383123
MCr50  =  49.946044
MCr52  =  51.940507
MCr53  =  52.940649
MCr54  =  53.938880
MMn55  =  54.938045
MFe54  =  53.939611
MFe56  =  55.934937
MFe57  =  56.935394
MFe58  =  57.933276
MNi58  =  57.935343
MNi60  =  59.930786
MNi61  =  60.931056
MNi62  =  61.928345
MNi64  =  63.927966
MZr90  =  89.904704
MZr91  =  90.905646
MZr92  =  91.905041
MZr94  =  93.906315
MZr96  =  95.908273
MMo92  =  91.906811
MMo94  =  93.905088
MMo95  =  94.905842
MMo96  =  95.904679
MMo97  =  96.906021
MMo98  =  97.905408
MMo100 =  99.90748
MAg107 = 106.905097
MAg109 = 108.904752
MCd106 = 105.90646
MCd108 = 107.90418
MCd110 = 109.903002
MCd111 = 110.904178
MCd112 = 111.902758
MCd113 = 112.904402
MCd114 = 113.903359
MCd116 = 115.904756
MIn113 = 112.904058
MIn115 = 114.903878
MSn112 = 111.904818
MSn114 = 113.902779
MSn115 = 114.903342
MSn116 = 115.901741
MSn117 = 116.902952
MSn118 = 117.901603
MSn119 = 118.903308
MSn120 = 119.902195
MSn122 = 121.903439
MSn124 = 123.905274
MU234  = 234.040952
MU235  = 235.043930
MU236  = 236.045568
MU238  = 238.050788


# naturally occurring abundances
aH1    = 0.99985
aH2    = 0.00015
aB10   = 0.199
aB11   = 0.801
aC12   = 0.9893
aC13   = 0.0107
aO16   = 0.99757
aO17   = 0.00038
aO18   = 0.00205
aN14   = 0.99636
aN15   = 0.00364
aSi28  = 0.92223
aSi29  = 0.04685
aSi30  = 0.03092
aAr36  = 0.003365
aAr38  = 0.000632
aAr40  = 0.996003
aCr50  = 0.04345
aCr52  = 0.83789
aCr53  = 0.09501
aCr54  = 0.02365
aFe54  = 0.05845
aFe56  = 0.91754
aFe57  = 0.02119
aFe58  = 0.00282
aNi58  = 0.680769
aNi60  = 0.262231
aNi61  = 0.011399
aNi62  = 0.036345
aNi64  = 0.009256
aZr90  = 0.5145
aZr91  = 0.1122
aZr92  = 0.1715
aZr94  = 0.1738
aZr96  = 0.0280
aMo92  = 0.1477
aMo94  = 0.0923
aMo95  = 0.1590
aMo96  = 0.1668
aMo97  = 0.0956
aMo98  = 0.2419
aMo100 = 0.0967
aAg107 = 0.51839
aAg109 = 0.48161
aCd106 = 0.0125
aCd108 = 0.0089
aCd110 = 0.1249
aCd111 = 0.1280
aCd112 = 0.2413
aCd113 = 0.1222
aCd114 = 0.2873
aCd116 = 0.0749
aIn113 = 0.0429
aIn115 = 0.9571
aSn112 = 0.0097
aSn114 = 0.0066
aSn115 = 0.0034
aSn116 = 0.1454
aSn117 = 0.0768
aSn118 = 0.2422
aSn119 = 0.0859
aSn120 = 0.3258
aSn122 = 0.0463
aSn124 = 0.0579

# naturally occurring molar masses
MH  = 1.00794
MB  = 10.811
MC  = 12.0107
MO  = 15.9994
MN  = 14.0067
MSi = 28.0855
MAr = 39.948
MCr = 51.9961
MFe = 55.845
MNi = 58.6934
MZr = 91.224
MMo = 95.94
MAg = 107.8682
MCd = 112.411
MIn = 114.818
MSn = 118.710

NA = 0.60221415 # avagadros number / 10^24

## Borated Water

class Water:
  def __init__(self):
    pass

  def get_water_mats(self,dens,ppm):
    

    # density at 2250 psia and T=560F NIST thermofluid properties
    self.rho = dens

    # ppm Boron (Figure 3.1 NDR)
    self.ppm = ppm

    # compute weight percent of natural boron in borated water
    self.wBph2o = self.ppm * 1e-6

    # compute borated water mass density
    self.rhoBW = dens/(1.0 - self.wBph2o)

    # molecular mass of water
    self.Mh2o = 2*MH + MO

    # compute number density
    self.Nh2o = self.rho * NA / self.Mh2o

    # compute number densities of elements
    self.NB = self.wBph2o * self.rhoBW * NA / MB
    self.NH = 2*self.Nh2o
    self.NO = self.Nh2o

    # compute isoptopic number densities
    self.NB10 = aB10 * self.NB
    self.NB11 = aB11 * self.NB
    self.NH1  = aH1  * self.NH
    self.NH2  = aH2  * self.NH
    self.NO16 = aO16 * self.NO
    self.NO17 = aO17 * self.NO
    self.NO18 = aO18 * self.NO

    # atom fractions
    self.aB10 = self.NB10/self.Nh2o
    self.aB11 = self.NB11/self.Nh2o
    self.aH   = self.NH/(self.NH + self.NO) * self.Nh2o/self.Nh2o # water is a molecule
    self.aO   = self.NO/(self.NH + self.NO) * self.Nh2o/self.Nh2o # water is a molecule
    self.aH1 = aH1 * self.aH
    self.aH2 = aH2 * self.aH
    self.O16 = aO16 * self.aO
    self.O17 = aO17 * self.aO
    self.O18 = aO18 * self.aO
  
  def print_latex(self):

    print 'Water Density: {0:8.6f}'.format(self.rhoBW)
    print 'B-10  & {0:10.4e} \\\\'.format(self.NB10)
    print 'B-11  & {0:10.4e} \\\\'.format(self.NB11)
    print 'H-1  & {0:10.4e} \\\\'.format(self.NH1)
    print 'H-2  & {0:10.4e} \\\\'.format(self.NH2)
    print 'O-16  & {0:10.4e} \\\\'.format(self.NO16)
    print 'O-17  & {0:10.4e} \\\\'.format(self.NO17+self.NO18)
