import h5py
import numpy


###############################################################################
# This file writes all of the materials data (multi-group nuclear
# cross-sections) for the OECD's C5G7 deterministic neutron transport
# benchmark problem to an HDF5 file. The script uses the h5py Python package
# to interact with the HDF5 file format. This may be a good example for those
# wishing ot write their nuclear data to an HDF5 file to import using the
# OpenMOC 'materialize' Python module.
###############################################################################


# Create the file to store C5G7 multi-groups cross-sections
f = h5py.File('c5g7-mgxs.h5')
f.attrs["# groups"] = 7

# Create a group to specify that MGXS are split by material (vs. cell)
material_group = f.create_group('material')


###############################################################################
################################      UO2      ################################
###############################################################################

# Create a subgroup for UO2 materials data
uo2 = material_group.create_group('UO2')

sigma_t = numpy.array([1.779490E-01, 3.298050E-01, 4.803880E-01,
                       5.543670E-01, 3.118010E-01, 3.951680E-01, 5.644060E-01])
sigma_s = numpy.array([1.275370E-01, 4.237800E-02, 9.437400E-06,
                       5.516300E-09, 0., 0., 0., 0., 3.244560E-01,
                       1.631400E-03, 3.142700E-09, 0., 0., 0., 0.,
                       0., 4.509400E-01, 2.679200E-03, 0., 0., 0.,
                       0., 0., 0., 4.525650E-01, 5.566400E-03, 0.,
                       0., 0., 0., 0., 1.252500E-04, 2.714010E-01,
                       1.025500E-02, 1.002100E-08, 0., 0., 0., 0.,
                       1.296800E-03, 2.658020E-01, 1.680900E-02,
                       0., 0., 0., 0., 0., 8.545800E-03, 2.730800E-01])

sigma_f = numpy.array([7.212060E-03, 8.193010E-04, 6.453200E-03,
                       1.856480E-02, 1.780840E-02, 8.303480E-02, 2.160040E-01])
nu_sigma_f = numpy.array([2.005998E-02, 2.027303E-03, 1.570599E-02,
                          4.518301E-02, 4.334208E-02, 2.020901E-01,
                          5.257105E-01])
chi = numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04,
                   1.17610E-07, 0., 0., 0.])

# Create datasets for each cross-section type
uo2.create_dataset('total', data=sigma_t)
uo2.create_dataset('scatter matrix', data=sigma_s)
uo2.create_dataset('fission', data=sigma_f)
uo2.create_dataset('nu-fission', data=nu_sigma_f)
uo2.create_dataset('chi', data=chi)



###############################################################################
##############################      MOX (4.3%)     ############################
###############################################################################

# Create a subgroup for MOX-4.3%  materials data
mox43 = material_group.create_group('MOX-4.3%')

sigma_t = numpy.array([1.787310E-01, 3.308490E-01, 4.837720E-01,
                       5.669220E-01, 4.262270E-01, 6.789970E-01, 6.828520E-01])
sigma_s = numpy.array([1.288760E-01, 4.141300E-02, 8.229000E-06,
                       5.040500E-09, 0., 0., 0., 0., 3.254520E-01, 1.639500E-03,
                       1.598200E-09, 0., 0., 0., 0., 0., 4.531880E-01,
                       2.614200E-03, 0., 0., 0., 0., 0., 0., 4.571730E-01,
                       5.539400E-03, 0., 0., 0., 0., 0., 1.604600E-04,
                       2.768140E-01, 9.312700E-03, 9.165600E-09, 0., 0., 0.,
                       0., 2.005100E-03, 2.529620E-01, 1.485000E-02, 0., 0.,
                       0., 0., 0., 8.494800E-03, 2.650070E-01])
sigma_f = numpy.array([7.62704E-03, 8.76898E-04, 5.69835E-03,
                       2.28872E-02, 1.07635E-02, 2.32757E-01, 2.48968E-01])
nu_sigma_f = numpy.array([2.175300E-02, 2.535103E-03, 1.626799E-02,
                          6.547410E-02, 3.072409E-02, 6.666510E-01,
                          7.139904E-01])
chi = numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07,
                          0., 0., 0.])

# Create datasets for each cross-section type
mox43.create_dataset('total', data=sigma_t)
mox43.create_dataset('scatter matrix', data=sigma_s)
mox43.create_dataset('fission', data=sigma_f)
mox43.create_dataset('nu-fission', data=nu_sigma_f)
mox43.create_dataset('chi', data=chi)


###############################################################################
##############################      MOX (7%)     ##############################
###############################################################################

# Create a subgroup for MOX-7% materials data
mox7 = material_group.create_group('MOX-7%')

sigma_t = numpy.array([1.813230E-01, 3.343680E-01, 4.937850E-01,
                       5.912160E-01, 4.741980E-01, 8.336010E-01, 8.536030E-01])
sigma_s = numpy.array([1.304570E-01, 4.179200E-02, 8.510500E-06,
                       5.132900E-09, 0., 0., 0., 0., 3.284280E-01, 1.643600E-03,
                       2.201700E-09, 0., 0., 0., 0., 0., 4.583710E-01,
                       2.533100E-03, 0., 0., 0., 0., 0., 0., 4.637090E-01,
                       5.476600E-03, 0., 0., 0., 0., 0., 1.761900E-04,
                       2.823130E-01, 8.728900E-03, 9.001600E-09, 0., 0., 0.,
                       0., 2.276000E-03, 2.497510E-01, 1.311400E-02, 0., 0.,
                       0., 0., 0., 8.864500E-03, 2.595290E-01])
sigma_f = numpy.array([8.25446E-03, 1.32565E-03, 8.42156E-03,
                       3.28730E-02, 1.59636E-02, 3.23794E-01, 3.62803E-01])
nu_sigma_f = numpy.array([2.381395E-02, 3.858689E-03, 2.413400E-02,
                          9.436622E-02, 4.576988E-02, 9.281814E-01,
                          1.043200E+00])
chi = numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07,
                   0., 0., 0.])

# Create datasets for each cross-section type
mox7.create_dataset('total', data=sigma_t)
mox7.create_dataset('scatter matrix', data=sigma_s)
mox7.create_dataset('fission', data=sigma_f)
mox7.create_dataset('nu-fission', data=nu_sigma_f)
mox7.create_dataset('chi', data=chi)


###############################################################################
##############################      MOX (8.7%)     ############################
###############################################################################

# Create a subgroup for MOX-8.7% materials data
mox87 = material_group.create_group('MOX-8.7%')

sigma_t = numpy.array([1.830450E-01, 3.367050E-01, 5.005070E-01,
                       6.061740E-01, 5.027540E-01, 9.210280E-01, 9.552310E-01])
sigma_s = numpy.array([1.315040E-01, 4.204600E-02, 8.697200E-06,
                    5.193800E-09, 0., 0., 0., 0., 3.304030E-01, 1.646300E-03,
                       2.600600E-09, 0., 0., 0., 0., 0., 4.617920E-01,
                       2.474900E-03, 0., 0., 0., 0., 0., 0., 4.680210E-01,
                       5.433000E-03, 0., 0., 0., 0., 0., 1.859700E-04,
                       2.857710E-01, 8.397300E-03, 8.928000E-09, 0., 0.,
                       0., 0., 2.391600E-03, 2.476140E-01, 1.232200E-02,
                       0., 0., 0., 0., 0., 8.968100E-03, 2.560930E-01])
sigma_f = numpy.array([8.67209E-03, 1.62426E-03, 1.02716E-02,
                       3.90447E-02, 1.92576E-02, 3.74888E-01, 4.30599E-01])
nu_sigma_f = numpy.array([2.518600E-02, 4.739509E-03, 2.947805E-02,
                          1.122500E-01, 5.530301E-02, 1.074999E+00,
                          1.239298E+00])
chi = numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07,
                   0., 0., 0.])

# Create datasets for each cross-section type
mox87.create_dataset('total', data=sigma_t)
mox87.create_dataset('scatter matrix', data=sigma_s)
mox87.create_dataset('fission', data=sigma_f)
mox87.create_dataset('nu-fission', data=nu_sigma_f)
mox87.create_dataset('chi', data=chi)


###############################################################################
############################      Fission Chamber     #########################
###############################################################################

# Create a subgroup for fission chamber materials data
fiss_chamber = material_group.create_group('Fission Chamber')

sigma_t = numpy.array([1.260320E-01, 2.931600E-01, 2.842500E-01,
                       2.810200E-01, 3.344600E-01, 5.656400E-01,
                       1.172140E+00])
sigma_s = numpy.array([6.616590E-02, 5.907000E-02, 2.833400E-04,
                       1.462200E-06, 2.064200E-08, 0., 0., 0.,
                       2.403770E-01, 5.243500E-02, 2.499000E-04,
                       1.923900E-05, 2.987500E-06, 4.214000E-07,
                       0., 0., 1.834250E-01, 9.228800E-02, 6.936500E-03,
                       1.079000E-03, 2.054300E-04, 0., 0., 0.,
                       7.907690E-02, 1.699900E-01, 2.586000E-02,
                       4.925600E-03, 0., 0., 0., 3.734000E-05,
                       9.975700E-02, 2.067900E-01, 2.447800E-02, 0.,
                       0., 0., 0., 9.174200E-04, 3.167740E-01,
                       2.387600E-01, 0., 0., 0., 0., 0., 4.979300E-02,
                       1.09910E+00])
sigma_f = numpy.array([4.79002E-09, 5.82564E-09, 4.63719E-07,
                       5.24406E-06, 1.45390E-07, 7.14972E-07, 2.08041E-06])
nu_sigma_f = numpy.array([1.323401E-08, 1.434500E-08, 1.128599E-06,
                          1.276299E-05, 3.538502E-07, 1.740099E-06,
                          5.063302E-06])
chi = numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04,
                   1.17610E-07, 0., 0., 0.])

# Create datasets for each cross-section type
fiss_chamber.create_dataset('total', data=sigma_t)
fiss_chamber.create_dataset('scatter matrix', data=sigma_s)
fiss_chamber.create_dataset('fission', data=sigma_f)
fiss_chamber.create_dataset('nu-fission', data=nu_sigma_f)
fiss_chamber.create_dataset('chi', data=chi)


###############################################################################
##############################      Guide Tube      ###########################
###############################################################################

# Create a subgroup for guide tube materials data
guide_tube = material_group.create_group('Guide Tube')

sigma_t = numpy.array([1.260320E-01, 2.931600E-01, 2.842400E-01,
                       2.809600E-01, 3.344400E-01, 5.656400E-01,
                       1.172150E+00])
sigma_s = numpy.array([6.616590E-02, 5.907000E-02, 2.833400E-04,
                       1.462200E-06, 2.064200E-08, 0., 0., 0., 2.403770E-01,
                       5.243500E-02, 2.499000E-04, 1.923900E-05,
                       2.987500E-06, 4.214000E-07, 0., 0., 1.832970E-01,
                       9.239700E-02, 6.944600E-03, 1.0803000E-03,
                       2.056700E-04, 0., 0., 0., 7.885110E-02, 1.701400E-01,
                       2.588100E-02, 4.929700E-03, 0., 0., 0., 3.733300E-05,
                       9.973720E-02, 2.067900E-01, 2.447800E-02, 0., 0., 0.,
                       0., 9.172600E-04, 3.167650E-01, 2.387700E-01, 0., 0.,
                       0., 0., 0., 4.979200E-02, 1.099120E+00])

sigma_f = numpy.zeros(7)
nu_sigma_f = numpy.zeros(7)
chi = numpy.zeros(7)

# Create datasets for each cross-section type
guide_tube.create_dataset('total', data=sigma_t)
guide_tube.create_dataset('scatter matrix', data=sigma_s)
guide_tube.create_dataset('fission', data=sigma_f)
guide_tube.create_dataset('nu-fission', data=nu_sigma_f)
guide_tube.create_dataset('chi', data=chi)


###############################################################################
################################      Water      ##############################
###############################################################################

# Create a subgroup for water materials data
water = material_group.create_group('Water')

sigma_t = numpy.array([1.592060E-01, 4.129700E-01, 5.903100E-01,
                       5.843500E-01, 7.180000E-01, 1.254450E+00, 2.650380E+00])
sigma_s = numpy.array([4.447770E-02, 1.134000E-01, 7.234700E-04,
                       3.749900E-06, 5.318400E-08, 0., 0., 0., 2.823340E-01,
                       1.299400E-01, 6.234000E-04, 4.800200E-05, 7.448600E-06,
                       1.045500E-06, 0., 0., 3.452560E-01, 2.245700E-01,
                       1.699900E-02, 2.644300E-03, 5.034400E-04, 0., 0., 0.,
                       9.102840E-02, 4.155100E-01, 6.373200E-02, 1.213900E-02,
                       0., 0., 0., 7.143700E-05, 1.391380E-01, 5.118200E-01,
                       6.122900E-02, 0., 0., 0., 0., 2.215700E-03, 6.999130E-01,
                       5.373200E-01, 0., 0., 0., 0., 0., 1.324400E-01,
                       2.480700E+00])
sigma_f = numpy.zeros(7)
nu_sigma_f = numpy.zeros(7)
chi = numpy.zeros(7)

# Create datasets for each cross-section type
water.create_dataset('total', data=sigma_t)
water.create_dataset('scatter matrix', data=sigma_s)
water.create_dataset('fission', data=sigma_f)
water.create_dataset('nu-fission', data=nu_sigma_f)
water.create_dataset('chi', data=chi)


###############################################################################
################################   Control Rod   ##############################
###############################################################################

# Create a subgroup for control rod materials data
control_rod = material_group.create_group('Control Rod')

sigma_t = numpy.array([2.16768E-01, 4.80098E-01, 8.86369E-01,
                       9.70009E-01, 9.10482E-01, 1.13775E+00,
                       1.84048E+00])
sigma_s = numpy.array([1.70563E-01, 4.44012E-02, 9.83670E-05,
                       1.27786E-07, 0., 0., 0., 0., 4.71050E-01,
                       6.85480E-04, 3.91395E-10, 0., 0.,
                       0., 0., 0., 8.01859E-01, 7.20132E-04,
                       0., 0., 0., 0., 0., 0.,
                       5.70752E-01, 1.46015E-03, 0., 0.,
                       0., 0., 0., 6.55562E-05, 2.07838E-01, 3.81486E-03,
                       3.69760E-09, 0., 0., 0., 0., 1.02427E-03, 2.02465E-01,
                       4.75290E-03, 0., 0., 0., 0., 0., 3.53043E-03,
                       6.58597E-01])
sigma_f = numpy.zeros(7)
nu_sigma_f = numpy.zeros(7)
chi = numpy.zeros(7)

# Create datasets for each cross-section type
control_rod.create_dataset('total', data=sigma_t)
control_rod.create_dataset('scatter matrix', data=sigma_s)
control_rod.create_dataset('fission', data=sigma_f)
control_rod.create_dataset('nu-fission', data=nu_sigma_f)
control_rod.create_dataset('chi', data=chi)


###############################################################################
################################      Clad       ##############################
###############################################################################

# Create a subgroup for Clad materials data
clad = material_group.create_group('Clad')

sigma_t = numpy.array([1.30060E-01, 3.05480E-01, 3.29910E-01, 
                       2.69700E-01, 2.72780E-01, 2.77940E-01,
                       2.95630E-01])
sigma_s = numpy.array([9.72490E-02, 3.25480E-02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
  0.00000E+00, 3.03980E-01, 7.72850E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
  0.00000E+00, 0.00000E+00, 3.24280E-01, 5.94050E-04, 0.00000E+00, 0.00000E+00, 0.00000E+00,
  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.63200E-01, 5.31350E-03, 0.00000E+00, 0.00000E+00,
  0.00000E+00, 0.00000E+00, 0.00000E+00, 2.12680E-03, 2.53950E-01, 1.39080E-02, 0.00000E+00,
  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.48500E-02, 2.41850E-01, 1.65340E-02,
  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.98990E-02, 2.57160E-01])
sigma_f = numpy.zeros(7)
nu_sigma_f = numpy.zeros(7)
chi = numpy.zeros(7)

# Create datasets for each cross-section type
clad.create_dataset('total', data=sigma_t)
clad.create_dataset('scatter matrix', data=sigma_s)
clad.create_dataset('fission', data=sigma_f)
clad.create_dataset('nu-fission', data=nu_sigma_f)
clad.create_dataset('chi', data=chi)

# Close the hdf5 data file
f.close()
