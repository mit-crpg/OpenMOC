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
f = h5py.File('single-group-materials.h5')
f.attrs["# groups"] = 1

# Create a group to specify that MGXS are split by material (vs. cell)
material_group = f.create_group('material')


###############################################################################
################################      UO2      ################################
###############################################################################

# Create a subgroup for UO2 materials data
'''
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
'''
###############################################################################
#                            Creating Materials
###############################################################################

# Create a 1 group fuel material
single = material_group.create_group('single')

sigma_t = numpy.array([3.7])

sigma_s = numpy.array([0])

sigma_f = numpy.array([1.7])

nu_sigma_f = numpy.array([2])

chi = numpy.array([1])

# Create datasets for each cross-section type
single.create_dataset('total', data=sigma_t)
single.create_dataset('scatter matrix', data=sigma_s)
single.create_dataset('fission', data=sigma_f)
single.create_dataset('nu-fission', data=nu_sigma_f)


# Close the hdf5 data file
f.close()
