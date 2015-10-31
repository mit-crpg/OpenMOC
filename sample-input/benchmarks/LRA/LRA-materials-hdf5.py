import h5py
import numpy


###############################################################################
# This file writes all of the materials data (multi-group nuclear 
# cross-sections) for the LRA diffusion
# benchmark problem to an HDF5 file. The script uses the h5py Python package
# to interact with the HDF5 file format. This may be a good example for those
# wishing ot write their nuclear data to an HDF5 file to import using the
# OpenMOC 'materialize' Python module.
###############################################################################


# Create the file to store LRA multi-groups cross-sections
f = h5py.File('LRA-mgxs.h5')
f.attrs["# groups"] = 2

# Create a group to specify that MGXS are split by material (vs. cell)
material_group = f.create_group('material')


###############################################################################
################################   region 1    ################################
###############################################################################

# Create a subgroup for region 1 materials data
region_1 = material_group.create_group('region_1')

sigma_t = numpy.array([0.2656, 1.5798])
sigma_s = numpy.array([0.2318925, 0.02533, 0.00, 1.4794789])
sigma_f = numpy.array([0.002, 0.05])
nu_sigma_f = numpy.array([0.004602, 0.1091])
chi = numpy.array([1.0, 0.0])

region_1.create_dataset('total', data=sigma_t)
region_1.create_dataset('scatter matrix', data=sigma_s)
region_1.create_dataset('fission', data=sigma_f)
region_1.create_dataset('nu-fission', data=nu_sigma_f)
region_1.create_dataset('chi', data=chi)


###############################################################################
################################   region 2    ################################
###############################################################################

# Create a subgroup for region 2 materials data
region_2 = material_group.create_group('region_2')

sigma_t = numpy.array([0.2629, 1.7525])
sigma_s = numpy.array([0.2279222, 0.02767, 0.00, 1.68201098])
sigma_f = numpy.array([0.002, 0.045])
nu_sigma_f = numpy.array([0.004609, 0.08675])
chi = numpy.array([1.0, 0.0])

region_2.create_dataset('total', data=sigma_t)
region_2.create_dataset('scatter matrix', data=sigma_s)
region_2.create_dataset('fission', data=sigma_f)
region_2.create_dataset('nu-fission', data=nu_sigma_f)
region_2.create_dataset('chi', data=chi)


###############################################################################
################################   region 3    ################################
###############################################################################

# Create a subgroup for region 3 materials data
region_3 = material_group.create_group('region_3')

sigma_t = numpy.array([0.2648, 1.5941])
sigma_s = numpy.array([0.2305021, 0.02617, 0.00, 1.51063909])
sigma_f = numpy.array([0.002, 0.045])
nu_sigma_f = numpy.array([0.004663, 0.1021])
chi = numpy.array([1.0, 0.0])

region_3.create_dataset('total', data=sigma_t)
region_3.create_dataset('scatter matrix', data=sigma_s)
region_3.create_dataset('fission', data=sigma_f)
region_3.create_dataset('nu-fission', data=nu_sigma_f)
region_3.create_dataset('chi', data=chi)


###############################################################################
################################   region 4    ################################
###############################################################################

# Create a subgroup for region 4 materials data
region_4 = material_group.create_group('region_4')

sigma_t = numpy.array([0.2648, 1.5941])
sigma_s = numpy.array([0.230462, 0.02617, 0.00, 1.520789])
sigma_f = numpy.array([0.002, 0.045])
nu_sigma_f = numpy.array([0.004663, 0.1021])
chi = numpy.array([1.0, 0.0])

region_4.create_dataset('total', data=sigma_t)
region_4.create_dataset('scatter matrix', data=sigma_s)
region_4.create_dataset('fission', data=sigma_f)
region_4.create_dataset('nu-fission', data=nu_sigma_f)
region_4.create_dataset('chi', data=chi)


###############################################################################
################################   region 5    ################################
###############################################################################

# Create a subgroup for region 5 materials data
region_5 = material_group.create_group('region_5')

sigma_t = numpy.array([0.2648, 1.5941])
sigma_s = numpy.array([0.230462, 0.02617, 0.00, 1.510672])
sigma_f = numpy.array([0.002, 0.045])
nu_sigma_f = numpy.array([0.004663, 0.1021])
chi = numpy.array([1.0, 0.0])

region_5.create_dataset('total', data=sigma_t)
region_5.create_dataset('scatter matrix', data=sigma_s)
region_5.create_dataset('fission', data=sigma_f)
region_5.create_dataset('nu-fission', data=nu_sigma_f)
region_5.create_dataset('chi', data=chi)


###############################################################################
################################   region 6    ################################
###############################################################################

# Create a subgroup for region 6 materials data
region_6 = material_group.create_group('region_6')

sigma_t = numpy.array([0.2652, 2.0938])
sigma_s = numpy.array([0.216931, 0.04754, 0.00, 2.074676])
sigma_f = numpy.array([0.0, 0.0])
nu_sigma_f = numpy.array([0.0, 0.0])
chi = numpy.array([1.0, 0.0])

region_6.create_dataset('total', data=sigma_t)
region_6.create_dataset('scatter matrix', data=sigma_s)
region_6.create_dataset('fission', data=sigma_f)
region_6.create_dataset('nu-fission', data=nu_sigma_f)
region_6.create_dataset('chi', data=chi)

# Close the hdf5 data file
f.close()
