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

# The neutron multiplication factor for all fissionable materials
nu = 2.43

###############################################################################
################################   region 1    ################################
###############################################################################

# Create a subgroup for region 1 materials data
region_1 = material_group.create_group('Region 1')

sigma_t = numpy.array([0.265604, 1.579779])
sigma_s = numpy.array([0.232022, 0.02533, 0.00, 1.479479])
sigma_f = numpy.array([0.004602, 0.1091]) / nu
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
region_2 = material_group.create_group('Region 2')

sigma_t = numpy.array([0.262881, 1.752541])
sigma_s = numpy.array([0.228030, 0.02767, 0.00, 1.682071])
sigma_f = numpy.array([0.004609, 0.08675]) / nu
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
region_3 = material_group.create_group('Region 3')

sigma_t = numpy.array([0.26476, 1.594134])
sigma_s = numpy.array([0.230588, 0.02617, 0.00, 1.510694])
sigma_f = numpy.array([0.004663, 0.1021]) / nu
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
region_4 = material_group.create_group('Region 4')

sigma_t = numpy.array([0.26476, 1.594134])
sigma_s = numpy.array([0.230588, 0.02617, 0.00, 1.52081])
sigma_f = numpy.array([0.004663, 0.1021]) / nu
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
region_5 = material_group.create_group('Region 5')

sigma_t = numpy.array([0.265182, 2.093802])
sigma_s = numpy.array([0.217039, 0.04754, 0.00, 2.074692])
sigma_f = numpy.array([0.0, 0.0])
nu_sigma_f = numpy.array([0.0, 0.0])
chi = numpy.array([1.0, 0.0])

region_5.create_dataset('total', data=sigma_t)
region_5.create_dataset('scatter matrix', data=sigma_s)
region_5.create_dataset('fission', data=sigma_f)
region_5.create_dataset('nu-fission', data=nu_sigma_f)
region_5.create_dataset('chi', data=chi)

# Close the hdf5 data file
f.close()
