#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
import openmoc
import numpy as np

"""Unit test for the exponential approximations"""

# Choose a large domain of tau values to show the range of validity
test_array = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100]
tol = 1e-8

expf1_p = np.zeros(len(test_array))
expf2_p = np.zeros(len(test_array))
exph_p = np.zeros(len(test_array))
expf1_r = np.zeros(len(test_array))
expf2_r = np.zeros(len(test_array))
exph_r = np.zeros(len(test_array))
expg_r = np.zeros(len(test_array))
expg2_r = np.zeros(len(test_array))
expv_r = np.zeros(len(test_array))

# Compute all exponentials
for i, x in enumerate(test_array):
    value = openmoc.expF1_poly(x)
    expf1_p[i] = value

    value = openmoc.expF2_poly(x)
    expf2_p[i] = value

    value = openmoc.expH_poly(x)
    exph_p[i] = value

    value = openmoc.expF1_fractional(x)
    expf1_r[i] = value

    value = openmoc.expF2_fractional(x)
    expf2_r[i] = value

    value = openmoc.expH_fractional(x)
    exph_r[i] = value

    value = openmoc.expG_fractional(x)
    expg_r[i] = value

    value = openmoc.expG2_fractional(x)
    expg2_r[i] = value

    value = openmoc.cram7(x)
    expv_r[i] = value

# Test polynomial forms
assert all(expf1_p - np.array([0.9999999950013811, 0.9999999500138135, 0.9999995001382836,
                   0.9999950013977791, 0.9999500154719884, 0.9995003041034797,
                   0.9950179432863955, 0.9516342044055093, 0.6321158759617262,
                   -1.2492783407705708, -230617368.72428462]) < tol)

assert all(expf2_p - np.array([1.6666567586545621e-09, 1.6666566836731985e-08,
                   1.6666559338597827e-07, 1.6666484357478334e-06,
                   1.6665734568489468e-05, 0.00016658238898725115,
                   0.0016583503730099513, 0.01585782946791598,
                   0.1036377806982629, -1.1392815325771322,
                   -2288519837.2429357]) < tol)

assert all(exph_p - np.array([0.4999999966669252, 0.4999999666692532, 0.4999996666926448,
                  0.49999666693768496, 0.49996667050055443, 0.49966681734657153,
                  0.4966793781803916, 0.46788544865980963, 0.2642408690992135,
                  0.7688656031312362, 1954689084.2831757]) < tol)

# Test rational fractions
assert all(expf1_r - np.array([0.9999999950003422, 0.9999999500034228, 0.9999995000343784,
                   0.9999950003587617, 0.9999500050851808, 0.9995002005718668,
                   0.9950169413610763, 0.9516272442309586, 0.6321211831479621,
                   0.09999547551656497, 0.010000005017389789]) < tol)

assert all(expf2_r - np.array([1.6666614617049907e-09, 1.6666613867190775e-08,
                   1.6666606368601697e-07, 1.6666531382932901e-06,
                   1.666578154844449e-05, 0.00016658285423037965,
                   0.0016583545640157496, 0.015857824983353568,
                   0.10363825395873519, 0.0800056000658619,
                   0.009799985072224623]) < tol)

assert all(exph_r - np.array([0.4999999966663039, 0.499999966663041, 0.49999966663052203,
                  0.49999666631649736, 0.49996666429260567, 0.4996667556587436,
                  0.4966787994230765, 0.46788257450436915, 0.2642404988745949,
                  0.009995036061745774, 0.00010000003576523983]) < tol)

assert all(expg_r - np.array([0.49999999833335007, 0.4999999833335013, 0.4999998333350512,
                  0.4999983333542612, 0.4999833339174937, 0.499833376655743,
                  0.4983375072555877, 0.48374187782582545, 0.3678794669042375,
                  0.09000046675901234, 0.009899999709267566]) < tol)

assert all(expg2_r - np.array([-8.335775827058522e-10, -8.335775300276534e-09,
                   -8.335770032458422e-08, -8.335717354455771e-07,
                   -8.335190592273496e-06, -8.329924754387514e-05,
                   -0.0008277444282564627, -0.007769958754069781,
                   -0.04061078253892451, 0.07865975028023744, 0.15676951276097334]) < tol)

assert all(expv_r - np.array([-1.0000001480267286e-08, -1.0000001930272869e-07,
                  -1.0000006430330166e-06, -1.0000051431051701e-05,
                  -0.00010000501453121614, -0.001000500315964985,
                  -0.010050169180882236, -0.10517105848100763,
                  -1.7195674765901083, 1.0604692107575897, 1.0000098580016492]) < tol)
