import numpy
#from numpy.linalg import inv
import openmoc

sigma_t = numpy.array([1.779490E-01, 3.298050E-01, 4.803880E-01,
                       5.543670E-01, 3.118010E-01, 3.951680E-01, 5.644060E-01])
sigma_s = numpy.array([[1.275370E-01, 4.237800E-02, 9.437400E-06, 5.516300E-09,
                          0., 0., 0.],
                        [0., 3.244560E-01,1.631400E-03, 3.142700E-09,
                          0., 0., 0.],
                        [0., 0., 4.509400E-01, 2.679200E-03, 0., 0., 0.], 
                        [0., 0., 0., 4.525650E-01, 5.566400E-03, 0., 0.], 
                        [0., 0., 0., 1.252500E-04, 2.714010E-01, 1.025500E-02,
                            1.002100E-08], 
                        [0., 0., 0., 0., 1.296800E-03, 2.658020E-01,
                                                                  1.680900E-02],
                        [0.,0., 0., 0., 0., 8.545800E-03, 2.730800E-01]])

sigma_f = numpy.array([7.212060E-03, 8.193010E-04, 6.453200E-03,
                       1.856480E-02, 1.780840E-02, 8.303480E-02, 2.160040E-01])
nu_sigma_f = numpy.array([2.005998E-02, 2.027303E-03, 1.570599E-02,
                          4.518301E-02, 4.334208E-02, 2.020901E-01,
                          5.257105E-01])
chi = numpy.array([5.87910E-01, 4.11760E-01, 3.39060E-04,
                   1.17610E-07, 0., 0., 0.])

a = numpy.zeros((7,7))
for i in range(0,7):
  for j in range(0,7):
    if i == j:
      a[i][i] = sigma_t[i] - sigma_s[i][i]
    else:
      a[i][j] = - sigma_s[j][i]
      print i
      print j
      print
#print a[6][5]

f = numpy.zeros((7,7))
for i in range(0,7):
  for j in range(0,7):
    f[i][j] = chi[i]*nu_sigma_f[j]

A = numpy.matrix(a)
F = numpy.matrix(f)
print
print a
print
print A
print
print f
print
print F
print
E = numpy.linalg.inv(A)*F
print
print numpy.linalg.eigvals(E)

print


