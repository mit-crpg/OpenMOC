from numpy.random import rand
from _sum_array import *

# Initialize a random NumPy array
input_array = rand(5)

# Sum the values in the random NumPy array
sum = sum_array(input_array)

print 'The sum of the array is %f' % (sum)
