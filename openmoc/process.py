##
# @file process.py
# @package openmoc.process
# @brief
# 
# @author William Boyd (wboyd@mit.edu)
# @date April 27, 2013

import numpy as np
import matplotlib.pyplot as plt
from openmoc import *
from log import *
import scipy.integrate as integrate
import numpy
import os


## A static variable for the output directory in which to save plots
subdirectory = "/plots/"


#
# @brief Returns an array of the center values for a tally's bins.
# @details A wrapper function to make it easier to access a tally's bin center
#          data array through SWIG. This would be invoked a PINSPEC Python
#          input file as follows:
#
# @code
#          bin_center_array = pinspec.process.getTallyCenters(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally bin centers
def strongScalingStudy(solver, num_threads=None, max_iters=25,
                     title='', filename=''):

    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    if (num_threads is None):
        num_threads = numpy.linspace(1,12,12)
    else:
        num_threads = numpy.array(num_threads)

    times = []
    timer = Timer()

    for threads in num_threads:
        solver.setNumThreads(int(threads))
        
        timer.resetTimer()
        timer.startTimer()
        solver.convergeSource(max_iters)
        timer.stopTimer()
        timer.recordSplit('' + str(threads) + ' threads')
        times.append(timer.getTime())


    timer.printSplits()

    # Plot Runtime
    if title == '':
        title1 = 'OpenMOC OpenMP Strong Scaling Time per Iteration'
    else:
        title1 = title + ' Runtime'

    if filename == '':
        filename1 = directory + 'strong-scaling-runtime.png'
    else:
        filename1 = directory + filename + '-runtime.png'

    times = numpy.array(times) / solver.getNumIterations()
    fig = plt.figure()
    plt.plot(num_threads, times, linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Runtime [sec]')
    plt.title(title1)
    plt.grid()
    plt.savefig(filename1)

    # Plot Speedup
    if title == '':
        title2 = 'OpenMOC OpenMP Strong Scaling Speedup'
    else:
        title2 = title + ' Speedup'

    if filename == '':
        filename2 = directory + 'strong-scaling-speedup.png'
    else:
        filename2 = directory + filename + '-speedup.png'

    speedup = times[0] / times
    fig = plt.figure()
    plt.plot(num_threads, speedup, linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Speedup')
    plt.title(title2)
    plt.grid()
    plt.savefig(filename2)

    # Plot parallel efficiency
    if title == '':
        title3 = 'OpenMOC OpenMP Strong Scaling Parallel Efficiency'
    else:
        title3 = title + ' Parallel Efficiency'

    if filename == '':
        filename3 = directory + 'strong-scaling-efficiency.png'
    else:
        filename3 = directory + filename + '-efficiency.png'
    
    efficiency = (times[0] * num_threads[0]) / (num_threads * times)
    fig = plt.figure()
    plt.plot(num_threads, efficiency, linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Efficiency')
    plt.title(title3)
    plt.grid()
    plt.savefig(filename3)


#
# @brief Returns an array of the center values for a tally's bins.
# @details A wrapper function to make it easier to access a tally's bin center
#          data array through SWIG. This would be invoked a PINSPEC Python
#          input file as follows:
#
# @code
#          bin_center_array = pinspec.process.getTallyCenters(tally)
# @endcode
#
# @param tally the tally of interest
# @return a numpy array with the tally bin centers
def weakScalingStudy(solver, track_generator,  num_threads=None, \
                         max_iters=25, title='', filename=''):

    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    if (num_threads is None):
        num_threads = numpy.linspace(1,12,12)
    else:
        num_threads = numpy.array(num_threads)

    times = []
    timer = Timer()

    num_azim = num_threads * 4

    for i in range(len(num_azim)):

        track_generator.setNumAzim(int(num_azim[i]))
        track_generator.generateTracks()
                                  
        solver.setNumThreads(int(num_threads[i]))
        
        timer.resetTimer()
        timer.startTimer()
        solver.convergeSource(max_iters)
        timer.stopTimer()
        times.append(timer.getTime())
        timer.recordSplit('' + str(num_threads[i]) + ' threads')


    timer.printSplits()

    # Plot Runtime
    if title == '':
        title1 = 'OpenMOC OpenMP Weak Scaling Time per Iteration'
    else:
        title1 = title + ' Runtime'

    if filename == '':
        filename1 = directory + 'weak-scaling-runtime.png'
    else:
        filename1 = directory + filename + '-runtime.png'

    times = numpy.array(times) / solver.getNumIterations()
    fig = plt.figure()
    plt.plot(num_threads, times, linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Runtime [sec]')
    plt.title(title1)
    plt.grid()
    plt.savefig(filename1)

    # Plot Speedup
    if title == '':
        title2 = 'OpenMOC OpenMP Weak Scaling Speedup'
    else:
        title2 = title + ' Speedup'

    if filename == '':
        filename2 = directory + 'weak-scaling-speedup.png'
    else:
        filename2 = directory + filename + '-speedup.png'

    speedup = num_threads * times[0] / times
    fig = plt.figure()
    plt.plot(num_threads, speedup, linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Speedup')
    plt.title(title2)
    plt.grid()
    plt.savefig(filename2)

    # Plot parallel efficiency
    if title == '':
        title3 = 'OpenMOC OpenMP Weak Scaling Parallel Efficiency'
    else:
        title3 = title + ' Parallel Efficiency'

    if filename == '':
        filename3 = directory + 'weak-scaling-efficiency.png'
    else:
        filename3 = directory + filename + '-efficiency.png'
    
    efficiency = (times[0] * num_threads[0]) / (num_threads * times)
    fig = plt.figure()
    plt.plot(num_threads, efficiency, linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Efficiency')
    plt.title(title3)
    plt.grid()
    plt.savefig(filename3)


def profile(solver, num_threads=4,  max_iters=25, title='', filename=''):

    global subdirectory

    directory = getOutputDirectory() + subdirectory

    # Make directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    solver.setNumThreads(int(num_threads))
    
    solver.convergeSource(max_iters)

    # The slices will be ordered and plotted counter-clockwise.
    labels = ['Computing FSR sources', 'Computing k-eff', 
              'Flattening FSR flux to zero', 'Normalizing fluxes',
              'Reducing FSR scalar fluxes', 'Transport sweep']

    times = []
    timer = Timer()

    for label in labels:
        times.append(timer.getSplit(label))

    timer.printSplits()
    times = numpy.array(times)
    tot_time = numpy.sum(times)
    fracs = times / tot_time

    # Plot a pie chart showing the break down in compute time
    if title == '':
        title = 'OpenMOC Compute Time Profile'
    else:
        title = title + ' Compute Time Profile'

    if filename == '':
        filename = directory + 'compute-time-profile.png'
    else:
        filename = directory + filename + '-profile.png'
    
    fig = plt.figure()
    plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
    plt.title(title)
    plt.savefig(filename)
