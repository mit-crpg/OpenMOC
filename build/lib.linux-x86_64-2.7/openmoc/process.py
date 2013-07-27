##
# @file process.py
# @package openmoc.process
# @brief
# 
# @author William Boyd (wboyd@mit.edu)
# @date April 27, 2013

import matplotlib.pyplot as plt
from log import *
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
def strongScalingStudy(geometry, num_azim=48, track_spacing=0.1, 
                       num_threads=None, max_iters=10, 
                       compiler='gnu', precision='double', 
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

    if compiler == 'all':
        compiler = ['gnu', 'intel']
    elif compiler == 'intel':
        compiler = ['intel']
    else:
        compiler = ['gnu']

    if precision == 'all':
        precision = ['single', 'double']
    elif precision == 'double':
        precision = ['double']
    else:
        precision = ['single']

    legend = []
    times = numpy.zeros((len(precision), len(compiler), num_threads.size))

    for i in range(len(precision)):
        for j in range(len(compiler)):
            
            fp = precision[i]
            cc = compiler[j]

            if fp == 'single' and cc == 'gnu':
                import gnu.single as openmoc
            elif fp == 'single' and cc == 'intel':
                import intel.single as openmoc
            elif fp == 'double' and cc == 'gnu':
                import gnu.double as openmoc
            elif fp == 'double' and cc == 'intel':
                import intel.double as openmoc

            track_generator = openmoc.TrackGenerator(geometry, num_azim, 
                                                     track_spacing)
            track_generator.generateTracks()

            py_printf('NORMAL', '# FSRs = %d, # tracks = %d, # segments = %d',
                      geometry.getNumFSRs(), track_generator.getNumTracks(),
                      track_generator.getNumSegments())

            solver = openmoc.Solver(geometry, track_generator)

            timer = openmoc.Timer()
            
            legend.append('' + cc + '-' + fp)

            for k in range(num_threads.size):
                
                solver.setNumThreads(int(num_threads[k]))

                timer.resetTimer()
                timer.startTimer()
                solver.convergeSource(max_iters)
                timer.stopTimer()
                timer.recordSplit('' + str(num_threads[k]) + ' threads' + ' ' 
                                  + cc + '-' + fp)
                times[i][j][k] = timer.getTime()
                timer.printSplits()

            timer.printSplits()

    times /= max_iters

    # Plot Runtime
    if title == '':
        title1 = 'OpenMOC OpenMP Strong Scaling Time per Iteration'
    else:
        title1 = title + ' Runtime'

    if filename == '':
        filename1 = directory + 'strong-scaling-runtime.png'
    else:
        filename1 = directory + filename + '-runtime.png'

    for i in range(len(precision)):
        for j in range(len(compiler)):
            plt.plot(num_threads, times[i][j], linewidth=3)

    fig = plt.figure()
    for i in range(len(precision)):
        for j in range(len(compiler)):
            plt.plot(num_threads, times[i][j], linewidth=3 )
    plt.xlabel('# threads')
    plt.ylabel('Runtime [sec]')
    plt.title(title1)
    if (len(legend) > 1):
        plt.legend(legend)
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

    speedup = times[:,:,0][:,:,numpy.newaxis] / times
    fig = plt.figure()
    for i in range(len(precision)):
        for j in range(len(compiler)):
            plt.plot(num_threads, speedup[i][j], linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Speedup')
    plt.title(title2)
    if (len(legend) > 1):
        plt.legend(legend)
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
    
    efficiency = (times[:,:,0][:,:,numpy.newaxis] * num_threads[0]) / \
                 (num_threads * times)
    fig = plt.figure()
    for i in range(len(precision)):
        for j in range(len(compiler)):
            plt.plot(num_threads, efficiency[i][j], linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Efficiency')
    plt.title(title3)
    if (len(legend) > 1):
        plt.legend(legend)
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
def weakScalingStudy(geometry, num_azim=4, track_spacing=0.1, 
                     num_threads=None, max_iters=10, 
                     compiler='gnu', precision='double', 
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

    if compiler == 'all':
        compiler = ['gnu', 'intel']
    elif compiler == 'intel':
        compiler = ['intel']
    else:
        compiler = ['gnu']

    if precision == 'all':
        precision = ['single', 'double']
    elif precision == 'double':
        precision = ['double']
    else:
        precision = ['single']

    legend = []
    times = numpy.zeros((len(precision), len(compiler), num_threads.size))
    num_azim *= num_threads

    for i in range(len(precision)):
        for j in range(len(compiler)):
            
            fp = precision[i]
            cc = compiler[j]

            if fp == 'single' and cc == 'gnu':
                import gnu.single as openmoc
            elif fp == 'single' and cc == 'intel':
                import intel.single as openmoc
            elif fp == 'double' and cc == 'gnu':
                import gnu.double as openmoc
            elif fp == 'double' and cc == 'intel':
                import intel.double as openmoc

            timer = openmoc.Timer()
            
            legend.append('' + cc + '-' + fp)

            for k in range(len(num_azim)):

                track_generator = openmoc.TrackGenerator(geometry, 
                                                         int(num_azim[k]),
                                                         track_spacing)
                track_generator.generateTracks()

                py_printf('NORMAL', '# FSRs = %d, # tracks = %d, # ' + \
                          'segments = %d', geometry.getNumFSRs(), 
                          track_generator.getNumTracks(),
                          track_generator.getNumSegments())

                solver = openmoc.Solver(geometry, track_generator)
                solver.setNumThreads(int(num_threads[k]))
        
                timer.resetTimer()
                timer.startTimer()
                solver.convergeSource(max_iters)
                timer.stopTimer()
                timer.recordSplit('' + str(num_threads[k]) + ' threads' + ' ' 
                                  + cc + '-' + fp)
                times[i][j][k] = timer.getTime()
                timer.printSplits()

            timer.printSplits()

    times /= max_iters

    # Plot Runtime
    if title == '':
        title1 = 'OpenMOC OpenMP Weak Scaling Time per Iteration'
    else:
        title1 = title + ' Runtime'

    if filename == '':
        filename1 = directory + 'weak-scaling-runtime.png'
    else:
        filename1 = directory + filename + '-runtime.png'

    fig = plt.figure()
    for i in range(len(precision)):
        for j in range(len(compiler)):
            plt.plot(num_threads, times[i][j], linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Runtime [sec]')
    plt.title(title1)
    if (len(legend) > 1):
        plt.legend(legend)
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

    speedup = num_threads * (times[:,:,0][:,:,numpy.newaxis] / times)
    fig = plt.figure()
    for i in range(len(precision)):
        for j in range(len(compiler)):
            plt.plot(num_threads, speedup[i][j], linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Speedup')
    plt.title(title2)
    if (len(legend) > 1):
        plt.legend(legend)
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
    
    efficiency = times[:,:,0][:,:,numpy.newaxis] / times
    fig = plt.figure()
    for i in range(len(precision)):
        for j in range(len(compiler)):
            plt.plot(num_threads, efficiency[i][j], linewidth=3)
    plt.xlabel('# threads')
    plt.ylabel('Efficiency')
    plt.title(title3)
    if (len(legend) > 1):
        plt.legend(legend)
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

    times = numpy.zeros(3)
    timer = Timer()

    for label in labels:
        if label == 'Transport sweep':
            times[0] = timer.getSplit(label)
        elif label == 'Computing FSR sources':
            times[1] = timer.getSplit(label)
        else:
            times[2] += timer.getSplit(label)

    timer.printSplits()
    tot_time = numpy.sum(times)
    fracs = times / tot_time
    labels = ['Transport Sweep', 'FSR Sources', 'Other']

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
