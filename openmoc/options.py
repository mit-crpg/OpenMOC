##
# @file options.py
# @package openmoc.options
# @brief
# 
# @author William Boyd (wboyd@mit.edu)
# @date July 24, 2013

import log
import getopt, sys

class options:

    # The default number of azimuthal angles
    num_azim = 4

    # The default track spacing
    track_spacing = 0.1

    # The default number of OpenMP threads
    num_omp_threads = 1

    # The default number of GPU threadblocks
    num_thread_blocks = 64

    # The default number of GPU threads per threadblock
    num_gpu_threads = 64

    # The default tolerance on the source distribution convergence
    tolerance = 1E-5

    # The default maximum number of source iterations
    max_iters = 1000

    def parseArguments(self):

        try:
            opts, args = getopt.getopt(sys.argv[1:], 
                                       'a:s:t:',
                                       ['num-azim=', 
                                        'track-spacing=', 
                                        'num-omp-threads=',
                                        'num-thread-blocks=', 
                                        'num-gpu-threads=', 
                                        'tolerance=', 
                                        'max-iters='])

        except getopt.GetoptError as err:
            log.py_printf('WARNING', str(err))
            pass

        # Parse the command line arguments - error checking will occur
        # at the setter method level in C++
        for opt, arg in opts:

            if opt in ('-a', '--num-azim'):
                self.num_azim = int(arg)

            elif opt in ('-s', '--track-spacing'):
                self.track_spacing = float(arg)

            elif opt in ('-t', '--num-omp-threads'):
                self.num_omp_threads = int(arg)

            elif opt in ('--num-thread-blocks'):
                self.num_thread_blocks = int(arg)

            elif opt in ('--num-gpu-threads'):
                self.num_gpu_threads = int(arg)

            elif opt in ('--tolerance'):
                self.tolerance = float(arg)

            elif opt in ('--max-iters'):
                self.max_iters = int(arg)
