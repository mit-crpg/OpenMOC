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

    # The default maximum number of source iterations
    max_iters = 1000

    # The default tolerance on the source distribution convergence
    tolerance = 1E-5

    # The default number of OpenMP threads
    num_omp_threads = 1

    # The default number of GPU threadblocks
    num_thread_blocks = 64

    # The default number of GPU threads per threadblock
    num_gpu_threads = 64


    def parseArguments(self):

        try:
            opts, args = getopt.getopt(sys.argv[1:], 
                                       'ha:s:i:c:t:b:g:',
                                       ['help',
                                        'num-azim=',
                                        'track-spacing=',
                                        'tolerance=',
                                        'max-iters=',
                                        'num-omp-threads=',
                                        'num-thread-blocks=', 
                                        'num-gpu-threads='])

        except getopt.GetoptError as err:
            log.py_printf('WARNING', str(err))
            pass

        # Parse the command line arguments - error checking will occur
        # at the setter method level in C++
        for opt, arg in opts:

            # Print a report of all supported runtime options and exit
            if opt in ('-h', '--help'):

                print '{:-^80}'.format('')
                print '{: ^80}'.format('OpenMOC v.0.1.1 runtime options')
                print '{:-^80}'.format('')
                print

                help = '\t{: <35}'.format('-h, --help')
                help += 'Report OpenMOC runtime options\n'
                print help

                num_azim = '\t{: <35}'.format('-a, --num-azim=<4>')
                num_azim += 'the number of azimuthal angles\n'
                print num_azim

                track_spacing = '\t{: <35}'.format('-s, --track-spacing=<0.1>')
                track_spacing += 'The track spacing [cm]\n'
                print track_spacing

                max_iters = '\t{: <35}'.format('-i, --max-iters=<1000>')
                max_iters += 'The max number of source iterations\n'
                print max_iters

                tolerance = '\t{: <35}'.format('-c, --tolerance=<1E-5>')
                tolerance += 'The source convergence tolerance\n'
                print tolerance

                num_omp_threads = '\t{: <35}'.format('-t, --num-omp-threads=<1>')
                num_omp_threads += 'The number of OpenMP threads\n'
                print num_omp_threads

                num_gpu_threadblocks = '\t{: <35}'.format('-b, --num-gpu-threadblocks=<64>')
                num_gpu_threadblocks += 'The number of GPU threadblocks\n'
                print num_gpu_threadblocks

                num_gpu_threads = '\t{: <35}'.format('-g, --num-gpu-threads=<64>')
                num_gpu_threads += 'The number of GPU threads per block\n'
                print num_gpu_threads

                sys.exit()

            elif opt in ('-a', '--num-azim'):
                self.num_azim = int(arg)

            elif opt in ('-s', '--track-spacing'):
                self.track_spacing = float(arg)

            elif opt in ('i', '--max-iters'):
                self.max_iters = int(arg)

            elif opt in ('-c, --tolerance'):
                self.tolerance = float(arg)

            elif opt in ('-t', '--num-omp-threads'):
                self.num_omp_threads = int(arg)

            elif opt in ('-b', '--num-thread-blocks'):
                self.num_thread_blocks = int(arg)

            elif opt in ('-g', '--num-gpu-threads'):
                self.num_gpu_threads = int(arg)

