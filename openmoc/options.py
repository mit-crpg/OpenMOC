import sys
import getopt
import multiprocessing

# For Python 2.X.X
if sys.version_info[0] == 2:
    from log import py_printf
# For Python 3.X.X
else:
    from openmoc.log import py_printf


class Options:
    """Command-line options for runtime configuration of OpenMOC.

    This class parses, interprets and encapsulates the runtime options for
    OpenMOC. This class can be instantiate in any OpenMOC Python script which
    imports the openmoc.process module. The user may request the value of any
    or all command line options and use them as they wish.

    Attributes
    ----------
    num_azim : Integral
        The number of azimuthal angles (default is 4)
    track_spacing : Real
        The track spacing in centimeters (default is 0.1 cm)
    max_iters : Integral
        The maximum number of source iterations (default is 1000)
    tolerance : Real
        The tolerance on convergence (default is 1E-5)
    num_omp_threads : Integral
        The number of OpenMP threads (default is the # of CPU cores)
    num_thread_blocks : Integral
        The number of GPU threadblocks (default is 64)
    num_gpu_threads : Integral
        The number of GPU threads per threadblock (default is 64)

    Examples
    --------
    An example of how the Options class would be used is given as follows:

        >>> import from openmoc.options import Options

        >>> # Create an Options object which will parse all arguments
        >>> # from the command line
        >>> options = Options()

        >>> # Retrieve useful command line options
        >>> num_azim = options.num_azim
        >>> ...

        >>> # Do something useful with command line arguments
        >>> ...

    """

    def __init__(self):
        """Initialize default values for each runtime parameter and parses the
        command line arguments assigns the appropriate value to each."""

        self._num_azim = 4
        self._track_spacing = 0.1
        self._max_iters = 1000
        self._tolerance = 1E-5
        self._num_omp_threads = multiprocessing.cpu_count()
        self._num_thread_blocks = 64
        self._num_gpu_threads = 64
        self.parseArguments()

    def parseArguments(self):
        """This method parses command line options and assigns the appropriate
        values to the corresponding class attributes."""

        try:
            opts, args = getopt.getopt(sys.argv[1:],
                                      'hfa:s:i:c:t:b:g:r:l:',
                                      ['help',
                                       'num-azim=',
                                       'track-spacing=',
                                       'tolerance=',
                                       'max-iters=',
                                       'num-omp-threads=',
                                       'num-thread-blocks=',
                                       'num-gpu-threads='])

        except getopt.GetoptError as err:
            py_printf('ERROR', str(err))

        # Parse the command line arguments - error checking will occur
        # at the setter method level in C++
        for opt, arg in opts:

            # Print a report of all supported runtime options and exit
            if opt in ('-h', '--help'):

                print('{:-^80}'.format(''))
                print('{: ^80}'.format('OpenMOC v.0.1.1 runtime options'))
                print('{:-^80}'.format(''))
                print('')

                help_msg = '\t{: <35}'.format('-h, --help')
                help_msg += 'Report OpenMOC runtime options\n'
                print(help_msg)

                num_azim = '\t{: <35}'.format('-a, --num-azim=<4>')
                num_azim += 'the number of azimuthal angles\n'
                print(num_azim)

                track_spacing = '\t{: <35}'.format('-s, --track-spacing=<0.1>')
                track_spacing += 'The track spacing [cm]\n'
                print(track_spacing)

                max_iters = '\t{: <35}'.format('-i, --max-iters=<1000>')
                max_iters += 'The max number of source iterations\n'
                print(max_iters)

                tolerance = '\t{: <35}'.format('-c, --tolerance=<1E-5>')
                tolerance += 'The source convergence tolerance\n'
                print(tolerance)

                num_omp_threads = '\t{: <35}'.format('-t, --num-omp-threads=<1>')
                num_omp_threads += 'The number of OpenMP threads\n'
                print(num_omp_threads)

                num_gpu_threadblocks = '\t{: <35}'.format('-b, ' + \
                                       '--num-gpu-threadblocks=<64>')
                num_gpu_threadblocks += 'The number of GPU threadblocks\n'
                print(num_gpu_threadblocks)

                num_gpu_threads = '\t{: <35}'.format('-g, --num-gpu-threads=<64>')
                num_gpu_threads += 'The number of GPU threads per block\n'
                print(num_gpu_threads)

                sys.exit()

            elif opt in ('-a', '--num-azim'):
                    self._num_azim = int(arg)
            elif opt in ('-s', '--track-spacing'):
                self._track_spacing = float(arg)
            elif opt in ('-i', '--max-iters'):
                self._max_iters = int(arg)
            elif opt in ('-c', '--tolerance'):
                self._tolerance = float(arg)
            elif opt in ('-t', '--num-omp-threads'):
                self._num_omp_threads = int(arg)
            elif opt in ('-b', '--num-thread-blocks'):
                self._num_thread_blocks = int(arg)
            elif opt in ('-g', '--num-gpu-threads'):
                self._num_gpu_threads = int(arg)

    @property
    def num_azim(self):
        return self._num_azim

    @property
    def track_spacing(self):
        return self._track_spacing

    @property
    def max_iters(self):
        return self._max_iters

    @property
    def tolerance(self):
        return self._tolerance

    @property
    def num_omp_threads(self):
        return self._num_omp_threads

    @property
    def num_thread_blocks(self):
        return self._num_thread_blocks

    @property
    def num_threads_per_block(self):
        return self._num_gpu_threads
