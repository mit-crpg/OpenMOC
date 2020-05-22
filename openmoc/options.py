import sys
import getopt
import multiprocessing
import os

# For Python 2.X.X
if sys.version_info[0] == 2:
    from log import py_printf
# For Python 3.X.X
else:
    from openmoc.log import py_printf


class Options(object):
    """Command-line options for runtime configuration of OpenMOC.

    This class parses, interprets and encapsulates the runtime options for
    OpenMOC. This class can be instantiate in any OpenMOC Python script which
    imports the openmoc.process module. The user may request the value of any
    or all command line options and use them as they wish.

    Attributes
    ----------
    num_azim : Integral
        The number of azimuthal angles (default is 4)
    num_polar : Integral
        The number of azimuthal angles (default is 6)
    azim_spacing : Real
        The azimuthal track spacing in centimeters (default is 0.1 cm)
    polar_spacing : Real
        The polar track spacing in centimeters (default is 0.1 cm)
    max_iters : Integral
        The maximum number of source iterations (default is 1000)
    tolerance : Real
        The tolerance on convergence (default is 1E-5)
    num_omp_threads : Integral
        The number of OpenMP threads (default is the # of CPU cores)
    num_thread_blocks : Integral
        The number of GPU threadblocks (default is 64)
    num_threads_per_block : Integral
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

        self._short_args = 'fha:s:p:l:i:c:t:b:g:'
        self._long_args = ['help',
                           'num-azim=',
                           'azim-spacing=',
                           'num-polar=',
                           'polar-spacing=',
                           'max-iters=',
                           'tolerance=',
                           'num-omp-threads=',
                           'num-thread-blocks=',
                           'num-threads-per-block=']

        self._num_azim = 4
        self._num_polar = 6
        self._azim_spacing = 0.1
        self._polar_spacing = 1.5
        self._max_iters = 1000
        self._tolerance = 1E-5
        try:
            self._num_omp_threads = int(os.environ["OMP_NUM_THREADS"])
        except KeyError:
            self._num_omp_threads = multiprocessing.cpu_count()
        self._num_thread_blocks = 64
        self._num_threads_per_block = 64

        self._opts = None
        self._args = None

        self.parseArguments()

    @property
    def short_args(self):
        return self._short_args

    @property
    def long_args(self):
        return self._long_args

    @property
    def opts(self):
        return self._opts

    @property
    def args(self):
        return self._args

    @property
    def num_azim(self):
        return self._num_azim

    @property
    def num_polar(self):
        return self._num_polar

    @property
    def azim_spacing(self):
        return self._azim_spacing

    @property
    def polar_spacing(self):
        return self._polar_spacing

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
        return self._num_threads_per_block

    def parseArguments(self):
        """This method parses command line options and assigns the appropriate
        values to the corresponding class attributes."""

        try:
            self._opts, self._args = \
                getopt.getopt(sys.argv[1:], self.short_args, self.long_args)
        except getopt.GetoptError as err:
            py_printf('ERROR', str(err))

        # Parse the command line arguments - error checking will occur
        # at the setter method level in C++
        for opt, arg in self.opts:

            # Print a report of all supported runtime options and exit
            if opt in ('-h', '--help'):

                print('{:-^80}'.format(''))
                print('{: ^80}'.format('OpenMOC v.0.4.0 runtime options'))
                print('{:-^80}'.format(''))
                print('')

                help_msg = '\t{: <35}'.format('-h, --help')
                help_msg += 'Report OpenMOC runtime options\n'
                print(help_msg)

                num_azim = '\t{: <35}'.format('-a, --num-azim=<4>')
                num_azim += 'the number of azimuthal angles\n'
                print(num_azim)

                num_polar = '\t{: <35}'.format('-p, --num-polar=<6>')
                num_polar += 'the number of polar angles\n'
                print(num_polar)

                azim_spacing = '\t{: <35}'.format('-s, --azim-spacing=<0.1>')
                azim_spacing += 'The azimuthal track spacing [cm]\n'
                print(azim_spacing)

                polar_spacing = '\t{: <35}'.format('-l, --polar-spacing=<1.5>')
                polar_spacing += 'The polar track spacing [cm]\n'
                print(polar_spacing)

                max_iters = '\t{: <35}'.format('-i, --max-iters=<1000>')
                max_iters += 'The max number of source iterations\n'
                print(max_iters)

                tolerance = '\t{: <35}'.format('-c, --tolerance=<1E-5>')
                tolerance += 'The source convergence tolerance\n'
                print(tolerance)

                num_omp_threads = '\t{: <35}'.format('-t, --num-omp-threads=<1>')
                num_omp_threads += 'The number of OpenMP threads\n'
                print(num_omp_threads)

                num_threadblocks = '\t{: <35}'.format('-b, ' + \
                                       '--num-thread-blocks=<64>')
                num_threadblocks += 'The number of GPU threadblocks\n'
                print(num_threadblocks)

                num_threads_per_block = \
                    '\t{: <35}'.format('-g, --num-threads-per-block=<64>')
                num_threads_per_block += 'The number of GPU threads per block\n'
                print(num_threads_per_block)

                sys.exit()

            elif opt in ('-a', '--num-azim'):
                self._num_azim = int(arg)
            elif opt in ('-p', '--num-polar'):
                self._num_polar = int(arg)
            elif opt in ('-s', '--azim-spacing'):
                self._azim_spacing = float(arg)
            elif opt in ('-l', '--polar-spacing'):
                self._polar_spacing = float(arg)
            elif opt in ('-i', '--max-iters'):
                self._max_iters = int(arg)
            elif opt in ('-c', '--tolerance'):
                self._tolerance = float(arg)
            elif opt in ('-t', '--num-omp-threads'):
                self._num_omp_threads = int(arg)
            elif opt in ('-b', '--num-thread-blocks'):
                self._num_thread_blocks = int(arg)
            elif opt in ('-g', '--num-threads-per-block'):
                self._num_threads_per_block = int(arg)
