##
# @file options.py
# @package openmoc.options
# @brief The options module provides the Options class to parse, interpret
#        and encapsulate command line arguments for OpenMOC runtime options.
# @author William Boyd (wboyd@mit.edu)
# @date July 24, 2013

import log
import getopt, sys


##
# @class Options process.py "openmoc/options.py"
# @brief Command-line options for runtime configuration of OpenMOC.
# @details This class parses, interprets and encapsulates the runtime options
#          for OpenMOC. This class can be instantiate in any OpenMOC Python
#          script which imports the openmoc.process module. The user may
#          request the value of any or all command line options and use them
#          as they wish. An example of how the Options class would be used
#          is given as follows:
#
# @code
#          import from openmoc.options import Options
#
#          # Create an Options object which will parse all arguments
#          # from the command line
#          options = Options()
#
#          # Retrieve useful command line options
#          num_azim = options.num_azim
#          ...
#
#          # Do something useful with command line arguments
#          ...
# @endcode
#
class Options:


  ##
  # @brief Options class constructor.
  # @details This method initializes an Options object with default values
  #          for each runtime parameter, and then parses the arguments
  #          passed in at runtime and assigns the appropriate value to each.
  #
  def __init__(self):

    ## The default number of azimuthal angles
    self._num_azim = 4

    ## The default track spacing
    self._track_spacing = 0.1

    ## The default maximum number of source iterations
    self._max_iters = 1000

    ## The default tolerance on the source distribution convergence
    self._tolerance = 1E-5

    ## The default number of OpenMP threads
    self._num_omp_threads = 1

    ## The default number of GPU threadblocks
    self._num_thread_blocks = 64

    ## The default number of GPU threads per threadblock
    self._num_gpu_threads = 64

    ## The default CMFD acceleration flag
    self._use_cmfd_acceleration = False

    ## The default CMFD relaxation factor
    self._cmfd_relax_factor = 0.6

    ## The default CMFD mesh level
    self._cmfd_mesh_level = -1

    # Parse in arguments from the command line
    self.parseArguments()


  ##
  # @brief This method parses command line options using the Python getopt
  #        module and assigns the appropriate values to the corresponding
  #        Options class attributes.
  # @param self the Options object pointer
  def parseArguments(self):
    try:
      opts, args = getopt.getopt(sys.argv[1:],
                                 'ha:s:i:c:t:b:g:r:f:l:',
                                 ['help',
                                  'num-azim=',
                                  'track-spacing=',
                                  'tolerance=',
                                  'max-iters=',
                                  'num-omp-threads=',
                                  'num-thread-blocks=',
                                  'num-gpu-threads=',
                                  'relax-factor=',
                                  'acceleration=',
                                  'mesh-level='])

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

        help_msg = '\t{: <35}'.format('-h, --help')
        help_msg = 'Report OpenMOC runtime options\n'
        print help_msg

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

        num_gpu_threadblocks = '\t{: <35}'.format('-b, ' + \
                               '--num-gpu-threadblocks=<64>')
        num_gpu_threadblocks += 'The number of GPU threadblocks\n'
        print num_gpu_threadblocks

        num_gpu_threads = '\t{: <35}'.format('-g, --num-gpu-threads=<64>')
        num_gpu_threads += 'The number of GPU threads per block\n'
        print num_gpu_threads

        relax_factor = '\t{: <35}'.format('-r, --relax-factor=<0.6>')
        relax_factor += 'The cmfd relaxation factor\n'
        print relax_factor

        acceleration = '\t{: <35}'.format('-f, --acceleration=<False>')
        relax_factor += 'The cmfd acceleration flag\n\n'

        mesh_level = '\t{: <35}'.format('-l, --mesh-level=<-1>')
        mesh_level += 'The mesh level\n'
        print mesh_level

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

      elif opt in ('-f', '--acceleration'):
        self._use_cmfd_acceleration = bool(arg)

      elif opt in ('-r', '--relax-factor'):
        self._cmfd_relax_factor = float(arg)

      elif opt in ('-l', '--mesh-level'):
        self._cmfd_mesh_level = int(arg)


  ##
  # @brief Returns the number of azimuthal angles.
  # @return the number of azimuthal angles
  def getNumAzimAngles(self):
    return self._num_azim


  ##
  # @brief Returns the track spacing [cm].
  # @return the track spacing [cm].
  def getTrackSpacing(self):
    return self._track_spacing


  ##
  # @brief Returns the maximum number of source iterations.
  # @return the maximum number of source iterations
  def getMaxIterations(self):
    return self._max_iters


  ##
  # @brief Returns the source convergence tolerance.
  # @return the source convergence tolerance
  def getTolerance(self):
    return self._tolerance


  ##
  # @brief Returns the number of OpenMP multi-core CPU threads.
  # @return the number of OpenMP threads
  def getNumThreads(self):
    return self._num_omp_threads


  ##
  # @brief Returns the number of CUDA thread blocks for a GPU.
  # @return the number of CUDA thread blocks
  def getNumThreadBlocks(self):
    return self._num_thread_blocks


  ##
  # @brief Returns the number of CUDA threads per block for a GPU.
  # @return the number of CUDA threads per block
  def getNumThreadsPerBlock(self):
    return self._num_gpu_threads


  ##
  # @brief Returns the whether or not to use CMFD acceleration.
  # @return use CMFD acceleration (true) or not (false)
  def useCmfdAccleration(self):
    return self._use_cmfd_acceleration


  ##
  # @brief Returns the CMFD relaxation factor.
  # @return the CMFD relaxation factor.
  def getCmfdRelaxationFactor(self):
    return self._cmfd_relax_factor


  ##
  # @brief Returns the number of CMFD multigrid mesh levels.
  # @return the number of CMFD multigrid mesh levels
  def getCmfdMeshLevel(self):
    return self._cmfd_mesh_level
