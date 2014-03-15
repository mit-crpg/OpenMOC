import signal, sys

# For Python 2.X.X
if (sys.version_info[0] == 2):
  import openmoc
  import _openmoc_cuda
  from openmoc_cuda_double import *
# For Python 3.X.X
else:
  import openmoc.openmoc as openmoc
  import _openmoc_cuda
  from openmoc.cuda.double.openmoc_cuda_double import *

# Tell Python to recognize CTRL+C and stop the C++ extension module
# when this is passed in from the keyboard
signal.signal(signal.SIGINT, signal.SIG_DFL)

Timer = openmoc.Timer
