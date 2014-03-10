import signal, sys

import _openmoc_gnu_single

# For Python 2.X.X
if (sys.version_info[0] == 2):
  import openmoc
  from openmoc_gnu_single import *
# For Python 3.X.X
else:
  import openmoc.openmoc as openmoc
  from openmoc.gnu.single.openmoc_gnu_single import *


# Tell Python to recognize CTRL+C and stop the C++ extension module
# when this is passed in from the keyboard
signal.signal(signal.SIGINT, signal.SIG_DFL)

openmoc.set_log_level(str(openmoc.get_log_level()))
openmoc.set_output_directory(openmoc.get_output_directory())
openmoc.set_log_filename(openmoc.get_log_filename())

Timer = openmoc.Timer
