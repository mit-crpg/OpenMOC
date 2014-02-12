import openmoc
import _openmoc_intel_double
from openmoc_intel_double import *
import signal

# Tell Python to recognize CTRL+C and stop the C++ extension module
# when this is passed in from the keyboard
signal.signal(signal.SIGINT, signal.SIG_DFL)

set_log_level(str(openmoc.get_log_level()))
set_output_directory(openmoc.get_output_directory())
set_log_filename(openmoc.get_log_filename())

Timer = openmoc.Timer
