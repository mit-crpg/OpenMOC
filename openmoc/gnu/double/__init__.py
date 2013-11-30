import openmoc
import _openmoc_gnu_double
from openmoc_gnu_double import *
import signal

# Tell Python to recognize CTRL+C and stop the C++ extension module
# when this is passed in from the keyboard
signal.signal(signal.SIGINT, signal.SIG_DFL)

openmoc.setLogLevel(str(openmoc.getLogLevel()))
openmoc.setOutputDirectory(openmoc.getOutputDirectory())
openmoc.setLogfileName(openmoc.getLogfileName())

Timer = openmoc.Timer
