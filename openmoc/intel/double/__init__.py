import openmoc
import _openmoc_intel_double
from openmoc_intel_double import *

setLogLevel(str(openmoc.getLogLevel()))
setOutputDirectory(openmoc.getOutputDirectory())
setLogfileName(openmoc.getLogfileName())

Timer = openmoc.Timer
