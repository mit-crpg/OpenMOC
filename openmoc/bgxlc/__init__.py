import openmoc
import _openmoc_bgxlc_double
from openmoc_bgxlc_double import *

setLogLevel(str(openmoc.getLogLevel()))
setOutputDirectory(openmoc.getOutputDirectory())
setLogfileName(openmoc.getLogfileName())

Timer = openmoc.Timer
