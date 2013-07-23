import openmoc
import _openmoc_bgxlc_single
from openmoc_bgxlc_single import *

setLogLevel(str(openmoc.getLogLevel()))
setOutputDirectory(openmoc.getOutputDirectory())
setLogfileName(openmoc.getLogfileName())

Timer = openmoc.Timer
