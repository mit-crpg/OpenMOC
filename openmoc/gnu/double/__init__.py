import openmoc
import _openmoc_gnu_double
from openmoc_gnu_double import *

setLogLevel(str(openmoc.getLogLevel()))
setOutputDirectory(openmoc.getOutputDirectory())
setLogfileName(openmoc.getLogfileName())

Timer = openmoc.Timer
