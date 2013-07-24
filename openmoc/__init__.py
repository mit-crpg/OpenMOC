import _openmoc
from openmoc import *
import options
import os
import random
import datetime

# Set a log file name using a date and time
now = datetime.datetime.now()
current_time = str(now.month) + '-' + str(now.day) + '-' + str(now.year) + '--' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second)
setLogfileName('log/openmoc-' + current_time + '.log');

Timer = Timer()

options = options.options()
options.parseArguments()
