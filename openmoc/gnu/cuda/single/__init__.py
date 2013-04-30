import _openmoc
from openmoc import *
import os
import random
import datetime

# Set a default logging level
log_setlevel(NORMAL)

# Set a log file name using a date and time
now = datetime.datetime.now()
current_time = str(now.month) + '-' + str(now.day) + '-' + str(now.year) + '--' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second)
setLogfileName('log/openmoc-' + current_time + '.log');

Timer = Timer()
