import _openmoc
from openmoc import *
import os
import random
import datetime
import signal

# Tell Python to recognize CTRL+C and stop the C++ extension module
# when this is passed in from the keyboard
signal.signal(signal.SIGINT, signal.SIG_DFL)

# Set a log file name using a date and time
now = datetime.datetime.now()
current_time = str(now.month) + '-' + str(now.day) + '-' + str(now.year) + '--'
current_time = current_time + str(now.hour) + ':' + str(now.minute)
current_time = current_time + ':' + str(now.second)
set_log_filename('log/openmoc-' + current_time + '.log');

Timer = Timer()
