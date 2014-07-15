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
current_time = str(now.month).zfill(2) + '-' + str(now.day).zfill(2) + '-' + str(now.year) + '--'
current_time = current_time + str(now.hour).zfill(2) + ':' + str(now.minute).zfill(2)
current_time = current_time + ':' + str(now.second).zfill(2)
set_log_filename('log/openmoc-' + current_time + '.log');

Timer = Timer()
