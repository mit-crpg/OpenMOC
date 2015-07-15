import sys, os
import random
import datetime
import signal

# For Python 2.X.X
if (sys.version_info[0] == 2):
  from openmoc import *
# For Python 3.X.X
else:
  from openmoc.openmoc import *

# Tell Python to recognize CTRL+C and stop the C++ extension module
# when this is passed in from the keyboard
signal.signal(signal.SIGINT, signal.SIG_DFL)

# Set a log file name using the date and time
now = datetime.datetime.now()
time = (now.month, now.day, now.year, now.hour, now.minute, now.second)
curr_time = '-'.join(map(lambda x: x.zfill(2), map(str, time)))
initialize_logger()
set_log_filename('openmoc-' + current_time + '.log');

# Create singleton Timer object for shared access throughout the module
Timer = Timer()
