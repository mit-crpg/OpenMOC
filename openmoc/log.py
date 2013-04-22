## 
# @file log.py
# @package openmoc.log 
# @brief Utility functions for writing log messages to the screen.
#
# @details This module includes a set of wrapper functions for
#          the logging routines provided by OpenMOC's C++ source 
#          code. These Python methods provide an interface for 
#          creating formatted log messages using level-based loggin
#          and to print them to the screen as well as a logfile.
# @author Samuel Shaner
# @date March 15, 2013


import numpy
from openmoc import *


## 
# @brief Function to print a log message to the screen
# @details This method is a wrapper to the log_printf C++ routine. It
#          allows for formatted messages to be printed to the screen 
#          in a similar fashion to the C/C++ printf method, but with
#          additional formatting provided by the OpenMOC logging utilities.
#          An example of how this might be used in a OpenMOC Python script
#          is as follows:
#
# @code
#          value1 = 25
#          value2 = 26.0
#          log.py_printf('NORMAL', 'My name is Will and I am %d going on'\
#                               ' %f years of age', value1, value2)
# @endcode
#
# @param level the logging level for this message
# @param my_str the string to print to the screen
# @param *args a variable length list of values for the message string
def py_printf(level, my_str, *args):
    if level == 'DEBUG':
        log_printf(DEBUG, my_str % args)
    elif level == 'INFO':
        log_printf(INFO, my_str % args)
    elif level == 'NORMAL':
        log_printf(NORMAL, my_str % args)
    elif level == 'SEPARATOR':
        log_printf(SEPARATOR, my_str % args)
    elif level == 'HEADER':
        log_printf(HEADER, my_str % args)
    elif level == 'TITLE':
        log_printf(TITLE, my_str % args)
    elif level == 'WARNING':
        log_printf(WARNING, my_str % args)
    elif level == 'CRITICAL':
        log_printf(CRITICAL, my_str % args)
    elif level == 'RESULT':
        log_printf(RESULT, my_str % args)
    elif level == 'UNITTEST':
        log_printf(UNITTEST, my_str % args)
    elif level == 'ERROR':
        log_printf(ERROR, my_str % args)


##
# @brief Assigns the lowest level logging message.
# @details Sets the lowest level logging message to print to the screen.
#          This controls the lowest level for both logging messages in the
#          C++ source code as well as the user's OpenMOC Python input file.
#          This function would be called at the beginning of the input file
#          as follows:
#
# @code
#          log.py_setlevel('INFO')
# @endcode
#
# @param level the minimum logging level ('DEBUG', 'INFO', etc)
def py_setlevel(level):
    
    if level == 'DEBUG':
        log_setlevel(DEBUG)
    elif level == 'INFO':
        log_setlevel(INFO)
    elif level == 'NORMAL':
        log_setlevel(NORMAL)
    elif level == 'SEPARATOR':
        log_setlevel(SEPARATOR)
    elif level == 'HEADER':
        log_setlevel(HEADER)
    elif level == 'TITLE':
        log_setlevel(TITLE)
    elif level == 'WARNING':
        log_setlevel(WARNING)
    elif level == 'CRITICAL':
        log_setlevel(CRITICAL)
    elif level == 'RESULT':
        log_setlevel(RESULT)
    elif level == 'UNITTEST':
        log_setlevel(UNITTEST)
    elif level == 'ERROR':
        log_setlevel(ERROR)
    else:
        py_printf('Cannot set log level to unsupported log level %s', str(level))


