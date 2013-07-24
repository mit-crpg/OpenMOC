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


import openmoc


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
        openmoc.log_printf(openmoc.DEBUG, my_str % args)
    elif level == 'INFO':
        openmoc.log_printf(openmoc.INFO, my_str % args)
    elif level == 'NORMAL':
        openmoc.log_printf(openmoc.NORMAL, my_str % args)
    elif level == 'SEPARATOR':
        openmoc.log_printf(openmoc.SEPARATOR, my_str % args)
    elif level == 'HEADER':
        openmoc.log_printf(openmoc.HEADER, my_str % args)
    elif level == 'TITLE':
        openmoc.log_printf(openmoc.TITLE, my_str % args)
    elif level == 'WARNING':
        openmoc.log_printf(openmoc.WARNING, my_str % args)
    elif level == 'CRITICAL':
        openmoc.log_printf(openmoc.CRITICAL, my_str % args)
    elif level == 'RESULT':
        openmoc.log_printf(openmoc.RESULT, my_str % args)
    elif level == 'UNITTEST':
        openmoc.log_printf(openmoc.UNITTEST, my_str % args)
    elif level == 'ERROR':
        openmoc.log_printf(openmoc.ERROR, my_str % args)


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
def setLogLevel(level):

    if level == 'DEBUG':
        openmoc.setLogLevel('DEBUG')
    elif level == 'INFO':
        openmoc.setLogLevel('INFO')
    elif level == 'NORMAL':
        openmoc.setLogLevel('NORMAL')
    elif level == 'SEPARATOR':
        openmoc.setLogLevel('SEPARATOR')
    elif level == 'HEADER':
        openmoc.setLogLevel('HEADER')
    elif level == 'TITLE':
        openmoc.setLogLevel('TITLE')
    elif level == 'WARNING':
        openmoc.setLogLevel('WARNING')
    elif level == 'CRITICAL':
        openmoc.setLogLevel('CRITICAL')
    elif level == 'RESULT':
        openmoc.setLogLevel('RESULT')
    elif level == 'UNITTEST':
        openmoc.setLogLevel('UNITTEST')
    elif level == 'ERROR':
        openmoc.setLogLevel('ERROR')
    else:
        py_printf('Cannot set log level to unsupported log level %s', 
                  str(level))
