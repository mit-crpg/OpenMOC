import openmoc


def py_printf(level, my_str, *args):
    """Print a log message to the screen and the log file.

    This method is a wrapper to the log_printf C++ routine. It allows for
    formatted messages to be printed to the screen in a similar fashion to the
    C/C++ printf method, but with additional formatting provided by the OpenMOC
    logging utilities.

    Parameters
    ----------
    level : str
        The logging level for this message (i.e., 'NORMAL')
    my_str : str
        The string to print to the screen
    *args : list
        A variable length list of values for the message string

    Examples
    --------
    An example of how this may be used in a OpenMOC Python script is as follows:

        >>> value1 = 25
        >>> value2 = 26.0
        >>> log.py_printf('NORMAL', 'My name is Will and I am %d going on ' \
                         '%f years of age', value1, value2)

    """

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
    elif level == 'ERROR':
        openmoc.log_printf(openmoc.ERROR, my_str % args)
    else:
        openmoc.log_printf(openmoc.ERROR, "Unknown message log level.")


def set_log_level(level):
    """Assign the lowest level logging message to be reported.

    Sets the lowest level logging message to print to the screen. This controls
    the lowest level for both logging messages in the C++ source code as well as
    the user's OpenMOC Python input file.

    Parameters
    ----------
    level : str
        The minimum logging level (i.e., 'DEBUG', 'INFO')

    Examples
    --------
    This routine may be called in an OpenMOC Python script as follows:

        >>> log.set_log_level('INFO')

    """

    if level == 'DEBUG':
        openmoc.set_log_level('DEBUG')
    elif level == 'INFO':
        openmoc.set_log_level('INFO')
    elif level == 'NORMAL':
        openmoc.set_log_level('NORMAL')
    elif level == 'SEPARATOR':
        openmoc.set_log_level('SEPARATOR')
    elif level == 'HEADER':
        openmoc.set_log_level('HEADER')
    elif level == 'TITLE':
        openmoc.set_log_level('TITLE')
    elif level == 'WARNING':
        openmoc.set_log_level('WARNING')
    elif level == 'CRITICAL':
        openmoc.set_log_level('CRITICAL')
    elif level == 'RESULT':
        openmoc.set_log_level('RESULT')
    elif level == 'ERROR':
        openmoc.set_log_level('ERROR')
    else:
        py_printf('Cannot set log level to unsupported level %s', str(level))
