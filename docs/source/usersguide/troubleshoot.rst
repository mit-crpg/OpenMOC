.. _usersguide_troubleshoot:

===============
Troubleshooting
===============

If you are experiencing problems trying to compile OpenMOC, first check if the
error you are receiving is among the following options.

If not, you may email the `mailing list`_, and a developer or another user may help you.


-------------------------
Problems with Compilation
-------------------------

NameError: ... is not defined
-----------------------------

In some cases, the symbolic links for the code are not established correctly and users may receive a `NameError`, and example of which is given below::

  Traceback (most recent call last):
    File "simple-lattice.py", line 1, in <module>
      from openmoc import *
    File "/home/wbinventor/.local/lib/python2.7/site-packages/openmoc/__init__.py", line 15, in <module>
      setLogfileName('log/openmoc-' + current_time + '.log');
  NameError: name 'setLogfileName' is not defined

This error is commonly seen by new users after just installing the code. As noted in the :ref:`Installation and Configuration <install>`, some Python distributions require that OpenMOC be installed twice the first time it is installed. The recommended solution in this case would be to execute the installation command twice in sequence::

  python setup.py install --user
  python setup.py install --user



.. _mailing list: https://groups.google.com/forum/?hl=en#!forum/openmoc-users
