.. _debugging:

=================
Debugging OpenMOC
=================

This section describes some recommended debugging techniques for developers working with OpenMOC source code. Although many novice code developers are accustomed to use **printf debugging** with the ``printf(...)`` C/C++ routing. The printf debugging methodology is **highly discouraged** since it leads to a very long debug iteration cycle between source code revision, compilation, and testing. Instead, developers are advised to use some form of a debug environment for C/C++ and/or Python.

In particular, this section describes **GNU debugger** ``gdb`` in detail with some basic examples of how to use it with the shared library C++ extension module that is built for OpenMOC. In addition, some debugging techniques for the Python source are overviewed.

----------------------
The GNU Debugger (GDB)
----------------------

The `GNU Project Debugger`_ (GDB) provides the free open source ``gdb`` program for debugging of programs written in C/C++ (and other languages). The GNU debugger is typically installed when the GNU compiler collection (*e.g.*, ``gcc``, ``g++``) is installed. However, if for some reason you do not have ``gdb``, you may install it in Ubuntu using the ``apt-get`` package manager as follows:

.. code-block:: none

    sudo apt-get install gdb

Likewise, the GNU debugger can be installed using MacPorts on machines running Mac OS X as follows:

.. code-block:: none
   
    sudo port install gdb

There are a number of different ways to use to use and interact with GDB. There are a variety of graphical user interfaces for GDB, including DDD_ and `Eclipse CDT`_ which may be helpful for advanced code developers. This section, however, will focus on using the command line interface to GDB via the console.


---------------------------
Debugging C/C++ Source Code
---------------------------

The majority of the compute-intensive code in OpenMOC is written in C/C++ as a compiled back-end shared library. This section overviews some of the key steps to using GDB to debug OpenMOC's C/C++ source code. By no means does this section cover all of GDB's capabilities. Advanced developers should consult one of the many outlets for GDB documentation available online, including the following

* `Official GDB Documentation`_
* `GNU GDB Debugger Command Cheat Sheet`_
* `The Art of Debugging`_
* `Debugging with GDB`_


Installing OpenMOC with Debug Symbols
-------------------------------------

In order to debug OpenMOC C/C++ source code, you must build and install the C/C++ extension module for Python with `debug symbols`_. This is easily done by appending the :option:`-g` command line option for the ``gcc`` compiler. The build system for OpenMOC can automatically append this flag for the compilation process with the :option:`--debug-mode` option:

.. code-block:: none

    python setup.py install --user --debug-mode


Starting GDB
------------

The GNU Debugger may be started from the console with the ``gdb`` command:

.. code-block:: none

    $ gdb
    GNU gdb (GDB) 7.5-ubuntu
    Copyright (C) 2012 Free Software Foundation, Inc.
		License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
    This is free software: you are free to change and redistribute it.
    There is NO WARRANTY, to the extent permitted by law.  Type "show copying"
    and "show warranty" for details.
    This GDB was configured as "x86_64-linux-gnu".
    For bug reporting instructions, please see:
    <http://www.gnu.org/software/gdb/bugs/>.
    (gdb)

If you are debugging code written for the GPU, the CUDA GNU Debugger should be used:

.. code-block:: none

    $ cuda-gdb
    NVIDIA (R) CUDA Debugger
    5.5 release
    Portions Copyright (C) 2007-2013 NVIDIA Corporation
    GNU gdb (GDB) 7.2
    Copyright (C) 2010 Free Software Foundation, Inc.
    License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
    This is free software: you are free to change and redistribute it.
    There is NO WARRANTY, to the extent permitted by law.  Type "show copying"
    and "show warranty" for details.
    This GDB was configured as "x86_64-unknown-linux-gnu".
    For bug reporting instructions, please see:
    <http://www.gnu.org/software/gdb/bugs/>.
    (cuda-gdb)

The CUDA GDB debugger oftentimes does not fail at the exact location of the error, so it is recommended to turn cuda memcheck on:

.. code-block:: none

    (cuda-gdb) set cuda memcheck on

These commands will leave you in an interactive session with GDB, complete with its own command line interface. From now on we will give all examples with gdb, but the same commands apply for cuda-gdb.


Initializing Python
-------------------

From within the interactive GDB session, you must inform GDB that you plan to use Python as your binary executable:

.. code-block:: none

    (gdb) file python
    Reading symbols from /usr/bin/python...Reading symbols from /usr/lib/debug/usr/bin/python2.7...done.


Running Python in GDB
---------------------

Next you can select a Python file to execute using the ``run`` command with the Python script as the argument. For example, to run the :file:`/OpenMOC/sample-input/simple-lattice/simple-lattice.py` from within GDB, simply execute the following:

.. code-block:: none

    (gdb) run simple-lattice.py -i 5
    Starting program: /usr/bin/python simple-lattice.py -i 5
    [Thread debugging using libthread_db enabled]
    Using host libthread_db library "/lib/x86_64-linux-gnu/libthread_db.so.1".
    [  NORMAL ]  Importing materials data from HDF5...
    [  NORMAL ]  Creating surfaces...
    [  NORMAL ]  Creating cells...
    [  NORMAL ]  Creating simple 4 x 4 lattice...
    [  NORMAL ]  Creating geometry...
    [  NORMAL ]  Number of flat source regions: 60
    [  NORMAL ]  Initializing the track generator...
    [  NORMAL ]  Returning from readTracksFromFile
    [  NORMAL ]  Computing azimuthal angles and track spacings...
    [  NORMAL ]  Generating track start and end points...
    [  NORMAL ]  Segmenting tracks...
    [  NORMAL ]  Initializing track boundary conditions...
    [  NORMAL ]  Converging the source...
    [  NORMAL ]  Iteration 0: 	k_eff = 1.000000	res = 0.000E+00
    [  NORMAL ]  Iteration 1: 	k_eff = 1.270088	res = 2.594E+02
    [  NORMAL ]  Iteration 2: 	k_eff = 1.290540	res = 4.797E-01
    [  NORMAL ]  Iteration 3: 	k_eff = 1.309195	res = 1.919E-01
    [  NORMAL ]  Iteration 4: 	k_eff = 1.318423	res = 1.680E-01
    [ WARNING ]  Unable to converge the source after 5 iterations
    [  TITLE  ]  *******************************************************************
    [  TITLE  ]                             TIMING REPORT                           
    [  TITLE  ]  *******************************************************************
    [  RESULT ]  Total time to solution...............................3.9251E-03 sec
    [  RESULT ]  Solution time per unknown............................9.3454E-06 sec
    [  RESULT ]  Solution time per iteration..........................7.8502E-04 sec
    [  RESULT ]  Integration time per segment integration.............2.0098E-08 sec
    [SEPARATOR]  -------------------------------------------------------------------
    [  RESULT ]             # tracks          # segments          # FSRs
    [SEPARATOR]  -------------------------------------------------------------------
    [  RESULT ]                 116               930               60               
    [SEPARATOR]  -------------------------------------------------------------------
    [  NORMAL ]  Plotting data...
    [  TITLE  ]  *******************************************************************
    [  TITLE  ]                               Finished                             
    [  TITLE  ]  *******************************************************************
    [Inferior 1 (process 25820) exited normally]
    (gdb) 

To obtain more information about program execution, GDB can be run in verbose mode using the :option:`-v` optional argument to the ``run`` command:

.. code-block:: none

    (gdb) run -v simple-lattice.py
    ...
    [Inferior 1 (process 25820) exited normally]
    (gdb) 


Set Breakpoints
---------------

A **breakpoint** is an intentional stopping or pausing place in a program for debugging purposes. A breakpoint can be set using ``gdb`` using the ``breakpoint`` or ``br`` commands. The ``br`` command can be set at a specific line number in a specific source file. For example, if we wanted to run an OpenMOC program until the ``Material::setNumEnergyGroups(...)`` routine was called, we could set a breakpoint to the first line in that routine and then execute :file:`simple-lattice.py` as follows: 

.. code-block:: none

    (gdb) br Material.cpp:349
    No source file named Material.cpp.
    Make breakpoint pending on future shared library load? (y or [n]) y
    
    Breakpoint 1 (Material.cpp:349) pending.
    (gdb) run simple-lattice.py 
    Starting program: /usr/bin/python simple-lattice.py
    [Thread debugging using libthread_db enabled]
    Using host libthread_db library "/lib/x86_64-linux-gnu/libthread_db.so.1".
    [  NORMAL ]  Importing materials data from HDF5...

    Breakpoint 1, Material::setNumEnergyGroups (this=0x1b6ef80, num_groups=7)
		at src/Material.cpp:349
    349	    if (num_groups < 0)
    (gdb) 

As shown in the snippet above, the code executes until the breakpoint is reached and short summary of the source code is printed to the console. Alternatively, we could have set the breakpoint using the name of the routine instead:

.. code-block:: none

    (gdb) br Material::setNumEnergyGroups
    Function "Material::setNumEnergyGroups" not defined.
    Make breakpoint pending on future shared library load? (y or [n]) y
    
    Breakpoint 1 (Material::setNumEnergyGroups) pending.
    (gdb) run simple-lattice.py 
    Starting program: /usr/bin/python simple-lattice.py
    [Thread debugging using libthread_db enabled]
    Using host libthread_db library "/lib/x86_64-linux-gnu/libthread_db.so.1".
    [  NORMAL ]  Importing materials data from HDF5...

    Breakpoint 1, Material::setNumEnergyGroups (this=0x1b6ef80, num_groups=7)
		at src/Material.cpp:347
    347	void Material::setNumEnergyGroups(const int num_groups) {
    (gdb) 

In each case, a breakpoint is set and the program is executed until that line is reached. The entire program state is stored and the execution is simply interrupted until further notice to GDB is given by the user.


Set Watchpoints
---------------

A **watchpoint** is a *conditional breakpoint*, or a breakpoint that is only reached when a certain condition is met. The condition may be the reading, writing, or modification of a specific location in memory. For example, if we wanted to watch the value of the ``Material::_num_groups`` private class attribute we could place a watchpoint on it. First, we might start gdb and place a breakpoint on the ``Material::setNumEnergyGroups(...)`` routine as shown in the preceding section. Then we could place a watchpoint as follows:

.. code-block:: none

    (gdb) watch _num_groups
    Watchpoint 2: _num_groups
    (gdb) continue
    Continuing.
    Watchpoint 2: _num_groups
    
    Old value = 0
    New value = 7
    Material::setNumEnergyGroups (this=0x1b6ef80, num_groups=7)
		at src/Material.cpp:358
    358	    if (_data_aligned) {
    (gdb) 

As illustrated, GDB stepped through the program until ``_num_groups`` was modified or used and reported its value to the console. Note that this snippet made use of the ``continue`` command which is covered in th next section. At this point, it suffices to say that ``continue`` resumes program execution from a breakpoint until it a new breakpoint or watchpoint is reached.

GDB provides a variety of method to `set watchpoints`_ during a program's execution. For example, instead of placing a watchpoint on ``_num_groups``, we could instead have placed a watchpoint on the condition that ``_num_groups`` > 0 as follows:

.. code-block:: none

    (gdb) watch _num_groups > 0
    Watchpoint 2: _num_groups > 0
    (gdb) continue
    Continuing.
    Watchpoint 2: _num_groups > 0
    
    Old value = false
    New value = true
    Material::setNumEnergyGroups (this=0x1b6ef80, num_groups=7)
		at src/Material.cpp:358
    358	    if (_data_aligned) {
    (gdb)

In this case, the result of the conditional is reported to the screen.


Step through the Program
------------------------

This section highlights a few of the key commands which may be used to control program execution using GDB.

* **Continue**

  The ``continue`` or ``c`` command is used to instruct GDB to continue program execution until the next breakpoint or watchpoint is reached (*i.e.*, useful for loops). An optional integer argument :option:`<number>` may be given to ``continue`` to instruct GDB to ignore the current breakpoint some number of times.

  .. code-block:: none

     (gdb) continue <number>


* **Step**

  The ``step`` or ``s`` command will step to the next line of code. If the next line of code is a function call, the ``step`` command **will** step into the function. An optional integer argument :option:`<number>` may be given to ``step`` to instruct GDB to step through some number of lines.

  .. code-block:: none
		  
     (gdb) step <number>


* **Next**

  The ``next`` or ``n`` command will step to the next line of code.  If the next line of code is a function call, the ``step`` command **will not** step into the function. An optional integer argument :option:`<number>` may be given to ``next`` to instruct GDB to step through some number of lines (without entering functions).

  .. code-block:: none

     (gdb) next <number>


* **Until**

  The ``until`` command will continue processing until reaching a specifed line number :option:`<number>`. This is akin to setting a breakpoint which is only used once and which is immediately deleted following its first use.

  .. code-block:: none
		  
      (gdb) until <number>


* **Where**

  The ``where`` command will show which line number you are at and which function you are in.

  .. code-block:: none

      (gdb) where


Examine Variables
-----------------

The ``print`` or ``p`` command may be used to examine variables within some scope of the code with GDB. For example, if you were interested in the value of :option:`variable`, you might set a breakpoint at the entrance point to the code region of interest, and step through the region while printing the value as follows:

.. code-block:: none

    (gdb) print variable


Report the Debugger State
-------------------------

The ``info`` or ``i`` command may be used to report debugger state information to the console. For example, to list all breakpoints - including the file and line numbers where each is set - the ``info`` command is used with the :option:`breakpoints` option:

.. code-block:: none

    (gdb) info breakpoints

To list breakpoint numbers only, ``info`` command is used with the :option:`break` option:

.. code-block:: none

    (gdb) info break

Likewise, to list all watchpoints, the ``info`` command is used with the :option:`watchpoints` option:

.. code-block:: none

    (gdb) info watchpoints


Disable Breakpoints
-------------------

The ``disable`` command is used to disable breakpoints with GDB. When a breakpoint is disabled, it is still retained by GDB but is not used during program execution. The ``enable`` command may be used to continue using the breakpoint again at a later time. The following illustrates how to cancel breakpoints 1, 3, 4, 5, and 6, and re-enable breakpoints 4 and 5:

.. code-block:: none

    (gdb) disable 1 3-6
    (gdb) enable 4-5


Printing the Stack
------------------

The ``backtrace`` or ``bt`` command may be used to show the trace of the function the program is currently in:

.. code-block:: none

    (gdb) bt

The ``backtrace`` command can be particularly useful when debugging `segmentation faults`_. In particular, GDB may be used to run the program until the segmentation fault is reached. At this point, the use of ``backtrace`` will print the function call stack, showing where the segmentation faul occurred.

Debugging with OpenMP threads
-----------------------------

Most of the time, OpenMOC will be run with more than one thread, and bugs might not manifest themselves on all threads. The ``info`` command can be run to check what each thread is doing.

.. code-block:: none

    (gdb) info thread

To switch between threads and print variables local to each thread, one can use the ``thread`` command with the desired thread number.

.. code-block:: none

    (gdb) thread 4

Debugging with multiple MPI processes
-------------------------------------

Some simulations may require OpenMOC to be run on more than one computing node. Debugging those situations with GDB is rather difficult, but not impossible. xterm can be used to spawn terminals which host each MPI process. On a cluster, a sshx connection and an interactive node reservation is required.

.. code-block:: none

    (gdb)  mpirun -np 2 xterm -e gdb python

Another handy command is to tell GDB to start each process, and to run backtrace when it crashes/reaches a breakpoint using "-ex".

.. code-block:: none

    (gdb)  mpirun -np 2 gdb -batch -ex "run assembly-case.py" -ex "bt" python

Stop Program Execution
----------------------

The ``kill`` command may be used to stop a program's execution while keeping the GDB process running:

.. code-block:: none

    (gdb) kill


Exiting GDB
-----------

The ``quit`` or ``q`` command may be used to exit the ``gdb`` debugger and return to the console:

.. code-block:: none

    (gdb) quit()


----------------------------
Debugging Python Source Code
----------------------------

There are a a number of resources which one may use to debug Python code. Many popular `Integrated Development Environments`_ (IDEs) for Python include interactive visual debugging support, including `PyCharm`_, `Eclipse PyDev`_, and `Wing IDE`_. It is **highly recommended** that code developers use one of these IDEs for Python development and debugging. The `PyCharm`_ IDE is especially recommended for OpenMOC users developing input and data processing modules in Python. Although PyCharm is a commercial product, a community version is provided for free with many of the most essential features including the following:

* Syntax highlighting
* Auto-indentation
* Code formatting
* Code completion
* Line/block commenting
* Refactoring
* Python interpreter
* Integrated debugger

In addition, advanced developers should consult one of the many online outlets for documentation on debugging Python programs, including the following:

* `Debugging in Python`_
* `Interactive Debugging in Python`_


.. _printf: http://www.cplusplus.com/reference/cstdio/printf/
.. _GNU Project Debugger: https://www.gnu.org/software/gdb/
.. _MacPorts: http://www.macports.org/
.. _debug symbols: http://en.wikipedia.org/wiki/Debug_symbol
.. _DDD: http://www.gnu.org/software/ddd/
.. _Eclipse CDT: http://www.eclipse.org/cdt/
.. _Official GDB Documentation: http://www.gnu.org/software/gdb/documentation/
.. _GNU GDB Debugger Command Cheat Sheet: http://www.yolinux.com/TUTORIALS/GDB-Commands.html
.. _The Art of Debugging: http://www.nostarch.com/debugging.htm
.. _Debugging with GDB: http://www.amazon.com/Debugging-GDB-The-Source-Level-Debugger/dp/1882114884
.. _set watchpoints: https://sourceware.org/gdb/onlinedocs/gdb/Set-Watchpoints.html
.. _segmentation faults: http://en.wikipedia.org/wiki/Segmentation_fault
.. _Integrated Development Environments: http://en.wikipedia.org/wiki/Integrated_development_environment
.. _Pycharm: http://www.jetbrains.com/pycharm/
.. _Eclipse PyDev: http://pydev.org/
.. _Wing IDE: https://wiki.python.org/moin/Wing%20IDE 
.. _Debugging in Python: http://pythonconquerstheuniverse.wordpress.com/2009/09/10/debugging-in-python/
.. _Interactive Debugging in Python: http://www.onlamp.com/pub/a/python/2005/09/01/debugger.html
