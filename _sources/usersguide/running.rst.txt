.. _running:

===============
Running OpenMOC
===============

Once you have created your OpenMOC Python script, it is relatively straightforward to run the code. Instructions for running directly in the shell, or on clusters using the TORQUE_ and Cobalt_ job schedulers are described in the following sections.


----------------------------
Running OpenMOC in the Shell
----------------------------

If you wish to run your code from the shell, simply type the following in the console::

    python your-script.py <option1> <option2> ...



------------------------
Submitting a TORQUE_ Job
------------------------

To schedule a job to run your code on a cluster using the TORQUE job scheduler, you will need to write a submission script. This script will need to specify the desired compute resources (number of nodes and processors per node), the maximum walltime, among `other flags`_. An example script requesting a single node with 12 cores for one hour is presented below:

.. code-block:: guess 

   #!/bin/sh
   ###################################################
   # Specify nodes, processors per node and maximum
   # running time
   ###################################################

   #PBS -l nodes=1:ppn=12
   #PBS -l walltime=01:00:00
   
   ###################################################
   # Enter directory and set PATH
   ###################################################

   cd $PBS_O_WORKDIR
   PATH=$PBS_O_PATH

   ###################################################
   # Run OpenMOC
   ###################################################
 
   python your-script.py <option1> <option2> ...


To submit this script (``my-job.pbs``) your job to TORQUE, simply use the following command in the console::

    qsub my-job.pbs


.. _TORQUE: http://www.adaptivecomputing.com/products/open-source/torque/
.. _Cobalt: https://www.alcf.anl.gov/user-guides/cobalt-job-control
.. _other flags: https://www.osc.edu/supercomputing/batch-processing-at-osc/pbs-directives-summary

------------------------
Submitting a SLURM Job
------------------------

To schedule a job to run your code on a cluster using the SLURM job scheduler, you will need to write a submission script. This script will need to specify the desired compute resources (number of nodes and processors per node), the maximum walltime, among other flags. An example script requesting a single node with 12 cores for one hour is presented below:

.. code-block:: guess 

   #!/bin/sh
   ###################################################
   # Specify nodes, processors per node, maximum
   # running time and memory
   ###################################################

   #SBATCH -N 1
   #SBATCH -n 12
   #SBATCH --time=1:00:00
   #SBATCH --mem 120000  #in MB

   ###################################################
   # Select partition and customize name of output
   ###################################################

   #SBATCH -p sched_mit_nse
   #SBATCH -o stderror-%j.out

   ###################################################
   # Run OpenMOC
   ###################################################
 
   python your-script.py <option1> <option2> ...


To submit this script (``my-job``) your job to SLURM, simply use the following command in the console:

    sbatch my-job

To run your script interactively instead, simply use the following command in the console:

    srun -p <partition_name> -I -N 1 --pty -t 0:60:00 --mem 124000 python your-script.py <option1> ...

------------------------
Submitting a Cobalt_ Job
------------------------

Cobalt is a job scheduler that is used on large scale machines at Argonne National Laboratory Leadership Computing Facility (ALCF_), namely the IBM BlueGene machines, such as Mira_. To schedule a job to run your code on a cluster using Cobalt, you can submit your job using :program:`qsub` with options to specify the desired compute resources (number of nodes and processors per node), the maximum walltime, etc.

An example command to request a single node for 20 minutes is presented below::

  qsub -A <your-alcf-account> -n 1 -t 20 --mode=c1 --env PYTHONPATH=<path-to-_openmoc.so>:<path-to-any-other-openmoc-shared-library-file> <path-to-python>/python your-script.py <option1> <option2> ... 

NOTE: You must specify the path to the location where OpenMOC was installed in the :envvar:`PYTHONPATH` environment variable. For OpenMOC ditributions built with the :option:`--user` option, this location will be `~/.local/lib/pythonX.X/site-packages/...`. In particular, the :envvar:`PYTHONPATH` given to :program:`qsub` must include the path to all OpenMOC shared libraries (files with a ``.so`` extension) needed for your script (*i.e.*, ``_openmoc.so`` for the ``openmoc`` Python module).


.. _ALCF: http://www.alcf.anl.gov/
.. _Mira: https://www.alcf.anl.gov/mira


.. _runtime_options:

---------------
Runtime Options
---------------

This section provides a brief overview of each of the runtime options that are supported by OpenMOC. The runtime arguments and the method to retrieve them in a Python script are elaborated upon further in :ref:`Simulation Parameters <simulation_params>`. Each of these can be passed into your Python script as follows::

    python your-script.py <option1> <option2> ...

.. option:: -h, --help

Reports all supported OpenMOC runtime options to the screen.


.. option:: -a, --num-azim=<4>

The number of azimuthal angles for ray tracing. The default is 4.


.. option:: -s, --azim-spacing=<0.1>

The track azimuthal spacing (in cm) for ray tracing. The default is 0.1 cm.


.. option:: -p, --num-polar=<6>

The number of polar angles for 3D ray tracing. The default is 3/6 (2D/3D).


.. option:: -l, --polar-spacing=<1.5>

The track spacing in the axial direction (in cm) for ray tracing. The default is 1.5 cm.


.. option:: -i, --max-iters=<1000>

The maximum number of source convergence iterations. The MOC solvers will execute as many iterations needed to converge the source, up to this limiting value. The default is 1000.


.. option:: -c, --tolerance=<1E-5>

The tolerance on the convergence. The default is 1E-5.


.. option:: -t, --num-omp-threads=<# CPU cores>

The number of OpenMP threads to use. This option only applies to scripts which use OpenMOC's :cpp:class:`CPUSolver`, or derived classes such as :cpp:class:`VectorizedSolver`. The default is the number of CPU cores available on the machine of interest.


.. option:: -b, --num-threadblocks=<64>

The number of CUDA threadblocks. This option only applies to scripts which use OpenMOC's :cpp:class:`GPUSolver` class. The default is 64 threadblocks.


.. option:: -g, --num-threads-per-block=<64>

The number of CUDA threads per threadblock. This option only applies to scripts which use OpenMOC's :cpp:class:`GPUSolver` class. This value must be a multiple of the number of threads in a CUDA warp_. At the time of this writing, nearly all NVIDIA GPUs have a warp size of 32, though this may change for future NVIDIA GPUs. If the value set for this option is not a multiple of 32, the CUDA source code will round up to the nearest multiple of 32 threads. The default is 64 threads. 

.. note:: If you are unsure what the warp size is for your GPU, you can use the ``openmoc.cuda`` module to find out. The following Python code will report the warp size for your GPU to the console:

.. code-block:: python

     import openmoc.cuda as cuda

     if cuda.machine_contains_gpu():
         num_threads = cuda.get_num_threads_per_warp()
         print 'This machines GPU contains %d threads per warp' % (num_threads)
  
     else:
         print 'This machine does not contain an NVIDIA CUDA-enabled GPU'


-------------------------------
Canceling an OpenMOC Simulation
-------------------------------

To cancel an OpenMOC job running in your shell, you can use the ``CTRL+C`` keyboard combination. This will kill the Python script as well as the underlying computation running in the C/C++/CUDA shared library.


.. _warp: http://www.pgroup.com/lit/articles/insider/v2n1a5.htm
