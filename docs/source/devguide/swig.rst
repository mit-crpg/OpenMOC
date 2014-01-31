.. _swig:

=============================================
Simplified Wrapper Interface Generator (SWIG)
=============================================


[THIS SITE IS UNDER CONSTRUCTION]


--------------
SWIG Wrap File
--------------

C/C++ Wrap File
---------------

Python Bindings File
--------------------


--------------------------
SWIG from the Command Line
--------------------------

SWIG is provided as the ``swig`` executable and is called on the command line to wrap C/C++ source code given a SWIG interface file (see :ref:`SWIG Interface File <swig_interface_file>`). The call to ``swig`` produces a SWIG wrap file, which is designated following the :option:`-o` argument. The language used for the input source files (*e.g.*, C/C++) as well as the language targeted for the bindings (*e.g.*, Python) must be specified as flags to the ``swig`` executable. An example of how one might issue a command from the shell to wrap the source code in the :file:`my_module.i` SWIG interface file could be the following::

  swig -python -c++ -o my_module_wrap.cpp my_module.i

------------------
Default Parameters
------------------

It is highly recommended that developers make use of `default parameters`_ for routines when possible. Default arguments in a C++ routine are wrapped by ``swig`` given the :option:`-keyword` to provide `keyword arguments`_ (also known as named parameters) in the Python binding for that routine. There are several benefits for defining default arguments in the C/C++ source code:

* **Readability** - Keyword arguments make code more readable, especially in example input files for users new to OpenMOC

* **Ordering** - Keyword arguments can be provided in any order, lessening the burden to the user to remember a specific ordering

* **Flexibility** - Function parameters with useful default values are not required at runtime, making Python scripts easier to comprehend

An example of function parameters with default values and the use of keyword arguments to override the default values in Python is given below:

.. code-block:: python

    # Define a function with two arguments with default values
    def multiverseHelloWorld(count=5, greeting='Hello'):

      for i in range(count):
        print '%s from World %d!' % (greeting, i)

    # Call routine and override default keyword arguments 
    # The keyword arguments can be provided in any order
    multiverseHelloWorld(greeting='Hola', count=7)

Likewise, an example of how to define default values for function parameters - which will be provided through the Python interface as `SWIG default arguments`_ - in C/C++ is given below:

.. code-block:: c

    /* Define a function prototype with two arguments with default values */
    void multiverseHelloWorld(int count=5, char* greeting="Hello");

    /* Function implementation doesn't include default values */
    void multiverseHelloWorld(int count, char* greeting) {

      for (int i=0; i < count; i++)
        printf("%s from World %d!", greeting, i)
    }



------------------
Exception Handling
------------------


--------------
NumPy Typemaps
--------------

It is often useful to input/return NumPy data structures to/from C/C++ routines. The `NumPy C API`_ makes this functionality possibility through **array conversions**. In addition, it is possible to automatically *embed* the NumPy C API directly into the source code with the use of `NumPy typemaps`_. Typemaps are a mechanism to match **function signatures** through a list of of function arguments. When SWIG finds a function which matches the typemap, it will target and subsequently modify the function to include the NumPy C API in order to input/output NumPy data arrays. Two types of parameters must be specified in the C/C++ function(s) of interest in order to match a NumPy typemap:

* **Array Pointer** - The data type and pointer to the array
* **Array Dimensions** - An integer for each array dimension

The :file:`numpy.i` interface file defines the typemaps and is shipped with OpenMOC in the :file:`/OpenMOC/openmoc` directory. In order to utilize NumPy typemaps, the following should be appended to the SWIG interface file used for the C/C++ extension module of interest:

.. code-block:: bash

    %include "numpy.i"

    %{
      #define SWIG_FILE_WITH_INIT
    %}

    %init %{
      import_array();
    %}

The following sections overview the basic steps to utilize NumPy typemaps to input NumPy data from Python into C/C++ routines, and return data from C/C++ routines as NumPy arrays.


Input NumPy Data Arrays
-----------------------

The :file:`numpy.i` interface file provides two particular typemaps for inputting a NumPy data array into a C/C++ routine. The :envvar:`IN_ARRAY*` defines an array which is passed into a routine but is not modified in-place by the C/C++ function and is not returned to the user. The :envvar:`INPLACE_ARRAY*` typemap defines arrays that are modified in-place. In each case, the :envvar:`*` represents the number of dimensions for the input array. For example, in order to input a 3D array to be modified in-place, one would use the :envvar:`INPLACE_ARRAY3` typemap. The array dimension(s) are included in each typemap through the use of the :envvar:`DIM*` parameter. 

The following is an example C/C++ in which which we wish to wrap some function ``sum_array(...)`` with SWIG and provide the capability to input a NumPy array as a function parameter. Note that the function prototype includes a first paramter for the pointer to the input double array and a second parameter for the length of the array (which together form the function signature). The function prototype is given below in the :file:`sum_array.h` file below:

.. code-block:: c

    /* Define function prototype to take in a NumPy array */
    double sum_array(double* input_array, int length);

One possible implementation of the ``sum_array(...)`` routine is given in the :file:`sum_array.c` file below:

.. code-block:: c

    /* Define function implementation */
    double sum_array(double* input_array, int length) {

      /* Initialize sum */
      double sum = 0.;

      /* Compute sum of array elements */
      for (int i=0; i < length; i++)
        sum += input_array[i];

      return sum;
    }

The following would be the required SWIG interface file to wrap :file:`sum_array.h` into the ``_sum_array`` C/C++ extension module for Python. The second-to-last line defines the NumPy typemap - the first tuple is a pair of the typemap (array type and dimension) while the second is the function signature to match using that typemap.

.. code-block:: bash

    %module sum_array

    %{
      #define SWIG_FILE_WITH_INIT
      #include "sum_array.h"
    %}

    %include "numpy.i"

    %init %{
      import_array();
    %}

    %apply (double* IN_ARRAY1, int DIM1) {(double* input_array, int length)};

    %include "sum_array.h"

After ``swig`` is used to generate the wrap file and it is compiled into the :file:`_sum_array.so` shared library, the module may be imported into Python and the routine used with a NumPy array as follows:

.. code-block:: python

    from numpy.random import rand 
    from _sum_array import *

    # Initialize a random NumPy array
    input_array = rand(5)

    # Sum the values in the random NumPy array
    sum = sum_array(input_array)

.. note:: More detailed information on :envvar:`IN_ARRAY` and :envvar:`INPLACE_ARRAY` typemaps is provided in the official `NumPy.i`_ documetation.


Return NumPy Data Arrays
------------------------

The :file:`numpy.i` interface file also provides two typemaps for returning a NumPy data array from a C/C++ routine. The :envvar:`ARGOUT_ARRAY*` used in situations where you would allocate an array on the heap and call the function to fill the array values. In Python, the arrays are allocated for you and returned as new array objects. As was the case for array input, the :envvar:`*` represents the number of dimensions for the input array. For example, in order to input a 3D array to be modified in-place, one would use the :envvar:`ARGOUT_ARRAY3` typemap. The array dimension(s) are included in each typemap through the use of the :envvar:`DIM*` parameter. 

The following is an example C/C++ in which which we wish to wrap some function ``get_rand_array(...)`` with SWIG and provide the capability to convert a C/C++ array into an output NumPy array. Based on the function signature, it would appear that the output array is input into the function and nothing is returned. Instead, SWIG will modify the source code with the NumPy C API such that a NumPy array is initialized and input as a C/C++ array and subsequently returned as a NumPy array.

The function prototype is given below in the :file:`get_rand_array.h` file below:

.. code-block:: c

    #include <stdlib.h>

    /* Define function prototype to take in a NumPy array */
    void get_rand_array(double* output_array, int length);

One possible implementation of the ``get_rand_array(...)`` routine is given in the :file:`get_rand_array.c` file below:

.. code-block:: c

    /* Define function implementation */
    double get_rand_array(double* output_array, int length) {

      /* Populate input NumPy array with random numbers */
      for (int i=0; i < length; i++)
        output_array[i] = rand();

      return;
    }

The following would be the required SWIG interface file to wrap :file:`get_rand_array.h` into the ``_get_rand_array`` C/C++ extension module for Python. The second-to-last line defines the NumPy typemap - the first tuple is a pair of the typemap (array type and dimension) while the second is the function signature to match using that typemap.

.. code-block:: bash

    %module get_rand_array

    %{
      #define SWIG_FILE_WITH_INIT
      #include "get_rand_array.h"
    %}

    %include "numpy.i"

    %init %{
      import_array();
    %}

    %apply (double* ARGOUT_ARRAY1, int DIM1) {(double* output_array, int length)};

    %include "get_rand_array.h"

After ``swig`` is used to generate a wrap file and it is compiled into the :file:`_get_rand_array.so` shared library, the module may be imported into Python and the routine used as follows:

.. code-block:: python

    from numpy.random import rand 
    from _get_rand_array import *

    # Sum the values in the random NumPy array
    length = 100
    output_array = get_rand_array(length)

.. note:: More detailed information on the :envvar:`ARGOUT_ARRAY` typemap is provided in the official `NumPy.i`_ documetation.


-------------
SWIG Typemaps
-------------

:file:`typemaps.i` which is included in the standard SWIG installation.
:option:`--no-numpy` flag for machines where NumPy is not available.


------
Macros
------

SWIG provides preprocessing_ capabilities for interface files. `Macro expansions`_ may be defined in the interface file using the traditional syntax for C/C++:

.. code-block:: c

    #ifndef PI
    #define PI 3.14159
    #endif

SWIG also includes special `SWIG macros`_ with more enhanced capabilities for interface files.


--------
Typedefs
--------

SWIG provides functionality to define typedefs_ in interface files. `SWIG typedefs`_ can be defined using the same syntax as in C/C++. As discussed in the SWIG online documentation, the typedef must be defined twice in the interface file for in order for it to be propagated to the generated wrapper file:

.. code-block:: bash

   %{

     /* Include in the generated wrapper file */
     typedef unsigned int size_t;

   %}

   /* Tell SWIG about it */
   typedef unsigned int size_t;




.. _SWIG Project: http://www.swig.org/
.. _NumPy typemaps: http://docs.scipy.org/doc/numpy/reference/swig.interface-file.html
.. _Numpy.i: http://docs.scipy.org/doc/numpy/reference/swig.interface-file.html
.. _NumPy C API: http://docs.scipy.org/doc/numpy/reference/c-api.html
.. _SWIG typemaps: http://www.swig.org/Doc1.3/Typemaps.html
.. _SWIG Basics Tutorial: http://www.swig.org/Doc1.3/SWIG.html#SWIG_nn20
.. _SWIG C++ Tutorial: http://www.swig.org/Doc1.3/SWIGPlus.html

.. _SWIG default arguments: http://www.swig.org/Doc1.3/SWIGPlus.html#SWIGPlus_default_args
.. _default parameters: http://www.learncpp.com/cpp-tutorial/77-default-parameters/
.. _keyword arguments: http://en.wikipedia.org/wiki/Named_parameter
.. _preprocessing: http://en.wikipedia.org/wiki/C_preprocessor
.. _Macro expansions: http://www.swig.org/Doc1.3/Preprocessor.html#Preprocessor_nn5
.. _SWIG macros: http://www.swig.org/Doc1.3/Preprocessor.html#Preprocessor_nn6
.. _typedefs: http://en.wikipedia.org/wiki/Typedef
.. _SWIG typedefs: http://www.swig.org/Doc1.3/SWIG.html#SWIG_nn20
