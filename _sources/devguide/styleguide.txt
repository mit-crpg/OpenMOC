.. _devguide_styleguide:

=======================
Style Guide for OpenMOC
=======================

In order to keep the OpenMOC code base consistent in style, this guide specifies
a number of rules which should be adhered to when modified existing code or
adding new code in OpenMC.

-------------
General Rules
-------------


------------------
Source Code Syntax
------------------

    * **Class Methods**


    * **Class Attributes**


    * **Function Names**


    * **Variable Names**


    * **Code Comments**

      OpenMOC uses doxygen_ for automated generation of Application Programming Interface (API_) documentation based upon code comments. Please adhere to the doxygen standard for code comments in both C/C++ and Python source code. In particular, 

      - C/C++ functions should be preceded by doxygen-style comments as illustrated in the following code snippet:

	.. code-block:: c

	   /**
            * @brief Single precision A*X plus Y
            * @details Multiplies a single precision vectors a
            *          and x and adds vector y. The output
            *          vector is stored in y.
            * @param n the size of the vectors x and y
            * @param a a single precision vector of length n
            * @param x a single precision vector of length n
            * @param y a single precision vector of length n
            */
	   int saxpy(int n, float* a, float* x, float* y) {
	       for (int i=0; i<n; i++)
	           y[i] = a[i] * x[i] + y[i];

               return;
	   }

      - C/C++ class definitions in a header file should be preceded by doxygen-style comments as illustrated in the following code snippet:

	.. code-block:: c

	   /**
            * @file Person.h
            * @brief The Person class.
            * @data January 23, 2014.
            * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
            */

	    ...
	    
	   /**
            * @class Person Person.h "openmoc/src/Person.h"
            * @brief The Person class represents a generic human being.
            * @details The Person class contains attributes representing
            *          a generic human being, such as name and age.
            */
	   class Person {

	       private:
	           /** The person's name */
	           char* _name;

		   /** The person's age */
	           int _age;
		   ...

	       public:
		   char* getName();
	           int getAge();
		   ...
	   }



      - Python functions should be preceded by doxygen-style comments as illustrated in the following code snippet:

      - Python classes should be preceded by doxygen-style comments as illustrated in the following code snippet:


    * **Line Length**

      For readability, source code in OpenMOC is limited to a maximum of 80 characters for each line. For your convenience in adhering to this policy, you mayupdate your text editor (gedit, emacs, vim, etc.) to display the right margin at column 80.


    * **Indentation**
      
      For readability, OpenMOC uses tabs composed of 4 white spaces per indentation level. For your convenience in adhering to this policy, you may update your text editor (gedit, emacs, vim, etc.) preferences to use a tab width of 4 spaces and to insert spaces instead of tabs. Emacs users should include the following line in their .emacs file:

      .. code-block:: common-lisp

	 (setq-default indent-tabs-mode nil)

      vim users should include the following line in their .vimrc file::

          set expandtab


    * **Whitespace**
      
      Use a single space between arguments to procedures.

      Avoid extraneous whitespace in the following situations:

      - In function calls::

	  Yes: myfunc(x, y(2), z)
	  No: myfunc ( x, y( 2 ), z )

      - In logical expressions, use one space around operators but nowhere else::

	  Yes: if(variable == 2) then
	  No: if ( variable==2 ) then


.. _doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _sphinx: http://sphinx-doc.org/
.. _api: http://en.wikipedia.org/wiki/Application_programming_interface
