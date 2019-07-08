.. _style_guide:

=======================
Style Guide for OpenMOC
=======================

In order to keep the OpenMOC code base consistent in style, this guide specifies
a number of rules which should be adhered to when modified existing code or
adding new code in OpenMC. This guide should mostly agree with the `Google style guides`_
for Python and C++.

-------------
General Rules
-------------

* **Compilation** - Make sure code can be compiled with most common compilers, especially gcc and the Intel C++ compiler icpc.

* **Dependencies** - Avoid using special extensions or libraries which add dependencies to the OpenMOC source.

* **Code Comments** - Always include comments to describe what your code is doing. Do not be afraid of using copious amounts of comments.


------------------
Source Code Syntax
------------------

The following describes some of the rules for source code syntax, including the class naming convention, rules for class methods and attributes as well as normal routines and variables.


Class Names
-----------
Class names are CamelCase_ starting with a **capital** letter and without **"_"** between words as illustrated in the C++ snippet below:

.. code-block:: c

   /* Yes */
   class MyClass { 
     ... 
   };
   
   /* No */
   class myClass { 
     ... 
   };

   /* No */
   class My_Class{ 
     ... 
   };


Class Methods
-------------
Class methods are camelCase_ starting with a **lowercase** letter and without **"_"** between words as illustrated in the Python snippet below:

.. code-block:: python

   class MyClass:
   
     # Yes
     def getterMethod(self):
       ...

     # No
     def Setter_Method(self):
       ...

     # No
     def logic_method(self):
       ...


Class Attributes
----------------
Class attributes are **private** and **lowercase** which start with **"_"** and with "_" between words as illustrated in the C++ snippet below:

.. code-block:: c

   class MyClass {

     private:

       /* Yes */
       int _first_attribute;

       /* No */
       double second_attribute;
       
       /* No */
       float _ThirdAttribute;

       ...

   };


Function Names
--------------
Functions (not class methods) are all **lowercase** with **"_"** between words as illustrated in the Python snippet below:

.. code-block:: python

   # Yes
   def my_function(a, b):
     ...

   # No
   def myFunction(a, b):
     ...

   # No
   def My_Function(a, b):
     ...


Variable Names
--------------
Temporary variables (e.g. variables defined within the scope of a function) are **lowercase** with **"_"** between words as illustrated in the C/C++ snippet below:

.. code-block:: c

   void my_function(int a, int b) {

     /* Yes */
     int first_variable;

     /* No */
     int secondVariable;

     /* No */
     int Second_Variable;

     ...

   }


.. _code_comments:

-------------
Code Comments
-------------

OpenMOC uses Doxygen_ for automated generation of Application Programming Interface (API_) documentation based upon code comments. Please adhere to the Doxygen standard for code comments in both C/C++ and Python source code. In particular, 

* **C/C++ functions** should be preceded by Doxygen-style comments as illustrated in the following code snippet:

.. code-block:: c

    /**
     * @brief Single precision A*X plus Y
     * @details Multiplies single precision vectors a and x and 
     *          adds vector y. The output vector is stored in y.
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


* **C/C++ class** definitions in a header file should be preceded by Doxygen-style comments as illustrated in the following code snippet:

.. code-block:: c

   /**
    * @file Person.h
    * @brief The Person class.
    * @date January 23, 2014.
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



* **Python functions** should be preceded by Doxygen-style comments as illustrated in the following code snippet:

.. code-block:: python

    ##
    # @brief Single precision A*X plus Y
    # @details Multiplies single precision vectors a and x and 
    #          adds vector y. The output vector is stored in y.
    # @param n the size of the vectors x and y
    # @param a a single precision vector of length n
    # @param x a single precision vector of length n
    # @param y a single precision vector of length n
    def saxpy(n, a, x, y):
      for i in range(n):
        y[i] = a[i] * x[i] + y[i]


* **Python classes** should be preceded by Doxygen-style comments as illustrated in the following code snippet:

.. code-block:: python

   ##
   # @file Person.py
   # @brief The Person class.
   # @date January 23, 2014.
   # @author William Boyd, MIT, Course 22 (wboyd@mit.edu)

   ...

   ##
   # @class Person Person.y "openmoc/src/Person.y"
   # @brief The Person class represents a generic human being.
   # @details The Person class contains attributes representing
   #          a generic human being, such as name and age.
   class Person:

     ##
     # @brief The default person constructor.
     # @return A handle to a person object.
     def __init__(self):
       _name = None
       _age = None

     ##
     # @brief Assigns a person an age.
     # @param age The person's age (years)
     def setAge(self, age):
       _age = age

     ##
     # @brief Assigns a person a name.
     # @param name The person's name (a string)
     def setName(self, name):
       _name = name

     ##
     # @brief Retrieves the person's age.
     # @return The person's age (years).
     def getAge(self):
       return _age

     ##
     # @brief Retrieves the person's name.
     # @return The person's name.
     def getName(self):
       return _name


-----------
Line Length
-----------

For readability, source code in OpenMOC is limited to a maximum of **80 characters** for each line. For your convenience in adhering to this policy, you mayupdate your text editor (gedit, emacs, vim, etc.) to display the right margin at column 80.


-----------
Indentation
-----------

For readability, OpenMOC uses tabs composed of **2 white spaces** per indentation level. For your convenience in adhering to this policy, you may update your text editor (gedit, emacs, vim, etc.) preferences to use a tab width of 2 spaces and to insert spaces instead of tabs. Emacs users should include the following line in their .emacs file:

.. code-block:: common-lisp

    (setq-default indent-tabs-mode nil)

vim users should include the following line in their .vimrc file::

  set expandtab


----------
Whitespace
----------
      
Use a **single space** between arguments to procedures.

Avoid extraneous whitespace in the following situations:

  * In function calls::

      Yes: myfunc(x, y(2), z)
      No: myfunc ( x, y( 2 ), z )

  * In logical expressions, use one space around operators but nowhere else::

      Yes: if(variable == 2) then
      No: if ( variable==2 ) then


.. _Doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _sphinx: http://sphinx-doc.org/
.. _api: http://en.wikipedia.org/wiki/Application_programming_interface
.. _CamelCase: http://en.wikipedia.org/wiki/CamelCase
.. _camelCase: http://en.wikipedia.org/wiki/CamelCase
.. _Google style guide: http://google.github.io/styleguide
