.. _documentation:

=============
Documentation
=============

The following describes the suggested steps that developers should take to edit, revise or add to OpenMOC's documentation. OpenMOC uses the Sphinx_ tool for code due to its growing popularity for open source codes. Sphinx uses reStructuredText_ (rst) files to encapuslate text-based content with semantic markup for figures, tables, URLs, and more. Sphinx provides a compiler which parses the reStructuredText source files to generate a target output in the form of HTML or a PDF. OpenMOC also uses Doxygen_ for automated generated of Application Programming Interface (API) documentation based on source code comments. The used of Doxygen-style comments is explained further in :ref:`Code Comments <code_comments>`. 


-------------------------
Preliminary Configuration
-------------------------

First, in order to build OpenMOC's documentation, you must have installed both Sphinx and Doxygen. Since both tools are widely used, they are available with most source package system. On Debian-based systems, such as Ubuntu, you may install both tools with the following command in the console:

.. code-block:: none

    sudo apt-get install python-sphinx doxygen

Likewise, on Mac OS with MacPorts_ you may install each with the following command:

.. code-block:: none

    sudo port install py27-sphinx doxygen

Next, it is important to note that GitHub_ reserves the use of the **gh-pages** branch to host a project website for each repository. OpenMOC makes uses of this feature for its code documentation website. However, the use of Sphinx coupled with Doxygen - which requires access to the source code on the ``master`` branch - presents a challenge to properly configure the ``gh-pages`` branch in GitHub. A nice workaround has been developed, however, based on `Rick Foos' blog post`_ blog post on the issue. The following steps will guide you through the workaround such that you will be able to easily edit and revise the documentation.

First, you need to create a :file:`/OpenMOC/sphinx` directory inside of the repository on your machine:

.. code-block:: none

    cd OpenMOC
    mkdir sphinx

Note that the :file:`sphinx` directory is intentionally in the :file:`.gitignore` file in the OpenMOC repository. As a result, the directory will not be tracked by Git in the traditional sense, which is what we want.

Next, enter the :file:`/OpenMOC/sphinx` directory and clone the OpenMOC repository there:

.. code-block:: none

    cd sphinx
    git clone https://github.com/mit-crpg/OpenMOC.git .

You have effectively created a `nested repository` of OpenMOC within the original local OpenMOC repository. Now, within the nested repository in the :file:`/OpenMOC/sphinx` directory, checkout the ``gh-pages`` branch, ignoring the warning message:

.. code-block:: none

    git checkout origin/gh-pages -b gh-pages

Now, your nested repository is on the ``gh-pages`` branch within the original repository on its original branch (*i.e.*, the ``master`` branch). Next, delete the ``master`` branch within the nested repository in :file:`/OpenMOC/sphinx`, ignoring the warning:

.. code-block:: none

    git branch -d master

Congratulations! You are now configured to edit and revise the documentation, as discussed in the following section.

.. note:: If you clone the repository to a different directory on the same machine, or to a new machine, you will have to repeat these steps in order to configure the newly cloned repository since the :file:`/OpenMOC/sphinx` directory will not be tracked within the :file:`master` branch.


-------------------------
Writing New Documentation
-------------------------

This tutorial will not present any of the semantics of how to write documentation for Sphinx as that is well covered by its own tutorials. Likewise, the same is true for Doxygen-style comments, although the basics are presented in :ref:`Code Comments <code_comments>`. Instead, this section presents the location of the source code and the basic console commands used to generate HTML documentation.

First, the :file:`/sphinx/docs/source` directory contains the reStructuredText source for the Sphinx documentation. In order to compile the Sphinx documentation into HTML, you can run the following from within the :file:`sphinx/docs` directory:

.. code-block:: none

    cd sphinx/docs
    make html

This command will use the Sphinx tool to generate HTML in the :file:`/sphinx/docs/build` directory. You may use a web browser, such as Mozilla Firefox, to view your newly generated documentation as follows:

.. code-block:: none

    cd sphinx/docs/build
    firefox index.html

Likewise, in order to build a revised version of the :ref:`OpenMOC API <api>` from Doxygen-style comments in the Python and C++ source code, you may use the following console command from within the :file:`sphinx/docs/doxygen` directory:

.. code-block:: none

    cd sphinx/docs/doxygen
    doxygen Doxyfile

This will build a new version of the API in the :file:`sphinx/doxygen/html` directory. You may use a web browser, such as Mozilla Firefox, to view your newly generated documentation as follows:

.. code-block:: none

    cd sphinx/doxygen/html
    firefox index.html

In order to **expose** your newly generated HTML documentation to GitHub such that it will present your website, you must copy the HTML into the main :file:`sphinx` directory:

.. code-block:: none

    cd sphinx
    cp -rf docs/build/html/. .

Finally, one should note that all of the commands needed to build the Sphinx documentation, build the Doxygen API, and copy the files to the outer directory are encapsulated in the :file:`sphinx/configure.sh` shell script. Hence, everything that was presented in this section may be executed using the following single command:

.. code-block:: none

    cd sphinx
    ./configure.sh


.. note:: Since the Doxygen API likely does not need to be re-built each time you compile the Sphinx source code, it may be easiest to debug your Sphinx source code using the individual commands presented at the beginning of this section, rather than using :file:`configure.sh`.

-------------------------------------
Version Control for New Documentation
-------------------------------------

The :file:`configure.sh` script automatically adds your newly generated HTML files and image files to Git for tracking of the ``gh-pages`` branch using the following console commands:

.. code-block:: none

    cd sphinx
    git add .
    git add -f _images/*
    git add -f _images/math*
    git add -f doxygen/html/* .

You may commit your changes locally for the ``gh-pages`` branch just as you normally would expect to using Git:

.. code-block:: none

    git commit -m 'my commit message here'

Likewise, to push the changes made to your personal forked repository of OpenMOC, use the following command in the console:

.. code-block:: none

    git push

In order to submit new and/or revised documentation to the main OpenMOC reposotiry, please use a `pull request`_ to the ``gh-pages`` branch of the OpenMOC repository on GitHub as described in :ref:`Development Workflow <work_flow>`.





.. _Sphinx: http://sphinx-doc.org/
.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _Doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _MacPorts: http://www.macports.org/
.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _Rick Foos' blog post: http://rickfoosusa.blogspot.com/2011/10/howto-use-doxygen-with-github.html
.. _pull request: https://help.github.com/articles/using-pull-requests
