.. _documentation:

=============
Documentation
=============

The following describes the suggested steps that developers should take to edit, revise or add to OpenMOC's documentation. When new features are added to OpenMOC, the documentation located in :file:`/docs` should be updated in the same pull request as the new feature. OpenMOC uses the Sphinx_ tool for documentation due to its growing popularity for open source codes. Sphinx uses reStructuredText_ (rst) files to encapuslate text-based content with semantic markup for figures, tables, URLs, and more. Sphinx provides a compiler which parses the reStructuredText source files to generate a target output in the form of HTML or a PDF. OpenMOC also uses Doxygen_ for automated generated of Application Programming Interface (API) documentation based on source code comments. The use of Doxygen-style comments is explained further in :ref:`Code Comments <code_comments>`.


-------------------------
Preliminary Configuration
-------------------------

First, in order to build OpenMOC's documentation, you must have installed both Sphinx and Doxygen. Since both tools are widely used, they are available with most source package system. On Debian-based systems, such as Ubuntu, you may install both tools with the following command in the console:

.. code-block:: none

    sudo apt-get install python-sphinx doxygen

Likewise, on Mac OS with MacPorts_ you may install each with the following command:

.. code-block:: none

    sudo port install py27-sphinx doxygen

These packages are all that is required to build and update the documentation in your pull request. Next, we will discuss how you can build and update the documentation when you modify the source code.

-------------------------
Writing New Documentation
-------------------------

This tutorial will not present any of the semantics of how to write documentation for Sphinx as that is well covered by its own tutorials. Likewise, the same is true for Doxygen-style comments, although the basics are presented in :ref:`Code Comments <code_comments>`. Instead, this section presents the location of the source code and the basic console commands used to generate HTML documentation.

First, the :file:`/docs/source` directory contains the reStructuredText source for the Sphinx documentation. In order to compile the Sphinx documentation into HTML, you can run the following from within the :file:`/docs` directory:

.. code-block:: none

    cd docs
    make html

This command will use the Sphinx tool to generate HTML in the :file:`/docs/build` directory. From within the :file:`/docs` directory, you may use a web browser, such as Mozilla Firefox, to view your newly generated documentation as follows:

.. code-block:: none

    firefox build/html/index.html

The reStructuredText source files that should be updated with each pull request are located in the :file:`/docs/source` directory. If developers wish to add new pictures to the documentation, image files should be added to the :file:`/docs/img` directory in either png or svg format. Developers can test out their changes to the documentation by recompiling the Sphinx documentation into HTML and viewing the updated HTML files.

Likewise, in order to build a revised version of the :ref:`OpenMOC API <api>` from Doxygen-style comments in the Python and C++ source code, you may use the following console commands from within the :file:`/docs` directory:

.. code-block:: none

    cd doxygen
    doxygen Doxyfile
    ./doxy2swig.py

This will build a new version of the API in the :file:`/docs/doxygen/html` directory. In addition, the :file:`docs/doxy2swig.py` script will parse the Doxygen-style comments and generate Python docstring-equivalents in an updated :file:`openmoc/docstring.i` file for SWIG to inject into the :code:`openmoc` Python module. From within the :file:`/docs/doxygen` directory, you may use a web browser, such as Mozilla Firefox, to view your newly generated documentation as follows:

.. code-block:: none

    firefox html/index.html

Note that this procedure updates the documentation shipped with the source code in the **develop** branch, but changes will not be reflected in the documentation on the website. GitHub_ reserves the use of the **gh-pages** branch to host a project website for each repository. The **gh-pages** branch contains the documentation source code for the latest public release of OpenMOC and thus will only be periodically updated with the documentation from the **develop** branch. The procedure for updating the website documentation will be discussed in the next section.

-------------------------------------
Version Control for New Documentation
-------------------------------------

With each new release, the website documentation will be updated by one of the core developers with the latest documentation from the **develop** branch. In order to **expose** newly generated HTML documentation to GitHub such that it will present on the website, the Sphinx and Doxygen build files must be copied to the **gh-pages** branch. Once the documentation is ready to be copied over to the **gh-pages** branch, the Sphinx and Doxygen HTML files need to be built:

.. code-block:: none

    cd docs
    make html
    cd doxygen
    doxygen Doxyfile
    cd ../..

Now that the HTML files have been built, you can checkout the **gh-pages** branch and then checkout a new branch to store the new version of the documentation:

.. code-block:: none

    git checkout gh-pages
    git checkout -b gh-pages-update

You can then copy the files from the :file:`/docs` directory to the root OpenMOC directory on your new branch:

.. code-block:: none

    cp -r docs/build/html/* .
    cp -r docs/doxygen/html doxygen/

Since the :file:`.gitignore` file does not allow image and HTML files to be tracked, you will need to force all image and HTML files to be added to tracking history. The :file:`configure.sh` script contains the commands to add all the necessary files. The :file:`configure.sh` can be executed and the changes committed and pushed to GitHub:

.. code-block:: none

    ./configure.sh
    git commit -am 'my commit message here'
    git push origin gh-pages-update

The final step in updating the website documentation is to `pull request`_ to the ``gh-pages`` branch of the OpenMOC repository on GitHub as described in :ref:`Development Workflow <work_flow>`.

.. note:: The diff for the pull request into **gh-pages** will likely be massive! The purpose of this pull request is not to conduct a code review, but rather to provide a transparent and publicly-available record of when the website documentation has been changed.

.. _Sphinx: http://sphinx-doc.org/
.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _Doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _MacPorts: http://www.macports.org/
.. _GitHub: https://github.com/mit-crpg/OpenMOC
.. _pull request: https://help.github.com/articles/using-pull-requests
