.. _usersguide_beginners:

=============================
A Beginner's Guide to OpenMOC
=============================

---------------------
What does OpenMOC do?
---------------------

OpenMOC is a neutron transport solver. It solves the neutron transport equation in a nuclear system, like
a fuel assembly or even a reactor core. This allows the reactor engineers to know where the power is
generated in the system, which is the starting point for solving other physics, like thermal hydraulics or
the fuel depletion.

-----------------
How does it work?
-----------------

In order to do anything, the code first needs to have a model of some problem of
interest. This could be a nuclear reactor or any other physical system with
fissioning material or a neutron source. You, as the code user, will need to
describe the model so that the code can do something with it. A basic model consists
of a few things:

- **Geometry** - A description of the geometry split into regions of homogeneous materials.
- **Materials** - A description of the nuclear cross-sections for each material
- **Parameters** - Various parameters for the numerical algorithm used in the simulation


-----------------------
What do I need to know?
-----------------------

If you are starting to work with OpenMOC, there are a few things you should be
familiar with. Whether you plan on working in Linux, Mac OS X, or Windows, you
should be comfortable working in a command line environment. There are many
resources online for learning command line environments. If you are using Linux
or Mac OS X (also Unix-derived), `this tutorial
<http://www.ee.surrey.ac.uk/Teaching/Unix/>`_ will help you get acquainted with
commonly-used commands. It is also helpful to be familiar with `Python
<http://www.python.org/>`_, as most of the post-processing utilities provided
with OpenMOC rely on it for data manipulation and results visualization.

OpenMOC uses a version control software called `git`_ to keep track of changes to
the code, document bugs and issues, and other development tasks. While you don't
necessarily have to have git installed in order to download and run OpenMOC, it
makes it much easier to receive updates if you do have it installed and have a
basic understanding of how it works. There are a list of good `git tutorials`_
at the git documentation website. The `OpenMOC source code`_ and documentation
are hosted at `GitHub`_. In order to receive updates to the code directly,
submit `bug reports`_, and perform other development tasks, you may want to sign
up for a free account on GitHub. Once you have an account, you can follow `these
instructions <http://help.github.com/set-up-git-redirect>`_ on how to set up
your computer for using GitHub.

If you are new to nuclear engineering, you may want to review the NRC's `Reactor
Concepts Manual`_. This manual describes the basics of nuclear power for
electricity generation, the fission process, and the overall systems in a
pressurized or boiling water reactor. Another resource that is a bit more
technical than the Reactor Concepts Manual but still at an elementary level is
the DOE Fundamentals Handbook on Nuclear Physics and Reactor Theory `Volume I`_
and `Volume II`_. You may also find it helpful to review the following terms:

- `Neutron cross-section`_
- `Effective multiplication factor`_
- `Neutron Flux`_

.. _git: http://git-scm.com/
.. _git tutorials: http://git-scm.com/documentation
.. _Reactor Concepts Manual: http://www.amazon.com/Electrical-Generation-Concepts-Technical-Training/dp/B009XXK564
.. _Volume I: http://energy.gov/sites/prod/files/2013/06/f2/h1019v1.pdf
.. _Volume II: http://energy.gov/sites/prod/files/2013/06/f2/h1019v2.pdf
.. _OpenMOC source code: https://github.com/mit-crpg/OpenMOC
.. _GitHub: https://github.com/
.. _bug reports: https://github.com/mit-crpg/OpenMOC/issues
.. _Neutron cross-section: http://en.wikipedia.org/wiki/Neutron_cross_section
.. _Effective multiplication factor: http://en.wikipedia.org/wiki/Effective_multiplication_factor
.. _Neutron Flux: http://en.wikipedia.org/wiki/Neutron_flux
