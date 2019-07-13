.. _work_flow:


====================
Development Workflow
====================

Anyone wishing to make contributions to OpenMOC should be fully acquainted and
comfortable working with git_ and GitHub_. The primary means of modifying and
making contributions to OpenMOC is through a GitHub `pull request`_. This is
what's known as a `fork and pull development model`_. The steps for this are as
follows:

* **Fork the Repository** -  Fork_ the main OpenMOC repository from `mit-crpg/OpenMOC`_. This will create a repository with the same name under your personal account.

* **Create a Git Branch** - `Create a git branch`_ to your local repository and make commits that you intend to merge back into the main OpenMOC repository.

* **Create a Pull Request** - Create a `pull request`_ from GitHub and select the branch to merge from your local repository with the appropriate branch in `mit-crpg/OpenMOC`_.

* **Merge the Pull Request** - The OpenMOC integration manager will review your pull request and make sure it compiles and runs correctly and conforms to the :ref:`style_guide`.

* **Code Integration** - After the pull request has been thoroughly vetted, the integration manager will merge it back into `mit-crpg/OpenMOC`_.

While the process above depends on the fork of the OpenMOC repository being
publicly available on GitHub, you may also wish to do development on a private
repository for research or commercial purposes. The proper way to do this is to
create a complete copy of the OpenMOC repository (not a fork from GitHub). The
private repository can then either be stored just locally or in conjunction with
a private repository on Github (this requires a `paid plan`_). If you want to
merge some changes you've made in your private repository back to
`mit-crpg/OpenMOC`_ repository, simply follow the steps above with an extra step
of pulling a branch from your private repository into your public fork.

======================
Continuous integration
======================

In order to guarantee the correctness of OpenMOC, every pull request to OpenMOC is automatically tested
by `Travis CI`_. The test suite makes sure that each feature of the code keeps producing the same results / are
not affected by the new code, within acceptable tolerances. In order for new features that you may develop to
also provide this guarantee, any code you add needs to be covered by existing tests, or more often by a new test.
New tests can either be classified as unit or regression tests, both types may and often should be used to test
new functionalities. If your Github repository is public, you can also use Travis CI to test your personal fork.


`Unit tests`_ are meant to test a single function or a single attribute of an object. They are usually short
and run very fast. Because of their limited scope, it is fairly easy to find the cause of a failed unit test.


`Regression tests`_ usually involve running OpenMOC entirely, from track generation to sweeping, and then looking
at the eigenvalue and fluxes to see if they changed. Using the provided input sets, they are fairly easy to write,
but it is harder to find the root of the problem when they fail.


`coveralls`_, a code coverage tool gives information on how often each line of the code is run by the test
suite. However, developers should make sure that every functionality is tested, not that every line of code
is run.

.. _git: http://git-scm.com/
.. _GitHub: https://github.com/
.. _pull request: https://help.github.com/articles/using-pull-requests
.. _fork and pull development model: https://help.github.com/articles/using-pull-requests
.. _Fork: https://help.github.com/articles/fork-a-repo
.. _Create a git branch: http://git-scm.com/book/en/Git-Branching-Basic-Branching-and-Merging
.. _mit-crpg/OpenMOC: https://github.com/mit-crpg/OpenMOC
.. _paid plan: https://github.com/plans
.. _Travis CI: https://travis-ci.org/mit-crpg/OpenMOC/
.. _Unit tests: https://en.wikipedia.org/wiki/Unit_testing
.. _Regression tests: https://en.wikipedia.org/wiki/Regression_testing
.. _coveralls: https://coveralls.io
