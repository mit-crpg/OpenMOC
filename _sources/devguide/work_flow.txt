.. _work_flow:


====================
Development Workflow
====================

Anyone wishing to make contributions to OpenMOC should be fully acquianted and
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

.. _git: http://git-scm.com/
.. _GitHub: https://github.com/
.. _pull request: https://help.github.com/articles/using-pull-requests
.. _fork and pull development model: https://help.github.com/articles/using-pull-requests
.. _Fork: https://help.github.com/articles/fork-a-repo
.. _Create a git branch: http://git-scm.com/book/en/Git-Branching-Basic-Branching-and-Merging
.. _mit-crpg/OpenMOC: https://github.com/mit-crpg/OpenMOC
.. _paid plan: https://github.com/plans
