.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where 
the path before the ``%`` sign denotes the directory in which the commands following the
 ``%`` should be entered.

Clone the Repository
====================

First of all, clone the `Atlite-agora repository <https://github.com/agoenergy/Atlite-agora>`_ 
using the version control system ``git``. The path to the directory into which the 
``git repository`` is cloned, must **not** have any spaces! If you do not have ``git`` 
installed, follow installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. code:: bash
    /some/other/path % cd /some/path/without/spaces
    /some/path/without/spaces % git clone https://github.com/agoenergy/Atlite-agora.git
.. _deps:

Install Python Dependencies
===============================

pypsa-agora relies on a set of other Python packages to function.
We recommend using the package manager and environment management system ``conda`` to 
install them. Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, 
which is a mini version of `Anaconda <https://www.anaconda.com/>`_ that includes only 
``conda`` and its dependencies or make sure ``conda`` is already installed on your system.
For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

Once you have installed conda we recommend you install Mamba to speed up the environment 
installation. Once installed open your terminal and run 

.. code:: bash

    .../pypsa-agora % conda install mamba -n base -c conda-forge

The python package requirements are curated in the `envs/environment.yaml <https://github.com/agoenergy/pypsa-agora/blob/main/env/environment.yaml>`_ file.
The environment can be installed and activated using

.. code:: bash

    .../pypsa-agora % mamba env create -f envs/environment.yaml

    .../pypsa-agora % conda activate pypsa-agora

Note that activation is local to the currently open shell!
After opening a new terminal window, one needs to reissue the second command!

.. warning::
    On Windows, new versions of ``ipopt`` have caused problems. Consider downgrading to version 3.11.1.

.. _defaultconfig:

Set Up the Default Configuration
================================

pypsa-agora has several configuration options that must be specified in a ``config.yaml`` file located in the root directory.
An example configuration ``config.default.yaml`` is maintained in the repository.
More details on the configuration options are in :ref:`config`.

Before first use, create a ``config.yaml`` by copying the example.

.. code:: bash

    .../pypsa-agora % cp config.default.yaml config.yaml

Users are advised to regularly check their own ``config.yaml`` against changes in the ``config.default.yaml``
when pulling a new version from the remote repository.
