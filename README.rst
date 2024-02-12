  .. SPDX-FileCopyrightText: 2024 The spatial_cost_optimisation Author

  .. SPDX-License-Identifier: CC-BY-4.0

======
Spatial cost optimisation
======

This is a GIS tool to evaluate land availability and execute cost optimisation on network planning (e.g. grid networks, road network, etc.). It consists of an GIS analysis to evaluate land availability, as well as a cost optimisation on network searching and growth using Greedy (Prim's) Algorithm. 

The GIS analaysis is based on the foundation of `Atlite <https://github.com/PyPSA/atlite>`_, a free software, xarray-based Python library for converting weather data (like wind speeds, solar influx) into energy systems data. 

Installation
============

To install you need a working installation running Python 3.6 or above
and we strongly recommend using either miniconda or anaconda for package
management. Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is a mini version of `Anaconda <https://www.anaconda.com/>`_ that includes only ``conda`` and its dependencies or make sure ``conda`` is already installed on your system. For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_. Once you have installed conda we recommend you install Mamba to speed up the environment installation. Once installed open your terminal and run 

.. code:: shell

         conda install mamba -n base -c conda-forge

After completion of installing Mamba, please clone the repository into your local path 

.. code:: shell

         git clone https://github.com/RichChang963/spatial_cost_optimisation.git

The python package requirements are curated in the ``envs/environment.yaml`` file. The environment can be installed and activated using

.. code:: shell

         mamba env create -f envs/environment.yaml

The environment can be installed and activated using

.. code:: shell

         conda activate spatial_cost_opt


Documentation
===============

Please check the ``documentation``.


Authors and Copyright
======================

Copyright (C) 2024 The spatial_cost_optimisation Author.


Licence
=======

This work is licensed under multiple licences:

-  All original source code is licensed under `MIT <https://github.com/RichChang963/spatial_cost_optimisation/blob/develop/LICENSES/MIT.txt>`_
-  Auxiliary code from SPHINX is licensed under `BSD-2-Clause <https://github.com/RichChang963/spatial_cost_optimisation/blob/develop/LICENSES/BSD-2-Clause.txt>`_.
-  The documentation is licensed under `CC-BY-4.0 <https://github.com/RichChang963/spatial_cost_optimisation/blob/develop/LICENSES/CC-BY-4.0.txt>`_.
-  Configuration and data files are mostly licensed under `CC0-1.0 <https://github.com/RichChang963/spatial_cost_optimisation/blob/develop/LICENSES/CC0-1.0.txt>`_.