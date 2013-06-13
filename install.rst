Installation
======================

System requirements
---------------------

* Python and its header files
* numpy
* python-lxml
* python-yaml
* spur.py

In the case of Ubuntu, the installation has be made as follows:

::

   % sudo apt-get install python-dev python-numpy python-lxml python-yaml python-paramiko python-pip
   % sudo pip install spur

spur.py is an interface for ssh. This is built on paramiko. So python-paramiko has to be installed. Ubuntu doesn't have spur.py package. Pip gives the easy way to install spur.py.

Tools
^^^^^^

crystalView
~~~~~~~~~~~~~

* mayavi2

::

   % sudo apt-get install mayavi2

Automation system
^^^^^^^^^^^^^^^^^

* Grid Engine (currently essential)
* VASP (for most first-principles calclulations)
* Phonopy (for phonon calculations)

Currently cogue automation system works only on a private cluster that
has the following systems:

   1. ``qsub`` and ``qstat`` of grid-engine queueing system has to
      work on the client node. Probably development for the other
      queueing system like ``torque`` is easy.
   2. VASP is used as the calculator.

For phonon calculation, ``phonopy`` is required to be installed and
correctly set-up. Please see `phonopy document
<http://phonopy.sf.net>`_.

Install
------------------------

Download Cogue using ``git`` command::

   git clone git://github.com/atztogo/cogue.git

Change directory to ``cogue`` and run ``setup.py`` script at current
directory::

   python setup.py install --home=.

Environment setup
-------------------

When using with VASP, the directory where the pseudopotentials are
stored has to be told to Cogue using the environment variable
``COGUE_POTCAR_PATH``. The details are given in the following sections.

.. toctree::
   :maxdepth: 1

   environment-setup
