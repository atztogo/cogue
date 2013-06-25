Installation
======================

System requirements
----------------------------

General system requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Python and its header files (python-dev)
* numpy (python-numpy)
* hdf5 interface (python-h5py)
* python-lxml
* python-yaml
* python-paramiko
* spur.py (https://github.com/mwilliamson/spur.py)

In the case of Ubuntu, the installation has be made as follows:

::

   % sudo apt-get install python-dev python-numpy python-lxml python-yaml python-paramiko python-h5py python-pip 
   % sudo pip install spur

spur.py is an interface for ssh. python-paramiko is used by spur.py.
Since Ubuntu doesn't contain spur.py deb-package, using pip is
recommended to install spur.py.

Automation system
^^^^^^^^^^^^^^^^^^^^^^

Cogue automation system controls computer clusters via queueing
system. It is supported that the following software are installed on
your computer cluster::

* Grid-engine queueing system 
* VASP (http://vasp.at)

PBE-like queueing system is not yet prepared.

For the automation of phonon related calculation, phonopy
(http://phonopy.sf.net) has to be installed on a local machine where
cogue is executed.

If you want to run everything on a local machine, see the following
contents for how to set-up the grid-engine queueing system:

.. toctree::
   :maxdepth: 1

   qsystem


Utility tools
^^^^^^^^^^^^^^^^^^^^

``crystalView`` requires mayavi (http://docs.enthought.com/mayavi/mayavi/).

::

   % sudo apt-get install mayavi2

Install
------------------------

Download Cogue from https://github.com/atztogo/cogue/tags .

Move into cogue directory and run ``setup.py`` script as follows::

   python setup.py install --home=.

Environment setup
-------------------

When using with VASP, location of pseudopotentials has to be told to
Cogue using the environment variable ``COGUE_POTCAR_PATH``. The
details are given in the following sections.

.. toctree::
   :maxdepth: 1

   environment-setup
