Installation
======================

.. contents::
   :depth: 2
   :local:

System requirements
----------------------------

General system requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Python and its header files (python-dev)
* numpy (python-numpy)
* hdf5 interface (python-h5py)
* python-yaml
* python-paramiko
* spur.py (https://github.com/mwilliamson/spur.py)
* phonopy (https://atztogo.github.io/phonopy/)
* phono3py (https://atztogo.github.io/phono3py/)

In the case of Ubuntu, the installation has be made as follows:

::

   % sudo apt-get install python-dev python-numpy python-yaml python-paramiko python-h5py python-pip
   % pip install spur --user
   * pip install phonopy --user
   * pip install phono3py --user

spur.py is an interface for ssh communication necessary to connect to
computer clusters. spur.py uses python-paramiko as its ssh kernel.

For the automation of phonon related calculations, phonopy
(https://atztogo.github.io/phonopy/) and phono3py
(https://atztogo.github.io/phono3py/) have to be installed on a local
machine where a cogue script is executed.


Automation system
^^^^^^^^^^^^^^^^^^

Cogue automation system controls computer clusters via queueing
system. It is supported that the following software are installed on
your computer cluster::

* Grid-engine queueing system
* VASP (http://vasp.at)

PBE-like queueing system is not yet prepared.

If you want to run everything on a local machine, see the following
contents for how to set-up the grid-engine queueing system:

.. toctree::
   :maxdepth: 1

   qsystem


Utility tools
^^^^^^^^^^^^^^

``crystalView`` requires mayavi (http://docs.enthought.com/mayavi/mayavi/).

::

   % sudo apt-get install mayavi2

Install
--------

Download Cogue from https://github.com/atztogo/cogue/tags .

Move into cogue directory and run ``setup.py`` script as follows::

   python setup.py install --home=.

PAW datasets for VASP
----------------------

When using with VASP, location of pseudopotentials has to be told to
Cogue using the environment variable ``COGUE_POTCAR_PATH``. The
details are given in the following sections.

For automatic generation of ``POTCAR`` files, directory of VASP
pseudopotentails has to be informed to cogue by setting environment
variable of ``COGUE_POTCAR_PATH``. This is specified by setting the
following line in your shell configuration file, e.g., ``.bashrc`` or ``.zshenv``.

::

   COGUE_POTCAR_PATH=$DIR_PREFIX/cogue_potcar_dir

where $DIR_PREFIX is the directory prefix, for instance home
directory, ``/home/your_unix_name``.

All the ``POTCAR`` files for every atom and electronic configuration
have to be extracted to plane text files. the filenames of the
extracted files are used to notify pseudopotential files used in
cogue through python dictionary.

For example, ``POTCAR`` files of GGA-PBE may be prepared by running
the following command in ``potpaw_PBE`` directory using zsh::

   for i in `/bin/ls -d *(/)`;do zcat $i/POTCAR.Z > cogue_potcar_dir/"$i_PBE";done

In this example, the PAW potential files are created in
``cogue_potcar_dir``. Then a dictionary is used to pass this filename
to a VASP task using the keyword ``pseudo_potential_map`` like::

   pseudo_potential_map = { 'Si': 'Si_PBE', 'O': 'O_PBE' }
