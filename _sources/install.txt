Installation
======================

System requirements
---------------------

Current Cogue works only on a private cluster that has the following
systems:

   1. Usual unix system, file system compatible of python ``os``
      module is required. 
   2. ``qsub`` and ``qstat`` of grid-engine queueing system has to
      work on the client node. Probably development for the other
      queueing system like ``torque`` is easy.
   3. VASP is used for the calculator.

The second condition can be critical for the use Cogue on
super-computers. Maybe in the future remote-job submitting should be
supported.

Install
------------------------

Firstly `spglib <http://spglib.sf.net>`_ has to be installed. Then
prepare ``settings.py`` on the same directory where ``setup.py``
exists. In ``settings.py``, complie optitions for spglib are written
as follows::

   library_dirs = [ Directory_of_libsymspg ]
   include_dirs = [ Directory_of_spglib_dot_h ]

Download Cogue using ``git`` command::

   git clone git://github.com/atztogo/cogue.git

 
Install by ``setup.py`` script at current directory::

   python setup.py install --home=.

Environment setup
-------------------

When using with VASP, the directory where the pseudopotentials are
stored has to be told to Cogue using the environment variable
``COGUE_POTCAR_PATH``. The details are given in the following sections.

.. toctree::
   :maxdepth: 2

   environment-setup
