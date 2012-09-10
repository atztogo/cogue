.. qsushi documentation master file, created by
   sphinx-quickstart on Tue Jan 10 14:32:09 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cogue
==================================

Cogue is a package of crystal simulation tools:

- Convenient tools to handle crystal structures
  - Crsytal format converter
- Wrappers of calculators and queueing systems for
  - `VASP5 <http://cms.mpi.univie.ac.at/vasp/vasp/>`_
  - `GridEngine <http://gridengine.org>`_
- Automation tools for various tasks

Tasks
------

- Total energy, eigenvalues, forces, stress
- Structure optimizations
- Elastic constants (finite difference approach in VASP)
- Phonon with `phonopy <http://phonopy.sf.net>`_ (finite difference approach)
- Bulk modulus (finite difference approach)

More tasks will be supported on demand.

Useful functions and classes
-----------------------------

- Crystal structure
- Calculator file parsers and writers
- Crystal symmetry handling with `spglib <http://spglib.sf.net>`_

Architecture of automation system
----------------------------------

Cogue is composed of three objects, task, taskset, and controller
(``autocalc``). Task corresponds to a calculation. Taskset is a set of
tasks and tasksets. A series of tasksets can be contained in a
task. Complicated process is designed by nesting tasks and tasksets.
Controller handles tasks and tasksets with the help of a queueing system.

Documentation
--------------

.. toctree::
   :maxdepth: 2

   contents



.. _example_rutile:

Example: calculation of bulk modulus
-------------------------------------

.. literalinclude:: ../examples/SnO2-bulkmodulus.py


About this code
------------------

The code is mainly written in Python.

- License: New BSD
- Contact: atz.togo@gmail.com
- Authour: Atsushi Togo

For the crystal structure comparison, Xtalcomp
(https://github.com/dlonie/XtalComp) by David Lonie is used.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


