.. qsushi documentation master file, created by
   sphinx-quickstart on Tue Jan 10 14:32:09 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cogue
==================================

Cogue is a package of crystal simulation tools:

- Automation tools for various :ref:`VASP tasks <VASP_tasks>` 
- Convenient tools to handle crystal structures, e.g., Crsytal format converter
- Wrappers of calculators and queueing systems for `VASP5
  <http://cms.mpi.univie.ac.at/vasp/vasp/>`_ and `GridEngine
  <http://gridengine.org>`_

Convenient command line tools
------------------------------

- Crystal format converters
- Transformation of crystal structure including supercell builder
- Symmetry finder
- Simple crystal viewer

Useful functions and classes
-----------------------------

- Crystal structure
- Calculator file parsers and writers
- Crystal symmetry handling with `spglib <http://spglib.sf.net>`_

Documentation
--------------

.. toctree::
   :maxdepth: 3

   contents



.. _example_rutile:

Example: calculation of bulk modulus
-------------------------------------

.. literalinclude:: ./examples/SnO2-bulkmodulus.py


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


