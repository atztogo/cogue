.. qsushi documentation master file, created by
   sphinx-quickstart on Tue Jan 10 14:32:09 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Cogue
==================================

Cogue is a package of crystal simulation tools constituting of
interfaces of calculator (VASP) and queueing system (GridEngine).

Documentation
--------------

.. toctree::
   :maxdepth: 2

   contents

- `Presentaion slide <https://github.com/atztogo/cogue/blob/gh-pages/cogue-presentation.pdf?raw=true>`_

.. _example_rutile:

Examples
---------

Bulk modulus task on a local machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ./examples/SnO2-bulkmodulus.py

Phonon task on a remote cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ./examples/ZnO-phonon-remote.py


About this code
------------------

The code is mainly written in Python.

- License: New BSD
- Contact: atz.togo@gmail.com
- Authour: Atsushi Togo

For the crystal structure comparison, Xtalcomp
(https://github.com/dlonie/XtalComp) by David Lonie is used.

