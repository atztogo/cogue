Example
========

``rutile.py`` in the ``example`` directory can be an introduction of
Cogue. Please see :ref:`example_rutile`. On the top of this example,
three lines of module imports are written::

   import cogue
   import cogue.calculator.vasp as vasp
   import cogue.qsystem.gridengine as ge

To import ``cogue``, a few method becomes ready to use as shownin
:py:mod:`cogue`. The methods :py:meth:`autocalc <cogue.autocalc>` and
:py:meth:`cell <cogue.cell>` return objects of
:py:class:`AutoCalc <cogue.controller.autocalc.AutoCalc>` and
:py:class:`Cell <cogue.crystal.cell.Cell>` classes,
respectively. These are designed to work for general calculators and queueing system. :py:mod:`vasp <cogue.calculator.vasp>` and 
:py:mod:`gridengine <cogue.qsystem.gridengine>` are necessary to run
the calculation on the specific system environment.

