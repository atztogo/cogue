Command line tools
====================

Symmetry finder
-----------------

``symPoscar``
^^^^^^^^^^^^^

Converter
----------

The following options are shared with these converters.

* Shift atomic positions (``--shift`` option)
* Transform crystal structure to a cell with shorter periodicity
  (``--tmat`` option).
* Transform primitive Rhombohedral cell to hexagonal Rhombohedral cell
  (``--r2h`` option)
* Transform crystal structure to, i.e., supercell (``--dim`` option).

The order that these options work is ``--shift``, ``--tmat``, ``--r2h``, ``--dim``.

::

   poscar2poscar -h
   Usage: poscar2poscar [options]
   
   Options:
     -h, --help          show this help message and exit
     --r2h               Transform primitive Rhombohedral to hexagonal
                         Rhombohedral. This has to be used exclusively to the
                         other options.
     --tmat=T_MAT        Multiply transformation matrix. Absolute value of
                         determinant has to be 1 or less than 1.
     --dim=S_MAT         Supercell matrix
     --shift=SHIFT       Origin shift
     -o OUTPUT_FILENAME  Output filename
     -v                  More information is output.



``poscar2poscar``
^^^^^^^^^^^^^^^^^^^

VASP POSCAR format is converted to POSCAR format. This is useful when
we want to transform cell.

``poscar2cif``
^^^^^^^^^^^^^^^^^^^

VASP POSCAR format is converted to CIF format with P1 symmetry.

``poscar2vsim``
^^^^^^^^^^^^^^^^^^^

VASP POSCAR format is converted to `v_sim <http://www-drfmc.cea.fr/sp2m/L_Sim/V_Sim/index.en.html>`_  format.

Crystal viewer
---------------

``crystalView``
^^^^^^^^^^^^^^^^^^^

Crystal structure is visualized. Mayavi2 is necessary to use this script.
