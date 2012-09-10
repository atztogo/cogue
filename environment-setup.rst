.. _environment_setups:

Environment setups
==================

Pseudopotentials for VASP
---------------------------

For automatic generation of ``POTCAR`` files, directory of VASP
pseudopotentails has to be informed to cogue by setting environment
variable of ``COGUE_POTCAR_PATH``. This is specified by setting the
following line in your shell configuration file, e.g., ``.bashrc`` or ``.zshenv``.

::

   COGUE_POTCAR_PATH=some_directory

All the ``POTCAR`` files for every atom and electronic configuration
have to be extracted to plane text files. the filenames of the
extracted files are used to notify pseudopotential files used in
cogue through python dictionary.

For example, ``POTCAR`` files of GGA-PBE may be prepared by running
the following command in ``potpaw_PBE`` directory using zsh::

   for i in `/bin/ls -d *(/)`;do zcat $i/POTCAR.Z > some_directory/"$i"_PBE;done

In this example, the PAW potential file of O named ``O_PBE`` is
created in ``some_directory``. Then a dictionary is used to pass this
filename to a VASP task using the keyword ``pseudo_potential_map`` like::

   pseudo_potential_map = { 'Si': 'Si_PBE', 'O': 'O_PBE' }

