#!/usr/bin/env python

import numpy as np

import cogue
import cogue.calculator.vasp as vasp
import cogue.qsystem.gridengine as ge

ps_map = {"Si": "Si_PBE"}

Si_str = """Si
1.
     5.4661639157319968    0.0000000000000000    0.0000000000000000
     0.0000000000000000    5.4661639157319968    0.0000000000000000
     0.0000000000000000    0.0000000000000000    5.4661639157319968
   Si
   8
Direct
  0.8750000000000000  0.8750000000000000  0.8750000000000000
  0.8750000000000000  0.3750000000000000  0.3750000000000000
  0.3750000000000000  0.8750000000000000  0.3750000000000000
  0.3750000000000000  0.3750000000000000  0.8750000000000000
  0.1250000000000000  0.1250000000000000  0.1250000000000000
  0.1250000000000000  0.6250000000000000  0.6250000000000000
  0.6250000000000000  0.1250000000000000  0.6250000000000000
  0.6250000000000000  0.6250000000000000  0.1250000000000000"""


def get_task(task_name):
    cell = vasp.parse_poscar(Si_str.split("\n"))

    # Vasp settings
    incar = vasp.incar()
    incar.set_structure_optimization()
    incar.set_nsw(20)
    incar.set_encut(300)

    incar_ph_rx = vasp.incar()
    incar_ph_rx.set_structure_optimization()
    incar_ph_rx.set_nsw(20)
    incar_ph_rx.set_encut(300)
    incar_ph_rx.set_isif(4)  # volume constant

    incar_ph_dsp = vasp.incar()
    incar_ph_dsp.set_electronic_structure()
    incar_ph_dsp.set_encut(300)

    # Queue
    job = ge.job(
        script="mpirun vasp5212mpi",
        shell="/bin/zsh",
        jobname=task_name,
        pe="mpi* 8",
        stdout="std.log",
        stderr="err.log",
    )

    # Task
    task = vasp.quasiharmonic_phonon(
        max_iteration=10,
        min_iteration=1,
        sampling_mesh=[20, 20, 20],
        supercell_matrix=np.diag([2, 2, 2]),
        cell=cell,
        pseudo_potential_map=ps_map,
        incar=[incar, incar_ph_rx, incar_ph_dsp],
        k_mesh=[[8, 8, 8], [8, 8, 8], [4, 4, 4]],
        job=[job, job.copy(task_name + "-ph_relax"), job.copy(task_name + "-ph")],
        traverse=False,
    )

    return task


# Automatic calculation
calc = cogue.autocalc(name=None, verbose=True)
task_name = "Si-QHA"
calc.append(task_name, get_task(task_name))  # More tasks can be appended.
calc.set_queue(ge.queue())
calc.run(check_period=10)
