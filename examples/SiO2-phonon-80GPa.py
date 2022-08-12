#!/usr/bin/env python

import numpy as np

import cogue
import cogue.calculator.vasp as vasp
import cogue.qsystem.gridengine as ge

task_name = "SiO2-80GPa"

# Crystal structure
symbols = ["Si"] * 2 + ["O"] * 4
lattice = [[4.65, 0, 0], [0, 4.75, 0], [0, 0, 3.25]]  # Orthorhombic
points = np.transpose(
    [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
)
cell = cogue.cell(lattice=lattice, points=points, symbols=symbols)

# Vasp settings
ps_map = {"Si": "Si_PBE", "O": "O_PBE"}
incar = vasp.incar()
incar.set_structure_optimization()
incar.set_nsw(20)
incar.set_pstress(800)

incar_phonon = vasp.incar()
incar_phonon.set_electronic_structure()

# Queue
job = ge.job(
    script="mpirun vasp5212mpi",
    shell="/bin/zsh",
    jobname=task_name,
    pe="mpi* 4",
    stdout="std.log",
    stderr="err.log",
)

# Task
task = vasp.phonon(
    max_iteration=10,
    min_iteration=1,
    supercell_matrix=np.diag([2, 2, 2]),
    cell=cell,
    pressure_target=80,
    pseudo_potential_map=ps_map,
    incar=[incar, incar_phonon],
    k_mesh=[[4, 4, 6], [2, 2, 3]],
    k_shift=[[0, 0, 0], [0, 0, 0.5]],
    job=job,
)

# Automatic calculation
calc = cogue.autocalc(name=task_name, verbose=True)
calc.append(task_name, task)  # More tasks can be appended.
calc.set_queue(ge.queue())
calc.run(check_period=10)
