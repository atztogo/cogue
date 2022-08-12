#!/usr/bin/env python

import numpy as np
import spur

import cogue
import cogue.calculator.vasp as vasp
import cogue.qsystem.gridengine as ge

poscar = """Zn O
   1.0
     3.2882532570702230    0.0000000000000000    0.0000000000000000
    -1.6441266285351115    2.8477108546997352    0.0000000000000000
     0.0000000000000000    0.0000000000000000    5.3061089751235242
   2   2
Direct
  0.3333333333333333  0.6666666666666667  0.9996795200513500
  0.6666666666666666  0.3333333333333333  0.4996795200513500
  0.3333333333333333  0.6666666666666667  0.3787634652515464
  0.6666666666666666  0.3333333333333333  0.8787634652515464
"""

task_name = "ZnO"

# Vasp settings
cell = vasp.parse_poscar(poscar.split("\n"))
ps_map = {"Zn": "Zn_PBE", "O": "O_PBE"}
incar = vasp.incar()
incar.set_structure_optimization()
incar.set_nsw(20)

incar_phonon = vasp.incar()
incar_phonon.set_electronic_structure()

# Grid engine job
job = ge.job(
    script="mpirun vasp5212mpi",
    shell="/bin/zsh",
    jobname=task_name,
    pe="mpi* 4",
    stdout="std.log",
    stderr="err.log",
)

# VASP phonon task
task = vasp.phonon(
    max_iteration=10,
    min_iteration=1,
    supercell_matrix=np.diag([2, 2, 2]),
    cell=cell,
    pseudo_potential_map=ps_map,
    incar=[incar, incar_phonon],
    k_mesh=[[6, 6, 4], [3, 3, 2]],
    k_shift=[[0.5, 0.5, 0], [0, 0, 0]],
    job=job,
)

# Automation system
calc = cogue.autocalc(name=task_name, verbose=True)

# Register task(s)
calc.append(task_name, task)  # More tasks can be appended.

# Register queue
shell = spur.SshShell(
    hostname="remotehost", missing_host_key=spur.ssh.MissingHostKey.accept
)
queue = ge.queue(ssh_shell=shell, temporary_dir="/home/bob/coguetmp")
calc.set_queue(queue)

# Run automation
calc.run(check_period=10)
