#!/usr/bin/env python

import numpy as np
import qsushi
import qsushi.calculator.vasp as vasp
import qsushi.qsystem.gridengine as ge

task_name = "sno2"

# Crystal structure
symbols = ['Sn'] * 2 + ['O'] * 4
lattice = [[4.75, 0, 0],
           [0, 4.75, 0],
           [0, 0, 3.25]]
points=np.transpose([[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.5],
                     [0.3, 0.3, 0.0],
                     [0.7, 0.7, 0.0],
                     [0.2, 0.8, 0.5],
                     [0.8, 0.2, 0.5]])
cell = qsushi.cell(lattice=lattice,
                   points=points,
                   symbols=symbols)

# Vasp settings
ps_map = {'Sn': 'Sn_PBE',
          'O': 'O_PBE'}
incar = vasp.incar()
incar.set_structure_optimization()
incar.set_encut(400)
incar.set_prec("Normal")

# Queue
job = ge.job(script="vasp5212serial",
             shell="/bin/zsh",
             jobname=task_name,
             stdout="std.log",
             stderr="err.log")

# Task
task = vasp.bulk_modulus(max_iteration=2,
                         cell=cell,
                         pseudo_potential_map=ps_map,
                         k_mesh=[4, 4, 6],
                         incar=incar,
                         job=job)

# Automatic calculation
calc = qsushi.autocalc()
calc.append(task_name, task) # More tasks can be appended.
calc.set_queue(ge.queue())
calc.run(check_period=5)
print "space group:", qsushi.symmetry(cell)['international']
print "status:", task.get_status()
# 201.411956183 GPa
print "bulk modulus:", task.get_bulk_modulus(), "GPa"
