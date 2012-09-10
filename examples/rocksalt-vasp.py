#!/usr/bin/env python

import os
import numpy as np
import time
import cogue
import cogue.calculator.vasp as vasp
import cogue.qsystem.gridengine as ge

a = 5.65
lattice = np.eye(3) * a
points = np.transpose([[0.0, 0.0, 0.0],
                       [0.0, 0.5, 0.5],
                       [0.5, 0.0, 0.5],
                       [0.5, 0.5, 0.0],
                       [0.5, 0.5, 0.5],
                       [0.5, 0.0, 0.0],
                       [0.0, 0.5, 0.0],
                       [0.0, 0.0, 0.5]])
symbols = ['Na'] * 4 + ['Cl'] * 4
cell = cogue.cell(lattice=lattice,
                  points=points,
                  symbols=symbols)

print cogue.symmetry(cell)['international']

ps_map = {'Na': 'Na_pv_PBE',
          'Cl': 'Cl_PBE'}
task_name = "nacl"
incar = vasp.incar()

# Structure optimization
incar.set_structure_optimization()
task = vasp.structure_optimization(max_iteration=3)

# One point calculation
# incar.set_electronic_structure()
# task = vasp.electronic_structure()

incar.set_prec("Normal")
incar.set_encut(400)
incar.set_lreal("Auto")
task.set_configurations(cell=cell,
                        pseudo_potential_map=ps_map,
                        k_mesh=[4, 4, 4],
                        incar=incar)
job = ge.job(script="vasp5212serial",
             shell="/bin/zsh",
             jobname=task_name,
             stdout="std.log",
             stderr="err.log")
task.set_job(job)

# Use autocalc as a controller
calc = cogue.autocalc()
calc.append(task_name, task) # More tasks can be appended.
calc.set_queue(ge.queue())
calc.run(check_period = 5)
print "space group:", cogue.symmetry(cell)['international']
print "status:", task.get_status()
