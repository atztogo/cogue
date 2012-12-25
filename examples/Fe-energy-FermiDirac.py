#!/usr/bin/env python

import numpy as np
import cogue
import cogue.calculator.vasp as vasp
import cogue.qsystem.gridengine as ge

task_name = "Fe"

# Crystal structure
symbols = ['Fe'] * 2
c = 2.8301794069183055
lattice = [[c, 0, 0],
           [0, c, 0],
           [0, 0, c]]
points=np.transpose([[0.0, 0.0, 0.0],
                     [0.5, 0.5, 0.5]])
cell = cogue.cell(lattice=lattice,
                  points=points,
                  symbols=symbols)

# Vasp settings
ps_map = {'Fe': 'Fe_PBE'}
incar = vasp.incar()
incar.set_electronic_structure()
incar.set_encut(300)
incar.set_ismear(-1)
incar.set_sigma(0.4)
incar.set_ispin(2)

# Queue
job = ge.job(script="mpirun vasp5212mpi",
             shell="/bin/zsh",
             jobname=task_name,
             pe="mpi* 4",
             stdout="std.log",
             stderr="err.log")

# Task
task = vasp.electronic_structure(cell=cell,
                                 pseudo_potential_map=ps_map,
                                 k_mesh=[8, 8, 8],
                                 incar=incar,
                                 job=job,
                                 traverse=False)

# Automatic calculation
calc = cogue.autocalc(name=task_name, verbose=True)
calc.append(task_name, task) # More tasks can be appended.
calc.set_queue(ge.queue())
calc.run(check_period=10)
