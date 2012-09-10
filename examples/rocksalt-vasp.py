#!/usr/bin/env python

import os
import numpy as np
import time
import qsushi
import qsushi.calculator.vasp as vasp
import qsushi.qsystem.gridengine as ge

a = 5.65
lattice = np.eye( 3 ) * a
points = np.array( [ [ 0.0, 0.0, 0.0 ],
                     [ 0.0, 0.5, 0.5 ],
                     [ 0.5, 0.0, 0.5 ],
                     [ 0.5, 0.5, 0.0 ],
                     [ 0.5, 0.5, 0.5 ],
                     [ 0.5, 0.0, 0.0 ],
                     [ 0.0, 0.5, 0.0 ],
                     [ 0.0, 0.0, 0.5 ] ], dtype=float ).T
symbols = ['Na'] * 4 + ['Cl'] * 4
cell = qsushi.cell( lattice = lattice,
                    points = points,
                    symbols = symbols )

print qsushi.symmetry( cell )['international']

ps_map = { 'Na': 'Na_pv_PBE',
           'Cl': 'Cl_PBE' }
task_name = "nacl"
incar = vasp.incar()

# Structure optimization
incar.set_structure_optimization()
task = vasp.structure_optimization( max_iteration = 3 )

# One point calculation
# incar.set_electronic_structure()
# task = vasp.electronic_structure()

incar.set_prec( "Normal" )
incar.set_encut( 400 )
incar.set_lreal( "Auto" )
task.set_configurations( cell = cell,
                         pseudo_potential_map = ps_map,
                         k_mesh = [ 4, 4, 4 ],
                         incar = incar )
job = ge.job( script = "vaspserial-togo",
              shell = "/bin/zsh",
              jobname = task_name,
              stdout = "std.log",
              stderr = "err.log" )
task.set_job( job )

# Use autocalc as a controller
calc = qsushi.autocalc()
calc.append( task_name, task ) # More tasks can be appended.
calc.set_queue( ge.queue() )
calc.run( check_period = 5 )
print "space group:", qsushi.symmetry( cell )['international']
print "status:", task.get_status()
