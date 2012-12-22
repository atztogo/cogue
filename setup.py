import os
from distutils.core import setup, Extension
# from numpy.distutils.misc_util import get_numpy_include_dirs
import numpy
include_dirs_numpy = [numpy.get_include()]

if os.path.exists("settings.py"):
    import settings
    
    try:
        print "Include directories", settings.include_dirs
        include_dirs = settings.include_dirs
    except NameError:
        include_dirs = []

    try:
        print "Library directories", settings.library_dirs
        library_dirs = settings.library_dirs
    except NameError:
        library_dirs = []

spglib = Extension('cogue._spglib',
                   include_dirs = ['ext/spglib'] + include_dirs_numpy,
                   # extra_compile_args=['-fopenmp'],
                   # extra_link_args=['-lgomp'],
                   sources = ['ext/_spglib.c',
                              'ext/spglib/refinement.c',
                              'ext/spglib/cell.c',
                              'ext/spglib/debug.c',
                              'ext/spglib/hall_symbol.c',
                              'ext/spglib/kpoint.c',
                              'ext/spglib/lattice.c',
                              'ext/spglib/mathfunc.c',
                              'ext/spglib/pointgroup.c',
                              'ext/spglib/primitive.c',
                              'ext/spglib/spacegroup.c',
                              'ext/spglib/spg_database.c',
                              'ext/spglib/spglib.c',
                              'ext/spglib/spin.c',
                              'ext/spglib/site_symmetry.c',
                              'ext/spglib/sitesym_database.c',
                              'ext/spglib/symmetry.c'] )

xtalcomp = Extension('cogue._xtalcomp',
                     libraries = ['stdc++'],
                     include_dirs = ['ext/xtalcomp'] + include_dirs_numpy,
                     sources = ['ext/_xtalcomp.c',
                                'ext/xtalcomp/xtalcomp_wrapper.cpp',
                                'ext/xtalcomp/xtalcomp.cpp',
                                'ext/xtalcomp/xctransform.cpp'])

setup(name = 'cogue',
      version = '0.1.0',
      description = 'Crystal simulation tools',
      author = 'Atsushi Togo',
      author_email = 'atz.togo@gmail.com',
      url = 'https://github.com/atztogo/cogue',
      packages = ['cogue',
                  'cogue.crystal',
                  'cogue.calculator',
                  'cogue.calculator.vasp',
                  'cogue.controller',
                  'cogue.electronic',
                  'cogue.interface',
                  'cogue.qsystem',
                  'cogue.task'],
      scripts = ['scripts/crystalView',
                 'scripts/poscar2poscar',
                 'scripts/poscar2cif',
                 'scripts/poscar2vsim',
                 'scripts/poscar2sketch',
                 'scripts/vasprun2poscars',
                 'scripts/symPoscar'],
      ext_modules = [spglib, xtalcomp])
