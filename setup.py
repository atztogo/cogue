import os
from distutils.core import Extension, setup

# from numpy.distutils.misc_util import get_numpy_include_dirs
import numpy

include_dirs_numpy = [numpy.get_include()]

if os.path.exists("settings.py"):
    import settings

    try:
        print("Include directories %s" % settings.include_dirs)
        include_dirs = settings.include_dirs
    except NameError:
        include_dirs = []

    try:
        print("Library directories %s" % settings.library_dirs)
        library_dirs = settings.library_dirs
    except NameError:
        library_dirs = []

xtalcomp = Extension(
    "cogue._xtalcomp",
    libraries=["stdc++"],
    include_dirs=["ext/xtalcomp"] + include_dirs_numpy,
    sources=[
        "ext/_xtalcomp.c",
        "ext/xtalcomp/xtalcomp_wrapper.cpp",
        "ext/xtalcomp/xtalcomp.cpp",
        "ext/xtalcomp/xctransform.cpp",
    ],
)

setup(
    name="cogue",
    version="0.1.0",
    description="Crystal simulation tools",
    author="Atsushi Togo",
    author_email="atz.togo@gmail.com",
    url="https://github.com/atztogo/cogue",
    packages=[
        "cogue",
        "cogue.crystal",
        "cogue.calculator",
        "cogue.controller",
        "cogue.electron",
        "cogue.interface",
        "cogue.phonon",
        "cogue.qsystem",
        "cogue.task",
        "cogue.visualizer",
    ],
    scripts=[
        "scripts/crystalView",
        "scripts/poscar2poscar",
        "scripts/poscar2cif",
        "scripts/poscar2vsim",
        "scripts/poscar2sketch",
        "scripts/vasprun2poscars",
        "scripts/symPoscar",
    ],
    ext_modules=[
        xtalcomp,
    ],
)
