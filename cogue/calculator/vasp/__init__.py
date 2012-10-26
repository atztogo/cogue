""" """

import os
import sys
import numpy as np
import cogue.crystal.vasp_io as vasp_io
from cogue.calculator.vasp.task import *

def incar(addgrid=None,
          ediff=None,
          ediffg=None,
          encut=None,
          ialgo=None,
          ibrion=None,
          isif=None,
          ismear=None,
          ispin=None,
          lcharg=None,
          lreal=None,
          lwave=None,
          nelmin=None,
          nsw=None,
          prec=None,
          pstress=None,
          sigma=None):
    """Returns Incar object"""

    return vasp_io.Incar(addgrid=addgrid,
                         ediff=ediff,
                         ediffg=ediffg,
                         encut=encut,
                         ialgo=ialgo,
                         ibrion=ibrion,
                         isif=isif,
                         ismear=ismear,
                         ispin=ispin,
                         lcharg=lcharg,
                         lreal=lreal,
                         lwave=lwave,
                         nelmin=nelmin,
                         nsw=nsw,
                         prec=prec,
                         pstress=pstress,
                         sigma=sigma)

def read_poscar(filename="POSCAR"):
    """Returns Cell object by reading POSCAR style file."""
    return vasp_io.read_poscar(filename)

def parse_poscar(lines):
    """Returns Cell object by parsing POSCAR lines list."""
    return vasp_io.parse_poscar(lines)

def write_poscar(cell, filename=None):
    """Write Cell object into POSCAR style file."""
    vasp_io.write_poscar(cell, filename)

def write_potcar(names, filename="POTCAR"):
    """Write POTCAR from filenames in the directory specified by cogue_POTCAR_PATH. """
    vasp_io.write_potcar(names, filename)


def electronic_structure(directory="electronic_structure",
                         name=None,
                         job=None,
                         traverse=False,
                         cell=None,
                         pseudo_potential_map=None,
                         k_mesh=None,
                         k_shift=None,
                         k_gamma=None,
                         k_length=None,
                         incar=None):

    es = ElectronicStructure(directory=directory,
                             name=name,
                             traverse=traverse)

    es.set_configurations(cell=cell,
                          pseudo_potential_map=pseudo_potential_map,
                          k_mesh=k_mesh,
                          k_shift=k_shift,
                          k_gamma=k_gamma,
                          k_length=k_length,
                          incar=incar)
    es.set_job(job)

    return es

def structure_optimization(directory="structure_optimization",
                           name=None,
                           job=None,
                           lattice_tolerance=0.1,
                           force_tolerance=1e-3,
                           pressure_target=0,
                           stress_tolerance=0.1, # GPa=kbar / 10
                           max_increase=1.5,
                           max_iteration=5,
                           min_iteration=1,
                           find_symmetry=True,
                           impose_symmetry=False,
                           symmetry_tolerance=0.1,
                           traverse=False,
                           cell=None,
                           pseudo_potential_map=None,
                           k_mesh=None,
                           k_shift=None,
                           k_gamma=None,
                           k_length=None,
                           incar=None):

    so = StructureOptimization(directory=directory,
                               name=name,
                               lattice_tolerance=lattice_tolerance,
                               force_tolerance=force_tolerance,
                               pressure_target=pressure_target,
                               stress_tolerance=stress_tolerance,
                               max_increase=max_increase,
                               max_iteration=max_iteration,
                               min_iteration=min_iteration,
                               find_symmetry=find_symmetry,
                               impose_symmetry=impose_symmetry,
                               symmetry_tolerance=symmetry_tolerance,
                               traverse=traverse)

    so.set_configurations(cell=cell,
                          pseudo_potential_map=pseudo_potential_map,
                          k_mesh=k_mesh,
                          k_shift=k_shift,
                          k_gamma=k_gamma,
                          k_length=k_length,
                          incar=incar)
    so.set_job(job)

    return so

def bulk_modulus(directory="bulk_modulus",
                 name=None,
                 job=None,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-3,
                 pressure_target=0,
                 stress_tolerance=10,
                 max_increase=1.5,
                 max_iteration=4,
                 min_iteration=1,
                 traverse=False,
                 is_cell_relaxed=False,
                 cell=None,
                 pseudo_potential_map=None,
                 k_mesh=None,
                 k_shift=None,
                 k_gamma=None,
                 k_length=None,
                 incar=None):

    bk = BulkModulus(directory=directory,
                     name=name,
                     lattice_tolerance=lattice_tolerance,
                     force_tolerance=force_tolerance,
                     pressure_target=pressure_target,
                     stress_tolerance=stress_tolerance,
                     max_increase=max_increase,
                     max_iteration=max_iteration,
                     min_iteration=min_iteration,
                     traverse=traverse,
                     is_cell_relaxed=is_cell_relaxed)

    bk.set_configurations(cell=cell,
                          pseudo_potential_map=pseudo_potential_map,
                          k_mesh=k_mesh,
                          k_shift=k_shift,
                          k_gamma=k_gamma,
                          k_length=k_length,
                          incar=incar)
    bk.set_job(job)

    return bk


def phonon(directory="phonon",
           name=None,
           job=None,
           supercell_matrix=np.eye(3, dtype=int),
           primitive_matrix=np.eye(3, dtype=int),
           distance=0.01,
           lattice_tolerance=0.1,
           force_tolerance=1e-3,
           pressure_target=0,
           stress_tolerance=10,
           max_increase=1.5,
           max_iteration=4,
           min_iteration=1,
           traverse=False,
           is_cell_relaxed=False,
           cell=None,
           pseudo_potential_map=None,
           k_mesh=None,
           k_shift=None,
           k_gamma=None,
           k_length=None,
           incar=None):

    ph = Phonon(directory=directory,
                name=name,
                supercell_matrix=supercell_matrix,
                primitive_matrix=primitive_matrix,
                distance=distance,
                lattice_tolerance=lattice_tolerance,
                force_tolerance=force_tolerance,
                pressure_target=pressure_target,
                stress_tolerance=stress_tolerance,
                max_increase=max_increase,
                max_iteration=max_iteration,
                min_iteration=min_iteration,
                traverse=traverse,
                is_cell_relaxed=is_cell_relaxed)

    ph.set_configurations(cell=cell,
                          pseudo_potential_map=pseudo_potential_map,
                          k_mesh=k_mesh,
                          k_shift=k_shift,
                          k_gamma=k_gamma,
                          k_length=k_length,
                          incar=incar)
    ph.set_job(job)

    return ph


def elastic_constants(directory="elastic_constants",
                      name=None,
                      job=None,
                      lattice_tolerance=0.1,
                      force_tolerance=1e-3,
                      pressure_target=0,
                      stress_tolerance=10,
                      max_increase=1.5,
                      max_iteration=4,
                      min_iteration=1,
                      traverse=False,
                      is_cell_relaxed=False,
                      cell=None,
                      pseudo_potential_map=None,
                      k_mesh=None,
                      k_shift=None,
                      k_gamma=None,
                      k_length=None,
                      incar=None):

    ec = ElasticConstants(directory=directory,
                          name=name,
                          lattice_tolerance=lattice_tolerance,
                          force_tolerance=force_tolerance,
                          pressure_target=pressure_target,
                          stress_tolerance=stress_tolerance,
                          max_increase=max_increase,
                          max_iteration=max_iteration,
                          min_iteration=min_iteration,
                          traverse=traverse,
                          is_cell_relaxed=is_cell_relaxed)

    ec.set_configurations(cell=cell,
                          pseudo_potential_map=pseudo_potential_map,
                          k_mesh=k_mesh,
                          k_shift=k_shift,
                          k_gamma=k_gamma,
                          k_length=k_length,
                          incar=incar)
    ec.set_job(job)

    return ec


def phonon_relax_element(directory="phonon_relax_element",
                         name=None,
                         job=None,
                         distance=0.01,
                         lattice_tolerance=0.1,
                         force_tolerance=1e-3,
                         pressure_target=0,
                         stress_tolerance=10,
                         max_increase=1.5,
                         max_iteration=4,
                         min_iteration=1,
                         symmetry_tolerance=0.1,
                         cutoff_eigenvalue=-0.02,
                         max_displacement=None,
                         num_sampling_points=60,
                         traverse=False,
                         cell=None,
                         pseudo_potential_map=None,
                         k_mesh=None,
                         k_shift=None,
                         k_gamma=None,
                         k_length=None,
                         incar=None):

    phre = PhononRelaxElement(directory=directory,
                              name=name,
                              distance=distance,
                              lattice_tolerance=lattice_tolerance,
                              force_tolerance=force_tolerance,
                              pressure_target=pressure_target,
                              stress_tolerance=stress_tolerance,
                              max_increase=max_increase,
                              max_iteration=max_iteration,
                              min_iteration=min_iteration,
                              symmetry_tolerance=symmetry_tolerance,
                              cutoff_eigenvalue=cutoff_eigenvalue,
                              max_displacement=max_displacement,
                              num_sampling_points=num_sampling_points,
                              traverse=traverse)

    phre.set_configurations(cell=cell,
                            pseudo_potential_map=pseudo_potential_map,
                            k_mesh=k_mesh,
                            k_shift=k_shift,
                            k_gamma=k_gamma,
                            k_length=k_length,
                            incar=incar)
    phre.set_job(job)

    return phre

def phonon_relax(directory="phonon_relax",
                 name=None,
                 job=None,
                 distance=0.01,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-3,
                 pressure_target=0,
                 stress_tolerance=10,
                 max_increase=1.5,
                 max_iteration=4,
                 min_iteration=1,
                 symmetry_tolerance=0.1,
                 restrict_offspring=False,
                 max_offspring=None,
                 cutoff_eigenvalue=-0.02,
                 max_displacement=None,
                 num_sampling_points=60,
                 traverse=False,
                 cell=None,
                 pseudo_potential_map=None,
                 k_mesh=None,
                 k_shift=None,
                 k_gamma=None,
                 k_length=None,
                 incar=None):

    phr = PhononRelax(directory=directory,
                      name=name,
                      distance=distance,
                      lattice_tolerance=lattice_tolerance,
                      force_tolerance=force_tolerance,
                      pressure_target=pressure_target,
                      stress_tolerance=stress_tolerance,
                      max_increase=max_increase,
                      max_iteration=max_iteration,
                      min_iteration=min_iteration,
                      symmetry_tolerance=symmetry_tolerance,
                      restrict_offspring=restrict_offspring,
                      max_offspring=max_offspring,
                      cutoff_eigenvalue=cutoff_eigenvalue,
                      max_displacement=max_displacement,
                      num_sampling_points=num_sampling_points,
                      traverse=traverse)

    phr.set_configurations(cell=cell,
                           pseudo_potential_map=pseudo_potential_map,
                           k_mesh=k_mesh,
                           k_shift=k_shift,
                           k_gamma=k_gamma,
                           k_length=k_length,
                           incar=incar)
    phr.set_job(job)

    return phr

