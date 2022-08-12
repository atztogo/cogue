__all__ = [
    "incar",
    "electronic_structure",
    "structure_optimization",
    "bulk_modulus",
    "phonon",
    "elastic_constants",
    "born_effective_charge",
    "mode_gruneisen",
    "quasiharmonic_phonon",
    "phonon_relax",
    "band_structure",
    "density_of_states",
]

import io
import numbers
import os
import shutil

import numpy as np

from cogue.crystal.cell import Cell
from cogue.crystal.converter import atoms2cell
from cogue.crystal.utility import klength2mesh
from cogue.interface.vasp_io import (
    Incar,
    Outcar,
    VaspCell,
    Vasprunxml,
    VasprunxmlExpat,
    change_point_order,
    get_atom_order_from_poscar_yaml,
    read_poscar,
    write_kpoints,
    write_potcar,
)
from cogue.task.band_structure import BandStructureBase
from cogue.task.born_effective_charge import BornEffectiveChargeBase
from cogue.task.bulk_modulus import BulkModulusBase
from cogue.task.density_of_states import DensityOfStatesBase
from cogue.task.elastic_constants import ElasticConstantsBase
from cogue.task.mode_gruneisen import ModeGruneisenBase
from cogue.task.oneshot_calculation import (
    BornEffectiveChargeElementBase,
    ElasticConstantsElementBase,
    ElectronicStructureBase,
    StructureOptimizationElementBase,
)
from cogue.task.phonon import PhononBase
from cogue.task.phonon_fc3 import PhononFC3Base
from cogue.task.phonon_relax import PhononRelaxBase, PhononRelaxElementBase
from cogue.task.quasiharmonic_phonon import QuasiHarmonicPhononBase
from cogue.task.structure_optimization import StructureOptimizationBase


def incar(
    addgrid=None,
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
    sigma=None,
):
    """Returns Incar object"""

    return Incar(
        addgrid=addgrid,
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
        sigma=sigma,
    )


def electronic_structure(
    directory="electronic_structure",
    name=None,
    job=None,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    es = ElectronicStructure(directory=directory, name=name, traverse=traverse)

    es.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    es.set_job(job)

    return es


def structure_optimization(
    directory="structure_optimization",
    name=None,
    job=None,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=0.1,  # GPa=kbar / 10
    max_increase=None,
    max_iteration=5,
    min_iteration=1,
    impose_symmetry=False,
    symmetry_tolerance=0.1,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    so = StructureOptimization(
        directory=directory,
        name=name,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        impose_symmetry=impose_symmetry,
        symmetry_tolerance=symmetry_tolerance,
        traverse=traverse,
    )

    so.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    so.set_job(job)

    return so


def bulk_modulus(
    directory="bulk_modulus",
    name=None,
    job=None,
    strains=None,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    is_cell_relaxed=False,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    bk = BulkModulus(
        directory=directory,
        name=name,
        strains=strains,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        is_cell_relaxed=is_cell_relaxed,
        traverse=traverse,
    )

    bk.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    bk.set_job(job)

    return bk


def band_structure(
    directory="band_structure",
    name=None,
    job=None,
    paths=None,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    is_cell_relaxed=False,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    bs = BandStructure(
        directory=directory,
        name=name,
        paths=paths,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        is_cell_relaxed=is_cell_relaxed,
        traverse=traverse,
    )

    bs.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    bs.set_job(job)

    return bs


def density_of_states(
    directory="density_of_states",
    name=None,
    job=None,
    is_partial_dos=False,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    is_cell_relaxed=False,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    dos = DensityOfStates(
        directory=directory,
        name=name,
        is_partial_dos=is_partial_dos,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        is_cell_relaxed=is_cell_relaxed,
        traverse=traverse,
    )

    dos.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    dos.set_job(job)

    return dos


def phonon(
    directory="phonon",
    name=None,
    job=None,
    supercell_matrix=None,
    primitive_matrix=None,
    nac=False,
    with_perfect=True,
    distance=0.01,
    displace_plusminus="auto",
    displace_diagonal=False,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    is_cell_relaxed=False,
    max_num_atoms=120,
    impose_symmetry=False,
    stop_condition=None,
    symmetry_tolerance=1e-5,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    ph = Phonon(
        directory=directory,
        name=name,
        supercell_matrix=supercell_matrix,
        primitive_matrix=primitive_matrix,
        nac=nac,
        with_perfect=with_perfect,
        distance=distance,
        displace_plusminus=displace_plusminus,
        displace_diagonal=displace_diagonal,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        is_cell_relaxed=is_cell_relaxed,
        max_num_atoms=max_num_atoms,
        impose_symmetry=impose_symmetry,
        stop_condition=stop_condition,
        symmetry_tolerance=symmetry_tolerance,
        traverse=traverse,
    )

    ph.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    ph.set_job(job)

    return ph


def phonon_fc3(
    directory="phonon_fc3",
    name=None,
    job=None,
    supercell_matrix=None,
    primitive_matrix=None,
    with_perfect=True,
    distance=0.03,
    is_diagonal=True,
    check_imaginary=True,
    cutoff_frequency=-0.5,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
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
    k_point=None,
    incar=None,
):

    ph3 = PhononFC3(
        directory=directory,
        name=name,
        supercell_matrix=supercell_matrix,
        primitive_matrix=primitive_matrix,
        with_perfect=with_perfect,
        distance=distance,
        is_diagonal=is_diagonal,
        check_imaginary=check_imaginary,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        traverse=traverse,
        is_cell_relaxed=is_cell_relaxed,
    )

    ph3.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    ph3.set_job(job)

    return ph3


def elastic_constants(
    directory="elastic_constants",
    name=None,
    job=None,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    is_cell_relaxed=False,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    ec = ElasticConstants(
        directory=directory,
        name=name,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        is_cell_relaxed=is_cell_relaxed,
        traverse=traverse,
    )

    ec.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    ec.set_job(job)

    return ec


def born_effective_charge(
    directory="born_effective_charge",
    name=None,
    job=None,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    is_cell_relaxed=False,
    symmetry_tolerance=1.0e-5,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    bec = BornEffectiveCharge(
        directory=directory,
        name=name,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        is_cell_relaxed=is_cell_relaxed,
        symmetry_tolerance=symmetry_tolerance,
        traverse=traverse,
    )

    bec.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    bec.set_job(job)

    return bec


def mode_gruneisen(
    directory="mode_gruneisen",
    name=None,
    job=None,
    delta_strain=0.001,
    strain=None,
    bias=None,
    supercell_matrix=None,
    primitive_matrix=None,
    distance=0.01,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=3,
    min_iteration=1,
    is_cell_relaxed=False,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    mg = ModeGruneisen(
        directory=directory,
        name=name,
        delta_strain=delta_strain,
        strain=strain,
        bias=bias,
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
        is_cell_relaxed=is_cell_relaxed,
        traverse=traverse,
    )

    mg.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    mg.set_job(job)

    return mg


def quasiharmonic_phonon(
    directory="quasiharmonic_phonon",
    name=None,
    job=None,
    strains=None,
    sampling_mesh=None,
    t_step=None,
    t_max=None,
    t_min=None,
    supercell_matrix=None,
    primitive_matrix=None,
    nac=False,
    distance=0.01,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=3,
    min_iteration=1,
    is_cell_relaxed=False,
    max_num_atoms=120,
    first_phonon_index=0,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    qh = QuasiHarmonicPhonon(
        directory=directory,
        name=name,
        strains=strains,
        sampling_mesh=sampling_mesh,
        t_step=t_step,
        t_max=t_max,
        t_min=t_min,
        supercell_matrix=supercell_matrix,
        primitive_matrix=primitive_matrix,
        nac=nac,
        distance=distance,
        lattice_tolerance=lattice_tolerance,
        force_tolerance=force_tolerance,
        pressure_target=pressure_target,
        stress_tolerance=stress_tolerance,
        max_increase=max_increase,
        max_iteration=max_iteration,
        min_iteration=min_iteration,
        is_cell_relaxed=is_cell_relaxed,
        max_num_atoms=max_num_atoms,
        first_phonon_index=first_phonon_index,
        traverse=traverse,
    )

    qh.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    qh.set_job(job)

    return qh


def phonon_relax_element(
    directory="phonon_relax_element",
    name=None,
    job=None,
    distance=0.01,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    symmetry_tolerance=0.1,
    cutoff_eigenvalue=-0.02,
    max_displacement=None,
    num_sampling_points=60,
    stop_condition=None,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    phre = PhononRelaxElement(
        directory=directory,
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
        stop_condition=stop_condition,
        traverse=traverse,
    )

    phre.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    phre.set_job(job)

    return phre


def phonon_relax(
    directory="phonon_relax",
    name=None,
    job=None,
    distance=0.01,
    lattice_tolerance=0.1,
    force_tolerance=1e-3,
    pressure_target=0,
    stress_tolerance=10,
    max_increase=None,
    max_iteration=4,
    min_iteration=1,
    symmetry_tolerance=0.1,
    restrict_offspring=False,
    max_offspring=None,
    cutoff_eigenvalue=-0.02,
    max_displacement=None,
    num_sampling_points=60,
    stop_condition=None,
    traverse=False,
    cell=None,
    pseudo_potential_map=None,
    k_mesh=None,
    k_shift=None,
    k_gamma=None,
    k_length=None,
    k_point=None,
    incar=None,
):

    phr = PhononRelax(
        directory=directory,
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
        stop_condition=stop_condition,
        traverse=traverse,
    )

    phr.set_configurations(
        cell=cell,
        pseudo_potential_map=pseudo_potential_map,
        k_mesh=k_mesh,
        k_shift=k_shift,
        k_gamma=k_gamma,
        k_length=k_length,
        k_point=k_point,
        incar=incar,
    )
    phr.set_job(job)

    return phr


class TaskVasp:
    def set_configurations(
        self,
        cell=None,
        pseudo_potential_map=None,
        k_mesh=None,
        k_shift=None,
        k_gamma=False,
        k_length=None,
        k_point=None,
        incar=None,
    ):
        if not cell:
            self._cell = None
        else:
            self._cell = cell.copy()

        if not pseudo_potential_map:
            self._pseudo_potential_map = None
        else:
            self._pseudo_potential_map = pseudo_potential_map.copy()

        if isinstance(k_mesh, np.ndarray):
            self._k_mesh = list(k_mesh)
        elif k_mesh:
            self._k_mesh = list(k_mesh)
        else:
            self._k_mesh = k_mesh

        if isinstance(k_shift, np.ndarray):
            self._k_shift = list(k_shift)
        elif k_shift:
            self._k_shift = list(k_shift)
        else:
            self._k_shift = k_shift

        if isinstance(k_gamma, tuple):
            self._k_gamma = list(k_gamma)
        else:
            self._k_gamma = k_gamma

        if isinstance(k_length, np.ndarray):
            self._k_length = list(k_length)
        elif isinstance(k_length, tuple):
            self._k_length = list(k_length)
        else:
            self._k_length = k_length

        if isinstance(k_point, np.ndarray):
            self._k_point = list(k_point)
        elif isinstance(k_point, tuple):
            self._k_point = list(k_point)
        else:
            self._k_point = k_point

        if isinstance(incar, tuple):
            self._incar = list(incar)
        else:
            self._incar = incar

    def set_copy_files(self, copy_files):
        self._copy_files = copy_files

    def _prepare(self):
        """
        Create input files for VASP

        We can suppose we are in the calculation directory.
        """
        if os.path.exists("vasprun.xml"):
            os.remove("vasprun.xml")
        if os.path.exists("CONTCAR"):
            os.remove("CONTCAR")

        self._vasp_cell = VaspCell(self._cell)
        self._vasp_cell.write(filename="POSCAR")
        self._vasp_cell.write_yaml(filename="POSCAR.yaml")

        ps_set = [self._pseudo_potential_map[x] for x in self._cell.get_symbols()]
        write_potcar(ps_set)

        if self._k_length:  # Overwrite k_mesh if k_length is given.
            k_mesh = klength2mesh(self._k_length, self._cell.lattice)
            k_gamma = True
            k_shift = [0.0, 0.0, 0.0]
        else:
            k_mesh = self._k_mesh
            k_gamma = self._k_gamma
            k_shift = self._k_shift

        write_kpoints(mesh=k_mesh, shift=k_shift, gamma=k_gamma, kpoint=self._k_point)
        self._incar.write()

        for (fsrc, fdst) in self._copy_files:
            shutil.copy(fsrc, fdst)

    def _choose_configuration(self, index=0):
        # incar
        if isinstance(self._incar, list):
            incar = self._incar[index].copy()
        else:
            incar = self._incar.copy()

        # k_mesh
        k_mesh = self._k_mesh
        if self._k_mesh:
            if None in self._k_mesh:
                k_mesh = self._k_mesh[index]
            elif np.array(self._k_mesh).ndim == 2:
                k_mesh = self._k_mesh[index]

        # k_shift
        k_shift = self._k_shift
        if self._k_shift:
            if None in self._k_shift:
                k_shift = self._k_shift[index]
            elif np.array(self._k_shift).ndim == 2:
                k_shift = self._k_shift[index]

        # k_gamma
        if isinstance(self._k_gamma, list):
            k_gamma = self._k_gamma[index]
        else:
            k_gamma = self._k_gamma

        # k_length
        if isinstance(self._k_length, list):
            k_length = self._k_length[index]
        else:
            k_length = self._k_length

        # k_point
        if self._k_point is None:
            k_point = None
        elif isinstance(self._k_point[0], numbers.Number):
            k_point = self._k_point
        elif isinstance(self._k_point, list):
            k_point = self._k_point[index]
        else:  # I don't know this case.
            k_point = self._k_point

        # job
        if isinstance(self._job, list):
            job = self._job[index]
        else:
            job = self._job

        return (
            job,
            incar,
            {
                "mesh": k_mesh,
                "shift": k_shift,
                "gamma": k_gamma,
                "length": k_length,
                "kpoint": k_point,
            },
        )

    def _get_equilibrium_task(
        self,
        index=0,
        cell=None,
        impose_symmetry=False,
        symmetry_tolerance=None,
        max_iteration=None,
        min_iteration=None,
        directory="equilibrium",
    ):
        if not cell:
            cell = self._cell
        job, incar, kpoints = self._choose_configuration(index=index)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        k_point = kpoints["kpoint"]

        if max_iteration is None:
            _max_iteration = self._max_iteration
        else:
            _max_iteration = max_iteration
        if min_iteration is None:
            _min_iteration = self._min_iteration
        else:
            _min_iteration = min_iteration

        if symmetry_tolerance:
            task = StructureOptimization(
                directory=directory,
                lattice_tolerance=self._lattice_tolerance,
                force_tolerance=self._force_tolerance,
                pressure_target=self._pressure_target,
                stress_tolerance=self._stress_tolerance,
                max_increase=self._max_increase,
                max_iteration=_max_iteration,
                min_iteration=_min_iteration,
                impose_symmetry=impose_symmetry,
                symmetry_tolerance=symmetry_tolerance,
                traverse=self._traverse,
            )
        else:  # Use default symmetry_tolerance
            task = StructureOptimization(
                directory=directory,
                lattice_tolerance=self._lattice_tolerance,
                force_tolerance=self._force_tolerance,
                pressure_target=self._pressure_target,
                stress_tolerance=self._stress_tolerance,
                max_increase=self._max_increase,
                max_iteration=_max_iteration,
                min_iteration=_min_iteration,
                impose_symmetry=impose_symmetry,
                traverse=self._traverse,
            )
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            k_point=k_point,
            incar=incar,
        )
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))
        return task

    def _get_charge_density_task(self, cell, index=1, directory="charge_density"):
        job, incar, kpoints = self._choose_configuration(index=index)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        incar.set_lcharg(True)
        incar.set_ibrion(-1)
        incar.set_nsw(None)
        incar.set_isif(None)
        incar.set_ediffg(None)

        task = ElectronicStructure(directory=directory, traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            incar=incar,
        )
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))

        return task


class ElectronicStructure(TaskVasp, ElectronicStructureBase):
    """ """

    def __init__(self, directory="electronic_structure", name=None, traverse=False):

        ElectronicStructureBase.__init__(
            self, directory=directory, name=name, traverse=traverse
        )

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._incar = None
        self._copy_files = []

    def _collect(self):
        """Collect information from output files of VASP.

        self._status of "done" or "terminate"  is stored.
        self._log: Terminate log is stored.

        """

        if os.path.exists("POSCAR.yaml"):
            atom_order = get_atom_order_from_poscar_yaml("POSCAR.yaml")
        else:
            atom_order = None

        if not os.path.exists("vasprun.xml"):
            self._log += "    vasprun.xml not exists.\n"
            self._status = "terminate"
        else:
            vxml = Vasprunxml("vasprun.xml")
            if (
                vxml.parse_calculation()
                and vxml.parse_eigenvalues()
                and vxml.parse_efermi()
                and vxml.parse_parameters()
            ):
                kpoints, weights = vxml.get_kpoints()
                if atom_order:
                    force_sets = vxml.get_forces()[:, atom_order, :]
                else:
                    force_sets = vxml.get_forces()
                self._properties = {
                    "stress": vxml.get_stress(),
                    "forces": force_sets,
                    "energies": vxml.get_energies()[:, 1],
                    "eigenvalues": vxml.get_eigenvalues(),
                    "occupancies": vxml.get_occupancies(),
                    "kpoints": kpoints,
                    "kpoint-weights": weights,
                    "fermi-energy": vxml.get_efermi(),
                    "nbands": vxml.get_nbands(),
                }
                self._status = "done"
            else:
                self._log += vxml.log
                self._log += "    Failed to parse vasprun.xml.\n"
                self._status = "terminate"


class StructureOptimizationElement(TaskVasp, StructureOptimizationElementBase):
    def __init__(
        self,
        directory="structure_optimization_element",
        name=None,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=0.1,  # GPa=kbar / 10
        max_increase=None,
        traverse=False,
    ):

        StructureOptimizationElementBase.__init__(
            self,
            directory=directory,
            name=name,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            traverse=traverse,
        )

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._incar = None
        self._copy_files = []
        self._atom_order = None

    def _collect(self):
        """Collect information from output files of VASP.

        self._current_cell: Final crystal structure of each relaxation steps.
        self._status: "next", "done", or "terminate" is stored.
        self._log: Logs

        """

        if os.path.exists("POSCAR.yaml"):
            self._atom_order = get_atom_order_from_poscar_yaml("POSCAR.yaml")
        else:
            self._atom_order = None

        if os.path.exists("POSCAR"):
            masses = self.get_current_cell().get_masses()
            cell = read_poscar("POSCAR")
            if self._atom_order:
                self._current_cell = change_point_order(cell, self._atom_order)
            else:
                self._current_cell = cell
            self._current_cell.set_masses(masses)

        if not os.path.exists("vasprun.xml"):
            self._log += "    vasprun.xml not exists.\n"
            self._status = "terminate"
        else:
            with io.open("vasprun.xml", "rb") as f:
                vxml = VasprunxmlExpat(f)
                is_success = vxml.parse()

            # [num_geomopt, 3, 3]
            lattice = np.array([lat.T for lat in vxml.get_lattice()])
            # [num_geomopt, 3, num_atoms]
            points = np.array([pos.T for pos in vxml.get_points()])
            forces = vxml.get_forces()  # [num_geomopt, num_atoms, 3]
            stress = vxml.get_stress()  # [num_geomopt, 3, 3]
            energies = vxml.get_energies()  # [num_geomopt, 3]

            max_iter = len(lattice)
            if len(points) < max_iter:
                max_iter = len(points)
            if len(forces) < max_iter:
                max_iter = len(forces)
            if len(stress) < max_iter:
                max_iter = len(stress)
            if len(energies) < max_iter:
                max_iter = len(energies)

            if is_success and max_iter > 0:
                self._stress = stress[max_iter - 1] / 10
                self._energy = energies[max_iter - 1, 1]
                if self._atom_order:
                    _points = points[max_iter - 1][:, self._atom_order]
                    self._forces = forces[max_iter - 1][self._atom_order, :]
                else:
                    _points = points[max_iter - 1]
                    self._forces = forces[max_iter - 1]
                self._judge(lattice[max_iter - 1], _points)
            elif (not is_success) and max_iter > 2:
                self._log += "    vasprun.xml is not cleanly closed.\n"
                self._stress = stress[max_iter - 3] / 10
                self._energy = energies[max_iter - 3, 1]
                if self._atom_order:
                    _points = points[max_iter - 3][:, self._atom_order]
                    self._forces = forces[max_iter - 3][self._atom_order, :]
                else:
                    _points = points[max_iter - 3]
                    self._forces = forces[max_iter - 3]
                self._judge(lattice[max_iter - 3], _points)
            else:
                self._log += vxml.log
                self._log += "    Failed to parse vasprun.xml.\n"
                self._current_cell = None
                self._status = "terminate"

    def _judge(self, lattice_last, points):
        lattice_init = self._current_cell.lattice
        vecs2_init = np.diag(np.dot(lattice_init.T, lattice_init))
        vecs2_last = np.diag(np.dot(lattice_last.T, lattice_last))
        d_vecs2_ratio = (vecs2_last - vecs2_init) / vecs2_init

        masses = self._current_cell.get_masses()
        cell = Cell(
            lattice=lattice_last,
            points=points,
            symbols=self._current_cell.get_symbols(),
            masses=masses,
        )

        if os.path.exists("CONTCAR"):
            try:
                self._current_cell = read_poscar("CONTCAR")
                if self._atom_order:
                    self._current_cell = change_point_order(cell, self._atom_order)
                else:
                    self._current_cell = cell
                self._current_cell.set_masses(masses)
            except:  # noqa E722
                self._current_cell = cell
            else:
                # Sometimes CONTCAR's structure becomes original
                # structure same as POSCAR. In this case, the third
                # relaxed structure from the last is used for the new
                # cell.
                if (abs(self._current_cell.lattice - lattice_init) < 1e-12).all():
                    self._current_cell = cell
        else:
            self._current_cell = cell

        # Non termination conditions
        if self._lattice_tolerance is not None:
            if (abs(d_vecs2_ratio) > self._lattice_tolerance**2).any():
                self._log += "    Lattice is not enough relaxed.\n"
                self._status = "next"

        if self._stress_tolerance is not None:
            if (
                abs(self._stress - np.eye(3) * self._pressure_target)
                > self._stress_tolerance
            ).any():
                self._log += "    Stress is not enough relaxed.\n"
                self._status = "next"

        if (abs(self._forces) > self._force_tolerance).any():
            self._log += "    Forces are not enough relaxed.\n"
            self._status = "next"

        if not self._status == "next":
            self._status = "done"

        # Termination condition
        vol_last = np.linalg.det(lattice_last)
        vol_init = np.linalg.det(lattice_init)
        if self._max_increase is not None:
            if vol_last > self._max_increase * vol_init:
                self._log += "    Too large volume expansion.\n"
                self._status = "terminate"


class StructureOptimization(TaskVasp, StructureOptimizationBase):
    def __init__(
        self,
        directory="structure_optimization",
        name=None,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=0.1,  # GPa=kbar / 10
        max_increase=None,
        max_iteration=5,
        min_iteration=1,
        impose_symmetry=False,
        symmetry_tolerance=0.1,
        traverse=False,
    ):

        StructureOptimizationBase.__init__(
            self,
            directory=directory,
            name=name,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            impose_symmetry=impose_symmetry,
            symmetry_tolerance=symmetry_tolerance,
            traverse=traverse,
        )

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._k_point = None
        self._incar = None

    def _get_next_task(self, cell):
        task = StructureOptimizationElement(
            directory="structopt-%d" % self._stage,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            traverse=self._traverse,
        )

        task.set_configurations(
            cell=cell.copy(),
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            k_point=self._k_point,
            incar=self._incar.copy(),
        )

        task.set_job(self._job.copy("%s-%s" % (self._job.get_jobname(), self._stage)))
        return task


class BulkModulus(TaskVasp, BulkModulusBase):
    """Task to calculate bulk modulus by VASP."""

    def __init__(
        self,
        directory="bulk_modulus",
        name=None,
        strains=None,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=3,
        min_iteration=1,
        is_cell_relaxed=False,
        traverse=False,
    ):

        BulkModulusBase.__init__(
            self,
            directory=directory,
            name=name,
            strains=strains,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            traverse=traverse,
        )

    def _get_bm_task(self, cell, directory):
        job, incar, kpoints = self._choose_configuration(index=1)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        incar.set_isif(4)

        task = StructureOptimization(
            directory=directory,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            max_iteration=1,
            min_iteration=1,
            traverse=self._traverse,
        )

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            incar=incar,
        )
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))

        return task


class BandStructure(TaskVasp, BandStructureBase):
    """Task to calculate band structure by VASP."""

    def __init__(
        self,
        directory="band_structure",
        name=None,
        paths=None,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=3,
        min_iteration=1,
        is_cell_relaxed=False,
        traverse=False,
    ):

        BandStructureBase.__init__(
            self,
            directory=directory,
            name=name,
            paths=paths,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            traverse=traverse,
        )

    def _get_band_point_tasks(self, cell, properties=None, chgcar_dir="charge_density"):
        directory = "path"
        tasks = []
        all_kpoints = []
        for kpoints in self._paths:
            all_kpoints += list(kpoints)
        for i, kpoint in enumerate(all_kpoints):
            job, incar, kpoints = self._choose_configuration(index=1)
            incar.set_icharg(11)
            incar.set_ibrion(-1)
            incar.set_nsw(None)
            incar.set_isif(None)
            incar.set_ediffg(None)
            if properties:
                if "nbands" in properties:
                    incar.set_nbands(2 * properties["nbands"])
            task = ElectronicStructure(
                directory=directory + "%04d" % i, traverse=self._traverse
            )
            task.set_configurations(
                cell=cell,
                pseudo_potential_map=self._pseudo_potential_map,
                k_point=kpoint,
                incar=incar,
            )
            task.set_job(job.copy("%s-%s%d" % (job.get_jobname(), directory, i)))
            task.set_copy_files([("../%s/CHGCAR" % chgcar_dir, "CHGCAR")])
            tasks.append(task)

        return tasks


class DensityOfStates(TaskVasp, DensityOfStatesBase):
    """Task to calculate density of states by VASP."""

    def __init__(
        self,
        directory="density_of_states",
        name=None,
        is_partial_dos=False,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=3,
        min_iteration=1,
        is_cell_relaxed=False,
        traverse=False,
    ):

        DensityOfStatesBase.__init__(
            self,
            directory=directory,
            name=name,
            is_partial_dos=is_partial_dos,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            traverse=traverse,
        )

    def _get_dos_task(self, cell, properties=None, chgcar_dir="charge_density"):
        directory = "dos"

        job, incar, kpoints = self._choose_configuration(index=2)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        k_point = kpoints["kpoint"]

        incar.set_icharg(11)
        incar.set_ibrion(-1)
        incar.set_nsw(None)
        incar.set_isif(None)
        incar.set_ediffg(None)
        if self._is_partial_dos:
            incar.set_lorbit(12)
        if properties:
            if "nbands" in properties:
                incar.set_nbands(2 * properties["nbands"])
        task = ElectronicStructure(directory=directory, traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            k_point=k_point,
            incar=incar,
        )
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))
        task.set_copy_files([("../%s/CHGCAR" % chgcar_dir, "CHGCAR")])

        return task


class TaskVaspPhonon:
    """Phonon calculation configuration class

    Normally the configurations are stored in a list. Each index is
    fixed for specific type of calculation in the task list.

    index:
        0: Structure optimization
        1: Displacement. For job configuration, if this is a list and k_mesh is
           Gamma-only, the second one is used.
        2: Born effective charge. If incar is given as a list, structure
           optimization is with fixed lattice is executed.

    """

    def _get_vasp_displacement_tasks(
        self, phonon, start=None, stop=None, digit_number=3
    ):
        incar = self._incar[1].copy()
        if start is None:
            istart = 0
        else:
            istart = start
        supercells = phonon.supercells_with_displacements
        if stop is None:
            disp_cells = supercells[istart:]
        else:
            disp_cells = supercells[istart:stop]

        tasks = []
        if start is None and self._with_perfect:
            tasks.append(
                self._get_disp_task(
                    atoms2cell(phonon.supercell), incar, 0, digit_number=digit_number
                )
            )

        for i, disp in enumerate(disp_cells):
            tasks.append(
                self._get_disp_task(
                    atoms2cell(disp), incar, i + 1 + istart, digit_number=digit_number
                )
            )
        return tasks

    def _get_disp_task(self, cell, incar, disp_number, digit_number=3):
        job, incar, kpoints = self._choose_configuration(index=1)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        k_point = kpoints["kpoint"]

        directory = ("disp-%0" + "%d" % digit_number + "d") % disp_number

        if k_length:
            k_mesh = klength2mesh(k_length, cell.lattice)
            k_gamma = True
            k_shift = [0.0, 0.0, 0.0]

        # For Gamma-only VASP, take secon element of job
        if isinstance(job, list):
            job_disp = job[0]
            if (np.array(k_mesh) == 1).all() and k_shift:
                if (np.abs(k_shift) < 1e-10).all():
                    job_disp = job[1]
            job = job_disp

        task = ElectronicStructure(directory=directory, traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_point=k_point,
            incar=incar,
        )
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))

        return task

    def _get_nac_task(self, is_cell_relaxed=True, directory="nac"):
        job, incar, kpoints = self._choose_configuration(index=2)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        k_point = kpoints["kpoint"]

        _is_cell_relaxed = is_cell_relaxed

        if isinstance(incar, list):
            _is_cell_relaxed = False
            incar[0].set_isif(2)
            incar[0].set_lepsilon(None)
            if incar[0].get_nsw() is None:
                incar[0].set_nsw(10)
                incar[0].set_ediffg(-1.0e-8)
        elif not is_cell_relaxed:
            _is_cell_relaxed = False
            incar_rx = incar.copy()
            incar_rx.set_ibrion(2)
            incar_rx.set_nsw(10)
            incar_rx.set_isif(2)
            incar_rx.set_ediffg(-1.0e-8)
            incar = [incar_rx, incar]

        if not _is_cell_relaxed:
            job_rx = job.copy(job.get_jobname() + "-rx")
            job = [job_rx, job]

        task = BornEffectiveCharge(
            directory="nac",
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=None,
            max_increase=self._max_increase,
            max_iteration=self._max_iteration,
            min_iteration=self._min_iteration,
            is_cell_relaxed=_is_cell_relaxed,
            traverse=self._traverse,
        )

        task.set_configurations(
            cell=self.get_cell(),
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            k_point=k_point,
            incar=incar,
        )

        task.set_job(job)

        return task


class TaskVaspQHA:
    """QHA calculation configuration class

    Normally the configurations are stored in a list. Each index is
    fixed for specific type of calculation in the task list.

    index:
        0: Structure optimization
        0-1: Mapped to Bulk modulus configuration
             ISIF=4 is forced for index=1.
        1-3: Mapped to 0-2 of Phonon configuration

    """

    def _get_phonon_task(
        self, cell, directory, stress_tolerance=None, is_cell_relaxed=False
    ):
        task = Phonon(
            directory=directory,
            supercell_matrix=self._supercell_matrix,
            primitive_matrix=self._primitive_matrix,
            nac=self._nac,
            distance=self._distance,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=self._max_increase,
            max_iteration=self._max_iteration,
            min_iteration=self._min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            traverse=self._traverse,
        )

        if isinstance(self._job, list):
            job = [
                j.copy("%s-%s" % (j.get_jobname(), directory)) for j in self._job[1:]
            ]
        else:
            job = self._job.copy("%s-%s" % (self._job.get_jobname(), directory))

        k_mesh = self._k_mesh
        if self._k_mesh:
            if None in self._k_mesh:
                k_mesh = self._k_mesh[1:]
            elif np.array(self._k_mesh).ndim > 1:
                k_mesh = self._k_mesh[1:]

        k_shift = self._k_shift
        if self._k_shift:
            if None in self._k_shift:
                k_shift = self._k_shift[1:]
            elif np.array(self._k_shift).ndim > 1:
                k_shift = self._k_shift[1:]

        if isinstance(self._k_gamma, list):
            k_gamma = self._k_gamma[1:]
        else:
            k_gamma = self._k_gamma

        if isinstance(self._k_length, list):
            k_length = self._k_length[1:]
        else:
            k_length = self._k_length

        if isinstance(self._k_point, list):
            k_point = self._k_point[1:]
        else:
            k_point = self._k_point

        if isinstance(self._incar, list):
            incar = [x.copy() for x in self._incar[1:]]
        else:
            incar = self._incar.copy()

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            k_point=k_point,
            incar=incar,
        )
        task.set_job(job)

        return task


class Phonon(TaskVasp, TaskVaspPhonon, PhononBase):
    def __init__(
        self,
        directory="phonon",
        name=None,
        supercell_matrix=None,
        primitive_matrix=None,
        nac=False,
        with_perfect=True,
        distance=0.01,
        displace_plusminus="auto",
        displace_diagonal=False,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=3,
        min_iteration=1,
        is_cell_relaxed=False,
        max_num_atoms=120,
        impose_symmetry=False,
        stop_condition=None,
        symmetry_tolerance=None,
        traverse=False,
    ):

        PhononBase.__init__(
            self,
            directory=directory,
            name=name,
            supercell_matrix=supercell_matrix,
            primitive_matrix=primitive_matrix,
            nac=nac,
            with_perfect=with_perfect,
            distance=distance,
            displace_plusminus=displace_plusminus,
            displace_diagonal=displace_diagonal,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            max_num_atoms=max_num_atoms,
            impose_symmetry=impose_symmetry,
            stop_condition=stop_condition,
            symmetry_tolerance=symmetry_tolerance,
            traverse=traverse,
        )

    def _get_displacement_tasks(self, start=None, stop=None):
        return self._get_vasp_displacement_tasks(self._phonon, start=start, stop=stop)


class PhononFC3(TaskVasp, TaskVaspPhonon, PhononFC3Base):
    def __init__(
        self,
        directory="phonon_fc3",
        name=None,
        supercell_matrix=None,
        primitive_matrix=None,
        with_perfect=True,
        distance=0.03,
        is_diagonal=True,
        check_imaginary=True,
        cutoff_frequency=-0.5,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=3,
        min_iteration=1,
        is_cell_relaxed=False,
        traverse=False,
    ):

        PhononFC3Base.__init__(
            self,
            directory=directory,
            name=name,
            supercell_matrix=supercell_matrix,
            primitive_matrix=primitive_matrix,
            with_perfect=with_perfect,
            distance=distance,
            is_diagonal=is_diagonal,
            check_imaginary=check_imaginary,
            cutoff_frequency=cutoff_frequency,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            traverse=traverse,
        )

    def _get_displacement_tasks(self, start=None, stop=None):
        return self._get_vasp_displacement_tasks(
            self._phonon_fc3, start=start, stop=stop, digit_number=5
        )


class ElasticConstants(TaskVasp, ElasticConstantsBase):
    def __init__(
        self,
        directory="elastic_constants",
        name=None,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=4,
        min_iteration=1,
        is_cell_relaxed=False,
        traverse=False,
    ):

        ElasticConstantsBase.__init__(
            self,
            directory=directory,
            name=name,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            traverse=traverse,
        )

    def _get_ec_task(self, cell):
        job, incar, kpoints = self._choose_configuration(index=1)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        incar.set_ibrion(6)
        incar.set_isif(3)
        incar.set_nsw(1)

        directory = "elastic_constants"
        task = ElasticConstantsElement(directory=directory, traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            incar=incar,
        )
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))
        return task


class ElasticConstantsElement(TaskVasp, ElasticConstantsElementBase):
    """ """

    def __init__(self, directory="elastic_constants", name=None, traverse=False):

        ElasticConstantsElementBase.__init__(
            self, directory=directory, name=name, traverse=traverse
        )

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._k_point = None
        self._incar = None
        self._copy_files = []

    def _collect(self):
        """Collect information from output files of VASP.

        self._status of "done" or "terminate"  is stored.
        self._log: Terminate log is stored.

        """

        if not os.path.exists("OUTCAR"):
            self._log += "    OUTCAR not exists.\n"
            self._status = "terminate"
        else:
            outcar = Outcar("OUTCAR")
            if outcar.parse_elastic_constants():
                self._elastic_constants = outcar.get_elastic_constants()
                self._status = "done"
            else:
                self._log += "    Failed to parse OUTCAR.\n"
                self._status = "terminate"


class BornEffectiveCharge(TaskVasp, BornEffectiveChargeBase):
    def __init__(
        self,
        directory="born_effective_charge",
        name=None,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=4,
        min_iteration=1,
        is_cell_relaxed=False,
        symmetry_tolerance=1.0e-5,
        traverse=False,
    ):

        BornEffectiveChargeBase.__init__(
            self,
            directory=directory,
            name=name,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            symmetry_tolerance=symmetry_tolerance,
            traverse=traverse,
        )

    def _get_bec_task(self, cell):
        job, incar, kpoints = self._choose_configuration(index=1)
        k_mesh = kpoints["mesh"]
        k_shift = kpoints["shift"]
        k_gamma = kpoints["gamma"]
        k_length = kpoints["length"]
        incar.set_ibrion(-1)
        incar.set_nsw(None)
        incar.set_lepsilon(True)

        directory = "born_effective_charge"
        task = BornEffectiveChargeElement(directory=directory, traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            incar=incar,
        )
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))
        return task


class BornEffectiveChargeElement(TaskVasp, BornEffectiveChargeElementBase):
    """ """

    def __init__(self, directory="born_effective_charge", name=None, traverse=False):

        BornEffectiveChargeElementBase.__init__(
            self, directory=directory, name=name, traverse=traverse
        )

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._k_point = None
        self._incar = None
        self._copy_files = []

    def _collect(self):
        """Collect information from vasprun.xml

        self._status of "done" or "terminate"  is stored.
        self._log: Terminate log is stored.

        """

        if not os.path.exists("vasprun.xml"):
            self._log += "    vasprun.xml not exists.\n"
            self._status = "terminate"
        else:
            with io.open("vasprun.xml", "rb") as f:
                vxml = VasprunxmlExpat(f)
                is_success = vxml.parse()
            if is_success:
                born = vxml.get_born()
                epsilon = vxml.get_epsilon()

            if is_success and born is not None and epsilon is not None:
                if os.path.exists("POSCAR.yaml"):
                    atom_order = get_atom_order_from_poscar_yaml("POSCAR.yaml")
                    self._born = born[atom_order]
                else:
                    self._born = born
                self._epsilon = epsilon
                self._status = "done"
            else:
                self._log += "    Failed to parse vasprun.xml for\n"
                self._log += "    Born effective charge and dielectric constant"
                self._log += ".\n"
                self._status = "terminate"


class ModeGruneisen(TaskVasp, TaskVaspQHA, ModeGruneisenBase):
    """Task to calculate mode Gruneisen parameters by VASP."""

    def __init__(
        self,
        directory="mode_gruneisen",
        name=None,
        delta_strain=0.001,
        strain=None,
        bias=None,
        supercell_matrix=None,
        primitive_matrix=None,
        distance=0.01,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=3,
        min_iteration=1,
        is_cell_relaxed=False,
        traverse=False,
    ):

        ModeGruneisenBase.__init__(
            self,
            directory=directory,
            name=name,
            delta_strain=delta_strain,
            strain=strain,
            bias=bias,
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
            is_cell_relaxed=is_cell_relaxed,
            traverse=traverse,
        )


class QuasiHarmonicPhonon(TaskVasp, TaskVaspQHA, QuasiHarmonicPhononBase):
    """Task to calculate quasi-harmonic phonons by VASP."""

    def __init__(
        self,
        directory="quasiharmonic_phonon",
        name=None,
        strains=None,
        sampling_mesh=None,
        t_step=None,
        t_max=None,
        t_min=None,
        supercell_matrix=None,
        primitive_matrix=None,
        nac=False,
        distance=0.01,
        lattice_tolerance=0.1,
        force_tolerance=1e-3,
        pressure_target=0,
        stress_tolerance=10,
        max_increase=None,
        max_iteration=3,
        min_iteration=1,
        is_cell_relaxed=False,
        max_num_atoms=120,
        first_phonon_index=0,
        traverse=False,
    ):

        QuasiHarmonicPhononBase.__init__(
            self,
            directory=directory,
            name=name,
            strains=strains,
            sampling_mesh=sampling_mesh,
            t_step=t_step,
            t_max=t_max,
            t_min=t_min,
            supercell_matrix=supercell_matrix,
            primitive_matrix=primitive_matrix,
            nac=nac,
            distance=distance,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            max_num_atoms=max_num_atoms,
            first_phonon_index=first_phonon_index,
            traverse=traverse,
        )

    def _get_bulk_modulus_task(
        self,
        cell,
        strains,
        is_cell_relaxed=True,
        max_iteration=3,
        min_iteration=1,
        directory="eos",
    ):
        task = BulkModulus(
            directory=directory,
            strains=strains,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            max_iteration=max_iteration,
            min_iteration=min_iteration,
            is_cell_relaxed=is_cell_relaxed,
            traverse=self._traverse,
        )

        if isinstance(self._job, list):
            job = [
                j.copy("%s-%s" % (j.get_jobname(), directory)) for j in self._job[:2]
            ]
        else:
            job = self._job.copy("%s-%s" % (self._job.get_jobname(), directory))

        k_mesh = self._k_mesh
        if self._k_mesh:
            if None in self._k_mesh:
                k_mesh = self._k_mesh[:2]
            elif np.array(self._k_mesh).ndim > 1:
                k_mesh = self._k_mesh[:2]

        k_shift = self._k_shift
        if self._k_shift:
            if None in self._k_shift:
                k_shift = self._k_shift[:2]
            elif np.array(self._k_shift).ndim > 1:
                k_shift = self._k_shift[:2]

        if isinstance(self._k_gamma, list):
            k_gamma = self._k_gamma[:2]
        else:
            k_gamma = self._k_gamma

        if isinstance(self._k_length, list):
            k_length = self._k_length[:2]
        else:
            k_length = self._k_length

        if isinstance(self._k_point, list):
            k_point = self._k_point[:2]
        else:
            k_point = self._k_point

        if isinstance(self._incar, list):
            incar = [x.copy() for x in self._incar[:2]]
        else:
            incar = self._incar.copy()

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            k_point=k_point,
            incar=incar,
        )
        task.set_job(job)

        return task


class PhononRelax(TaskVasp, PhononRelaxBase):
    def __init__(
        self,
        directory="phonon_relax",
        name=None,
        ancestral_cells={},
        distance=0.01,
        lattice_tolerance=0.1,
        force_tolerance=1e-5,
        pressure_target=0,
        stress_tolerance=1,
        max_increase=None,
        max_iteration=10,
        min_iteration=5,
        symmetry_tolerance=0.1,
        restrict_offspring=False,
        max_offspring=None,
        cutoff_eigenvalue=-0.02,
        max_displacement=None,
        num_sampling_points=60,
        stop_condition=None,
        traverse=False,
    ):

        PhononRelaxBase.__init__(
            self,
            directory=directory,
            name=name,
            ancestral_cells=ancestral_cells,
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
            stop_condition=stop_condition,
            traverse=traverse,
        )

    def _get_phonon_relax_element_task(self, cell):
        task = PhononRelaxElement(
            directory="phonon_relax_element",
            ancestral_cells=self._ancestral_cells,
            tid_parent=self._tid,
            distance=self._distance,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            max_iteration=self._max_iteration,
            min_iteration=self._min_iteration,
            symmetry_tolerance=self._symmetry_tolerance,
            cutoff_eigenvalue=self._cutoff_eigenvalue,
            max_displacement=self._max_displacement,
            num_sampling_points=self._num_sampling_points,
            stop_condition=self._stop_condition,
            traverse=self._traverse,
        )

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            k_point=self._k_point,
            incar=self._incar,
        )
        task.set_job(self._job)

        return task

    def _get_phonon_relax_task(self, cell, ancestral_cells, directory):
        task = PhononRelax(
            directory=directory,
            ancestral_cells=ancestral_cells,
            distance=self._distance,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            max_iteration=self._max_iteration,
            min_iteration=self._min_iteration,
            symmetry_tolerance=self._symmetry_tolerance,
            restrict_offspring=self._restrict_offspring,
            max_offspring=self._max_offspring,
            cutoff_eigenvalue=self._cutoff_eigenvalue,
            max_displacement=self._max_displacement,
            num_sampling_points=self._num_sampling_points,
            stop_condition=self._stop_condition,
            traverse=self._traverse,
        )

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            k_point=self._k_point,
            incar=self._incar,
        )
        task.set_job(self._job)

        return task


class PhononRelaxElement(TaskVasp, PhononRelaxElementBase):
    def __init__(
        self,
        directory="phonon_relax_element",
        name=None,
        ancestral_cells={},
        tid_parent=None,
        distance=0.01,
        lattice_tolerance=0.1,
        force_tolerance=1e-5,
        pressure_target=0,
        stress_tolerance=1,
        max_increase=None,
        max_iteration=10,
        min_iteration=5,
        symmetry_tolerance=0.1,
        cutoff_eigenvalue=-0.02,
        max_displacement=None,
        num_sampling_points=60,
        stop_condition=None,
        traverse=False,
    ):

        PhononRelaxElementBase.__init__(
            self,
            directory=directory,
            name=name,
            ancestral_cells=ancestral_cells,
            tid_parent=tid_parent,
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
            stop_condition=stop_condition,
            traverse=traverse,
        )

    def _get_phonon_task(self, cell, supercell_matrix, directory):
        task = Phonon(
            directory=directory,
            supercell_matrix=supercell_matrix,
            distance=self._distance,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            max_iteration=1,
            min_iteration=1,
            is_cell_relaxed=False,
            stop_condition=self._stop_condition,
            traverse=self._traverse,
        )

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            k_point=self._k_point,
            incar=self._incar[1:3],
        )
        task.set_job(self._job)

        return task
