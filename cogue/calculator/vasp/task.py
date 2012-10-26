import os
import numpy as np
from cogue.crystal.cell import Cell
from cogue.crystal.converter import atoms2cell
from cogue.task.oneshot_calculation import *
from cogue.task.structure_optimization import *
from cogue.task.bulk_modulus import *
from cogue.task.phonon import *
from cogue.task.phonon_relax import *
from cogue.task.elastic_constants import *
from cogue.crystal.vasp_io import *

def klength2mesh(k_length, lattice):
    """Convert length to mesh in k-point sampling

    This conversion follows VASP manual.

    """
    rec_lattice = np.linalg.inv(lattice).T
    rec_lat_lengths = np.sqrt(np.diagonal(np.dot(rec_lattice.T, rec_lattice)))
    k_mesh = (rec_lat_lengths * k_length + 0.5).astype(int)
    return np.maximum(k_mesh, [1, 1, 1])

class TaskVasp:
    def set_configurations(self,
                           cell=None,
                           pseudo_potential_map=None,
                           k_mesh=None,
                           k_shift=None,
                           k_gamma=False,
                           k_length=None,
                           incar=None):

        if not cell:
            self._cell = None
        else:
            self._cell = cell.copy()

        if not pseudo_potential_map:
            self._pseudo_potential_map = None
        else:
            self._pseudo_potential_map = pseudo_potential_map.copy()

        self._k_mesh = k_mesh
        self._k_gamma = k_gamma
        self._k_shift = k_shift
        self._k_length = k_length
        self._incar = incar

    def _prepare(self):
        """
        Create input files for VASP
        
        We can suppose we are in the calculation directory.
        """
        if os.path.exists("vasprun.xml"):
            os.remove("vasprun.xml")
        if os.path.exists("CONTCAR"):
            os.remove("CONTCAR")
        write_poscar(self._cell, "POSCAR")
        ps_set = [self._pseudo_potential_map[x]
                  for x in self._cell.get_symbols()]
        write_potcar(ps_set)

        if self._k_length: # Overwrite k_mesh if k_length is given.
            k_mesh = klength2mesh(self._k_length, self._cell.get_lattice())
            k_gamma = True
            k_shift = [0.0, 0.0, 0.0]
        else:
            k_mesh = self._k_mesh
            k_gamma = self._k_gamma
            k_shift = self._k_shift
            
        write_kpoints(mesh=k_mesh,
                      shift=k_shift,
                      gamma=k_gamma)
        self._incar.write()

    def _choose_configuration(self, index=0):
        if (isinstance(self._incar, list) or
            isinstance(self._incar, tuple)):
            incar = self._incar[index].copy()
        else:
            incar = self._incar.copy()

        if np.array(self._k_mesh).ndim == 2:
            k_mesh = self._k_mesh[index]
        else:
            k_mesh = self._k_mesh

        if self._k_shift == None:
            k_shift = None
        elif np.array(self._k_shift).ndim == 2:
            k_shift = self._k_shift[index]
        else:
            k_shift = self._k_shift

        if not self._k_gamma:
            k_gamma = None
        elif (isinstance(self._k_gamma, list) or 
              isinstance(self._k_gamma, tuple)):
            k_gamma = self._k_gamma[index]
        else:
            k_gamma = self._k_gamma

        if not self._k_length:
            k_length = None
        elif (isinstance(self._k_length, list) or 
              isinstance(self._k_length, tuple)):
            k_length = self._k_length[index]
        else:
            k_length = self._k_length

        if (isinstance(self._job, list) or
            isinstance(self._job, tuple)):
            job = self._job[index]
        else:
            job = self._job

        return job, incar, k_mesh, k_shift, k_gamma, k_length

    def _get_equilibrium_task(self,
                              index=0,
                              cell=None,
                              impose_symmetry=False,
                              symmetry_tolerance=None,
                              directory="equilibrium"):
        if not cell:
            cell = self._cell
        job, incar, k_mesh, k_shift, k_gamma, k_length = \
            self._choose_configuration(index=index)

        if symmetry_tolerance:
            task = StructureOptimization(
                directory=directory,
                lattice_tolerance=self._lattice_tolerance,
                force_tolerance=self._force_tolerance,
                pressure_target=self._pressure_target,
                stress_tolerance=self._stress_tolerance,
                max_increase=self._max_increase,
                max_iteration=self._max_iteration,
                min_iteration=self._min_iteration,
                impose_symmetry=impose_symmetry,
                symmetry_tolerance=symmetry_tolerance,
                traverse=self._traverse)
        else: # Use default symmetry_tolerance
            task = StructureOptimization(
                directory=directory,
                lattice_tolerance=self._lattice_tolerance,
                force_tolerance=self._force_tolerance,
                pressure_target=self._pressure_target,
                stress_tolerance=self._stress_tolerance,
                max_increase=self._max_increase,
                max_iteration=self._max_iteration,
                min_iteration=self._min_iteration,
                impose_symmetry=impose_symmetry,
                traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            incar=incar)
        task.set_job(job.copy("%s-%s" %
                              (job.get_jobname(), directory)))
        return task
    
class ElectronicStructure(TaskVasp, ElectronicStructureBase):
    """ """
    def __init__(self,
                 directory="electronic_structure",
                 name=None,
                 traverse=False):

        ElectronicStructureBase.__init__(
            self,
            directory=directory,
            name=name,
            traverse=traverse)

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._incar = None
        
    def _collect(self):
        """Collect information from output files of VASP.

        self._status of "done" or "terminate"  is stored.
        self._log: Terminate log is stored.

        """
        if not os.path.exists("vasprun.xml"):
            self._log += "vasprun.xml not exists.\n"
            self._status = "terminate"
        else:
            vxml = Vasprunxml("vasprun.xml")
            if vxml.parse_calculation():
                vxml.parse_eigenvalues()
                kpoints, weights = vxml.get_kpoints()
                self._properties = {'stress': vxml.get_stress(),
                                    'forces': vxml.get_forces(),
                                    'energies': vxml.get_energies()[:,1],
                                    'eigenvalues': vxml.get_eigenvalues(),
                                    'kpoints': kpoints,
                                    'kpoint-weights': weights}
                self._status = "done"
            else:
                self._log += "Failed to parse vasprun.xml.\n"
                self._status = "terminate"

class StructureOptimizationElement(TaskVasp,
                                   StructureOptimizationElementBase):
    def __init__(self,
                 directory="structure_optimization_element",
                 name=None,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-3,
                 pressure_target=0,
                 stress_tolerance=0.1, # GPa=kbar / 10
                 max_increase=1.5,
                 traverse=False):

        StructureOptimizationElementBase.__init__(
            self,
            directory=directory,
            name=name,
            lattice_tolerance=lattice_tolerance,
            force_tolerance=force_tolerance,
            pressure_target=pressure_target,
            stress_tolerance=stress_tolerance,
            max_increase=max_increase,
            traverse=traverse)
        
        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._incar = None

    def _collect(self):
        """Collect information from output files of VASP.

        self._current_cell: Final crystal structure of each relaxation steps.
        self._status: "next", "done", or "terminate" is stored.
        self._log: Terminate log is stored.

        """
        self._log = ""
        self._status = "terminate"

        if os.path.exists("POSCAR"):
            self._current_cell = read_poscar("POSCAR")
    
        if not os.path.exists("vasprun.xml"):
            self._log += "vasprun.xml not exists.\n"
        else:
            vxml = VasprunxmlExpat("vasprun.xml")
            is_success = vxml.parse()
            
            lattice = vxml.get_lattice()
            points = vxml.get_points()
            forces = vxml.get_forces()
            stress = vxml.get_stress()
            energies = vxml.get_energies()

            if is_success and len(lattice) > 0:
                self._stress = stress[-1] / 10
                self._forces = forces[-1]
                self._energy = energies[-1,1]
                self._judge(lattice[-1], points[-1])
            elif (not is_success) and len(lattice) > 2:
                self._log += "vasprun.xml is not cleanly closed.\n"
                self._stress = stress[-3] / 10
                self._forces = forces[-3]
                if len(energies) > 3:
                    self._energy = energies[-3,1]
                self._judge(lattice[-3], points[-3])
            else:
                self._log += "Failed to parse vasprun.xml.\n"
                self._current_cell = None

    def _judge(self, lattice_last, points):
        lattice_init = self._current_cell.get_lattice()
        vecs2_init = np.diag(np.dot(lattice_init.T, lattice_init))
        vecs2_last = np.diag(np.dot(lattice_last.T, lattice_last))
        d_vecs2_ratio = (vecs2_last - vecs2_init) / vecs2_init

        cell = Cell(lattice=lattice_last,
                    points=points,
                    symbols=self._current_cell.get_symbols())
        if os.path.exists("CONTCAR"):
            try:
                self._current_cell = read_poscar("CONTCAR")
            except:
                self._current_cell = cell
            else:
                # Sometimes CONTCAR's structure becomes original
                # structure same as POSCAR. In this case, the third
                # relaxed structure from the last is used for the new
                # cell.
                if (abs(self._current_cell.get_lattice() - lattice_init)
                    < 1e-12).all():
                    self._current_cell = cell
        else:
            self._current_cell = cell

        if np.linalg.det(lattice_last) > \
                self._max_increase * np.linalg.det(lattice_init):
            self._log += "Too large volume expansion.\n"
        else:
            if (abs(d_vecs2_ratio) > self._lattice_tolerance ** 2).any():
                self._log += "Lattice is not enough relaxed.\n"
                self._status = "next"
            if (abs(self._stress - np.eye(3) * self._pressure_target)
                > self._stress_tolerance).any():
                self._log += "Stress is not enough relaxed.\n"
                self._status = "next"
            if (abs(self._forces) > self._force_tolerance).any():
                self._log += "Forces are not enough relaxed.\n"
                self._status = "next"
            if not self._status == "next":
                self._status = "done"


class StructureOptimization(TaskVasp, StructureOptimizationBase):
    def __init__(self,
                 directory="structure_optimization",
                 name=None,
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
                 traverse=False):

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
            find_symmetry=find_symmetry,
            impose_symmetry=impose_symmetry,
            symmetry_tolerance=symmetry_tolerance,
            traverse=traverse)

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._incar = None

    def _get_next_task(self, cell):
        task = StructureOptimizationElement(
            directory="structopt-%d" % self._stage,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            traverse=self._traverse)
        
        task.set_configurations(
            cell=cell.copy(),
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            incar=self._incar.copy())

        task.set_job(self._job.copy(
                "%s-%s" % (self._job.get_jobname(), self._stage)))
        return task

class BulkModulus(TaskVasp, BulkModulusBase):
    """Task to calculate bulk modulus by VASP."""
    
    def __init__(self,
                 directory="bulk_modulus",
                 name=None,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-3,
                 pressure_target=0,
                 stress_tolerance=10,
                 max_increase=1.5,
                 max_iteration=3,
                 min_iteration=1,
                 traverse=False,
                 is_cell_relaxed=False):

        BulkModulusBase.__init__(
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
            traverse=traverse,
            is_cell_relaxed=is_cell_relaxed)

    def _get_plus_minus_tasks(self, cell):
        lattice = cell.get_lattice()

        cell_plus = cell.copy()
        cell_plus.set_lattice(lattice * 1.01 ** (1.0/3))
        plus = self._get_bm_task(cell_plus, "plus")

        cell_minus = cell.copy()
        cell_minus.set_lattice(lattice * 0.99 ** (1.0/3))
        minus = self._get_bm_task(cell_minus, "minus")

        return plus, minus

    def _get_bm_task(self, cell, directory):
        job, incar, k_mesh, k_shift, k_gamma, k_length = \
            self._choose_configuration(index=1)
        incar.set_isif(4)

        task=StructureOptimization(
            directory=directory,
            lattice_tolerance=self._lattice_tolerance,
            force_tolerance=self._force_tolerance,
            pressure_target=self._pressure_target,
            stress_tolerance=self._stress_tolerance,
            max_increase=self._max_increase,
            max_iteration=1,
            min_iteration=1,
            traverse=self._traverse)

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            incar=incar)
        task.set_job(job.copy("%s-%s" %
                              (job.get_jobname(), directory)))
        return task


class Phonon(TaskVasp, PhononBase):
    def __init__(self,
                 directory="phonon",
                 name=None,
                 supercell_matrix=np.eye(3, dtype=int),
                 primitive_matrix=np.eye(3, dtype=int),
                 distance=0.01,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-3,
                 pressure_target=0,
                 stress_tolerance=10,
                 max_increase=1.5,
                 max_iteration=3,
                 min_iteration=1,
                 traverse=False,
                 is_cell_relaxed=False):

        PhononBase.__init__(
            self,
            directory=directory,
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

    def _get_displacement_tasks(self):
        incar = self._incar[1].copy()

        # Perfect supercell
        supercell = self._phonon.get_supercell()
        tasks = [self._get_disp_task(atoms2cell(supercell),
                                     incar,
                                     "perfect")]
        # Displacements
        for i, disp in enumerate(
            self._phonon.get_supercells_with_displacements()):
            tasks.append(self._get_disp_task(atoms2cell(disp),
                                             incar,
                                             "disp-%03d" % (i+1)))
        return tasks

    def _get_disp_task(self, cell, incar, directory):
        job, incar, k_mesh, k_shift, k_gamma, k_length = \
            self._choose_configuration(index=1)

        if k_length:
            k_mesh = klength2mesh(k_length, cell.get_lattice())
            k_gamma = True
            k_shift = [0.0, 0.0, 0.0]
        
        # For Gamma-only VASP, take third element of self._job
        if isinstance(self._job, list) or isinstance(self._job, tuple):
            if len(self._job) > 2:
                if (np.array(k_mesh) == 1).all():
                    if k_shift:
                        if (np.abs(k_shift) < 1e-10).all():
                            job = self._job[2]
                    else:
                        job = self._job[2]

        task = ElectronicStructure(directory=directory,
                                   traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            incar=incar)
        task.set_job(job.copy("%s-%s" % (job.get_jobname(), directory)))

        return task
            
class ElasticConstants(TaskVasp, ElasticConstantsBase):
    def __init__(self,
                 directory="elastic_constants",
                 name=None,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-3,
                 pressure_target=0,
                 stress_tolerance=10,
                 max_increase=1.5,
                 max_iteration=3,
                 traverse=False,
                 is_cell_relaxed=False):
        
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
            traverse=traverse,
            is_cell_relaxed=is_cell_relaxed)
    
    def _get_ec_task(self, cell):
        job, incar, k_mesh, k_shift, k_gamma, k_length = \
            self._choose_configuration(index=1)
        incar.set_ibrion(6)
        incar.set_isif(3)
        incar.set_nsw(1)

        directory = "elastic_constants"
        task = ElasticConstantsElement(directory=directory,
                                       traverse=self._traverse)
        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=k_mesh,
            k_shift=k_shift,
            k_gamma=k_gamma,
            k_length=k_length,
            incar=incar)
        task.set_job(job.copy("%s-%s" %
                              (job.get_jobname(), directory)))
        return task
    
class ElasticConstantsElement(TaskVasp, ElasticConstantsElementBase):
    """ """
    def __init__(self,
                 directory="elastic_constants",
                 name=None,
                 traverse=False):

        ElasticConstantsElementBase.__init__(
            self,
            directory=directory,
            name=name,
            traverse=traverse)

        self._pseudo_potential_map = None
        self._k_mesh = None
        self._k_shift = None
        self._k_gamma = None
        self._incar = None
        
    def _collect(self):
        """Collect information from output files of VASP.

        self._status of "done" or "terminate"  is stored.
        self._log: Terminate log is stored.

        """
        if not os.path.exists("OUTCAR"):
            self._log += "OUTCAR not exists.\n"
            self._status = "terminate"
        else:
            outcar = Outcar("OUTCAR")
            elastic_constants = outcar.get_elastic_constants()
            if outcar.parse_elastic_constants():
                self._elastic_constants = outcar.get_elastic_constants()
                self._status = "done"
            else:
                self._log += "Failed to parse OUTCAR.\n"
                self._status = "terminate"

class PhononRelax(TaskVasp, PhononRelaxBase):
    def __init__(self,
                 directory="phonon_relax",
                 name=None,
                 ancestral_cells={},
                 distance=0.01,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-5,
                 pressure_target=0,
                 stress_tolerance=1,
                 max_increase=1.5,
                 max_iteration=10,
                 min_iteration=5,
                 symmetry_tolerance=0.1,
                 restrict_offspring=False,
                 max_offspring=None,
                 cutoff_eigenvalue=-0.02,
                 max_displacement=None,
                 num_sampling_points=60,
                 traverse=False):

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
            traverse=traverse)

    def _get_phonon_relax_element_task(self, cell):
        task = PhononRelaxElement(directory="phonon_relax_element",
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
                                  traverse=self._traverse)

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            incar=self._incar)
        task.set_job(self._job)

        return task

    def _get_phonon_relax_task(self, cell, ancestral_cells, directory):
        task = PhononRelax(directory=directory,
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
                           traverse=self._traverse)

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            incar=self._incar)
        task.set_job(self._job)

        return task


class PhononRelaxElement(TaskVasp, PhononRelaxElementBase):
    def __init__(self,
                 directory="phonon_relax_element",
                 name=None,
                 ancestral_cells={},
                 tid_parent=None,
                 distance=0.01,
                 lattice_tolerance=0.1,
                 force_tolerance=1e-5,
                 pressure_target=0,
                 stress_tolerance=1,
                 max_increase=1.5,
                 max_iteration=10,
                 min_iteration=5,
                 symmetry_tolerance=0.1,
                 cutoff_eigenvalue=-0.02,
                 max_displacement=None,
                 num_sampling_points=60,
                 traverse=False):

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
            traverse=traverse)

    def _get_phonon_task(self, cell, supercell_matrix, directory):
        task = Phonon(directory=directory,
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
                      traverse=self._traverse)

        task.set_configurations(
            cell=cell,
            pseudo_potential_map=self._pseudo_potential_map,
            k_mesh=self._k_mesh,
            k_shift=self._k_shift,
            k_gamma=self._k_gamma,
            k_length=self._k_length,
            incar=self._incar[1:3])
        task.set_job(self._job)

        return task

