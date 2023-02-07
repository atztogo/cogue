import numbers
import os
import sys
import xml.etree.cElementTree as etree

import numpy as np
from phonopy.interface.vasp import VasprunxmlExpat as PhonopyVasprunExpat

from cogue.crystal.atom import atomic_symbols, atomic_weights
from cogue.crystal.cell import Cell


class VaspCell(Cell):
    def __init__(self, cell, is_vasp4=False, comment=None):
        self._cell = cell
        self._is_vasp4 = is_vasp4

        self._compressed_symbols = None
        self._set_compressed_symbols()
        self._atom_order = None
        self._set_atom_order()
        self._set_vasp_cell()

        self._comment = None
        self._set_comment(comment)

        self._poscar_lines = None
        self.create_poscar_lines()

        self._poscar_yaml_lines = None
        self.create_poscar_yaml_lines()

    def write(self, filename="POSCAR"):
        with open(filename, "w") as w:
            for line in self._poscar_lines[:-1]:
                w.write(line + "\n")
            w.write(self._poscar_lines[-1])

    def write_yaml(self, filename="POSCAR.yaml"):
        with open(filename, "w") as w:
            for line in self._poscar_yaml_lines[:-1]:
                w.write(line + "\n")
            w.write(self._poscar_yaml_lines[-1])

    def get_poscar_lines(self):
        return self._poscar_lines

    def get_poscar_yaml_lines(self):
        return self._poscar_yaml_lines

    def get_atom_order(self):
        return self._atom_order

    def get_comment(self):
        return self._comment

    def get_compressed_symbols(self):
        return self._compressed_symbols

    def create_poscar_lines(self):
        lines = []
        lines.append(self._comment)
        lines.append("   1.0")

        for v in self._lattice.T:
            lines.append(" %22.16f%22.16f%22.16f" % tuple(v))

        if not self._is_vasp4:
            lines.append(" " + " ".join(self._compressed_symbols))

        lines.append(
            " "
            + " ".join(
                ["%3d" % self._symbols.count(s) for s in self._compressed_symbols]
            )
        )

        lines.append("Direct")

        for v in self._points.T:
            v16 = (v - np.rint(v)).round(decimals=16)
            for i in range(3):
                if v16[i] < 0:
                    v16[i] += 1
            lines.append(" %20.16f%20.16f%20.16f" % tuple(v16))

        self._poscar_lines = lines

    def create_poscar_yaml_lines(self):
        lines = self._cell.get_yaml_lines()
        lines.append("poscar_order:")
        for i in self._atom_order:
            lines.append("- %d" % (i + 1))

        self._poscar_yaml_lines = lines

    def _set_compressed_symbols(self):
        symbols = self._cell.get_symbols()
        self._compressed_symbols = []
        for s in symbols:
            if s not in self._compressed_symbols:
                self._compressed_symbols.append(s)

    def _set_atom_order(self):
        symbols = self._cell.get_symbols()
        self._atom_order = []
        for cs in self._compressed_symbols:
            for i, s in enumerate(symbols):
                if s == cs:
                    self._atom_order.append(i)

    def _set_vasp_cell(self):
        symbols = [self._cell.get_symbols()[i] for i in self._atom_order]
        if self._cell.get_magnetic_moments() is not None:
            magmoms = self._cell.get_magnetic_moments()[self._atom_order]
        else:
            magmoms = None

        Cell.__init__(
            self,
            lattice=self._cell.lattice,
            points=(self._cell.get_points().T)[self._atom_order].T,
            symbols=symbols,
            magmoms=magmoms,
            masses=self._cell.get_masses()[self._atom_order],
        )

    def _set_comment(self, comment):
        if self._is_vasp4 or comment is None:
            self._comment = " " + " ".join(self._compressed_symbols)
        else:
            self._comment = comment.strip()


def parse_poscar(lines):
    isvasp5 = False
    symbols = []

    for w in lines[0].split():
        if w.strip() in atomic_symbols:
            symbols.append(w)

    for w in lines[5].split():
        if not w.isdigit():
            isvasp5 = True
            break

    if isvasp5:
        symbols_vasp5 = []
        for w in lines[5].split():
            if w.strip() in atomic_symbols:
                symbols_vasp5.append(w)
        lines.pop(5)

    # Read lines in VASP 4 style
    scale = float(lines[1])
    lattice = np.zeros((3, 3), dtype=float)
    lattice[:, 0] = [float(x) for x in lines[2].split()]
    lattice[:, 1] = [float(x) for x in lines[3].split()]
    lattice[:, 2] = [float(x) for x in lines[4].split()]

    lattice *= scale
    num_atoms = np.array([int(x) for x in lines[5].split()])

    # Check if symbols in VASP 5 style is valid
    if isvasp5:
        if len(num_atoms) == len(symbols_vasp5):
            symbols = symbols_vasp5

    # Assume Direct
    if lines[6][0] in "CcKk":
        print("Cartesian is not supported.")
        raise RuntimeError

    points = np.zeros((3, num_atoms.sum()), dtype=float)
    for i in range(num_atoms.sum()):
        points[:, i] = [float(x) for x in lines[i + 7].split()[:3]]

    # Expand symbols
    symbols_expanded = []
    if len(symbols) != len(num_atoms):
        sys.stderr.write("Chemical symbols are not specified correctly.\n")
        symbols = [atomic_weights[i + 1][0] for i in range(len(num_atoms))]

    for i, n in enumerate(num_atoms):
        symbols_expanded += [symbols[i]] * n

    return Cell(lattice=lattice, points=points, symbols=symbols_expanded)


def read_poscar(filename="POSCAR"):
    with open(filename) as f:
        cell = parse_poscar(f.readlines())
        return cell
    return None


def read_poscar_yaml(filename="POSCAR.yaml"):
    try:
        import yaml
    except ImportError:
        print("You need to install python-yaml.")
        sys.exit(1)

    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader

    with open(filename) as f:
        data = yaml.load(f, Loader=Loader)
        lattice = np.transpose(data["lattice"])
        points = np.transpose([x["coordinates"] for x in data["points"]])
        symbols = []
        masses = []
        for point in data["points"]:
            if "mass" in point:
                masses.append(point["mass"])
            if "symbol" in point:
                symbols.append(point["symbol"])
        if len(masses) != len(data["points"]):
            masses = None
        if len(symbols) != len(data["points"]):
            symbols = None

        poscar_order = data["poscar_order"]

        return (
            Cell(lattice=lattice, points=points, symbols=symbols, masses=masses),
            poscar_order,
        )


def change_point_order(cell, atom_order):
    symbols = [cell.get_symbols()[i] for i in atom_order]
    if cell.get_magnetic_moments() is not None:
        magmoms = cell.get_magnetic_moments()[atom_order]
    else:
        magmoms = None

    return Cell(
        lattice=cell.lattice,
        points=(cell.get_points().T)[atom_order].T,
        symbols=symbols,
        magmoms=magmoms,
        masses=cell.get_masses()[atom_order],
    )


def get_atom_order_from_poscar_yaml(filename):
    try:
        import yaml
    except ImportError:
        print("You need to install python-yaml.")
        sys.exit(1)

    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader

    with open(filename) as f:
        data = yaml.load(f, Loader=Loader)
        poscar_order = data["poscar_order"]
        inverse_order = {(j - 1): i for i, j in enumerate(poscar_order)}
        return [inverse_order[i] for i in range(len(inverse_order))]


def write_poscar(cell, filename=None, is_vasp4=False, comment=None):
    vasp_cell = VaspCell(cell, is_vasp4=is_vasp4, comment=comment)

    if filename is None:
        poscar_lines = vasp_cell.get_poscar_lines()
        return "\n".join(poscar_lines)
    else:
        vasp_cell.write(filename=filename)


def write_poscar_yaml(cell, filename=None):
    vasp_cell = VaspCell(cell)

    if filename is None:
        poscar_yaml_lines = vasp_cell.get_poscar_yaml_lines()
        return "\n".join(poscar_yaml_lines)
    else:
        vasp_cell.write_yaml(filename=filename)


def write_potcar(names, filename="POTCAR"):
    if "COGUE_POTCAR_PATH" in os.environ:
        potcarpath = os.environ["COGUE_POTCAR_PATH"]
    else:
        print("COGUE_POTCAR_PATH is not set correctly.")
        return False

    with open(filename, "w") as w:
        for i, s in enumerate(names):
            if i == 0 or not s == names[i - 1]:
                with open("%s/%s" % (potcarpath, s)) as f:
                    for line in f:
                        w.write(line)


def get_enmax_from_potcar(names):
    if "COGUE_POTCAR_PATH" in os.environ:
        potcarpath = os.environ["COGUE_POTCAR_PATH"]
    else:
        print("COGUE_POTCAR_PATH is not set correctly.")
        return False

    enmax = []
    for i, s in enumerate(names):
        if i == 0 or not s == names[i - 1]:
            with open("%s/%s" % (potcarpath, s)) as f:
                for line in f:
                    if "ENMAX" in line:
                        enmax.append(float(line[11:20]))
    return enmax


class Incar:
    def __init__(
        self,
        addgrid=None,
        aggac=None,
        algo=None,
        ediff=None,
        ediffg=None,
        emax=None,
        emin=None,
        encut=None,
        gga=None,
        ialgo=None,
        ibrion=None,
        icharg=None,
        isif=None,
        ismear=None,
        ispin=None,
        isym=None,
        ivdw=None,
        kpar=None,
        lcharg=None,
        lepsilon=None,
        lorbit=None,
        lreal=None,
        luse_vdw=None,
        lwave=None,
        magmom=None,
        nbands=None,
        nedos=None,
        nelm=None,
        nelmin=None,
        npar=None,
        nsw=None,
        prec=None,
        pstress=None,
        sigma=None,
        symprec=None,
    ):
        self._tagnames = {
            "addgrid": "ADDGRID",
            "aggac": "AGGAC",
            "algo": "ALGO",
            "ediff": "EDIFF",
            "ediffg": "EDIFFG",
            "emax": "EMAX",
            "emin": "EMIN",
            "encut": "ENCUT",
            "gga": "GGA",
            "ialgo": "IALGO",
            "ibrion": "IBRION",
            "icharg": "ICHARG",
            "isif": "ISIF",
            "ismear": "ISMEAR",
            "ispin": "ISPIN",
            "isym": "ISYM",
            "ivdw": "IVDW",
            "kpar": "KPAR",
            "lcharg": "LCHARG",
            "lepsilon": "LEPSILON",
            "lorbit": "LORBIT",
            "lreal": "LREAL",
            "luse_vdw": "LUSE_VDW",
            "lwave": "LWAVE",
            "magmom": "MAGMOM",
            "nbands": "NBANDS",
            "nedos": "NEDOS",
            "nelm": "NELM",
            "nelmin": "NELMIN",
            "npar": "NPAR",
            "nsw": "NSW",
            "prec": "PREC",
            "pstress": "PSTRESS",
            "sigma": "SIGMA",
            "symprec": "SYMPREC",
        }

        self._tagvals = {
            "addgrid": addgrid,
            "aggac": aggac,
            "algo": algo,
            "ediff": ediff,
            "ediffg": ediffg,
            "emax": emax,
            "emin": emin,
            "encut": encut,
            "gga": gga,
            "ialgo": ialgo,
            "ibrion": ibrion,
            "icharg": icharg,
            "isif": isif,
            "ismear": ismear,
            "ispin": ispin,
            "isym": isym,
            "ivdw": ivdw,
            "kpar": kpar,
            "lcharg": lcharg,
            "lepsilon": lepsilon,
            "lorbit": lorbit,
            "lreal": lreal,
            "luse_vdw": luse_vdw,
            "lwave": lwave,
            "magmom": magmom,
            "nbands": nbands,
            "nedos": nedos,
            "nelm": nelm,
            "nelmin": nelmin,
            "npar": npar,
            "nsw": nsw,
            "prec": prec,
            "pstress": pstress,
            "sigma": sigma,
            "symprec": symprec,
        }

        self._tagorder = [
            "prec",
            "gga",
            "ivdw",
            "luse_vdw",
            "aggac",
            "ibrion",
            "nsw",
            "algo",
            "nelm",
            "nelmin",
            "isif",
            "encut",
            "ediff",
            "ediffg",
            "icharg",
            "ispin",
            "lorbit",
            "magmom",
            "ismear",
            "sigma",
            "nbands",
            "emin",
            "emax",
            "nedos",
            "pstress",
            "ialgo",
            "lreal",
            "addgrid",
            "lwave",
            "lcharg",
            "lepsilon",
            "npar",
            "kpar",
            "isym",
            "symprec",
        ]

    def clear(self):
        for k in self._tagvals.keys():
            self._tagvals[k] = None

    def set_tag(self, k, v):
        if k in self._tagvals:
            self._tagvals[k] = v
        else:
            print("Key %s is not available." % k)

    def set_addgrid(self, x):
        self._tagvals["addgrid"] = x

    def get_addgrid(self):
        return self._tagvals["addgrid"]

    def set_aggac(self, x):
        self._tagvals["aggac"] = x

    def get_aggac(self):
        return self._tagvals["aggac"]

    def set_algo(self, x):
        if "ialgo" in self._tagvals:
            self._tagvals.pop("ialgo")
        self._tagvals["algo"] = x

    def get_algo(self):
        return self._tagvals["algo"]

    def set_ediff(self, x):
        self._tagvals["ediff"] = x

    def get_ediff(self):
        return self._tagvals["ediff"]

    def set_ediffg(self, x):
        self._tagvals["ediffg"] = x

    def get_ediffg(self):
        return self._tagvals["ediffg"]

    def set_emax(self, x):
        self._tagvals["emax"] = x

    def get_emax(self):
        return self._tagvals["emax"]

    def set_emin(self, x):
        self._tagvals["emin"] = x

    def get_emin(self):
        return self._tagvals["emin"]

    def set_encut(self, x):
        self._tagvals["encut"] = x

    def get_encut(self):
        return self._tagvals["encut"]

    def set_gga(self, x):
        self._tagvals["gga"] = x

    def get_gga(self):
        return self._tagvals["gga"]

    def set_ialgo(self, x):
        if "algo" in self._tagvals:
            self._tagvals.pop("algo")
        self._tagvals["ialgo"] = x

    def get_ialgo(self):
        return self._tagvals["ialgo"]

    def set_ibrion(self, x):
        self._tagvals["ibrion"] = x

    def get_ibrion(self):
        return self._tagvals["ibrion"]

    def set_icharg(self, x):
        self._tagvals["icharg"] = x

    def get_icharg(self):
        return self._tagvals["icharg"]

    def set_isif(self, x):
        self._tagvals["isif"] = x

    def get_isif(self):
        return self._tagvals["isif"]

    def set_ismear(self, x):
        self._tagvals["ismear"] = x

    def get_ismear(self):
        return self._tagvals["ismear"]

    def set_ispin(self, x):
        self._tagvals["ispin"] = x

    def get_ispin(self):
        return self._tagvals["ispin"]

    def set_isym(self, x):
        self._tagvals["isym"] = x

    def get_isym(self):
        return self._tagvals["isym"]

    def set_ivdw(self, x):
        self._tagvals["ivdw"] = x

    def get_ivdw(self):
        return self._tagvals["ivdw"]

    def set_kpar(self, x):
        self._tagvals["kpar"] = x

    def get_kpar(self):
        return self._tagvals["kpar"]

    def set_lcharg(self, x):
        self._tagvals["lcharg"] = x

    def get_lcharg(self):
        return self._tagvals["lcharg"]

    def set_lepsilon(self, x):
        self._tagvals["lepsilon"] = x

    def get_lepsilon(self):
        return self._tagvals["lepsilon"]

    def set_lorbit(self, x):
        self._tagvals["lorbit"] = x

    def get_lorbit(self):
        return self._tagvals["lorbit"]

    def set_lreal(self, x):
        self._tagvals["lreal"] = x

    def get_lreal(self):
        return self._tagvals["lreal"]

    def set_luse_vdw(self, x):
        self._tagvals["luse_vdw"] = x

    def get_luse_vdw(self):
        return self._tagvals["luse_vdw"]

    def set_lwave(self, x):
        self._tagvals["lwave"] = x

    def get_lwave(self):
        return self._tagvals["lwave"]

    def set_magmom(self, x):
        self._tagvals["magmom"] = x

    def get_magmom(self):
        return self._tagvals["magmom"]

    def set_nbands(self, x):
        self._tagvals["nbands"] = x

    def get_nbands(self):
        return self._tagvals["nbands"]

    def set_nedos(self, x):
        self._tagvals["nedos"] = x

    def get_nedos(self):
        return self._tagvals["nedos"]

    def set_nelm(self, x):
        self._tagvals["nelm"] = x

    def get_nelm(self):
        return self._tagvals["nelm"]

    def set_nelmin(self, x):
        self._tagvals["nelmin"] = x

    def get_nelmin(self):
        return self._tagvals["nelmin"]

    def set_npar(self, x):
        self._tagvals["npar"] = x

    def get_npar(self):
        return self._tagvals["npar"]

    def set_nsw(self, x):
        self._tagvals["nsw"] = x

    def get_nsw(self):
        return self._tagvals["nsw"]

    def set_prec(self, x):
        self._tagvals["prec"] = x

    def get_prec(self):
        return self._tagvals["prec"]

    def set_pstress(self, x):
        self._tagvals["pstress"] = x

    def get_pstress(self):
        return self._tagvals["pstress"]

    def set_sigma(self, x):
        self._tagvals["sigma"] = x

    def get_sigma(self):
        return self._tagvals["sigma"]

    def set_symprec(self, x):
        self._tagvals["symprec"] = x

    def get_symprec(self):
        return self._tagvals["symprec"]

    def set_electronic_structure(self):
        tags = self._tagvals
        tags["prec"] = "Accurate"
        tags["ibrion"] = -1
        tags["nelmin"] = 5
        tags["encut"] = 500
        tags["ediff"] = 1.0e-08
        tags["ismear"] = 0
        tags["sigma"] = 0.01
        tags["ialgo"] = 38
        tags["lreal"] = False
        tags["addgrid"] = True
        tags["lwave"] = False
        tags["lcharg"] = False

    def set_structure_optimization(self):
        tags = self._tagvals
        tags["prec"] = "Accurate"
        tags["ibrion"] = 2
        tags["nsw"] = 10
        tags["nelmin"] = 5
        tags["isif"] = 3
        tags["encut"] = 500
        tags["ediff"] = 1.0e-08
        tags["ediffg"] = -1.0e-08
        tags["ismear"] = 0
        tags["sigma"] = 0.01
        tags["ialgo"] = 38
        tags["lreal"] = False
        tags["addgrid"] = True
        tags["lwave"] = False
        tags["lcharg"] = False

    def copy(self):
        incar = Incar()
        for k in self._tagvals:
            incar.set_tag(k, self._tagvals[k])
        return incar

    def write(self, filename="INCAR"):
        names = self._tagnames

        with open(filename, "w") as w:
            for k in self._tagorder:
                v = self._tagvals[k]
                if isinstance(v, bool):
                    if v:
                        w.write("%10s = .TRUE.\n" % (names[k]))
                    else:
                        w.write("%10s = .FALSE.\n" % (names[k]))
                elif isinstance(v, int):
                    w.write("%10s = %d\n" % (names[k], v))
                elif isinstance(v, float):
                    if v < 1:
                        w.write("%10s = %e\n" % (names[k], v))
                    else:
                        w.write("%10s = %f\n" % (names[k], v))
                elif isinstance(v, str):
                    w.write("%10s = %s\n" % (names[k], v))


def write_kpoints(
    filename="KPOINTS", mesh=None, shift=None, gamma=False, length=None, kpoint=None
):
    with open(filename, "w") as w:
        if length:
            w.write("Automatic mesh\n")
            w.write("0\n")
            w.write("Auto\n")
            w.write("%4d\n" % length)
        elif kpoint is not None:
            if isinstance(kpoint[0], numbers.Number):
                w.write("Explicit k-point\n")
                w.write("1\n")
                w.write("Reciprocal\n")
                w.write("%10.7f %10.7f %10.7f  1\n" % tuple(kpoint))
            else:
                w.write("Automatic mesh\n")
                w.write("0\n")
                w.write("Reciprocal\n")
                for k in kpoint:
                    w.write("%10.7f %10.7f %10.7f\n" % tuple(k))
        elif mesh is not None:
            w.write("Automatic mesh\n")
            w.write("0\n")
            if gamma:
                w.write("Gamma\n")
            else:
                w.write("Monkhorst-pack\n")
            w.write(" %5d %5d %5d\n" % tuple(mesh))
            if shift is None:
                w.write("     0.    0.    0.\n")
            else:
                w.write(" %5.3f %5.3f %5.3f\n" % tuple(shift))


class Outcar:
    def __init__(self, filename="OUTCAR"):
        self._filename = filename
        self._elastic_constants = None

    def get_elastic_constants(self):
        return self._elastic_constants

    def parse_elastic_constants(self):
        with open(self._filename) as outcar:
            hooked = False
            for line in outcar:
                if line.strip() == "TOTAL ELASTIC MODULI (kBar)":
                    hooked = True
                    break

            if hooked:
                next(outcar)
                next(outcar)
                ec = []
                for i in range(6):
                    pos = 8
                    line = next(outcar)
                    for j in range(6):
                        try:
                            elem = float(line[pos : (pos + 12)])
                        except ValueError:
                            return False

                        ec.append(elem)
                        pos += 12

                self._elastic_constants = np.array(
                    np.reshape(ec, (6, 6)), dtype="double", order="C"
                )
                return True
            else:
                return False


class Vasprunxml(object):
    def __init__(self, filename="vasprun.xml"):
        self._filename = filename
        self._forces = None
        self._stress = None
        self._lattice = None
        self._points = None
        self._energies = None
        self._eigenvalues_spin1 = None
        self._eigenvalues_spin2 = None
        self._occupancies_spin1 = None
        self._occupancies_spin2 = None
        self._kpoints = None
        self._kpoint_weights = None
        self._born_charges = None
        self._epsilon = None
        self._nbands = None
        self._efermi = None

        self._log = ""

    def get_forces(self):
        return self._forces

    def get_stress(self):
        return self._stress

    def get_lattice(self):
        return self._lattice

    def get_points(self):
        return self._points

    def get_energies(self):
        return self._energies

    def get_eigenvalues(self):
        # Degeneracy of electrons
        # Spin components 1 and 2 are stored in tuple as
        # (spin1, spin2) or (spin1,) if no spin polarized.
        if self._eigenvalues_spin2 is None:
            return (self._eigenvalues_spin1,)
        else:
            return (self._eigenvalues_spin1, self._eigenvalues_spin2)

    def get_occupancies(self):
        # Degeneracy of electrons
        # Spin components 1 and 2 are stored in tuple as
        # (spin1, spin2) or (spin1,) if no spin polarized.
        if self._occupancies_spin2 is None:
            return (self._occupancies_spin1,)
        else:
            return (self._occupancies_spin1, self._occupancies_spin2)

    def get_kpoints(self):
        return self._kpoints, self._kpoint_weights

    def get_born_charges(self):
        return self._born_charges

    def get_epsilon(self):
        return self._epsilon

    def get_efermi(self):
        return self._efermi

    def get_nbands(self):
        return self._nbands

    @property
    def log(self):
        log = self._log
        self._log = ""
        return log

    def parse_calculation(self):
        forces = []
        stress = []
        lattice = []
        points = []
        energies = []
        born_charges = []
        epsilon = []

        try:
            for event, element in etree.iterparse(self._filename):
                if element.tag != "calculation":
                    continue

                for varray in element.findall("./varray"):
                    self._parse_forces_and_stress(varray, forces, stress)

                for varray in element.findall("./structure/varray"):
                    self._parse_points(varray, points)

                for varray in element.findall("./structure/crystal/varray"):
                    self._parse_lattice(varray, lattice)

                for energy in element.findall("./energy"):
                    self._parse_energies(energy, energies)

                for array in element.findall("./array"):
                    if array.attrib["name"] == "born_charges":
                        self._parse_born_charges(array, born_charges)

                for varray in element.findall("./varray"):
                    if varray.attrib["name"] == "epsilon":
                        self._parse_vectors(varray, epsilon)

            self._forces = np.array(forces)
            self._stress = np.array(stress)
            self._lattice = np.array(lattice)
            self._points = np.array(points)
            self._energies = np.array(energies)
            if born_charges:
                self._born_charges = np.array(born_charges)
            if epsilon:
                self._epsilon = np.array(epsilon)

            return True

        except:  # noqa E722
            self._log += "    [Vasprunxml] Failed parse_calculation\n"
            return False

    def parse_parameters(self):
        try:
            for event, element in etree.iterparse(self._filename):
                if element.tag != "parameters":
                    continue

                for separator in element.findall("./separator"):
                    if separator.attrib["name"] == "electronic":
                        for i in separator.findall("./i"):
                            if i.attrib["name"] == "NBANDS":
                                self._nbands = int(i.text)
            return True
        except:  # noqa E722
            self._log += "    [Vasprunxml] Failed parse_parameters\n"
            return False

    def parse_efermi(self):
        try:
            for event, element in etree.iterparse(self._filename):
                if element.tag != "dos":
                    continue

                for i in element.findall("./i"):
                    if i.attrib["name"] == "efermi":
                        efermi = float(i.text)

            self._efermi = efermi
            return True
        except:  # noqa E722
            self._log += "    [Vasprunxml] Failed parse_efermi\n"
            return False

    def parse_eigenvalues(self):
        spin1 = []
        spin2 = []
        occ1 = []
        occ2 = []
        try:
            for event, element in etree.iterparse(self._filename):
                if element.tag != "eigenvalues":
                    continue

                for array in element.findall("./array/set/set"):
                    if array.attrib["comment"] == "spin 1":
                        self._parse_eigenvalues_spin(array, spin1, occ1)

                    if array.attrib["comment"] == "spin 2":
                        self._parse_eigenvalues_spin(array, spin2, occ2)

            if spin1:
                self._eigenvalues_spin1 = np.array(spin1)
                self._occupancies_spin1 = np.array(occ1)
            if spin2:
                self._eigenvalues_spin2 = np.array(spin2)
                self._occupancies_spin2 = np.array(occ2)

            return self._parse_kpoints()

        except:  # noqa E722
            self._log += "    [Vasprunxml] Failed parse_eigenvalues\n"
            return False

    def _parse_kpoints(self):
        try:
            for event, element in etree.iterparse(self._filename):
                if element.tag != "kpoints":
                    continue

                kpoints = []
                weights = []
                for varray in element.findall("./varray"):
                    if varray.attrib["name"] == "kpointlist":
                        for v in varray.findall("./v"):
                            kpoints.append([float(x) for x in v.text.split()])

                    if varray.attrib["name"] == "weights":
                        for v in varray.findall("./v"):
                            weights.append(float(v.text))

            self._kpoints = np.array(kpoints)
            self._kpoint_weights = np.array(weights)
            return True

        except:  # noqa E722
            self._log += "    [Vasprunxml] Failed parse_kpoints\n"
            return False

    def _parse_eigenvalues_spin(self, array, eigenvals, occupancies):
        for kset in array.findall("./set"):
            eigs = []
            occs = []
            for r in kset.findall("./r"):
                vals = r.text.split()
                eigs.append(float(vals[0]))
                occs.append(float(vals[1]))
            eigenvals.append(eigs)
            occupancies.append(occs)

    def _parse_forces_and_stress(self, varray, forces, stress):
        # force
        if varray.attrib["name"] == "forces":
            forces_geomopt = []
            for v in varray.findall("./v"):
                forces_geomopt.append([float(x) for x in v.text.strip().split()])
            forces.append(forces_geomopt)

        # stress
        if varray.attrib["name"] == "stress":
            stress_geomopt = []
            for v in varray.findall("./v"):
                stress_geomopt.append([float(x) for x in v.text.strip().split()])
            stress.append(stress_geomopt)

    def _parse_points(self, varray, points):
        # points
        if varray.attrib["name"] == "positions":
            points_geomopt = []
            for v in varray.findall("./v"):
                points_geomopt.append([float(x) for x in v.text.strip().split()])
            points.append(np.transpose(points_geomopt))

    def _parse_lattice(self, varray, lattice):
        if varray.attrib["name"] == "basis":
            lattice_geomopt = []
            for v in varray.findall("./v"):
                lattice_geomopt.append(
                    np.transpose([float(x) for x in v.text.strip().split()])
                )
            lattice.append(np.transpose(lattice_geomopt))

    def _parse_energies(self, energy, energies):
        energies_geomopt = []
        for v in energy.findall("./i"):
            energies_geomopt.append([float(x) for x in v.text.strip().split()])
        energies.append(energies_geomopt)

    def _parse_born_charges(self, array, born_charges):
        for ion_set in array.findall("./set"):
            tensor = []
            for v in ion_set.findall("./v"):
                tensor.append([float(x) for x in v.text.strip().split()])
            born_charges.append(tensor)

    def _parse_vectors(self, varray, vectors):
        for v in varray.findall("./v"):
            vectors.append([float(x) for x in v.text.strip().split()])


class VasprunxmlExpat(PhonopyVasprunExpat):
    def __init__(self, fileptr):
        """Parsing vasprun.xml by Expat

        Args:
           fileptr: binary stream. Considering compatibility between python2.7
               and 3.x, it is prepared as follows:

               import io
               io.open(filename, "rb")

        """

        PhonopyVasprunExpat.__init__(self, fileptr)
        self._log = ""

    @property
    def log(self):
        log = self._log
        self._log = ""
        return log

    def get_cells(self):
        cells = []
        if len(self._all_points) == len(self._all_lattice):
            for p, l in zip(self._all_points, self._all_lattice):
                cells.append(Cell(lattice=l.T, points=p.T, symbols=self._symbols))
        return cells
