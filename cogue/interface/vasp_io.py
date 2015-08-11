import os
import sys
import numpy as np
import xml.parsers.expat
import xml.etree.cElementTree as etree
from cogue.crystal.atom import atomic_symbols, atomic_weights
from cogue.crystal.cell import Cell
import StringIO

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
    if lines[6][0] in 'CcKk':
        print "Cartesian is not supported"
        raise

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
    
    return Cell(lattice=lattice,
                points=points,
                symbols=symbols_expanded)
    

def read_poscar(filename="POSCAR"):
    return parse_poscar(open(filename).readlines())

def write_poscar(cell, filename=None, vasp4=False, comment=None):
    if filename is None:
        w = StringIO.StringIO()
    else:
        w = open(filename, 'w')
    
    symbols = cell.get_symbols()
    symbols_compressed = []
    for s in symbols:
        if not s in symbols_compressed:
            symbols_compressed.append(s)
    
    if vasp4 or comment is None:
        for s in symbols_compressed:
            w.write(" %s" % s)
        w.write("\n")
    else:
        w.write(comment.strip() + "\n")

    w.write("   1.0\n")

    for v in cell.get_lattice().T:
        w.write(" %22.16f%22.16f%22.16f\n" % tuple(v))

    if not vasp4:
        for s in symbols_compressed:
            w.write(" %s" % s)
        w.write("\n")
    
    for s in symbols_compressed:
        w.write(" %3d" % symbols.count(s))
    w.write("\n")

    w.write("Direct\n")

    points = cell.get_points().T
    count_atoms = 0
    num_atoms = len(points)
    for sc in symbols_compressed:
        for s, v in zip(symbols, points):
            if s == sc:
                v16 = v.round(decimals=16)
                for i in range(3):
                    if not v16[i] < 1:
                        v16[i] -= 1
                w.write(" %20.16f%20.16f%20.16f" % tuple(v16))
                count_atoms += 1
                if count_atoms < num_atoms:
                    w.write("\n")
                
                
    
    if filename is None:
        poscar_string = w.getvalue()
        w.close()
        return poscar_string
    else:
        w.close()

def write_potcar(names, filename="POTCAR"):
    if 'COGUE_POTCAR_PATH' in os.environ:
        potcarpath = os.environ['COGUE_POTCAR_PATH']
    else:
        print "COGUE_POTCAR_PATH is not set correctly."
        return False
    
    w = open(filename, 'w')
    for i, s in enumerate(names):
        if i == 0 or not s == names[i - 1]:
            for line in open("%s/%s" % (potcarpath, s)):
                w.write(line)
    w.close()

def get_enmax_from_potcar(names):
    if 'COGUE_POTCAR_PATH' in os.environ:
        potcarpath = os.environ['COGUE_POTCAR_PATH']
    else:
        print "COGUE_POTCAR_PATH is not set correctly."
        return False
    
    enmax = []
    for i, s in enumerate(names):
        if i == 0 or not s == names[i - 1]:
            for line in open("%s/%s" % (potcarpath, s)):
                if 'ENMAX' in line:
                    enmax.append(float(line[11:20]))
    return enmax

class Incar:
    def __init__(self,
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
                 ivdw=None,
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
                 symprec=None):

        self._tagnames = {
            'addgrid' : "ADDGRID",
            'aggac'   : "AGGAC",
            'algo'    : "ALGO",
            'ediff'   : "EDIFF",
            'ediffg'  : "EDIFFG",
            'emax'    : "EMAX",
            'emin'    : "EMIN",
            'encut'   : "ENCUT",
            'gga'     : "GGA",
            'ialgo'   : "IALGO",
            'ibrion'  : "IBRION",
            'icharg'  : "ICHARG",
            'isif'    : "ISIF",
            'ismear'  : "ISMEAR",
            'ispin'   : "ISPIN",
            'ivdw'    : "IVDW",
            'lcharg'  : "LCHARG",
            'lepsilon': "LEPSILON",
            'lorbit'  : "LORBIT",
            'lreal'   : "LREAL",
            'luse_vdw': "LUSE_VDW",
            'lwave'   : "LWAVE",
            'magmom'  : "MAGMOM",
            'nbands'  : "NBANDS",
            'nedos'   : "NEDOS",
            'nelm'    : "NELM",
            'nelmin'  : "NELMIN",
            'npar'    : "NPAR",
            'nsw'     : "NSW",
            'prec'    : "PREC",
            'pstress' : "PSTRESS",
            'sigma'   : "SIGMA",
            'symprec' : "SYMPREC"}

        self._tagvals = {
            'addgrid' : addgrid,
            'aggac'   : aggac,
            'algo'    : algo,
            'ediff'   : ediff, 
            'ediffg'  : ediffg,
            'emax'    : emax,
            'emin'    : emin,
            'encut'   : encut,
            'gga'     : gga,
            'ialgo'   : ialgo,
            'ibrion'  : ibrion,
            'icharg'  : icharg,
            'isif'    : isif,
            'ismear'  : ismear,
            'ispin'   : ispin,
            'ivdw'    : ivdw,
            'lcharg'  : lcharg,
            'lepsilon': lepsilon,
            'lorbit'  : lorbit,
            'lreal'   : lreal,
            'luse_vdw': luse_vdw,
            'lwave'   : lwave,
            'magmom'  : magmom,
            'nbands'  : nbands,
            'nedos'   : nedos,
            'nelm'    : nelm,
            'nelmin'  : nelmin,
            'npar'    : npar,
            'nsw'     : nsw,
            'prec'    : prec,
            'pstress' : pstress,
            'sigma'   : sigma,
            'symprec' : symprec}

        self._tagorder = ['prec',
                          'gga',
                          'ivdw',
                          'luse_vdw',
                          'aggac',
                          'ibrion',
                          'nsw',
                          'algo',
                          'nelm',
                          'nelmin',
                          'isif',
                          'encut',
                          'ediff',
                          'ediffg',
                          'icharg',
                          'ispin',
                          'lorbit',
                          'magmom',
                          'ismear',
                          'sigma',
                          'nbands',
                          'emin',
                          'emax',
                          'nedos',
                          'pstress',
                          'ialgo',
                          'lreal',
                          'addgrid',
                          'lwave',
                          'lcharg',
                          'lepsilon',
                          'npar',
                          'symprec']

    def clear(self):
        for k in self._tagvals.keys():
            self._tagvals[k] = None

    def set_tag(self, k, v):
        if k in self._tagvals:
            self._tagvals[k] = v
        else:
            print "Key %s is not available."

    def set_addgrid(self, x):
        self._tagvals['addgrid'] = x

    def get_addgrid(self):
        return self._tagvals['addgrid']

    def set_aggac(self, x):
        self._tagvals['aggac'] = x

    def get_aggac(self):
        return self._tagvals['aggac']

    def set_algo(self, x):
        self._tagvals['algo'] = x

    def get_algo(self):
        return self._tagvals['algo']

    def set_ediff(self, x):
        self._tagvals['ediff'] = x

    def get_ediff(self):
        return self._tagvals['ediff']

    def set_ediffg(self, x):
        self._tagvals['ediffg'] = x

    def get_ediffg(self):
        return self._tagvals['ediffg']

    def set_emax(self, x):
        self._tagvals['emax'] = x

    def get_emax(self):
        return self._tagvals['emax']

    def set_emin(self, x):
        self._tagvals['emin'] = x

    def get_emin(self):
        return self._tagvals['emin']

    def set_encut(self, x):
        self._tagvals['encut'] = x

    def get_encut(self):
        return self._tagvals['encut']

    def set_gga(self, x):
        self._tagvals['gga'] = x

    def get_gga(self):
        return self._tagvals['gga']

    def set_ialgo(self, x):
        self._tagvals['ialgo'] = x

    def get_ialgo(self):
        return self._tagvals['ialgo']

    def set_ibrion(self, x):
        self._tagvals['ibrion'] = x

    def get_ibrion(self):
        return self._tagvals['ibrion']

    def set_icharg(self, x):
        self._tagvals['icharg'] = x

    def get_icharg(self):
        return self._tagvals['icharg']

    def set_isif(self, x):
        self._tagvals['isif'] = x

    def get_isif(self):
        return self._tagvals['isif']

    def set_ismear(self, x):
        self._tagvals['ismear'] = x

    def get_ismear(self):
        return self._tagvals['ismear']

    def set_ispin(self, x):
        self._tagvals['ispin'] = x

    def get_ispin(self):
        return self._tagvals['ispin']

    def set_ivdw(self, x):
        self._tagvals['ivdw'] = x

    def get_ivdw(self):
        return self._tagvals['ivdw']

    def set_lcharg(self, x):
        self._tagvals['lcharg'] = x

    def get_lcharg(self):
        return self._tagvals['lcharg']

    def set_lepsilon(self, x):
        self._tagvals['lepsilon'] = x

    def get_lepsilon(self):
        return self._tagvals['lepsilon']

    def set_lorbit(self, x):
        self._tagvals['lorbit'] = x

    def get_lorbit(self):
        return self._tagvals['lorbit']

    def set_lreal(self, x):
        self._tagvals['lreal'] = x

    def get_lreal(self):
        return self._tagvals['lreal']

    def set_luse_vdw(self, x):
        self._tagvals['luse_vdw'] = x

    def get_luse_vdw(self):
        return self._tagvals['luse_vdw']

    def set_lwave(self, x):
        self._tagvals['lwave'] = x

    def get_lwave(self):
        return self._tagvals['lwave']

    def set_magmom(self, x):
        self._tagvals['magmom'] = x

    def get_magmom(self):
        return self._tagvals['magmom']

    def set_nbands(self, x):
        self._tagvals['nbands'] = x

    def get_nbands(self):
        return self._tagvals['nbands']

    def set_nedos(self, x):
        self._tagvals['nedos'] = x

    def get_nedos(self):
        return self._tagvals['nedos']

    def set_nelm(self, x):
        self._tagvals['nelm'] = x

    def get_nelm(self):
        return self._tagvals['nelm']

    def set_nelmin(self, x):
        self._tagvals['nelmin'] = x

    def get_nelmin(self):
        return self._tagvals['nelmin']

    def set_npar(self, x):
        self._tagvals['npar'] = x

    def get_npar(self):
        return self._tagvals['npar']

    def set_nsw(self, x):
        self._tagvals['nsw'] = x

    def get_nsw(self):
        return self._tagvals['nsw']

    def set_prec(self, x):
        self._tagvals['prec'] = x

    def get_prec(self):
        return self._tagvals['prec']

    def set_pstress(self, x):
        self._tagvals['pstress'] = x

    def get_pstress(self):
        return self._tagvals['pstress']

    def set_sigma(self, x):
        self._tagvals['sigma'] = x

    def get_sigma(self):
        return self._tagvals['sigma']

    def set_symprec(self, x):
        self._tagvals['symprec'] = x

    def get_symprec(self):
        return self._tagvals['symprec']

    def set_electronic_structure(self):
        tags = self._tagvals
        tags['prec']    = "Accurate"
        tags['ibrion']  = -1
        tags['nelmin']  = 5
        tags['encut']   = 500
        tags['ediff']   = 1.0e-08
        tags['ismear']  = 0
        tags['sigma']   = 0.01
        tags['ialgo']   = 38
        tags['lreal']   = False
        tags['addgrid'] = True
        tags['lwave']   = False
        tags['lcharg']  = False

    def set_structure_optimization(self):
        tags = self._tagvals
        tags['prec']    = "Accurate"
        tags['ibrion']  = 2
        tags['nsw']     = 10
        tags['nelmin']  = 5
        tags['isif']    = 3
        tags['encut']   = 500
        tags['ediff']   = 1.0e-08
        tags['ediffg']  = -1.0e-08
        tags['ismear']  = 0
        tags['sigma']   = 0.01
        tags['ialgo']   = 38
        tags['lreal']   = False
        tags['addgrid'] = True
        tags['lwave']   = False
        tags['lcharg']  = False

    def copy(self):
        incar = Incar()
        for k, v in self._tagvals.iteritems():
            incar.set_tag(k, v)
        return incar

    def write(self, filename="INCAR"):
        names = self._tagnames

        w = open(filename, 'w')
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

        w.close()

def write_kpoints(filename="KPOINTS",
                  mesh=None,
                  shift=None,
                  gamma=False,
                  length=None,
                  kpoint=None):

    w = open(filename, 'w')
    if length:
        w.write("Automatic mesh\n")
        w.write("0\n")
        w.write("Auto\n")
        w.write("%4d\n" % length)
    elif kpoint is not None:
        w.write("Explicit k-point\n")
        w.write("1\n")
        w.write("Reciprocal\n")
        w.write("%10.7f %10.7f %10.7f  1\n" % tuple(kpoint))
    elif mesh is not None:
        w.write("Automatic mesh\n")
        w.write("0\n")
        if gamma:
            w.write("Gamma\n")
        else:
            w.write("Monkhorst-pack\n")
        w.write(" %5d %5d %5d\n" % tuple(mesh))
        if shift == None:
            w.write("     0.    0.    0.\n")
        else:
            w.write(" %5.3f %5.3f %5.3f\n" % tuple(shift))
        
    w.close()

class Outcar:
    def __init__(self, filename="OUTCAR"):
        self._filename = filename
        self._elastic_constants = None

    def get_elastic_constants(self):
        return self._elastic_constants

    def parse_elastic_constants(self):
        outcar = open(self._filename)
        hooked = False
        for line in outcar:
            if line.strip() == 'TOTAL ELASTIC MODULI (kBar)':
                hooked = True
                break

        if hooked:
            outcar.next()
            outcar.next()
            ec = []
            for i in range(6):
                pos = 8
                line = outcar.next()
                for j in range(6):
                    ec.append(float(line[pos:(pos+12)]))
                    pos += 12
    
            self._elastic_constants = np.reshape(ec, (6, 6))
            return True
        else:
            return False

 
class Vasprunxml:
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
        if self._eigenvalues_spin2 == None:
            return (self._eigenvalues_spin1,)
        else:
            return (self._eigenvalues_spin1, self._eigenvalues_spin2)

    def get_occupancies(self):
        # Degeneracy of electrons
        # Spin components 1 and 2 are stored in tuple as
        # (spin1, spin2) or (spin1,) if no spin polarized.
        if self._occupancies_spin2 == None:
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

                if element.tag != 'calculation':
                    continue
    
                for varray in element.findall('./varray'):
                    self._parse_forces_and_stress(varray, forces, stress)
        
                for varray in element.findall('./structure/varray'):
                    self._parse_points(varray, points)
                    
                for varray in element.findall('./structure/crystal/varray'):
                    self._parse_lattice(varray, lattice)
    
                for energy in element.findall('./energy'):
                    self._parse_energies(energy, energies)

                for array in element.findall('./array'):
                    if array.attrib['name'] == 'born_charges':
                        self._parse_born_charges(array, born_charges)

                for varray in element.findall('./varray'):
                    if varray.attrib['name'] == 'epsilon':
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

        except:
            return False

    def parse_parameters(self):
        try:
            for event, element in etree.iterparse(self._filename):

                if element.tag != 'parameters':
                    continue
                
                for separator in element.findall('./separator'):
                    if separator.attrib['name'] == 'electronic':
                        for i in separator.findall('./i'):
                            if i.attrib['name'] == 'NBANDS':
                                self._nbands = int(i.text)
            return True
        except:
            return False
            
    def parse_efermi(self):
        try:
            for event, element in etree.iterparse(self._filename):

                if element.tag != 'dos':
                    continue
                
                for i in element.findall('./i'):
                    if i.attrib['name'] == 'efermi':
                        efermi = float(i.text)
                        
            self._efermi = efermi
            return True
        except:
            return False
            
    def parse_eigenvalues(self):
        spin1 = []
        spin2 = []
        occ1 = []
        occ2 = []
        try:
            for event, element in etree.iterparse(self._filename):

                if element.tag != 'eigenvalues':
                    continue
                
                for array in element.findall('./array/set/set'):

                    if array.attrib['comment'] == 'spin 1':
                        self._parse_eigenvalues_spin(array, spin1, occ1)

                    if array.attrib['comment'] == 'spin 2':
                        self._parse_eigenvalues_spin(array, spin2, occ2)

            if spin1:
                self._eigenvalues_spin1 = np.array(spin1)
                self._occupancies_spin1 = np.array(occ1)
            if spin2:
                self._eigenvalues_spin2 = np.array(spin2)
                self._occupancies_spin2 = np.array(occ2)

            return self._parse_kpoints()
        
        except:
            return False

    def _parse_kpoints(self):
        try:
            for event, element in etree.iterparse(self._filename):

                if element.tag != 'kpoints':
                    continue
                
                kpoints = []
                weights = []
                for varray in element.findall('./varray'):
                    if varray.attrib['name'] == 'kpointlist':
                        for v in varray.findall('./v'):
                            kpoints.append(
                                [float(x) for x in v.text.split()])

                    if varray.attrib['name'] == 'weights':
                        for v in varray.findall('./v'):
                            weights.append(float(v.text))

            self._kpoints = np.array(kpoints)
            self._kpoint_weights = np.array(weights)
            return True

        except:
            return False

    def _parse_eigenvalues_spin(self, array, eigenvals, occupancies):
        for kset in array.findall('./set'):
            eigs = []
            occs = []
            for r in kset.findall('./r'):
                vals = r.text.split()
                eigs.append(float(vals[0]))
                occs.append(float(vals[1]))
            eigenvals.append(eigs)
            occupancies.append(occs)

    def _parse_forces_and_stress(self, varray, forces, stress):
        # force
        if varray.attrib['name'] == 'forces':
            forces_geomopt = []
            for v in varray.findall('./v'):
                forces_geomopt.append(
                    [float(x) for x in v.text.strip().split()])
            forces.append(forces_geomopt)

        # stress
        if varray.attrib['name'] == 'stress':
            stress_geomopt = []
            for v in varray.findall('./v'):
                stress_geomopt.append(
                    [float(x) for x in v.text.strip().split()])
            stress.append(stress_geomopt)
        
    def _parse_points(self, varray, points):
        # points
        if varray.attrib['name'] == 'positions':
            points_geomopt = []
            for v in varray.findall('./v'):
                points_geomopt.append(
                    [float(x) for x in v.text.strip().split()])
            points.append(np.transpose(points_geomopt))

    def _parse_lattice(self, varray, lattice):
        if varray.attrib['name'] == 'basis':
            lattice_geomopt = []
            for v in varray.findall('./v'):
                lattice_geomopt.append(
                    np.transpose(
                        [float(x) for x in v.text.strip().split()]))
            lattice.append(np.transpose(lattice_geomopt))

    def _parse_energies(self, energy, energies):
        energies_geomopt = []
        for v in energy.findall('./i'):
            energies_geomopt.append(
                [float(x) for x in v.text.strip().split()])
        energies.append(energies_geomopt)

    def _parse_born_charges(self, array, born_charges):
        for ion_set in array.findall('./set'):
            tensor = []
            for v in ion_set.findall('./v'):
                tensor.append([float(x) for x in v.text.strip().split()])
            born_charges.append(tensor)

    def _parse_vectors(self, varray, vectors):
        for v in varray.findall('./v'):
            vectors.append([float(x) for x in v.text.strip().split()])

class VasprunxmlExpat:
    def __init__(self, filename):
        self._filename = filename
        
        self._is_forces = False
        self._is_stress = False
        self._is_positions = False
        self._is_symbols = False
        self._is_basis = False
        self._is_energy = False

        self._is_v = False
        self._is_i = False
        self._is_rc = False
        self._is_c = False

        self._is_scstep = False
        self._is_structure = False

        self._all_forces = []
        self._all_stress = []
        self._all_points = []
        self._all_lattice = []
        self._symbols = []
        self._all_energies = []
        self._forces = None
        self._stress = None
        self._points = None        
        self._lattice = None
        self._energies = None

        self._p = xml.parsers.expat.ParserCreate()
        self._p.buffer_text = True
        self._p.StartElementHandler = self._start_element
        self._p.EndElementHandler = self._end_element
        self._p.CharacterDataHandler = self._char_data
    
    def parse(self):
        try:
            self._p.ParseFile(open(self._filename))
        except:
            return False
        else:
            return True

    def get_forces(self):
        return np.array(self._all_forces)

    def get_stress(self):
        return np.array(self._all_stress)

    def get_points(self):
        return np.array(self._all_points)

    def get_lattice(self):
        return np.array(self._all_lattice)
    
    def get_symbols(self):
        return self._symbols

    def get_cells(self):
        cells = []
        if len(self._all_points) == len(self._all_lattice):
            for p, l in zip(self._all_points, self._all_lattice):
                cells.append(Cell(lattice=l,
                                  points=p,
                                  symbols=self._symbols))
        return cells

    def get_energies(self):
        return np.array(self._all_energies)

    def _start_element(self, name, attrs):
        # Used not to collect energies in <scstep>
        if name == 'scstep': 
            self._is_scstep = True

        # Used not to collect basis and positions in
        # <structure name="initialpos" >
        # <structure name="finalpos" >
        if name == 'structure': 
            if 'name' in attrs.keys():
                self._is_structure = True

        if (self._is_forces or 
            self._is_stress or 
            self._is_positions or 
            self._is_basis):
            if name == 'v':
                self._is_v = True

        if name == 'varray':
            if 'name' in attrs.keys():
                if attrs['name'] == 'forces':
                    self._is_forces = True
                    self._forces = []
    
                if attrs['name'] == 'stress':
                    self._is_stress = True
                    self._stress = []
    
                if not self._is_structure:
                    if attrs['name'] == 'positions':
                        self._is_positions = True
                        self._points = []
        
                    if attrs['name'] == 'basis':
                        self._is_basis = True
                        self._lattice = []

        if self._is_energy and name == 'i':
            self._is_i = True

        if name == 'energy' and (not self._is_scstep):
            self._is_energy = True
            self._energies = []

        if self._is_symbols and name == 'rc':
            self._is_rc = True

        if self._is_symbols and self._is_rc and name == 'c':
            self._is_c = True
            
        if name == 'array':
            if 'name' in attrs.keys():
                if attrs['name'] == 'atoms':
                    self._is_symbols = True
            
                
    def _end_element(self, name):
        if name == 'scstep':
            self._is_scstep = False

        if name == 'structure' and self._is_structure:
            self._is_structure = False
        
        if name == 'varray':
            if self._is_forces:
                self._is_forces = False
                self._all_forces.append(self._forces)

            if self._is_stress:
                self._is_stress = False
                self._all_stress.append(self._stress)

            if self._is_positions:
                self._is_positions = False
                self._all_points.append(np.transpose(self._points))

            if self._is_basis:
                self._is_basis = False
                self._all_lattice.append(np.transpose(self._lattice))

        if name == 'array':
            if self._is_symbols:
                self._is_symbols = False
                

        if name == 'energy' and (not self._is_scstep):
            self._is_energy = False
            self._all_energies.append(self._energies)

        if name == 'v':
            self._is_v = False

        if name == 'i':
            self._is_i = False

        if name == 'rc':
            self._is_rc = False
            if self._is_symbols:
                self._symbols.pop(-1)

        if name == 'c':
            self._is_c = False
    
    def _char_data(self, data):
        if self._is_v:
            if self._is_forces:
                self._forces.append(
                    [float(x) for x in data.split()])

            if self._is_stress:
                self._stress.append(
                    [float(x) for x in data.split()])

            if self._is_positions:
                self._points.append(
                    [float(x) for x in data.split()])

            if self._is_basis:
                self._lattice.append(
                    [float(x) for x in data.split()])

        if self._is_i:
            if self._is_energy:
                self._energies.append(float(data.strip()))

        if self._is_c:
            if self._is_symbols:
                self._symbols.append(str(data.strip()))





if __name__ == '__main__':
    vxml = VasprunxmlExpat(sys.argv[1])
    vxml.parse()
    
    print "VasprunxmlExpat"
    print "Forces:"
    print vxml.get_forces()
    print "Stress:"
    print vxml.get_stress()
    print "Atomic points:"
    print vxml.get_points()
    print "Lattice:"
    print vxml.get_lattice()
    print "Energy:"
    print vxml.get_energies()

    print "Vasprunxml"
    vxml = Vasprunxml(sys.argv[1])
    vxml.parse_calculation()
    vxml.parse_eigenvalues()
    print "Forces:"
    print vxml.get_forces()
    print "Stress:"
    print vxml.get_stress()
    print "Lattice:"
    print vxml.get_lattice()
    print "Atomic points:"
    print vxml.get_points()
    print "Energy:"
    print vxml.get_energies()
    print "Eigenvalues:"
    print vxml.get_eigenvalues()
    print "Kpoints:"
    print vxml.get_kpoints()
    print "Occupancies:"
    print vxml.get_occupancies()
    print "Born charges"
    print vxml.get_born_charges()
    print "Epsilon"
    print vxml.get_epsilon()
