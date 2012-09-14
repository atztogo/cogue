import sys
import numpy as np
from cogue.crystal.cell import Cell
from cogue.crystal.symmetry import get_symmetry_dataset, get_crystallographic_cell

#
# Utility functions
#
def frac2val(string):
    num, denom = [float(x) for x in string.split('/')]
    return num / denom

def get_angles(lattice):
    a, b, c = get_cell_parameters(lattice)
    alpha = np.arccos(np.vdot(lattice[:,1], lattice[:,2]) / b / c)
    beta  = np.arccos(np.vdot(lattice[:,2], lattice[:,0]) / c / a)
    gamma = np.arccos(np.vdot(lattice[:,0], lattice[:,1]) / a / b)
    return alpha / np.pi * 180, beta / np.pi * 180, gamma / np.pi * 180

def lattice2cartesian(a, b, c, alpha, beta, gamma):
    """
    The conversion refers the wikipedia,
    http://en.wikipedia.org/wiki/Fractional_coordinates
    """
    cg = np.cos(gamma / 180 * np.pi)
    cb = np.cos(beta / 180 * np.pi)
    ca = np.cos(alpha / 180 * np.pi)
    sg = np.sin(gamma / 180 * np.pi)
    L = np.zeros((3, 3))
    L[0, 0] = a
    L[0, 1] = b * cg
    L[0, 2] = c * cb
    L[1, 1] = b * sg
    L[1, 2] = c * (ca - cb * cg) / sg
    L[2, 2] = c * np.sqrt(1 - ca ** 2 - cb ** 2 - cg ** 2 +
                          2 * ca * cb * cg) / sg
    return L

def get_cell_parameters(lattice):
    return np.sqrt(np.dot(lattice.T, lattice).diagonal())

def get_oriented_lattice(lattice):
    a, b, c = get_cell_parameters(lattice)
    alpha, beta, gamma = get_angles(lattice)
    alpha *= np.pi / 180
    beta *= np.pi / 180
    gamma *= np.pi / 180
    a1 = a
    a2 = 0.0
    a3 = 0.0
    b1 = np.cos(gamma)
    b2 = np.sin(gamma)
    b3 = 0.0
    c1 = np.cos(beta)
    c2 = (2 * np.cos(alpha) + b1**2 + b2**2 - 2 * b1 * c1 - 1) / (2 * b2)
    c3 = np.sqrt(1 - c1**2 - c2**2)
    lattice = np.zeros((3, 3), dtype=float)
    lattice[0, 0] = a
    lattice[:,1] = np.array([b1, b2, b3]) * b
    lattice[:,2] = np.array([c1, c2, c3]) * c
    return lattice

def get_primitive(cell, tolerance=1e-5):
    # spglib returns R-centred lattice for Rhombohedrals
    brv_cell = get_crystallographic_cell(cell, tolerance)
    sym_dataset = get_symmetry_dataset(brv_cell)
    spg_symbol = sym_dataset['international'][0]
    if spg_symbol == 'F':
        brv_cell = fc2prim(brv_cell)
    elif spg_symbol == 'I':
        brv_cell = bc2prim(brv_cell)
    elif spg_symbol == 'A':
        brv_cell = abc2prim(brv_cell)
    elif spg_symbol == 'B':
        brv_cell = bbc2prim(brv_cell)
    elif spg_symbol == 'C':
        brv_cell = cbc2prim(brv_cell)
    return brv_cell

def rhomb2hex(cell):
    shift = [[0, 0, 0],
             [1, 0, 0],
             [0, 0,-1]] 
    tmat = [[ 1, 0, 1],
            [-1, 1, 1],
            [ 0,-1, 1]]
    natom = len(cell.get_symbols())
    points = cell.get_points()
    lattice = get_oriented_lattice(np.dot(cell.get_lattice(), tmat))
    triple_pos = np.zeros((3, natom * 3), dtype=float)
    for i in range(3):
        triple_pos[:,i*natom:(i+1)*natom] = (points.T + shift[i]).T
    
    points_hex = np.dot(np.linalg.inv(tmat), triple_pos)
    points_hex -= np.floor(points_hex)
    
    return Cell(lattice=lattice,
                points=points_hex,
                symbols=cell.get_symbols() * 3)

def reduce_points(tmat, cell, tolerance=1e-5):
    points_prim = []
    symbols_prim = []
    symbols = cell.get_symbols()
    for i, p in enumerate(
        np.dot(np.linalg.inv(tmat), cell.get_points()).T):
        is_different = True
        for pp in points_prim:
            diff = pp - p
            if (abs(diff - diff.round()) < tolerance).all():
                is_different = False
                break
        if is_different:
            points_prim.append(p - np.floor(p))
            symbols_prim.append(symbols[i])

    return Cell(lattice=np.dot(cell.get_lattice(), tmat),
                points=np.transpose(points_prim),
                symbols=symbols_prim)

def fc2prim(cell):
    tmat = [[ 0.0, 0.5, 0.5],
            [ 0.5, 0.0, 0.5],
            [ 0.5, 0.5, 0.0]]
    return reduce_points(tmat, cell)

def bc2prim(cell):
    tmat = [[-0.5, 0.5, 0.5],
            [ 0.5,-0.5, 0.5],
            [ 0.5, 0.5,-0.5]]
    return reduce_points(tmat, cell)

def cbc2prim(cell):
    tmat = [[ 0.5, 0.5, 0.0],
            [-0.5, 0.5, 0.0],
            [ 0.0, 0.0, 1.0]]
    return reduce_points(tmat, cell)

def abc2prim(cell):
    tmat = [[ 0.0, 0.5, 0.5],
            [ 0.0,-0.5, 0.5],
            [ 1.0, 0.0, 0.0]]
    return reduce_points(tmat, cell)

def bbc2prim(cell):
    tmat = [[ 0.5, 0.0, 0.5],
            [ 0.5, 0.0,-0.5],
            [ 0.0, 1.0, 0.0]]
    return reduce_points(tmat, cell)

def atoms2cell(phonopy_cell):
    return Cell(lattice=phonopy_cell.get_cell().T,
                points=phonopy_cell.get_scaled_positions().T,
                symbols=phonopy_cell.symbols)

#######################
# Writers and readers #
#######################

#
# yaml
#
def read_yaml(filename):
    import yaml
    data = yaml.load(open(filename))

    if 'status' in data:
        if (not data['status'] == 'terminate' and
            'lattice' in data and
            'points' in data and
            'symbols' in data):
            lattice = np.transpose(data['lattice'])
            points = np.transpose(data['points'])
            symbols = data['symbols']
            return Cell(lattice = lattice,
                        points = points,
                        symbols = symbols)
    return None

#
# Cif
#
def write_cif_P1(cell, filename=None):
    a, b, c = get_cell_parameters(cell.get_lattice())
    alpha, beta, gamma = get_angles(cell.get_lattice())
    
    cif = """data_cogue_crystal_converter

_symmetry_space_group_name_H-M     'P 1'
_symmetry_Int_Tables_number        1

_cell_length_a                     %.5f
_cell_length_b                     %.5f
_cell_length_c                     %.5f
_cell_angle_alpha                  %.5f
_cell_angle_beta                   %.5f
_cell_angle_gamma                  %.5f
_cell_volume                       %.5f
_cell_formula_units_Z              1

loop_
_space_group_symop_operation_xyz
x,y,z

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy\n""" % (a, b, c, alpha, beta, gamma, cell.get_volume())

    symbols = []
    for s, p in zip(cell.get_symbols(), cell.get_points().T):
        symbols.append(s)
        cif += "%-7s%2s %10.5f%10.5f%10.5f   1.00000\n" % (s + "%d" % symbols.count(s), s, p[0], p[1], p[2])
        
    if filename:
        w = open(filename, 'w')
        w.write(cif)
        w.close()
    else:
        print cif,

def read_cif(str_cif):
    """
    A sample of cif file generated by openbabel
    
    data_I
    _chemical_name_common 'filename.cif'
    _cell_length_a 4.0021
    _cell_length_b 13.8984
    _cell_length_c 3.6111
    _cell_angle_alpha 90
    _cell_angle_beta 90
    _cell_angle_gamma 90
    _space_group_name_H-M_alt 'C m m m'
    _space_group_name_Hall '-C 2 2'
    loop_
        _symmetry_equiv_pos_as_xyz
        'x,y,z'
        '-x,-y,z'
        '-x,y,-z'
        'x,-y,-z'
        '-x,-y,-z'
        'x,y,-z'
        'x,-y,z'
        '-x,y,z'
        '1/2+x,1/2+y,z'
        '1/2-x,1/2-y,z'
        '1/2-x,1/2+y,-z'
        '1/2+x,1/2-y,-z'
        '1/2-x,1/2-y,-z'
        '1/2+x,1/2+y,-z'
        '1/2+x,1/2-y,z'
        '1/2-x,1/2+y,z'
    loop_
        _atom_site_type_symbol
        _atom_site_label
        _atom_site_Cartn_x
        _atom_site_Cartn_y
        _atom_site_Cartn_z
        Er   Er1    -0.00000    5.00759    1.80555
        Ni   Ni2    -0.00000    2.77968    0.00000
        Pb   Pb3     0.00000    0.00000    0.00000
    """

    atom_site_order = []
    str_atom_site = ""
    str_symmetry_equiv_pos_as_xyz = ""
    loop_mode = False
    for line in str_cif.splitlines():
        if not line.strip():
            loop_mode = False
            continue

        if loop_mode:
            if '_symmetry_equiv_pos_as_xyz' in line:
                loop_mode = 'symmetry_equiv_pos_as_xyz'
                continue
            elif '_atom_site_type_symbol' in line:
                atom_site_order.append('type_symbol')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_label' in line:
                atom_site_order.append('label')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_U_iso_or_equiv' in line:
                atom_site_order.append('U_iso_or_equiv')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_thermal_displace_type' in line:
                atom_site_order.append('thermal_displace_type')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_occupancy' in line:
                atom_site_order.append('occupancy')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_Cartn_x' in line:
                atom_site_order.append('Cartn_x')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_Cartn_y' in line:
                atom_site_order.append('Cartn_y')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_Cartn_z' in line:
                atom_site_order.append('Cartn_z')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_fract_x' in line:
                atom_site_order.append('fract_x')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_fract_y' in line:
                atom_site_order.append('fract_y')
                loop_mode = 'atom_site'
                continue
            elif '_atom_site_fract_z' in line:
                atom_site_order.append('fract_z')
                loop_mode = 'atom_site'
                continue


        if 'loop_' in line:
            loop_mode = True
        elif '_cell_length_a' in line:
            a = float(line.split()[1])
            loop_mode = False
            continue
        elif '_cell_length_b' in line:
            b = float(line.split()[1])
            loop_mode = False
            continue
        elif '_cell_length_c' in line:
            c = float(line.split()[1])
            loop_mode = False
            continue
        elif '_cell_angle_alpha' in line:
            alpha = float(line.split()[1])
            loop_mode = False
            continue
        elif '_cell_angle_beta' in line:
            beta = float(line.split()[1])
            loop_mode = False
            continue
        elif '_cell_angle_gamma' in line:
            gamma = float(line.split()[1])
            loop_mode = False
            continue

        if loop_mode == 'symmetry_equiv_pos_as_xyz':
            str_symmetry_equiv_pos_as_xyz += line + '\n'
        elif loop_mode == 'atom_site':
            str_atom_site += line  + '\n'

    return get_cell_from_cif_info(str_symmetry_equiv_pos_as_xyz,
                                  str_atom_site,
                                  atom_site_order,
                                  lattice2cartesian(a, b, c,
                                                    alpha, beta, gamma))
            
def get_cell_from_cif_info(str_sym, str_atom, order, lattice):
    symbols, points = read_atom_site_lines(str_atom, order, lattice)
    operations = expand_symmetry_equiv_pos_as_xyz(str_sym)
    points = np.vstack((np.array(points), np.ones(len(symbols))))
    all_points = np.dot(operations[0], points)
    for opn in operations[1:]:
        all_points = np.hstack((all_points, np.dot(opn, points)))
    symbols = symbols * len(operations)
    red_symbols, red_points = remove_overlapping_points(all_points,
                                                        symbols)

    return Cell(lattice=lattice,
                points=red_points,
                symbols=red_symbols)

def remove_overlapping_points(points, symbols, symprec=1e-2):
    red_pos = []
    red_sym = []

    for sym, pos in zip(symbols, points.T):
        ok = True
        for rpos in red_pos:
            diff = pos - rpos
            if (abs(diff.round() - diff) < symprec).all():
                ok = False
                break
        if ok:
            red_pos.append(pos)
            red_sym.append(sym)

    return red_sym, (np.array(red_pos) - np.floor(red_pos)).T

def read_atom_site_lines(str_lines, atom_site_order, lattice):
    ind_symbol = atom_site_order.index("type_symbol")
    if "Cartn_x" in atom_site_order:
        index_x = atom_site_order.index("Cartn_x")
        index_y = atom_site_order.index("Cartn_y")
        index_z = atom_site_order.index("Cartn_z")
    if "fract_x" in atom_site_order:
        index_x = atom_site_order.index("fract_x")
        index_y = atom_site_order.index("fract_y")
        index_z = atom_site_order.index("fract_z")
    points = []
    symbols = []
    for line in str_lines.splitlines():
        vals = line.split()
        symbols.append(vals[ind_symbol])
        x = float(vals[index_x])
        y = float(vals[index_y])
        z = float(vals[index_z])
        points.append([x, y, z])

    points = np.transpose(points)
    if "Cartn_x" in atom_site_order:
        points = np.dot(np.linalg.inv(lattice), points)
        
    return symbols, points

def expand_symmetry_equiv_pos_as_xyz(str_lines):
    operations = []
    for line in str_lines.splitlines():
        if line.strip():
            xyz = [x for x in line.translate(None, '\' ').split(',')]
            operation = np.zeros((3, 4))
            for i, symbol in enumerate(xyz):
                parts = split_xyz_symbol(symbol)
                for string in parts:
                    if '/' in string:
                        operation[i, 3] += frac2val(string)
                    if 'x' in string:
                        operation[i, 0] += int(string.replace('x', '1'))
                    if 'y' in string:
                        operation[i, 1] += int(string.replace('y', '1'))
                    if 'z' in string:
                        operation[i, 2] += int(string.replace('z', '1'))

            operations.append(operation)

    return np.array(operations)
                
def split_xyz_symbol(symbol):
    def plusminus(char):
        if char == '+' or char == '-':
            return True
        else:
            return False
    def xyz(char):
        if char == 'x' or char == 'y' or char == 'z':
            return True
        else:
            False
    def slash(char):
        if char == '/':
            return True
        else:
            False

    parts = []
    string = ''
    for char in symbol:
        if plusminus(char):
            if string:
                parts.append(string)
            string = char
            continue
        if xyz(char):
            string += char
            parts.append(string)
            string = ''
            continue
        if slash(char):
            string += char
            continue
        if char.isdigit():
            string += char
            continue
    if string:
        parts.append(string)

    return parts

#
# V_sim ascii
#    
def write_v_sim(cell, filename=None):
    lat = get_oriented_lattice(cell.get_lattice())
    text  = "# cogue generated file\n"
    # text += "%15.9f%15.9f%15.9f\n" % tuple(get_cell_parameters(lattice))
    # text += "%15.9f%15.9f%15.9f\n" % tuple(get_angles(lattice))
    text += "%15.9f%15.9f%15.9f\n" % (lat[0,0], lat[0,1], lat[1,1])
    text += "%15.9f%15.9f%15.9f\n" % (lat[0,2], lat[1,2], lat[2,2])
    text += "#keyword: reduced\n"
    # text += "#keyword: angdeg, reduced\n"
    for s, p in zip(cell.get_symbols(), cell.get_points().T):
        text += "%15.9f%15.9f%15.9f %2s\n" % (p[0], p[1], p[2], s)

    if filename:
        w = open(filename, 'w')
        w.write(text)
        w.close()
    else:
        print text,
                                                                        
if __name__ == '__main__':
    pass
