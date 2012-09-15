import sys
import numpy as np

def get_options(parser=None):
    if parser:
        (options, args) = parser.parse_args()
    else:
        (options, args) = get_parser().parse_args()

    return options, args

def get_parser():
    from optparse import OptionParser

    parser = OptionParser()
    parser.set_defaults(is_r2h=False,
                        is_bravais=False,
                        is_verbose=False,
                        output_filename=None,
                        s_mat=None,
                        t_mat=False,
                        shift=None)
    parser.add_option("--r2h",
                      dest="is_r2h",
                      action="store_true",
                      help="Transform primitive Rhombohedral to hexagonal Rhombohedral. This has to be used exclusively to the other options.")
    parser.add_option("--bravais",
                      dest="is_bravais",
                      action="store_true",
                      help="Transform to cell with Bravais lattice.")
    parser.add_option("--tmat",
                      dest="t_mat",
                      action="store",
                      type="string",                      
                      help="Multiply transformation matrix. Absolute value of determinant has to be 1 or less than 1.")
    parser.add_option("--dim",
                      dest="s_mat",
                      action="store",
                      type="string",                      
                      help="Supercell matrix")
    parser.add_option("--shift",
                      dest="shift",
                      action="store",
                      type="string",                      
                      help="Origin shift")
    parser.add_option("-o",
                      dest="output_filename",
                      action="store",
                      type="string",                      
                      help="Output filename")
    parser.add_option("-v",
                      dest="is_verbose",
                      action="store_true",
                      help="More information is output.")
    return parser

def get_lines(filenames):
    import fileinput
    file_obj = fileinput.input(filenames)
    lines = []
    filelines = []
    for line in file_obj:
        if file_obj.filelineno() == 1:
            if filelines:
                lines.append(filelines)
            filelines = []
        filelines.append(line)
    lines.append(filelines)
    return lines

def get_matrix(mat):
    if len(mat) == 3:
        return np.diag(mat)
    elif len(mat) == 9:
        mat.shape = (3, 3)
        return mat
    else:
        return False

def get_tmat_cell(cell, options):
    t_mat = np.array([float(x) for x in options.t_mat.split()])
    t_mat = get_matrix(t_mat)
    if t_mat is False:
        sys.stderr.write("Transformation matrix is not correctly set.\n")
        return False
    else:
        if options.is_verbose:
            print "Transform cell using transformation matrix:"
            print t_mat
        return reduce_points(t_mat, cell)

def get_smat_cell(cell, options):
    from phonopy.structure.atoms import Atoms
    from phonopy.structure.cells import get_supercell
    from cogue.crystal.converter import atoms2cell

    s_mat = np.array([int(x) for x in options.s_mat.split()])
    s_mat = get_matrix(s_mat)
    if s_mat is False:
        sys.stderr.write("Supercell matrix is not correctly set.\n")
        return False
    else:
        if options.is_verbose:
            print "Transform cell using supercell matrix:"
            print s_mat
        cell_phonopy = Atoms(cell=cell.get_lattice().T.copy(),
                             scaled_positions=cell.get_points().T.copy(),
                             symbols=cell.get_symbols()[:],
                             pbc=True)
        scell_phonopy = get_supercell(cell_phonopy, list(s_mat))
        return atoms2cell(scell_phonopy)

def transform_cell(cell, options):
    if options.is_r2h:
        cell = rhomb2hex(cell)
    else:
        if options.t_mat:
            cell = get_tmat_cell(cell, options)
        if options.s_mat:
            cell = get_smat_cell(cell, options)
    return cell

def write_cells(write_func, cells,
                input_filenames=None, output_filename=None):
    for i, cell in enumerate(cells):
        if len(cells) > 1:
            if output_filename:
                write_func(cell, output_filename + "%d" % (i + 1))
            else:
                print "-" * len(input_filenames[i])
                print input_filenames[i]
                print "-" * len(input_filenames[i])
                write_func(cell)
        else:
            if output_filename:
                write_func(cell, output_filename)
            else:
                write_func(cell)
    
