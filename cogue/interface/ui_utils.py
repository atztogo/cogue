"""UI utiles."""

import os
import sys

import numpy as np

from cogue.crystal.converter import reduce_points
from cogue.crystal.supercell import get_supercell
from cogue.crystal.utility import frac2val

r2h_observe_R3 = [[-1, 1, 1], [0, -1, 1], [1, 0, 1]]

h2r_observe_R3 = [
    [-1.0 / 3, -1.0 / 3, 2.0 / 3],
    [1.0 / 3, -2.0 / 3, 1.0 / 3],
    [1.0 / 3, 1.0 / 3, 1.0 / 3],
]

r2h_reverse_R3 = [[1, -1, 1], [0, 1, 1], [-1, 0, 1]]

h2r_reverse_R3 = [
    [1.0 / 3, 1.0 / 3, -2.0 / 3],
    [-1.0 / 3, 2.0 / 3, -1.0 / 3],
    [1.0 / 3, 1.0 / 3, 1.0 / 3],
]

r2h = r2h_observe_R3


def get_options(parser=None):
    """Get parser options."""
    if parser:
        (options, args) = parser.parse_args()
    else:
        (options, args) = get_parser().parse_args()

    return options, args


def get_parser():
    """Get OptionParser object."""
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Cogue tools.")
    parser.add_argument(
        "--bravais",
        dest="is_bravais",
        default=False,
        action="store_true",
        help="Transform to cell with Bravais lattice.",
    )
    parser.add_argument(
        "--dim",
        dest="s_mat",
        default=None,
        help="Supercell matrix",
    )
    parser.add_argument(
        "--noaxes",
        dest="with_axes",
        default=True,
        action="store_false",
        help=(
            "Transform primitive Rhombohedral to hexagonal Rhombohedral."
            "This has to be used exclusively to the other options."
        ),
    )
    parser.add_argument(
        "-o",
        dest="output_filename",
        default=None,
        help="Output filename",
    )
    parser.add_argument(
        "--r2h",
        dest="is_r2h",
        default=False,
        action="store_true",
        help=(
            "Transform primitive Rhombohedral to hexagonal Rhombohedral."
            "This has to be used exclusively to the other options."
        ),
    )
    parser.add_argument(
        "--shift",
        dest="shift",
        default=None,
        help="Origin shift",
    )
    parser.add_argument(
        "--tmat",
        dest="t_mat",
        default=None,
        help=(
            "Multiply transformation matrix. Absolute value of determinant "
            "has to be 1 or less than 1."
        ),
    )
    parser.add_argument(
        "-v",
        dest="is_verbose",
        default=False,
        action="store_true",
        help="More information is output.",
    )
    parser.add_argument("filenames", nargs="*", help="File names")
    return parser


def get_lines(filenames):
    """Return lines."""
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


def set_shift(cell, options):
    """Add shift to points."""
    shift = np.array([float(x) for x in options.shift.split()])
    if len(shift) == 3:
        points = cell.get_points()
        points += shift.reshape(3, 1)
        cell.set_points(points)
    else:
        sys.stderr.write("Atomic position shift is not correctly set.\n")


def transform_cell(cell_orig, options, is_shift=True):
    """Transform cell."""
    cell = cell_orig.copy()
    if options.shift and is_shift:
        set_shift(cell, options)
    if options.t_mat:
        cell = _get_tmat_cell(cell, options)
    if options.is_r2h:
        if options.is_verbose:
            print(
                "Transform cell by transformation matrix of rhombohedral "
                "to hexagonal:"
            )
            print(np.array(r2h))
        cell = get_supercell(cell, r2h)
    if options.s_mat:
        cell = _get_smat_cell(cell, options)

    return cell


def write_cells(write_func, cells, input_filenames=None, output_filename=None):
    """Write cells."""
    for i, cell in enumerate(cells):
        if len(cells) > 1:
            if output_filename:
                root, ext = os.path.splitext(output_filename)
                write_func(cell, root + "-%03d" % (i + 1) + ext)
            else:
                print("-" * len(input_filenames[i]))
                print("%s" % input_filenames[i])
                print("-" * len(input_filenames[i]))
                print(write_func(cell))
        else:
            if output_filename:
                write_func(cell, output_filename)
            else:
                sys.stdout.write(write_func(cell))


def _get_matrix(mat):
    if len(mat) == 3:
        return np.diag(mat)
    elif len(mat) == 9:
        mat.shape = (3, 3)
        return mat
    else:
        return False


def _get_tmat_cell(cell, options):
    t_mat = np.array([frac2val(x) for x in options.t_mat.split()])
    t_mat = _get_matrix(t_mat)
    if t_mat is False:
        sys.stderr.write("Transformation matrix is not correctly set.\n")
        return False
    else:
        if options.is_verbose:
            print("Transform cell using transformation matrix:")
            print(t_mat)
        return reduce_points(t_mat, cell)


def _get_smat_cell(cell, options):
    s_mat = np.array([int(x) for x in options.s_mat.split()])
    s_mat = _get_matrix(s_mat)
    if s_mat is False:
        sys.stderr.write("Supercell matrix is not correctly set.\n")
        return False
    else:
        if options.is_verbose:
            print("Transform cell using supercell matrix:")
            print(s_mat)
        return get_supercell(cell, s_mat)
