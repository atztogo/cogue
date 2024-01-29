"""Sketch interface."""

import numpy as np

from cogue.crystal.atom import covalent_radii


def get_tex_template(width, height):
    """Return latex header."""
    return r"""\documentclass{minimal}
\usepackage{tikz}
\usepackage[paperwidth=%fcm,paperheight=%fcm,hmargin=0cm,vmargin=0cm]{geometry}
\begin{document}
\include{filename}
\end{document}""" % (
        width,
        height,
    )


class Sketch:
    """Sketch class."""

    def __init__(self, lattice, points, symbols, scale=1.0):
        """Init method."""
        self._scale = scale
        self._lattice = lattice * scale
        self._points = points
        self._shift = -self._lattice.sum(axis=1) / 2
        self._symbols = symbols
        self._text = None
        self._picturebox = None

    def write(self, filename):
        """Write sketch file."""
        f = open(filename, "w")
        f.write(self._text)
        f.close()

    def get_text(self):
        """Return sketch lines."""
        return self._text

    def set_text(self, eye=None, look_at=None, up_vector=None, transform=None):
        """Set text."""
        # Sort atoms by symbols
        indices_of_symbols = {}
        for i, s in enumerate(self._symbols):
            if s not in indices_of_symbols:
                indices_of_symbols[s] = []
            indices_of_symbols[s].append(i)

        # Define atoms
        text = self._atoms_text(indices_of_symbols)

        # Define lattice
        text += self._lattice_text()

        # Show
        if transform:
            text += "put { %s } { {atoms}{lattice} }\n" % transform
        elif eye is not None and look_at is not None:
            # View
            text += "def eye (%f,%f,%f)\n" % tuple(eye)
            text += "def look_at (%f,%f,%f)\n" % tuple(look_at)
            if up_vector is None:
                text += "put { view((eye), (look_at)) } { {atoms}{lattice} }\n"
            else:
                text += "def up_vector [%f,%f,%f]\n" % tuple(up_vector)
                text += "put { view((eye), (look_at), [up_vector]) } "
                text += "{ {atoms}{lattice} }\n"
        else:
            text += "{atoms}{lattice}\n"

        # Footer
        text += "global {\n"
        text += "  language tikz\n"
        if self._picturebox is not None:
            text += "  picturebox(%f,%f)(%f,%f)\n" % (
                tuple(self._picturebox[0]) + tuple(self._picturebox[1])
            )
        text += "}\n"

        self._text = text

    def set_picturebox(self, p1, p2):
        """Set picturebox."""
        self._picturebox = [p1, p2]

    def _lattice_text(self):
        text = ""
        text += "def O (%f,%f,%f)\n" % tuple(self._shift)
        text += "def A (%f,%f,%f)\n" % tuple(self._lattice[:, 0] + self._shift)
        text += "def B (%f,%f,%f)\n" % tuple(self._lattice[:, 1] + self._shift)
        text += "def C (%f,%f,%f)\n" % tuple(self._lattice[:, 2] + self._shift)
        text += "def AB (%f,%f,%f)\n" % tuple(
            self._lattice[:, 0] + self._lattice[:, 1] + self._shift
        )
        text += "def BC (%f,%f,%f)\n" % tuple(
            self._lattice[:, 1] + self._lattice[:, 2] + self._shift
        )
        text += "def CA (%f,%f,%f)\n" % tuple(
            self._lattice[:, 2] + self._lattice[:, 0] + self._shift
        )
        text += "def ABC (%f,%f,%f)\n" % tuple(self._lattice.sum(axis=1) + self._shift)

        text += "def bottom polygon[cull=false, fill opacity=0, line width=2](O)(A)(AB)(B)\n"  # noqa E501
        text += "def top polygon[cull=false, fill opacity=0, line width=2](C)(CA)(ABC)(BC)\n"  # noqa E501
        text += "def frontA polygon[cull=false, fill opacity=0, line width=2](O)(B)(BC)(C)\n"  # noqa E501
        text += "def frontB polygon[cull=false, fill opacity=0, line width=2](O)(A)(CA)(C)\n"  # noqa E501
        text += "def backA polygon[cull=false, fill opacity=0, line width=2](A)(AB)(ABC)(CA)\n"  # noqa E501
        text += "def backB polygon[cull=false, fill opacity=0, line width=2](B)(AB)(ABC)(BC)\n"  # noqa E501

        text += "def lattice { {bottom}{top}{frontA}{frontB}{backA}{backB} }\n"

        return text

    def _atoms_text(self, indices_of_symbols):
        text = ""

        positions = np.dot(self._lattice, self._points)
        positions += self._shift.reshape(3, 1)
        for s, indices in indices_of_symbols.iteritems():
            for i in indices:
                p = positions[:, i]
                text += "def %s%d (%f, %f, %f)\n" % ((s, i) + tuple(p))

        for s, indices in indices_of_symbols.iteritems():
            color = np.array(atomic_colors[s], dtype=float) / 256
            radius = covalent_radii[s] * self._scale
            text += (
                r"special|\definecolor{%s_color}{rgb}{%f,%f,%f}|"
                "[lay=under]\n" % ((s,) + tuple(color))
            )
            text += "def %s_atoms dots[ball color=%s_color, " "dotsize=%f] " % (
                s,
                s,
                radius,
            )
            for i in indices:
                text += "(%s%i)" % (s, i)
            text += "\n"
        text += "def atoms { "
        for s in indices_of_symbols:
            text += "{%s_atoms}" % s
        text += " }\n"

        return text

    def _arrows_text(self):
        pass


class SketchCell(Sketch):
    """SketchCell class."""

    def __init__(self, cell, scale=1.0):
        """Init method."""
        Sketch.__init__(
            self, cell.lattice, cell.get_points(), cell.get_symbols(), scale=scale
        )


# this colors were copied from jmol-colors.
# http://jmol.sourceforge.net/jscolors/
atomic_colors = {
    "H": (255, 255, 255),
    "He": (217, 255, 255),
    "Li": (204, 128, 255),
    "Be": (194, 255, 0),
    "B": (255, 181, 181),
    "C": (144, 144, 144),
    "N": (48, 80, 248),
    "O": (255, 13, 13),
    "F": (144, 224, 80),
    "Ne": (179, 227, 245),
    "Na": (171, 92, 242),
    "Mg": (138, 255, 0),
    "Al": (191, 166, 166),
    "Si": (240, 200, 160),
    "P": (255, 128, 0),
    "S": (255, 255, 48),
    "Cl": (31, 240, 31),
    "Ar": (128, 209, 227),
    "K": (143, 64, 212),
    "Ca": (61, 255, 0),
    "Sc": (230, 230, 230),
    "Ti": (191, 194, 199),
    "V": (166, 166, 171),
    "Cr": (138, 153, 199),
    "Mn": (156, 122, 199),
    "Fe": (224, 102, 51),
    "Co": (240, 144, 160),
    "Ni": (80, 208, 80),
    "Cu": (200, 128, 51),
    "Zn": (125, 128, 176),
    "Ga": (194, 143, 143),
    "Ge": (102, 143, 143),
    "As": (189, 128, 227),
    "Se": (255, 161, 0),
    "Br": (166, 41, 41),
    "Kr": (92, 184, 209),
    "Rb": (112, 46, 176),
    "Sr": (0, 255, 0),
    "Y": (148, 255, 255),
    "Zr": (148, 224, 224),
    "Nb": (115, 194, 201),
    "Mo": (84, 181, 181),
    "Tc": (59, 158, 158),
    "Ru": (36, 143, 143),
    "Rh": (10, 125, 140),
    "Pd": (0, 105, 133),
    "Ag": (192, 192, 192),
    "Cd": (255, 217, 143),
    "In": (166, 117, 115),
    "Sn": (102, 128, 128),
    "Sb": (158, 99, 181),
    "Te": (212, 122, 0),
    "I": (148, 0, 148),
    "Xe": (66, 158, 176),
    "Cs": (87, 23, 143),
    "Ba": (0, 201, 0),
    "La": (112, 212, 255),
    "Ce": (255, 255, 199),
    "Pr": (217, 255, 199),
    "Nd": (199, 255, 199),
    "Pm": (163, 255, 199),
    "Sm": (143, 255, 199),
    "Eu": (97, 255, 199),
    "Gd": (69, 255, 199),
    "Tb": (48, 255, 199),
    "Dy": (31, 255, 199),
    "Ho": (0, 255, 156),
    "Er": (0, 230, 117),
    "Tm": (0, 212, 82),
    "Yb": (0, 191, 56),
    "Lu": (0, 171, 36),
    "Hf": (77, 194, 255),
    "Ta": (77, 166, 255),
    "W": (33, 148, 214),
    "Re": (38, 125, 171),
    "Os": (38, 102, 150),
    "Ir": (23, 84, 135),
    "Pt": (208, 208, 224),
    "Au": (255, 209, 35),
    "Hg": (184, 184, 208),
    "Tl": (166, 84, 77),
    "Pb": (87, 89, 97),
    "Bi": (158, 79, 181),
    "Po": (171, 92, 0),
    "At": (117, 79, 69),
    "Rn": (66, 130, 150),
    "Fr": (66, 0, 102),
    "Ra": (0, 125, 0),
    "Ac": (112, 171, 250),
    "Th": (0, 186, 255),
    "Pa": (0, 161, 255),
    "U": (0, 143, 255),
    "Np": (0, 128, 255),
    "Pu": (0, 107, 255),
    "Am": (84, 92, 242),
    "Cm": (120, 92, 227),
    "Bk": (138, 79, 227),
    "Cf": (161, 54, 212),
    "Es": (179, 31, 212),
    "Fm": (179, 31, 186),
    "Md": (179, 13, 166),
    "No": (189, 13, 135),
    "Lr": (199, 0, 102),
    "Rf": (204, 0, 89),
    "Db": (209, 0, 79),
    "Sg": (217, 0, 69),
    "Bh": (224, 0, 56),
    "Hs": (230, 0, 46),
    "Mt": (235, 0, 38),
}

if __name__ == "__main__":
    pass
