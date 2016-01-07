import numpy as np
from mayavi import mlab
from cogue.crystal.atom import atomic_jmol_colors, covalent_radii
from cogue.crystal.utility import get_lattice_parameters


def set_figure():
    mlab.figure(bgcolor=(1, 1, 1))

def show():
    mlab.show()

def savefig(filename, size=None):
    mlab.savefig(filename, size=size)

def plot_cell(cell, margin=1e-5, color=(1, 0, 0)):
    _plot_lattice(cell.get_lattice(), color=color)
    _plot_axes(cell.get_lattice(), color=color)
    _plot_atoms(cell, margin=margin)

def _line_plot(m, n, pt, color=None):
    mlab.plot3d([pt[m][0], pt[n][0]],
                [pt[m][1], pt[n][1]],
                [pt[m][2], pt[n][2]],
                tube_radius=0.015, opacity=1, color=color)

def _plot_lattice(lattice, color=None):
    lat = lattice.T
    origin = np.zeros(3)
    pt = [origin,
          lat[0],
          lat[1],
          lat[2],
          lat[0] + lat[1],
          lat[1] + lat[2],
          lat[2] + lat[0],
          lat[0] + lat[1] + lat[2]]

    pt = np.array(pt)
    
    _line_plot(0, 1, pt, color)
    _line_plot(0, 2, pt, color)
    _line_plot(0, 3, pt, color)
    _line_plot(4, 7, pt, color)
    _line_plot(5, 7, pt, color)
    _line_plot(6, 7, pt, color)
    _line_plot(1, 4, pt, color)
    _line_plot(2, 5, pt, color)
    _line_plot(3, 6, pt, color)
    _line_plot(1, 6, pt, color)
    _line_plot(2, 4, pt, color)
    _line_plot(3, 5, pt, color)

def _plot_axes(lattice, color=(1, 0, 0)):
    lat = np.transpose([x/np.linalg.norm(x) for x in lattice.T])
    mlab.quiver3d([0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0],
                  lat[0],
                  lat[1],
                  lat[2],
                  color=color,
                  line_width=3,
                  scale_factor=1)

    for c, v in zip(('a','b','c'), (lat * 1.3).T):
        t = mlab.text3d(v[0] + 0.15, v[1], v[2], c, color=color, scale=0.3)
        # Workaround the bug
        #   https://github.com/enthought/mayavi/issues/169
        t.vector_text.update() 

def _plot_lattice_points(lattice, dim):
    lat_points = []
    for i in range(-dim[0], dim[0] + 1):
        for j in range(-dim[1], dim[1] + 1):
            for k in range(-dim[2], dim[2] + 1):
                lat_points.append([i, j, k])

    lp = np.dot(lattice, np.transpose(lat_points))
    mlab.points3d(lp[0], lp[1], lp[2],
                  scale_factor=0.2, opacity=0.2, color=(0,0,0))              

def _plot_atoms(cell, margin=1e-5, shift=[0,0,0], atom_scale=0.4):
    points = cell.get_points()
    points += np.reshape(shift, (3, 1))
    points -= np.floor(points)
    points, symbols = _get_points_with_margin(cell, margin)
    
    xs, ys, zs = np.dot(cell.get_lattice(), points)
    for x, y, z, s in zip(xs, ys, zs, symbols):
        color = tuple(np.array(atomic_jmol_colors[s], dtype=float) / 256)
        mlab.points3d(x, y, z,
                      resolution=16,
                      scale_factor=covalent_radii[s],
                      color=color)

def _get_points_with_margin(cell, margin=1e-5):
    abc = get_lattice_parameters(cell.get_lattice())
    points = cell.get_points()
    points_new = []
    symbols_new = []
    for p, s in zip(points.T, cell.get_symbols()):
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    p_inspect = p + np.array([i, j, k])
                    if ((p_inspect > 0 - margin / abc).all() and
                        (p_inspect < 1 + margin / abc).all()):
                        points_new.append(p_inspect)
                        symbols_new.append(s)
    return np.transpose(points_new), symbols_new
