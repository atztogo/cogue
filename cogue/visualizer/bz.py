import numpy as np
from mayavi import mlab
from mayavi.sources.vtk_data_source import VTKDataSource
from scipy.spatial import Voronoi
from tvtk.api import tvtk


class VisualizeBrillouinZone:
    def __init__(self, lattice, magnitude=1):
        """
        Parameters
        ----------
        lattice : ndarray
            Reciprocal basis vectors as row vectors.

        """

        self._magnitude = magnitude
        self._lattice = lattice * self._magnitude

    def run(self, with_axis=True):
        if self._set_Brillouin_zone():
            self._draw_voronoi()
            if with_axis:
                self._plot_axes()
            return True
        else:
            return False

    def get_voronoi_vertices(self):
        return self._vertices

    def get_BZ_faces(self):
        return self._faces

    def get_lattice(self):
        return self._lattice

    def _plot_axes(self, color=(1, 0, 0)):
        shortest = min(np.linalg.norm(self._lattice, axis=1))
        lat = np.array([x / np.linalg.norm(x) for x in self._lattice])
        lat *= shortest
        mlab.quiver3d(
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            lat[0],
            lat[1],
            lat[2],
            color=color,
            line_width=3,
            scale_factor=1,
        )

        scale = shortest / 8
        for c, v in zip(("a", "b", "c"), lat * (1 + scale / 2)):
            x, y, z = v
            x += scale / 2
            y -= scale / 4
            z -= scale / 4
            mlab.text3d(x, y, z, c, color=color, scale=scale)

    def _set_Brillouin_zone(self):
        points = np.dot(np.array(list(np.ndindex(5, 5, 5))) - [2, 2, 2], self._lattice)
        voronoi = Voronoi(points)
        voronoi_cells = [cell for cell in voronoi.regions if cell and (-1 not in cell)]

        if len(voronoi_cells) == 0:
            print("BZ is not unique.")
            return False

        norm2s = np.sum(voronoi.vertices**2, axis=1)
        BZ_cell = voronoi_cells[
            np.argmin([np.sum(norm2s[vcell]) for vcell in voronoi_cells])
        ]
        faces = [edge for edge in voronoi.ridge_vertices if -1 not in edge]
        BZ_faces = [face for face in faces if np.all([x in BZ_cell for x in face])]

        unique_points = np.unique([v for face in BZ_faces for v in face])
        unique_vertices = np.array(voronoi.vertices[unique_points])
        indices = {up: i for i, up in enumerate(unique_points)}
        new_BZ_faces = [[indices[x] for x in one_face] for one_face in BZ_faces]

        self._vertices = unique_vertices
        self._faces = new_BZ_faces

        return True

    def _draw_voronoi(self):
        faces = self._polydata()
        src = VTKDataSource(data=faces)
        mlab.pipeline.surface(src, opacity=0.5, representation="wireframe")

    def _polydata(self):
        pointArr = self._vertices
        faceArr = self._faces
        faces = tvtk.PolyData()
        faces.points = pointArr
        faces.polys = faceArr
        faces.point_data.scalars = [1] * len(pointArr)
        faces.point_data.scalars.name = "Height"
        return faces
