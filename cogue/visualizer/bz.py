from tvtk.api import tvtk
from mayavi.sources.vtk_data_source import VTKDataSource
from mayavi.modules.surface import Surface
from mayavi import mlab

from scipy.spatial import Voronoi
import numpy as np
import sys

def plot_axes(lattice, color=(1, 0, 0)):
    lat = np.transpose([x / np.linalg.norm(x) for x in lattice])
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
        mlab.text3d(v[0] + 0.15, v[1], v[2], c, color=color, scale=0.3)

def draw_voronoi(vertices, faces):
    data = polydata(vertices, faces)
    src = VTKDataSource(data=data)
    mlab.pipeline.surface(src, opacity=0.5, representation='wireframe')

def get_Brillouin_zone(points):
    voronoi = Voronoi(points)
    voronoi_cells = [cell for cell in voronoi.regions
                     if cell and (-1 not in cell)]

    if len(voronoi_cells) == 0:
        print "BZ is not unique."
        sys.exit(0)

    norm2s = np.sum(voronoi.vertices ** 2, axis=1)
    BZ_cell = voronoi_cells[
        np.argmin([np.sum(norm2s[vcell]) for vcell in voronoi_cells])]
    faces = [edge for edge in voronoi.ridge_vertices if -1 not in edge]
    BZ_faces = [face for face in faces if np.all([x in BZ_cell for x in face])]

    return voronoi.vertices, BZ_faces

def polydata(vertices, faces):
    pointArr = vertices
    faceArr = faces
    faces = tvtk.PolyData()
    faces.points=pointArr
    faces.polys=faceArr
    faces.point_data.scalars = [1] * len(pointArr)
    faces.point_data.scalars.name = 'Height'

    return faces

a, b, c = np.reshape([float(x) for x in sys.argv[1].split()], (3, 3))
points = np.dot(np.array(list(np.ndindex(5, 5, 5))) - [2, 2, 2], [a, b, c])
vertices, faces = get_Brillouin_zone(points)
draw_voronoi(vertices, faces)
plot_axes([a, b, c])

mlab.show()
