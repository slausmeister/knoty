import trimesh
from shapely.geometry import Point, Polygon
import numpy as np
import sympy as sp


def export_to_stl(knot, filename, resolution=200, domain=[0, 1], circular_resolution=10, radius=0.1):
    try:
        len(knot[0])
    except:
        linspace = np.linspace(domain[0], domain[1], resolution)
        x, y, z = [], [], []
        t = sp.symbols('t')  # Define the symbolic variable
        for tval in linspace:
            x.append(knot(t)[0].subs(t, tval).evalf())
            y.append(knot(t)[1].subs(t, tval).evalf())
            z.append(knot(t)[2].subs(t, tval).evalf())
        vertices = np.vstack([x, y, z]).T
    else:
        linspace = np.linspace(domain[0], domain[1], resolution)
        t = sp.symbols('t')  # Define the symbolic variable
        for i in range(0, len(knot[0])):
            x, y, z = [], [], []
            for tval in linspace:
                x.append(knot[0][i].subs(t, tval).evalf())
                y.append(knot[1][i].subs(t, tval).evalf())
                z.append(knot[2][i].subs(t, tval).evalf())
            vertices = np.vstack([x, y, z]).T

    # Create tube

    print(vertices)

    if vertices.shape[1] != 3:
        raise ValueError("Vertices must be a 3D path")

    circle = Point(0, 0).buffer(1)  # unit circle centered at origin
    polygon = Polygon(circle.exterior.coords)  # convert to polygon

    # Sweep the polygon along the path
    mesh = trimesh.creation.sweep_polygon(polygon, vertices)

    # Export the mesh to an STL file
    mesh.export(filename)

    print('Exported to', filename)
