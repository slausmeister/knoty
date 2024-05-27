import numpy as np
import sympy as sp
import k3d

def _get_basis(matrix, x_val, y_val, z_val):
    """
    This function takes in a one form and a point in R3 and calculates the basis of the kernel

    :matrix: One form represented as a sympy matrix
    :x_val, y_val, z_val: Coordinates in R3
    :return: Returns basis vectors of kernel
    """
    # Define symbolic variables x, y, z
    x, y, z = sp.symbols('x y z')

    # Convert the matrix elements to SymPy expressions
    matrix_sympy = [sp.sympify(element) for element in matrix]

    # For a matrix [a, b, c], the kernel is solved by ax + by + cz = 0
    a, b, c = matrix_sympy

    # Constructing the basis vectors
    basis1 = np.array([1, 0, -a/c], dtype=object) if c != 0 else np.array([1, 0, 0], dtype=object)
    basis2 = np.array([0, 1, -b/c], dtype=object) if c != 0 else np.array([0 , 1, 0], dtype=object)

    # Function to substitute values into a symbolic expression or return the value if it's not symbolic
    def substitute_if_symbolic(expr, substitutions):
        return expr.subs(substitutions) if isinstance(expr, sp.Expr) else expr

    # Substituting the values of x, y, z
    substitutions = {x: x_val, y: y_val, z: z_val}
    basis1_evaluated = np.array([substitute_if_symbolic(el, substitutions) for el in basis1], dtype=float)
    basis2_evaluated = np.array([substitute_if_symbolic(el, substitutions) for el in basis2], dtype=float)

    return basis1_evaluated, basis2_evaluated

def _plot_plane(x, y, z, form, alpha, surfcolor='blue', plot_radius=0.1):
    """
    Plots a surface in R3 at a specific coordinate (x, y, z) spanned by two vectors (v1, v2) using K3D-Jupyter.

    Parameters:
    - x, y, z: The coordinates where the surface originates.
    - v1, v2: The vectors spanning the surface.
    - plot_radius: The half-length of the sides of the surface.
    - height_limit: The maximum height variation allowed from the center.
    - surfcolor: The color of the surface (in hexadecimal).
    - alpha: The transparency of the surface.
    """

    v1, v2 = _get_basis(form, x, y, z)

    # Create a grid on the plane
    u, v = np.meshgrid(np.linspace(-plot_radius, plot_radius, 2), np.linspace(-plot_radius, plot_radius, 2))
    plane_x = x + u * v1[0] + v * v2[0]
    plane_y = y + u * v1[1] + v * v2[1]
    plane_z = z + u * v1[2] + v * v2[2]

    # Clip the z-values to be within the height_limit

    # Flatten the x, y, z coordinates for K3D
    vertices = np.vstack([plane_x.ravel(), plane_y.ravel(), plane_z.ravel()]).T.astype(np.float32)

    # Generate indices for the triangles
    indices = []
    for i in range(vertices.shape[0] - np.sqrt(vertices.shape[0]).astype(int) - 1):
        if (i + 1) % np.sqrt(vertices.shape[0]).astype(int) != 0:
            indices.append([i, i + 1, i + np.sqrt(vertices.shape[0]).astype(int)])
            indices.append([i + 1, i + 1 + np.sqrt(vertices.shape[0]).astype(int), i + np.sqrt(vertices.shape[0]).astype(int)])
    indices = np.array(indices).astype(np.uint32)

    # Plot the surface

    tmp = k3d.mesh(vertices, indices, opacity=alpha)
    return tmp

def plot_contact_structure(plot, form = [0, 'x', 1], grid_size = 1.2, step = 0.5, alpha = 0.1, size = None):
    """
    Generates a grid and calles the plot_plane() on each gridpoint. It therefore initialtes the visualization of the CS

    :ax: A 3D matplotlib plot# x, y, z: Coordinates in R3
    :form: A one form. Per default it chooses the standard CS on R3
    :grid_size: Sets the extend of the grid where the CS is vizualized
    :step: Sets the step size of the grid
    :alp: Sets the transperancy of the planes indicating the CS
    """
    if size == None:
        size = grid_size / 10
    height_limit=0.3
    for x in np.arange(-grid_size, grid_size + step, step):
        for y in np.arange(-grid_size, grid_size + step, step):
            plot += _plot_plane(x, y, 0, form, alpha, plot_radius=size)
    return plot


def plot_knot(plot, knot, resolution=200, domain=[0, 1]):
    """
    Takes in a 3D matplotlib plot and a sympy equation of a knot and draws the knot onto the plot.

    :plot: 3D matplotlib plot
    :knot: Function that takes a sympy Symbol and returns a tuple of sympy expressions (x, y, z) or a list of such functions
    :resolution: Sets the number of points on which the knot equation is evaluated. Default value is 200
    :domain: The domain over which to evaluate the knot function. Default is [0, 1]
    """

    linspace = np.linspace(domain[0], domain[1], resolution)
    t = sp.symbols('t')  # Define the symbolic variable

    if callable(knot):  # Check if knot is a function
        x, y, z = [], [], []
        for tval in linspace:
            x_val, y_val, z_val = knot(t)
            x.append(x_val.subs(t, tval).evalf())
            y.append(y_val.subs(t, tval).evalf())
            z.append(z_val.subs(t, tval).evalf())
        vertices = np.vstack([x, y, z]).T
        plot += k3d.line(vertices)

    else:  # Assume knot is a list of functions
        for i in range(len(knot[0])):
            x, y, z = [], [], []
            for tval in linspace:
                x.append(knot[0][i].subs(t, tval).evalf())
                y.append(knot[1][i].subs(t, tval).evalf())
                z.append(knot[2][i].subs(t, tval).evalf())
            vertices = np.vstack([x, y, z]).T
            plot += k3d.line(vertices)

    # Plot the knot
    return plot



def plot_planes_along_knot(plot, knot, n=50, form = [0, 'x', 1], size = 0.1, height_limit = 0.3, alpha = 0.5, color = 'red', domain = [0,1] ):
    """
    Takes in a 3D matplotlib plot and a sympy equation of a knot and draws contact planes along the knot

    :ax: 3D matplotlib plot
    :knot: sympy exquation of knot
    :n: Sets the number of contact planes drawn. Default is 50
    :form: Sets the contact form that isto be visualized. Default is the standard CS on R3
    """

    # Differentiate between parameterized and drawn knots
    try:
        len(knot[0])
    except:
        t = sp.symbols('t')  # Define the symbolic variable
        linspace = np.linspace(domain[0], domain[1], n)
        x, y, z = [], [], []
        for tval in linspace:
            x.append(knot(t)[0].subs(t, tval).evalf())
            y.append(knot(t)[1].subs(t, tval).evalf())
            z.append(knot(t)[2].subs(t, tval).evalf())
        # Call plot_plane for these points
        for xi, yi, zi in zip(x,y,z):
            plot += _plot_plane(xi, yi, zi, form, alpha, surfcolor=color, plot_radius=size)
    else:
        linspace = np.linspace(domain[0], domain[1], n)
        t = sp.symbols('t')  # Define the symbolic variable
        for i in range(0,len(knot[0])):
            x, y, z = [], [], []
            for tval in linspace:
                x.append(knot[0][i].subs(t, tval).evalf())
                y.append(knot[1][i].subs(t, tval).evalf())
                z.append(knot[2][i].subs(t, tval).evalf())
            # Call plot_plane for these points
            for xi, yi, zi in zip(x,y,z):
                plot += _plot_plane(xi, yi, zi, form, alpha, surfcolor=color, plot_radius=size)

    return plot
