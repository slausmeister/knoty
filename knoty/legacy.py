import numpy as np
import sympy as sp
import matplotlib



def plot_plane_legacy(ax, x, y, z, form, size, height_limit, alpha, surfcolor='blue'):
    """
    This function draws a hyperplane at a point as to indicate the contact structure

    :ax: A 2D matplotlib plot# x, y, z: Coordinates in R3
    :form: A one form
    :size: Defines the extend of the plane indicating the contact structure
    :height_limit: Limits the vertical extend of the planes, as to make the picture less convoluted
    :surfacecolor: Sets the color of the plane
    :alpha: Sets the transperancy of the plane
    """
    v0, v2 = get_basis(form, x, y, z)

    # Create a grid on the plane
    u, v = np.meshgrid(np.linspace(-size, size, 9), np.linspace(-size, size, 10))
    plane_x = x + u * v0[0] + v * v2[0]
    plane_y = y + u * v0[1] + v * v2[1]
    plane_z = z + u * v0[2] + v * v2[2]


    plane_x = np.array(plane_x, dtype = float)
    plane_y = np.array(plane_y, dtype = float)
    plane_z = np.array(plane_z, dtype = float)


    # Clipping the z-values to be within the height_limit
    plane_z = np.clip(plane_z, z - height_limit/1, z + height_limit/2)
    

    # Plot the plane
    ax.plot_surface(plane_x, plane_y, plane_z, color=surfcolor, alpha=alpha)


def plot_knot_legacy(ax, knot, resolution=200, domain = [0,2 * np.pi]):
    """
    Takes in a 3D matplotlib plot and a sympy equation of a knot and draws the knot onto the plot.
    
    :ax: 3D matplotlib plot
    :knot: Function that takes a sympy Symbol and returns a tuple of sympy expressions (x, y, z)
    :resolution: Sets the number of points on which the knot equation is evaluated. Default value is 200
    """
    linspace = np.linspace(domain[0], domain[1], resolution)
    x, y, z = [], [], []
    
    # Differentiate between parameterized and drawn knots
    try:
        knot[0].free_symbols
    except:
        t = sp.symbols('t')  # Define the symbolic variable
        for tval in linspace:
            x.append(knot(t)[0].subs(t, tval).evalf())
            y.append(knot(t)[1].subs(t, tval).evalf())
            z.append(knot(t)[2].subs(t, tval).evalf())
    else:
        t = sp.symbols('t')  # Define the symbolic variable
        for tval in linspace:
            x.append(knot[0].subs(t, tval).evalf())
            y.append(knot[1].subs(t, tval).evalf())
            z.append(knot[2].subs(t, tval).evalf())
    
    # Plot the knot
    ax.plot(x, y, z, color='b')

def plot_planes_along_knot_legacy(ax, knot, n=50, form = [0, 'x', 1], size = 0.1, height_limit = 0.3, alpha = 0.5, color = 'red', domain = [0,2 * np.pi] ):
    """
    Takes in a 3D matplotlib plot and a sympy equation of a knot and draws contact planes along the knot

    :ax: 3D matplotlib plot
    :knot: sympy exquation of knot
    :n: Sets the number of contact planes drawn. Default is 50
    :form: Sets the contact form that isto be visualized. Default is the standard CS on R3
    """
    linspace = np.linspace(domain[0], domain[1], n)
    x, y, z = [], [], []

    # Differentiate between parameterized and drawn knots
    try:
        knot[0].free_symbols
    except:
        t = sp.symbols('t')  # Define the symbolic variable
        for tval in linspace:
            x.append(knot(t)[0].subs(t, tval).evalf())
            y.append(knot(t)[1].subs(t, tval).evalf())
            z.append(knot(t)[2].subs(t, tval).evalf())
    else:
        t = sp.symbols('t')  # Define the symbolic variable
        for tval in linspace:
            x.append(knot[0].subs(t, tval).evalf())
            y.append(knot[1].subs(t, tval).evalf())
            z.append(knot[2].subs(t, tval).evalf())

    # Call plot_plane for these points
    for xi, yi, zi in zip(x,y,z):
        plot_plane_legacy(ax, xi, yi, zi, form, size, height_limit, alpha, surfcolor=color)


def calculate_tb(knot):
    """
    WIP function. Takes in a knot and calculated its TB invariant. Currently broken.
    """
    # Define the symbols
    t1 = sp.symbols('t1', real=True)
    t2 = sp.symbols('t2', real=True)
    
    # Extract the y and z components of the knot function
    _, y1, z1 = knot(t1)
    _, y2, z2 = knot(t2)

    z2 = z2 + 0.1


    system = [z2-z1, y2-y1]


    
    # Solve the equations
    solutions = sp.solve(system)

    print(solutions)

    return solutions

def ot_in_cartesian(x_val, y_val, z_val):
    # Define symbolic variables
    x, y, z = sp.symbols('x y z')

    # Express r and theta in terms of x, y, z
    r = sp.sqrt(x**2 + y**2)
    theta = sp.atan2(y, x)

    # Calculate the derivatives of theta with respect to x and y
    dtheta_dx = sp.diff(theta, x)
    dtheta_dy = sp.diff(theta, y)

    # Components of the differential form in cylindrical coordinates
    dx_component = r * sp.sin(r) * dtheta_dx
    dy_component = r * sp.sin(r) * dtheta_dy
    dz_component = sp.cos(r)

    # Substitute specific values and evaluate each component
    dx_component_substituted = dx_component.subs({x: x_val, y: y_val, z: z_val}).evalf()
    dy_component_substituted = dy_component.subs({x: x_val, y: y_val, z: z_val}).evalf()
    dz_component_substituted = dz_component.subs({x: x_val, y: y_val, z: z_val}).evalf()

    return [dx_component_substituted, dy_component_substituted, dz_component_substituted]

def plot_contact_structure_cylindrical(ax, radius=1.2, radius_step=0.2, angle_step = np.pi / 4, alpha = 1):
    for r in np.arange(0.2, radius + radius_step, radius_step):
        for theta in np.arange(0, 2 * np.pi, angle_step):
            x = r * np.cos(theta)
            y = r * np.sin(theta)

            # Call the function with Cartesian coordinates
            ot_values = ot_in_cartesian(x, y, 0)

            size = 0.06
            height_limit = 0.3

            plot_plane_legacy(ax, x, y, 0, ot_values, size, height_limit, alpha)