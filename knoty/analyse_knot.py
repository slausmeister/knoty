import sympy as sp


def find_vanishing_cusps(knot):
    """
    WIP function. Takes in a Knot and calculates the coordinates of all cusps. Currently working.
    """
    t = sp.symbols('t')
    x, y, z = knot(t)

    # Differentiate y and z with respect to t
    dy_dt = sp.diff(y, t)
    dz_dt = sp.diff(z, t)

    # Solve for points where derivatives vanish
    vanishing_points_y = [sp.N(sol) for sol in sp.solve(dy_dt, t)]
    vanishing_points_z = [sp.N(sol) for sol in sp.solve(dz_dt, t)]

    # Find intersections considering the tolerance
    intersection = set(vanishing_points_y) & set(vanishing_points_z)

    x = [knot(t)[0].evalf() for t in intersection]
    y = [knot(t)[1].evalf() for t in intersection]
    z = [knot(t)[2].evalf() for t in intersection]

    return x, y, z
