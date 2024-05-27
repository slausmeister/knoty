import sympy as sp
import numpy as np


def unknot1(t):
    return 3 * sp.sin(2 * np.pi *t) * sp.cos(2 * np.pi *t), sp.cos(2 * np.pi *t), sp.sin(2 * np.pi *t)**3

def unknot2(t):
    return sp.cos(2 * np.pi *t), sp.sin(2 * np.pi *2*t), 2/3 * sp.sin(2 * np.pi *t)*sp.cos(2 * np.pi *2*t) -4/3 * sp.sin(2 * np.pi *2*t)*sp.cos(2 * np.pi *t)



# Lenny Ng's knot equations

def LNLegTrefoil(t):
    LegTref1_x = 3 * t
    LegTref1_y = 3/32 * (4 + 5 * t * (sp.sqrt(1 - t) - sp.sqrt(1 + t)) - 
                         3 * t**3 * (sp.sqrt(1 - t) - sp.sqrt(1 + t)) + 
                         2 * t**2 * (-2 + sp.sqrt(1 - t) + sp.sqrt(1 + t))) * \
                (-4 + ((1 - (1 - t)**(sp.Rational(3, 2))) * (2 + 3 * t - t**3) + (2 - 3 * t + t**3) * 
                (-1 + (1 + t)**(sp.Rational(3, 2))))**2)
    LegTref1_z = 1/16 * ((1 - (1 - t)**(sp.Rational(3, 2))) * (2 + 3 * t - t**3) + (2 - 3 * t + t**3) * 
                (-1 + (1 + t)**(sp.Rational(3, 2)))) * \
                (-12 + ((1 - (1 - t)**(sp.Rational(3, 2))) * (2 + 3 * t - t**3) + (2 - 3 * t + t**3) * 
                (-1 + (1 + t)**(sp.Rational(3, 2))))**2)

    LegTref3_y = -(9/4) * t * sp.sqrt(1 - t**2)
    LegTref3_z = 1 + 9/4 * (1 - t**2)**(sp.Rational(3, 2))

    x = sp.Piecewise((LegTref1_x.subs(t, t - 1), t <= 2), 
                     (LegTref1_x.subs(t, 3 - t), (t > 2) & (t <= 4)), 
                     (LegTref1_x.subs(t, t - 5), (t > 4) & (t <= 6)), 
                     (LegTref1_x.subs(t, 7 - t), t > 6))
    y = sp.Piecewise((LegTref1_y.subs(t, t - 1), t <= 2), 
                     (LegTref3_y.subs(t, 3 - t), (t > 2) & (t <= 4)), 
                     (-LegTref1_y.subs(t, t - 5), (t > 4) & (t <= 6)), 
                     (-LegTref3_y.subs(t, 7 - t), t > 6))
    z = sp.Piecewise((LegTref1_z.subs(t, t - 1), t <= 2), 
                     (LegTref3_z.subs(t, 3 - t), (t > 2) & (t <= 4)), 
                     (-LegTref1_z.subs(t, t - 5), (t > 4) & (t <= 6)), 
                     (-LegTref3_z.subs(t, 7 - t), t > 6))

    return x, y, z

def LNLegLHTrefoil(t):
    LegLHtref1_x = t
    LegLHtref1_y = 3/2 * (sp.sqrt(t) - 4 * (-1 + sp.sqrt(1 - t)) * t + (-4 + 11 * sp.sqrt(1 - t)) * t**2 - 7 * t**(sp.Rational(5, 2)) - 6 * sp.sqrt(1 - t) * t**3 + 6 * t**(sp.Rational(7, 2)))
    LegLHtref1_z = (1 - (1 - t)**(sp.Rational(3, 2))) * (3 * t**2 - 2 * t**3) + t**(sp.Rational(3, 2)) * (1 - 3 * t**2 + 2 * t**3)

    LegLHtref3_x = 2 - t
    LegLHtref3_y = -1.125 + 8.25 * sp.sqrt(3 - t) + 2.25 * sp.sqrt(-1 + t) + (1.5 - 14.625 * sp.sqrt(3 - t) - 8.625 * sp.sqrt(-1 + t)) * t + (-0.375 + 7.5 * sp.sqrt(3 - t) + 6 * sp.sqrt(-1 + t)) * t**2 + (-1.125 * sp.sqrt(3 - t) - 1.125 * sp.sqrt(-1 + t)) * t**3
    LegLHtref3_z = 1/4 * (-(0.5 + (3 - t)**(sp.Rational(3, 2))) * (-4 + t) * (-1 + t)**2 + (-3 + t)**2 * t * (1 - sp.sqrt(-1 + t) + sp.sqrt(-1 + t) * t))

    LegLHtref5_x = -4 + t
    LegLHtref5_y = 5.625 + 76.5 * sp.sqrt(5 - t) + 52.5 * sp.sqrt(-3 + t) + (-3 - 58.125 * sp.sqrt(5 - t) - 46.125 * sp.sqrt(-3 + t)) * t + (0.375 + 14.25 * sp.sqrt(5 - t) + 12.75 * sp.sqrt(-3 + t)) * t**2 + (-1.125 * sp.sqrt(5 - t) - 1.125 * sp.sqrt(-3 + t)) * t**3
    LegLHtref5_z = 1/4 * ((5 - t)**(sp.Rational(3, 2)) * (-6 + t) * (-3 + t)**2 + (0.5 - (-3 + t)**(sp.Rational(3, 2))) * (-5 + t)**2 * (-2 + t))

    LegLHtref6_x = 6 - t
    LegLHtref6_y = 3/2 * (120 - 1045 * sp.sqrt(6 - t) + 924 * sp.sqrt(-5 + t) + (-44 + 564 * sp.sqrt(6 - t) - 520 * sp.sqrt(-5 + t)) * t + (4 - 101 * sp.sqrt(6 - t) + 97 * sp.sqrt(-5 + t)) * t**2 + 6 * (sp.sqrt(6 - t) - sp.sqrt(-5 + t)) * t**3)
    LegLHtref6_z = -(1 - (6 - t)**(sp.Rational(3, 2))) * (-5 + t)**2 * (-13 + 2 * t) + (-6 + t)**2 * (-5 + t)**(sp.Rational(3, 2)) * (-9 + 2 * t)

    LegLHtref8_x = -6 + t
    LegLHtref8_y = 18 - 409.5 * sp.sqrt(8 - t) - 336 * sp.sqrt(-6 + t) + (-5.25 + 174 * sp.sqrt(8 - t) + 153 * sp.sqrt(-6 + t)) * t + (0.375 - 24.375 * sp.sqrt(8 - t) - 22.875 * sp.sqrt(-6 + t)) * t**2 + (1.125 * sp.sqrt(8 - t) + 1.125 * sp.sqrt(-6 + t)) * t**3
    LegLHtref8_z = 1/4 * (-(0.5 + (8 - t)**(sp.Rational(3, 2))) * (-9 + t) * (-6 + t)**2 + (1 + (-6 + t)**(sp.Rational(3, 2))) * (-8 + t)**2 * (-5 + t))

    LegLHtref10_x = 10 - t
    LegLHtref10_y = -30 - 864 * sp.sqrt(10 - t) - 742.5 * sp.sqrt(-8 + t) + (6.75 + 285 * sp.sqrt(10 - t) + 258 * sp.sqrt(-8 + t)) * t + (-0.375 - 31.125 * sp.sqrt(10 - t) - 29.625 * sp.sqrt(-8 + t)) * t**2 + (1.125 * sp.sqrt(10 - t) + 1.125 * sp.sqrt(-8 + t)) * t**3
    LegLHtref10_z = 1/4 * ((10 - t)**(sp.Rational(3, 2)) * (-11 + t) * (-8 + t)**2 + (0.5 - (-8 + t)**(sp.Rational(3, 2))) * (-10 + t)**2 * (-7 + t))

    x = sp.Piecewise((LegLHtref1_x, t <= 1), (LegLHtref3_x, t <= 3), (LegLHtref5_x, t <= 5), 
                     (LegLHtref6_x, t <= 6), (LegLHtref8_x, t <= 8), (LegLHtref10_x, True))
    y = sp.Piecewise((LegLHtref1_y, t <= 1), (LegLHtref3_y, t <= 3), (LegLHtref5_y, t <= 5),
                     (LegLHtref6_y, t <= 6), (LegLHtref8_y, t <= 8), (LegLHtref10_y, True))
    z = sp.Piecewise((LegLHtref1_z, t <= 1), (LegLHtref3_z, t <= 3), (LegLHtref5_z, t <= 5),
                     (LegLHtref6_z, t <= 6), (LegLHtref8_z, t <= 8), (LegLHtref10_z, True))

    return x, y, z




def LNChekanovA(t):
    Leg52A1_x = t
    Leg52A1_y = 3/2 * (sp.sqrt(t) - 4 * (-1 + sp.sqrt(1 - t)) * t + (-4 + 11 * sp.sqrt(1 - t)) * t**2 - 7 * t**(sp.Rational(5, 2)) - 6 * sp.sqrt(1 - t) * t**3 + 6 * t**(sp.Rational(7, 2)))
    Leg52A1_z = (1 - (1 - t)**(sp.Rational(3, 2))) * (3 - 2 * t) * t**2 + (-1 + t)**2 * t**(sp.Rational(3, 2)) * (1 + 2 * t)

    Leg52A2_x = 2 - t
    Leg52A2_y = 3/2 * (8 - 21 * sp.sqrt(2 - t) + 12 * sp.sqrt(-1 + t) + 4 * (-3 + 11 * sp.sqrt(2 - t) - 8 * sp.sqrt(-1 + t)) * t + (4 - 29 * sp.sqrt(2 - t) + 25 * sp.sqrt(-1 + t)) * t**2 + 6 * (sp.sqrt(2 - t) - sp.sqrt(-1 + t)) * t**3)
    Leg52A2_z = -(2 - (2 - t)**(sp.Rational(3, 2))) * (-1 + t)**2 * (-5 + 2 * t) + (1 + (-1 + t)**(sp.Rational(3, 2))) * (-2 + t)**2 * (-1 + 2 * t)

    Leg52A4_x = -2 + t
    Leg52A4_y = 3/8 * (-3 + t) * (4 * (7 * sp.sqrt(4 - t) + 4 * sp.sqrt(-2 + t)) - 4 * (5 * sp.sqrt(4 - t) + 4 * sp.sqrt(-2 + t)) * t + 3 * (sp.sqrt(4 - t) + sp.sqrt(-2 + t)) * t**2)
    Leg52A4_z = 1/4 * (-(2 + (4 - t)**(sp.Rational(3, 2))) * (-5 + t) * (-2 + t)**2 + (2 + (-2 + t)**(sp.Rational(3, 2))) * (-4 + t)**2 * (-1 + t))

    Leg52A6_x = 6 - t
    Leg52A6_y = -3/8 * (48 - 400 * sp.sqrt(6 - t) + 300 * sp.sqrt(-4 + t) + 20 * (-1 + 12 * sp.sqrt(6 - t) - 10 * sp.sqrt(-4 + t)) * t + (2 - 47 * sp.sqrt(6 - t) + 43 * sp.sqrt(-4 + t)) * t**2 + 3 * (sp.sqrt(6 - t) - sp.sqrt(-4 + t)) * t**3)
    Leg52A6_z = (1 + (6 - t)**(sp.Rational(3, 2))) * (1 - 3/4 * (6 - t)**2 + 1/4 * (6 - t)**3) + (3/4 * (6 - t)**2 - 1/4 * (6 - t)**3) * (2 - (-4 + t)**(sp.Rational(3, 2)))

    Leg52A9_x = -6 + t
    Leg52A9_y = -1/12 * (-15 + 2 * t) * (60 * sp.sqrt(9 - t) + 45 * sp.sqrt(-6 + t) - 2 * (8 * sp.sqrt(9 - t) + 7 * sp.sqrt(-6 + t)) * t + (sp.sqrt(9 - t) + sp.sqrt(-6 + t)) * t**2)
    Leg52A9_z = 1/54 * (54 + (9 - t)**(sp.Rational(3, 2)) * (-6 + t)**2 * (-21 + 2 * t) - (-9 + t)**2 * (-6 + t)**(sp.Rational(3, 2)) * (-9 + 2 * t))

    Leg52A11_x = 12 - t
    Leg52A11_y = 3/8 * (198 - 3150 * sp.sqrt(11 - t) + 2750 * sp.sqrt(-9 + t) + 5 * (-8 + 187 * sp.sqrt(11 - t) - 171 * sp.sqrt(-9 + t)) * t + (2 - 92 * sp.sqrt(11 - t) + 88 * sp.sqrt(-9 + t)) * t**2 + 3 * (sp.sqrt(11 - t) - sp.sqrt(-9 + t)) * t**3)
    Leg52A11_z = 1/4 * (-(2 - (11 - t)**(sp.Rational(3, 2))) * (-12 + t) * (-9 + t)**2 + (-11 + t)**2 * (-8 + t) * (1 - 9 * sp.sqrt(-9 + t) + sp.sqrt(-9 + t) * t))

    Leg52A13_x = -10 + t
    Leg52A13_y = 3/8 * (-12 + t) * (451 * sp.sqrt(13 - t) + 403 * sp.sqrt(-11 + t) - 2 * (37 * sp.sqrt(13 - t) + 35 * sp.sqrt(-11 + t)) * t + 3 * (sp.sqrt(13 - t) + sp.sqrt(-11 + t)) * t**2)
    Leg52A13_z = 1/4 * (-(2 + (13 - t)**(sp.Rational(3, 2))) * (-14 + t) * (-11 + t)**2 + (2 + (-11 + t)**(sp.Rational(3, 2))) * (-13 + t)**2 * (-10 + t))

    Leg52A14_x = 16 - t
    Leg52A14_y = -3/2 * (728 - 15093 * sp.sqrt(14 - t) + 14364 * sp.sqrt(-13 + t) + 4 * (-27 + 833 * sp.sqrt(14 - t) - 806 * sp.sqrt(-13 + t)) * t + (4 - 245 * sp.sqrt(14 - t) + 241 * sp.sqrt(-13 + t)) * t**2 + 6 * (sp.sqrt(14 - t) - sp.sqrt(-13 + t)) * t**3)
    Leg52A14_z = -(1 + (14 - t)**(sp.Rational(3, 2))) * (-13 + t)**2 * (-29 + 2 * t) + (2 - (-13 + t)**(sp.Rational(3, 2))) * (-14 + t)**2 * (-25 + 2 * t)

    Leg52A15_x = -12 + t
    Leg52A15_y = 3/2 * (sp.sqrt(15 - t) * (-14 + t)**2 * (-31 + 2 * t) - (-15 + t)**2 * sp.sqrt(-14 + t) * (-27 + 2 * t) + 4 * (1 - (-14 + t)**(sp.Rational(3, 2))) * (210 - 29 * t + t**2) - 4 * (15 - t)**(sp.Rational(3, 2)) * (210 - 29 * t + t**2))
    Leg52A15_z = -(15 - t)**(sp.Rational(3, 2)) * (-14 + t)**2 * (-31 + 2 * t) + (1 - (-14 + t)**(sp.Rational(3, 2))) * (-15 + t)**2 * (-27 + 2 * t)

    Leg52A18_x = 18 - t
    Leg52A18_y = 1/12 * (-33 + 2 * t) * (3 * (95 * sp.sqrt(18 - t) + 84 * sp.sqrt(-15 + t)) - 2 * (17 * sp.sqrt(18 - t) + 16 * sp.sqrt(-15 + t)) * t + (sp.sqrt(18 - t) + sp.sqrt(-15 + t)) * t**2)
    Leg52A18_z = 1/54 * ((18 - t)**(sp.Rational(3, 2)) * (-15 + t)**2 * (-39 + 2 * t) - (-18 + t)**2 * (-15 + t)**(sp.Rational(3, 2)) * (-27 + 2 * t))

    x = sp.Piecewise(
        (Leg52A1_x, t <= 1),
        (Leg52A2_x, t <= 2),
        (Leg52A4_x, t <= 4),
        (Leg52A6_x, t <= 6),
        (Leg52A9_x, t <= 9),
        (Leg52A11_x, t <= 11),
        (Leg52A13_x, t <= 13),
        (Leg52A14_x, t <= 14),
        (Leg52A15_x, t <= 15),
        (Leg52A18_x, True)
    )

    y = sp.Piecewise(
        (Leg52A1_y, t <= 1),
        (Leg52A2_y, t <= 2),
        (Leg52A4_y, t <= 4),
        (Leg52A6_y, t <= 6),
        (Leg52A9_y, t <= 9),
        (Leg52A11_y, t <= 11),
        (Leg52A13_y, t <= 13),
        (Leg52A14_y, t <= 14),
        (Leg52A15_y, t <= 15),
        (Leg52A18_y, True)
    )

    z = sp.Piecewise(
        (Leg52A1_z, t <= 1),
        (Leg52A2_z, t <= 2),
        (Leg52A4_z, t <= 4),
        (Leg52A6_z, t <= 6),
        (Leg52A9_z, t <= 9),
        (Leg52A11_z, t <= 11),
        (Leg52A13_z, t <= 13),
        (Leg52A14_z, t <= 14),
        (Leg52A15_z, t <= 15),
        (Leg52A18_z, True)
    )

    return x, y, z


def LNChekanovB(t):
    Leg52B1_x = t
    Leg52B1_y = 3/2 * (sp.sqrt(t) - 4 * (-1 + sp.sqrt(1 - t)) * t + (-4 + 11 * sp.sqrt(1 - t)) * t**2 - 7 * t**(sp.Rational(5, 2)) - 6 * sp.sqrt(1 - t) * t**3 + 6 * t**(sp.Rational(7, 2)))
    Leg52B1_z = (1 - (1 - t)**(sp.Rational(3, 2))) * (3 - 2 * t) * t**2 + (-1 + t)**2 * t**(sp.Rational(3, 2)) * (1 + 2 * t)

    Leg52B2_x = 2 - t
    Leg52B2_y = 3/2 * (8 - 21 * sp.sqrt(2 - t) + 12 * sp.sqrt(-1 + t) + 4 * (-3 + 11 * sp.sqrt(2 - t) - 8 * sp.sqrt(-1 + t)) * t + (4 - 29 * sp.sqrt(2 - t) + 25 * sp.sqrt(-1 + t)) * t**2 + 6 * (sp.sqrt(2 - t) - sp.sqrt(-1 + t)) * t**3)
    Leg52B2_z = -(2 - (2 - t)**(sp.Rational(3, 2))) * (-1 + t)**2 * (-5 + 2 * t) + (1 + (-1 + t)**(sp.Rational(3, 2))) * (-2 + t)**2 * (-1 + 2 * t)

    Leg52B4_x = -2 + t
    Leg52B4_y = 3/8 * (-3 + t) * (4 * (7 * sp.sqrt(4 - t) + 4 * sp.sqrt(-2 + t)) - 4 * (5 * sp.sqrt(4 - t) + 4 * sp.sqrt(-2 + t)) * t + 3 * (sp.sqrt(4 - t) + sp.sqrt(-2 + t)) * t**2)
    Leg52B4_z = 1/4 * (-(2 + (4 - t)**(sp.Rational(3, 2))) * (-5 + t) * (-2 + t)**2 + (2 + (-2 + t)**(sp.Rational(3, 2))) * (-4 + t)**2 * (-1 + t))

    Leg52B6_x = 6 - t
    Leg52B6_y = -3/8 * (48 - 400 * sp.sqrt(6 - t) + 300 * sp.sqrt(-4 + t) + 20 * (-1 + 12 * sp.sqrt(6 - t) - 10 * sp.sqrt(-4 + t)) * t + (2 - 47 * sp.sqrt(6 - t) + 43 * sp.sqrt(-4 + t)) * t**2 + 3 * (sp.sqrt(6 - t) - sp.sqrt(-4 + t)) * t**3)
    Leg52B6_z = (1 + (6 - t)**(sp.Rational(3, 2))) * (1 - 3/4 * (6 - t)**2 + 1/4 * (6 - t)**3) + (3/4 * (6 - t)**2 - 1/4 * (6 - t)**3) * (2 - (-4 + t)**(sp.Rational(3, 2)))

    Leg52B7_x = -6 + t
    Leg52B7_y = 3/2 * (sp.sqrt(7 - t) * (-6 + t)**2 * (-15 + 2 * t) - (-7 + t)**2 * sp.sqrt(-6 + t) * (-11 + 2 * t) + 4 * (1 - (-6 + t)**(sp.Rational(3, 2))) * (42 - 13 * t + t**2) - 4 * (7 - t)**(sp.Rational(3, 2)) * (42 - 13 * t + t**2))
    Leg52B7_z = -(7 - t)**(sp.Rational(3, 2)) * (-6 + t)**2 * (-15 + 2 * t) + (1 - (-6 + t)**(sp.Rational(3, 2))) * (-7 + t)**2 * (-11 + 2 * t)

    Leg52B8_x = 8 - t
    Leg52B8_y = 3/2 * (-sp.sqrt(8 - t) * (-7 + t)**2 * (-17 + 2 * t) + (-8 + t)**2 * sp.sqrt(-7 + t) * (-13 + 2 * t) + 4 * (-1 + (8 - t)**(sp.Rational(3, 2))) * (56 - 15 * t + t**2) + 4 * (-7 + t)**(sp.Rational(3, 2)) * (56 - 15 * t + t**2))
    Leg52B8_z = -(-1 + (8 - t)**(sp.Rational(3, 2))) * (-7 + t)**2 * (-17 + 2 * t) - (-8 + t)**2 * (-7 + t)**(sp.Rational(3, 2)) * (-13 + 2 * t)

    Leg52B11_x = -8 + t
    Leg52B11_y = 1/18 * (-704 + 5472 * sp.sqrt(11 - t) + 4389 * sp.sqrt(-8 + t) - 4 * (-38 + 429 * sp.sqrt(11 - t) + 372 * sp.sqrt(-8 + t)) * t + (-8 + 177 * sp.sqrt(11 - t) + 165 * sp.sqrt(-8 + t)) * t**2 - 6 * (sp.sqrt(11 - t) + sp.sqrt(-8 + t)) * t**3)
    Leg52B11_z = 1/27 * (-(1 - (11 - t)**(sp.Rational(3, 2))) * (-8 + t)**2 * (-25 + 2 * t) + (-1 - (-8 + t)**(sp.Rational(3, 2))) * (-11 + t)**2 * (-13 + 2 * t))

    Leg52B13_x = 14 - t
    Leg52B13_y = 3/8 * (286 - 5412 * sp.sqrt(13 - t) + 4836 * sp.sqrt(-11 + t) + (-48 + 1339 * sp.sqrt(13 - t) - 1243 * sp.sqrt(-11 + t)) * t + (2 - 110 * sp.sqrt(13 - t) + 106 * sp.sqrt(-11 + t)) * t**2 + 3 * (sp.sqrt(13 - t) - sp.sqrt(-11 + t)) * t**3)
    Leg52B13_z = 1/4 * (-(2 - (13 - t)**(sp.Rational(3, 2))) * (-14 + t) * (-11 + t)**2 + (1 + (-11 + t)**(sp.Rational(3, 2))) * (-13 + t)**2 * (-10 + t))

    Leg52B15_x = -12 + t
    Leg52B15_y = 3/8 * (-14 + t) * (611 * sp.sqrt(15 - t) + 555 * sp.sqrt(-13 + t) - 2 * (43 * sp.sqrt(15 - t) + 41 * sp.sqrt(-13 + t)) * t + 3 * (sp.sqrt(15 - t) + sp.sqrt(-13 + t)) * t**2)
    Leg52B15_z = 1/4 * (-(2 + (15 - t)**(sp.Rational(3, 2))) * (-16 + t) * (-13 + t)**2 + (2 + (-13 + t)**(sp.Rational(3, 2))) * (-15 + t)**2 * (-12 + t))

    Leg52B18_x = 18 - t
    Leg52B18_y = 1/18 * (-27 * (80 + 1045 * sp.sqrt(18 - t) + 924 * sp.sqrt(-15 + t)) + 12 * (22 + 423 * sp.sqrt(18 - t) + 390 * sp.sqrt(-15 + t)) * t - (8 + 303 * sp.sqrt(18 - t) + 291 * sp.sqrt(-15 + t)) * t**2 + 6 * (sp.sqrt(18 - t) + sp.sqrt(-15 + t)) * t**3)
    Leg52B18_z = 1/27 * ((18 - t)**(sp.Rational(3, 2)) * (-15 + t)**2 * (-39 + 2 * t) + (2 - (-15 + t)**(sp.Rational(3, 2))) * (-18 + t)**2 * (-27 + 2 * t))

    x = sp.Piecewise(
        (Leg52B1_x, t <= 1),
        (Leg52B2_x, t <= 2),
        (Leg52B4_x, t <= 4),
        (Leg52B6_x, t <= 6),
        (Leg52B7_x, t <= 7),
        (Leg52B8_x, t <= 8),
        (Leg52B11_x, t <= 11),
        (Leg52B13_x, t <= 13),
        (Leg52B15_x, t <= 15),
        (Leg52B18_x, True)
    )

    y = sp.Piecewise(
        (Leg52B1_y, t <= 1),
        (Leg52B2_y, t <= 2),
        (Leg52B4_y, t <= 4),
        (Leg52B6_y, t <= 6),
        (Leg52B7_y, t <= 7),
        (Leg52B8_y, t <= 8),
        (Leg52B11_y, t <= 11),
        (Leg52B13_y, t <= 13),
        (Leg52B15_y, t <= 15),
        (Leg52B18_y, True)
    )

    z = sp.Piecewise(
        (Leg52B1_z, t <= 1),
        (Leg52B2_z, t <= 2),
        (Leg52B4_z, t <= 4),
        (Leg52B6_z, t <= 6),
        (Leg52B7_z, t <= 7),
        (Leg52B8_z, t <= 8),
        (Leg52B11_z, t <= 11),
        (Leg52B13_z, t <= 13),
        (Leg52B15_z, t <= 15),
        (Leg52B18_z, True)
    )

    return x, y, z
