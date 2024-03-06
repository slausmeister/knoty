The goal of this Python Module is to provide an interactive framework for the visualization of Legendre knots. It generates interactive 3D plots of Legendre knots and contact structures.

# Installation
To install the knoty Module simply copy the "knoty" folder into the directory containing your Python file. Make sure, that all dependencies are installed.

## Dependencies
- Matplotlib
- Numpy
- Sympy
- SvgPathTools


# Basic Usage

First we load the dependencies and the module

```python
import knoty
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
```

As the goal of the package is to visualize structures, we initialize a 3D plot

```python
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

window_size = 1.5

ax.set_xlim([-window_size, window_size])
ax.set_ylim([-window_size, window_size])
ax.set_zlim([-window_size, window_size])
```

## Visualizing Contact Structures

We can specify a Contact Form in Cartesian coordinates using the notation

```python
form = [0, 'x', 1]
```

where the first entry describes the coefficient of the dx, the second entry the coefficient of dy and the third entry the coefficient of dz. We can also specify the three coordinates as coefficients using upticks. Therefore the expression above defines the form

$$xdy + dz$$

i.e. the standard contact form. We can now plot the contact structure via the knoty.plot_contact_structure() function. The function needs the 3D plot as an input. An example is given by the following:

```python
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

window_size = 1.5

ax.set_xlim([-window_size, window_size])
ax.set_ylim([-window_size, window_size])
ax.set_zlim([-window_size, window_size])

form = [0, 'x', 1]

knoty.plot_contact_structure(ax, form=form, alpha = 1)
```
![Alt text](/ReadmeImages/1.png?raw=true)
The argument form specifies the form. If no form is given, then the function defaults to the standard contact form. The argument alpha sets the transparency of the contact planes and defaults to 0.5.

## Visualizing Legendre Knots

There are two ways to specify Legendre knots in knoty. Either trough an explicit equation or though a drawn front projection. For drawn front projections please see "Specifying knots through SVGs".

To define a knot through an equation we first need to define a function for the knot as follows

```python
def knot(t):
    return 3 * sp.sin(t) * sp.cos(t), sp.cos(t), sp.sin(t)**3
```

Make sure that the return contains three coordinates and that the components are defined in a Sympy compatible manner. Here for example we choose sp.sin instead of np.sin to get the sympy version of the sine function.

We can now plot the knot through the knoty.plot_knot() function. This function need the 3D plot and the function defining the knot. An example is given by the following:

```python
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
knoty.plot_knot(ax, knot)
```
![Alt text](/ReadmeImages/3.png?raw=true)
We can also plot contact planes along the knot through the knoty.plot_planes_along_knot() function. It again takes in the 3D plot and the knot equation. Additionally we can layer everything onto a single plot. For example we could visualize the standard contact structure, plot the knot and indicate the contact planes in one 3D plot through the following code:

```python
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
knoty.plot_contact_structure(ax)
knoty.plot_knot(ax, knot)
knoty.plot_planes_along_knot(ax, knot)
```
![Alt text](/ReadmeImages/2.png?raw=true)
## Specifying knots through SVGs

A full example can be seen in the jupyter notebook attached. Example SVGs and a template can be found in Examples.

# TODO

- Full documentation for SVG import
- Calculating TB
- Fixing scaling issues in SVG import
- Fixing Jumping issues in sp.Piecewise in SVG import
- Publishing to pip
- Full in code documentation
