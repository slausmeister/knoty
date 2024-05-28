import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from svgpathtools import svg2paths, CubicBezier, Line, Path

def _distance(point1, point2):
    return np.sqrt((point1.real - point2.real) ** 2 + (point1.imag - point2.imag) ** 2)

def _process_segment(segment):
    t = sp.symbols("t")
    if isinstance(segment, CubicBezier):
        P0_y, P0_z = segment.start.real, -segment.start.imag
        P1_y, P1_z = segment.control1.real, -segment.control1.imag
        P2_y, P2_z = segment.control2.real, -segment.control2.imag
        P3_y, P3_z = segment.end.real, -segment.end.imag

        y_t = P0_y * (1 - t)**3 + 3 * P1_y * t * (1 - t)**2 + 3 * P2_y * t**2 * (1 - t) + P3_y * t**3
        z_t = P0_z * (1 - t)**3 + 3 * P1_z * t * (1 - t)**2 + 3 * P2_z * t**2 * (1 - t) + P3_z * t**3
        dy_t = -3 * P0_y * (1 - t)**2 + 3 * P1_y * (1 - t)**2 -6 * P1_y * t * (1 - t) + 6 * P2_y * t * (1 - t)  - 3 * P2_y * t**2 + 3 * P3_y * t**2
        dz_t = -3 * P0_z * (1 - t)**2 + 3 * P1_z * (1 - t)**2 -6 * P1_z * t * (1 - t) + 6 * P2_z * t * (1 - t)  - 3 * P2_z * t**2 + 3 * P3_z * t**2

    elif isinstance(segment, Line):
        start_y, start_z = segment.start.real, -segment.start.imag
        end_y, end_z = segment.end.real, -segment.end.imag

        y_t = (1 - t) * start_y + t * end_y
        z_t = (1 - t) * start_z + t * end_z
        dy_t = end_y - start_y
        dz_t = end_z - start_z
    else:
        raise ValueError("Unsupported segment type")

    x_t = -dz_t / dy_t
    return x_t, y_t, z_t

def _extract_path_data(paths):
    first_points = []
    last_points = []
    num_of_pieces = 0
    max_left = np.inf
    max_right = -np.inf
    max_top = -np.inf
    max_bottom = np.inf

    for path in paths:
        first_points.append(path[0].start)
        last_points.append(path[-1].end)
        for segment in path:
            if segment.start.real < max_left: max_left = segment.start.real
            if segment.control1.real < max_left: max_left = segment.control1.real
            if segment.control2.real < max_left: max_left = segment.control2.real
            if segment.end.real < max_left: max_left = segment.end.real

            if segment.start.real > max_right: max_right = segment.start.real
            if segment.control1.real > max_right: max_right = segment.control1.real
            if segment.control2.real > max_right: max_right = segment.control2.real
            if segment.end.real > max_right: max_right = segment.end.real

            if segment.start.imag < max_bottom: max_bottom = segment.start.imag
            if segment.control1.imag < max_bottom: max_bottom = segment.control1.imag
            if segment.control2.imag < max_bottom: max_bottom = segment.control2.imag
            if segment.end.imag < max_bottom: max_bottom = segment.end.imag

            if segment.start.imag > max_top: max_top = segment.start.imag
            if segment.control1.imag > max_top: max_top = segment.control1.imag
            if segment.control2.imag > max_top: max_top = segment.control2.imag
            if segment.end.imag > max_top: max_top = segment.end.imag

            num_of_pieces += 1

    return first_points, last_points, num_of_pieces, max_left, max_right, max_top, max_bottom

def _apply_scaling(paths, scaling_factor, max_left, max_right, max_top, max_bottom):
    scalar = scaling_factor * 1 / max([max_bottom, max_top, max_left, max_right])
    hori_correction = 0.5 * (max_right - max_left)
    vert_correction = 0.5 * (max_top - max_bottom)

    for path in paths:
        for i in range(len(path)):
            segment = path[i]
            start = complex(scalar * (segment.start.real - hori_correction), scalar * (segment.start.imag - vert_correction))
            control1 = complex(scalar * (segment.control1.real - hori_correction), scalar * (segment.control1.imag - vert_correction))
            control2 = complex(scalar * (segment.control2.real - hori_correction), scalar * (segment.control2.imag - vert_correction))
            end = complex(scalar * (segment.end.real - hori_correction), scalar * (segment.end.imag - vert_correction))

            path[i] = CubicBezier(start, control1, control2, end)

def _apply_correction(first_points, last_points, threshold):
    for i in range(len(first_points)):
        for j in range(len(first_points)):
            if i == j:
                continue

            if _distance(first_points[i], first_points[j]) < threshold and _distance(first_points[i], first_points[j]) > 1e-6:
                avg_point = (first_points[i] + first_points[j]) / 2
                first_points[i], first_points[j] = avg_point, avg_point

            if _distance(last_points[i], last_points[j]) < threshold and _distance(last_points[i], last_points[j]) > 1e-6:
                avg_point = (last_points[i] + last_points[j]) / 2
                last_points[i], last_points[j] = avg_point, avg_point

            if _distance(first_points[i], last_points[j]) < threshold and _distance(first_points[i], last_points[j]) > 1e-6:
                avg_point = (first_points[i] + last_points[j]) / 2
                first_points[i], last_points[j] = avg_point, avg_point

            if _distance(last_points[i], first_points[j]) < threshold and _distance(last_points[i], first_points[j]) > 1e-6:
                avg_point = (last_points[i] + first_points[j]) / 2
                last_points[i], first_points[j] = avg_point, avg_point


def _create_2d_plot(y, z):
    t_vals = np.linspace(0, 1, 100)
    
    plt.figure(figsize=(10, 6))
    
    # Plot each segment individually
    for y_t, z_t in zip(y, z):
        y_vals = [float(y_t.subs('t', t_val)) for t_val in t_vals]
        z_vals = [float(z_t.subs('t', t_val)) for t_val in t_vals]
        plt.plot(y_vals, z_vals, color='black')

    plt.xlabel('y')
    plt.ylabel('z')
    plt.title('Imported front')
    plt.grid(True)
    plt.show()

def _update_paths_with_adjusted_points(paths, first_points, last_points):
    adjusted_paths = []
    path_index = 0

    for path in paths:
        new_path = []
        for i, segment in enumerate(path):
            is_first_segment = i == 0
            is_last_segment = i == len(path) - 1

            if is_first_segment:
                if isinstance(segment, CubicBezier):
                    new_control1 = complex(segment.control1.real, first_points[path_index].imag)
                    new_segment = CubicBezier(start=first_points[path_index], control1=new_control1, control2=segment.control2, end=segment.end if not is_last_segment else last_points[path_index])
                elif isinstance(segment, Line):
                    new_segment = Line(start=first_points[path_index], end=segment.end if not is_last_segment else last_points[path_index])
                new_path.append(new_segment)

            if is_last_segment and not is_first_segment:
                if isinstance(segment, CubicBezier):
                    new_control2 = complex(segment.control2.real, last_points[path_index].imag)
                    new_segment = CubicBezier(start=segment.start, control1=segment.control1, control2=new_control2, end=last_points[path_index])
                elif isinstance(segment, Line):
                    new_segment = Line(start=segment.start, end=last_points[path_index])
                new_path.append(new_segment)

            if not is_first_segment and not is_last_segment:
                new_path.append(segment)

        adjusted_paths.append(Path(*new_path))
        path_index += 1

    return adjusted_paths

def import_from_svg(svg_file, correction=True, threshold=0.1, scaling=False, scaling_factor=1, show_front=True):

    # We import the control points from the SVG file.
    paths, _ = svg2paths(svg_file)
    paths = [path for path in paths if path]

    # We extract the beginning and end of each path, the number of pieces in each path, and the maximum and minimum
    # values of y and z coordinates in the paths. The latter may be used for scaling.
    first_points, last_points, num_of_pieces, max_left, max_right, max_top, max_bottom = _extract_path_data(paths)

    # We scale the paths to ensure nicer knots. WIP
    if scaling:
        _apply_scaling(paths, scaling_factor, max_left, max_right, max_top, max_bottom)

    # We apply corrections to the paths to ensure that the paths are connected.
    if correction:
        _apply_correction(first_points, last_points, threshold)
        adjusted_paths = _update_paths_with_adjusted_points(paths, first_points, last_points)
    else:
        adjusted_paths = paths


    # We assemble the knot as a piecewise function of t.
    x, y, z = [], [], []
    for path in adjusted_paths:
        for segment in path:
            # We generate the Bezier equations from the control points.
            x_t, y_t, z_t = _process_segment(segment)
            y.append(y_t)
            z.append(z_t)
            x.append(x_t)

    # We plot the front of the knot. 
    if show_front:
        _create_2d_plot(y, z)

    return x, y, z