import numpy as np
import math
import sympy as sp
from svgpathtools import svg2paths, CubicBezier, QuadraticBezier, Line, Path
import matplotlib.pyplot as plt


def distance(point1, point2):
    return np.sqrt((point1.real - point2.real) ** 2 + (point1.imag - point2.imag) ** 2)


# Function to process each segment
def process_segment(segment, num_of_pieces, current_piece):
    t = sp.symbols("t")
    if isinstance(segment, CubicBezier):
        # Extract real and imaginary parts of control points
        P0_y, P0_z = segment.start.real, -segment.start.imag
        P1_y, P1_z = segment.control1.real, -segment.control1.imag
        P2_y, P2_z = segment.control2.real, -segment.control2.imag
        P3_y, P3_z = segment.end.real, -segment.end.imag

        # Parametric equations for y(t) and z(t)
        y_t = P0_y * (1 - t)**3 + 3 * P1_y * t * (1 - t)**2 + 3 * P2_y * t**2 * (1 - t) + P3_y * t**3
        z_t = P0_z * (1 - t)**3 + 3 * P1_z * t * (1 - t)**2 + 3 * P2_z * t**2 * (1 - t) + P3_z * t**3
        dy_t = -3 * P0_y * (1 - t)**2 + 3 * P1_y * (1 - t)**2 -6 * P1_y * t * (1 - t) + 6 * P2_y * t * (1 - t)  - 3 * P2_y * t**2 + 3 * P3_y * t**2
        dz_t = -3 * P0_z * (1 - t)**2 + 3 * P1_z * (1 - t)**2 -6 * P1_z * t * (1 - t) + 6 * P2_z * t * (1 - t)  - 3 * P2_z * t**2 + 3 * P3_z * t**2

        y_values = [P0_y, P1_y, P2_y, P3_y]
        z_values = [P0_z, P1_z, P2_z, P3_z]

    elif isinstance(segment, Line):
        # Linear interpolation for y(t) and z(t)
        start_y, start_z = segment.start.real, -segment.start.imag
        end_y, end_z = segment.end.real, -segment.end.imag
        
        y_t = (1 - t) * start_y + t * end_y
        z_t = (1 - t) * start_z + t * end_z
        dy_t = end_y - start_y
        dz_t = end_z - start_z
    else:
        raise ValueError("Unsupported segment type")

    x_t = -dz_t / dy_t

    x_t = x_t.subs(t, num_of_pieces * t - current_piece)
    y_t = y_t.subs(t, num_of_pieces * t - current_piece)
    z_t = z_t.subs(t, num_of_pieces * t - current_piece)

    # Debug
    # linspace = np.linspace(0, 1/num_of_pieces, 100)
    # y_values_func = [y_t.subs(t, tval).evalf() for tval in linspace]
    # z_values_func = [z_t.subs(t, tval).evalf() for tval in linspace]
    # plt.figure(figsize=(8, 5))  # Optional: Adjusts the size of the figure
    # plt.plot(y_values_func, z_values_func, label='Cubic Bezier Curve')
    # plt.show()
    return x_t, y_t, z_t

def import_from_svg(svg_file, correction=True, threshold=0.1):
    paths, _ = svg2paths(svg_file)

    paths = [path for path in paths if path]

    first_points = []
    last_points = []
    curve_data = []
    num_of_pieces = 0
    t = sp.symbols('t')

    # Extract first and last points of each path, including placeholders for empty paths
    for path in paths:
        first_points.append(path[0].start)
        last_points.append(path[-1].end)
        for segment in path:
            num_of_pieces += 1

    if correction:
        # Compare and adjust similar start and end points between paths
        for i in range(len(first_points)):
            for j in range(len(first_points)):
                if i == j:
                    continue  # Skip comparing the same path or with empty paths

                # Adjust start points with other start points
                if distance(first_points[i], first_points[j]) < threshold and distance(first_points[i], first_points[j]) > 1e-6:
                    avg_point = (first_points[i] + first_points[j]) / 2
                    print("Correcting start point", first_points[i].conjugate(), " and start point", first_points[j].conjugate(), " to ", avg_point.conjugate())
                    first_points[i], first_points[j] = avg_point, avg_point

                # Adjust end points with other end points
                if distance(last_points[i], last_points[j]) < threshold and distance(last_points[i], last_points[j]) > 1e-6:
                    avg_point = (last_points[i] + last_points[j]) / 2
                    print("Correcting end point", last_points[i].conjugate(), " and end point", last_points[j].conjugate(), " to ", avg_point.conjugate())
                    last_points[i], last_points[j] = avg_point, avg_point

                # Adjust start points with other end points
                if distance(first_points[i], last_points[j]) < threshold and distance(first_points[i], last_points[j]) > 1e-6:
                    avg_point = (first_points[i] + last_points[j]) / 2
                    print("Correcting start point", first_points[i].conjugate(), " and end point", last_points[j].conjugate(), " to ", avg_point.conjugate())
                    first_points[i], last_points[j] = avg_point, avg_point

                # Adjust end points with other start points
                if distance(last_points[i], first_points[j]) < threshold and distance(last_points[i], first_points[j]) > 1e-6:
                    avg_point = (last_points[i] + first_points[j]) / 2
                    print("Correcting end point", last_points[i].conjugate(), " and start point", first_points[j].conjugate(), " to ", avg_point.conjugate())
                    last_points[i], first_points[j] = avg_point, avg_point

        # Update paths with adjusted points and control points
        adjusted_paths = []
        path_index = 0  # Track index for first_points and last_points
        for path in paths:
            new_path = []
            for i, segment in enumerate(path):
                is_first_segment = i == 0
                is_last_segment = i == len(path) - 1
        
                if is_first_segment:
                    # Adjusting the first segment of each path
                    if isinstance(segment, CubicBezier):
                        new_control1 = complex(segment.control1.real, first_points[path_index].imag)
                        # print("First adjusted: ", first_points[path_index])
                        new_segment = CubicBezier(start=first_points[path_index], 
                                                  control1=new_control1, 
                                                  control2=segment.control2, 
                                                  end=segment.end if not is_last_segment else last_points[path_index])
                    elif isinstance(segment, Line):
                        new_segment = Line(start=first_points[path_index], end=segment.end if not is_last_segment else last_points[path_index])
                    new_path.append(new_segment)
        
                if is_last_segment and not is_first_segment:
                    # Adjusting the last segment of each path, but only if it's not also the first segment
                    if isinstance(segment, CubicBezier):
                        new_control2 = complex(segment.control2.real, last_points[path_index].imag)
                        # print("Last adjusted: ", last_points[path_index])
                        new_segment = CubicBezier(start=segment.start, 
                                                  control1=segment.control1, 
                                                  control2=new_control2, 
                                                  end=last_points[path_index])  # Adjusted end point
                    elif isinstance(segment, Line):
                        new_segment = Line(start=segment.start, end=last_points[path_index])  # Adjusted end point
                    new_path.append(new_segment)
        
                if not is_first_segment and not is_last_segment:
                    # For all other segments
                    new_path.append(segment)
        
            adjusted_paths.append(Path(*new_path))
            path_index += 1
        
    else:
        adjusted_paths = paths

    # Clear curve_data before processing adjusted paths
    curve_data.clear()
    # Process adjusted paths

    x = []
    y = []
    z = []
    current_piece = 0
    for path in paths:
        for segment in path:
            x_t, y_t, z_t = process_segment(segment, num_of_pieces, current_piece)
            y.append(y_t)
            z.append(z_t)
            x.append(x_t)
            current_piece += 1

    tmp_x = []
    tmp_y = []
    tmp_z = []

    for piece in range(0, num_of_pieces):
        # Use sympy's And to define the condition for each piece
        condition = sp.And(t > piece/num_of_pieces, t < (piece + 1) / num_of_pieces)
        tmp_x.append((x[piece], condition))
        tmp_y.append((y[piece], condition))
        tmp_z.append((z[piece], condition))

    # Now directly create the Piecewise function
    piecewise_x = sp.Piecewise(*tmp_x)
    piecewise_y = sp.Piecewise(*tmp_y)
    piecewise_z = sp.Piecewise(*tmp_z)

    return piecewise_x, piecewise_y, piecewise_z