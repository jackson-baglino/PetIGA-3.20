#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def is_valid(candidate, circles, tol=1e-6):
    """Return True if candidate circle (x,y,r) does not overlap any in circles (touching allowed)."""
    x_c, y_c, r_c = candidate
    for (x, y, r) in circles:
        if np.hypot(x_c - x, y_c - y) < (r_c + r - tol):
            return False
    return True

def is_valid_inclusion(candidate, circles, gap=5e-6, tol=1e-6):
    """
    Return True if candidate air inclusion (x,y,r) is not touching or overlapping any
    inclusion in circles, enforcing a minimum gap.
    """
    x_c, y_c, r_c = candidate
    for (x, y, r) in circles:
        if np.hypot(x_c - x, y_c - y) < (r_c + r + gap - tol):
            return False
    return True

def generate_packing(num_circles, region, size_range, max_attempts, validate_func=is_valid):
    """
    Generate a circle packing in a rectangular region.
    region: (xmin, xmax, ymin, ymax)
    size_range: (r_min, r_max)
    validate_func: function to validate a candidate circle (default is is_valid)
    Returns a list of circles as tuples (x, y, r).
    """
    xmin, xmax, ymin, ymax = region
    circles = []
    # Place the first circle randomly inside the region
    for _ in range(int(max_attempts)):
        r = np.random.uniform(*size_range)
        x = np.random.uniform(xmin + r, xmax - r)
        y = np.random.uniform(ymin + r, ymax - r)
        circles.append((x, y, r))
        break

    iterations = 0
    while len(circles) < num_circles:
        iterations += 1
        if iterations > max_attempts:
            print("Max attempts reached in generate_packing, stopping.")
            break
        base = circles[np.random.randint(len(circles))]
        r_new = np.random.uniform(*size_range)
        angle = np.random.uniform(0, 2*np.pi)
        x_new = base[0] + (base[2] + r_new) * np.cos(angle)
        y_new = base[1] + (base[2] + r_new) * np.sin(angle)
        candidate = (x_new, y_new, r_new)
        # Ensure candidate is within the region
        if (x_new - r_new < xmin) or (x_new + r_new > xmax) or (y_new - r_new < ymin) or (y_new + r_new > ymax):
            continue
        if validate_func(candidate, circles):
            circles.append(candidate)
        else:
            # Alternatively, try a random candidate in the region
            x_rand = np.random.uniform(xmin + r_new, xmax - r_new)
            y_rand = np.random.uniform(ymin + r_new, ymax - r_new)
            candidate = (x_rand, y_rand, r_new)
            if validate_func(candidate, circles):
                circles.append(candidate)
    return circles

def generate_ice_grains(Lx, Ly, num_grains, size_range, max_attempts):
    """
    Generate densely packed ice grains above the midline.
    About 50% of the candidate grains are forced to have their bottom edge
    exactly touching the solid ice surface (y = Ly/2), ensuring that grains
    near the surface are actually in contact.
    """
    y_surface = Ly / 2.0
    grains = []
    attempts = 0
    while len(grains) < num_grains and attempts < max_attempts:
        attempts += 1
        # With 50% probability (or if no grain exists yet), force candidate to be on the surface.
        if np.random.rand() < 0.5 or len(grains) == 0:
            r = np.random.uniform(*size_range)
            x = np.random.uniform(r, Lx - r)
            y = y_surface + r  # Bottom edge exactly touches surface (y - r = y_surface)
            candidate = (x, y, r)
        else:
            # Otherwise, generate candidate based on an existing grain
            base = grains[np.random.randint(len(grains))]
            r = np.random.uniform(*size_range)
            angle = np.random.uniform(0, 2*np.pi)
            x = base[0] + (base[2] + r) * np.cos(angle)
            y = base[1] + (base[2] + r) * np.sin(angle)
            # If candidate would fall below the surface, force it to touch the surface.
            if y - r < y_surface:
                y = y_surface + r
            candidate = (x, y, r)
        # Check domain constraints
        if candidate[0] - candidate[2] < 0 or candidate[0] + candidate[2] > Lx:
            continue
        if candidate[1] - candidate[2] < y_surface or candidate[1] + candidate[2] > Ly:
            continue
        if is_valid(candidate, grains):
            grains.append(candidate)
    if len(grains) < num_grains:
        print(f"Warning: Only generated {len(grains)} out of {num_grains} ice grains after {attempts} attempts.")
    return grains

def generate_air_inclusions(Lx, Ly, num_inclusions, size_range, max_attempts):
    """
    Generate air inclusions in the bottom half of the domain.
    Ensures that inclusions do not touch or overlap.
    """
    region = (0, Lx, 0, Ly/2.0)
    return generate_packing(num_inclusions, region, size_range, max_attempts, validate_func=is_valid_inclusion)

def generate_nested_ice_inclusion(inclusion, num_nested, nested_size_range, max_attempts):
    """
    Given an air inclusion (x, y, R) that crosses the interface,
    generate nested ice grains along the bottom surface of the inclusion.
    """
    x0, y0, R = inclusion
    nested = []
    delta = 0.1 * R  # vertical band near the bottom of the inclusion
    for _ in range(num_nested):
        for _ in range(max_attempts):
            r = np.random.uniform(*nested_size_range)
            candidate_y = np.random.uniform(y0 - R + r, y0 - R + r + delta)
            max_dx = np.sqrt(max(0, (R - r)**2 - (candidate_y - y0)**2))
            candidate_x = np.random.uniform(x0 - max_dx, x0 + max_dx)
            candidate = (candidate_x, candidate_y, r)
            # Accept candidate if entirely inside inclusion and non-overlapping
            if np.hypot(candidate_x - x0, candidate_y - y0) + r <= R and is_valid(candidate, nested):
                nested.append(candidate)
                break
    return nested

def plot_circles(ice_grains, air_inclusions, nested_ice, Lx, Ly, title="Circle Packing"):
    fig, ax = plt.subplots()
    # Draw solid ice background (lower half)
    ax.add_patch(plt.Rectangle((0, 0), Lx, Ly/2, color='#d0f0ff', zorder=0))
    # Plot air inclusions first
    for (x, y, r) in air_inclusions:
        ax.add_patch(plt.Circle((x, y), r, edgecolor='black', facecolor='white', lw=0.5, zorder=1))
    # Plot ice grains and nested ice grains on top
    for (x, y, r) in ice_grains + nested_ice:
        ax.add_patch(plt.Circle((x, y), r, edgecolor='black', facecolor='#d0f0ff', lw=0.5, zorder=2))
    ax.set_aspect('equal')
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.axhline(Ly/2, color='gray', linestyle='--', lw=1)
    ax.set_title(title)
    output_path = f"./outputs/initial_domain_img/{title.replace(' ', '_').lower()}.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot as {output_path}")

def write_circle_data(filename, ice_grains, air_inclusions, nested_ice):
    """
    Write circle data (center and radius) to a CSV file.
    Each line contains: type,x,y,r.
    """
    with open(filename, "w") as f:
        # f.write("type,x,y,r\n")
        for (x, y, r) in ice_grains:
            f.write(f"{x} {y} {Lx/3.0} {r}\n")
        for (x, y, r) in air_inclusions:
            f.write(f"{x} {y} {Lx/3.0} {r}\n")
        for (x, y, r) in nested_ice:
            f.write(f"{x} {y} {Lx/3.0} {r}\n")
    print(f"Circle data written to {filename}")

if __name__ == '__main__':
    # Domain parameters
    Lx, Ly = 1.0e-3, 1.0e-3

    # Parameters for ice grains (top half)
    num_ice_grains = 300
    size_range_ice = (20e-6, 75e-6)

    # Parameters for air inclusions (bottom half)
    num_air_inclusions = 17
    size_range_air = (25e-6, 50e-6)

    max_attempts = int(1e7)

    # Generate ice grains in the top half
    ice_grains = generate_ice_grains(Lx, Ly, num_ice_grains, size_range_ice, max_attempts)

    # Generate air inclusions in the bottom half (ensuring they do not touch)
    air_inclusions = generate_air_inclusions(Lx, Ly, num_air_inclusions, size_range_air, max_attempts)

    # Generate nested ice grains for any air inclusion that crosses the interface
    nested_ice = []
    for inclusion in air_inclusions:
        x, y, r = inclusion
        if y > Ly/2 and y - r < Ly/2:  # Inclusion crosses midline and is centered above it
            nested = generate_nested_ice_inclusion(inclusion, num_nested=3, nested_size_range=(0.2, 0.4), max_attempts=max_attempts)
            nested_ice.extend(nested)

    # Plot the resulting circle packing
    plot_circles(ice_grains, air_inclusions, nested_ice, Lx, Ly, title="Ice Grains Pack Atop Solid Ice With Air Inclusions")
    
    # Write circle data to file
    write_circle_data("./outputs/initial_domain_img/circle_data.csv", ice_grains, air_inclusions, nested_ice)