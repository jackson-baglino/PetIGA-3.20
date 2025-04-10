#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(21)

def is_valid_gap(candidate, circles, gap=1e-4):
    """
    Check if the candidate circle (x, y, r) is separated from all circles in
    'circles' by at least the sum of the two radii plus a given gap.
    """
    x_c, y_c, r_c = candidate
    for (x, y, r) in circles:
        if np.hypot(x_c - x, y_c - y) < (r_c + r + gap):
            return False
    return True

def generate_circle_packing(region, num_circles, size_range, Lx, Ly, max_attempts=100000, gap=1e-6, extra_circles=None):
    """
    Generate a packing of circles whose centers lie in the specified 'region'
    (a tuple: (xmin, xmax, ymin, ymax) for allowed centers) while ensuring:
      1. The entire circle is contained in the domain [0, Lx] x [0, Ly].
      2. The candidate does not overlap any previously placed circles in this region
         or in extra_circles (e.g. from the opposite half).
    
    Parameters:
      region       : tuple (xmin, xmax, ymin, ymax) for allowed center positions.
      num_circles  : desired number of circles to generate.
      size_range   : tuple (r_min, r_max) for circle radii.
      Lx, Ly       : overall domain dimensions.
      max_attempts : maximum number of candidate attempts.
      gap          : required gap between circles.
      extra_circles: list of circles (tuples) that should not be overlapped.
    
    Returns:
      A list of circles, each a tuple (x, y, r).
    """
    if extra_circles is None:
        extra_circles = []
    
    xmin, xmax, ymin, ymax = region
    circles = []
    attempts = 0

    while len(circles) < num_circles and attempts < max_attempts:
        attempts += 1
        # Choose a random radius for the candidate circle.
        r_new = np.random.uniform(*size_range)
        # When no circles have been placed yet in this region,
        # generate a candidate uniformly from the allowed center region,
        # but also ensuring the circle remains within the domain.
        if not circles:
            x = np.random.uniform(max(xmin, r_new), min(xmax, Lx - r_new))
            y = np.random.uniform(max(ymin, r_new), min(ymax, Ly - r_new))
            candidate = (x, y, r_new)
            if (x - r_new < 0) or (x + r_new > Lx) or (y - r_new < 0) or (y + r_new > Ly):
                continue
            if is_valid_gap(candidate, circles + extra_circles, gap):
                circles.append(candidate)
        else:
            base = circles[np.random.randint(len(circles))]
            angle = np.random.uniform(0, 2*np.pi)
            x_new = base[0] + (base[2] + r_new + gap) * np.cos(angle)
            y_new = base[1] + (base[2] + r_new + gap) * np.sin(angle)
            candidate = (x_new, y_new, r_new)
            # Check two conditions:
            # 1. The circle is completely inside the domain.
            # 2. The candidate's center lies within the allowed region.
            if (x_new - r_new < 0) or (x_new + r_new > Lx) or \
               (y_new - r_new < 0) or (y_new + r_new > Ly):
                valid_candidate = False
            elif (x_new < xmin) or (x_new > xmax) or (y_new < ymin) or (y_new > ymax):
                valid_candidate = False
            else:
                valid_candidate = True

            if valid_candidate and is_valid_gap(candidate, circles + extra_circles, gap):
                circles.append(candidate)
            else:
                # Fallback: attempt a candidate with a uniformly sampled center.
                x_rand = np.random.uniform(max(xmin, r_new), min(xmax, Lx - r_new))
                y_rand = np.random.uniform(max(ymin, r_new), min(ymax, Ly - r_new))
                candidate = (x_rand, y_rand, r_new)
                if (x_rand - r_new < 0) or (x_rand + r_new > Lx) or \
                   (y_rand - r_new < 0) or (y_rand + r_new > Ly):
                    continue
                if is_valid_gap(candidate, circles + extra_circles, gap):
                    circles.append(candidate)
    
    if len(circles) < num_circles:
        print(f"Warning: Only generated {len(circles)} circles after {attempts} attempts.")
    else:
        print(f"Successfully generated {len(circles)} circles after {attempts} attempts.")
    return circles

def save_vectorized_and_binary_images(top_circles, bottom_circles, Lx, Ly, resolution=1000):
    """
    Save two outputs:
      1. A vectorized image (SVG) of the domain with filled white circles on a black background.
      2. A binary raster image (PNG) where pixels inside any circle are 1 (white) and outside are 0 (black).
    """
    all_circles = top_circles + bottom_circles
    
    # ----- Vectorized Image (SVG) -----
    # Create a figure with black background.
    fig, ax = plt.subplots(figsize=(6,6))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    
    # Draw filled circles (white) for each circle.
    for (x, y, r) in all_circles:
        circle = plt.Circle((x, y), r, facecolor='white', edgecolor='none')
        ax.add_patch(circle)
    
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout(pad=0)
    svg_filename = "./outputs/RubyIceGrains.svg"
    plt.savefig(svg_filename, format="svg", bbox_inches='tight')
    print(f"Vectorized image saved as: {svg_filename}")
    plt.close(fig)
    
    # ----- Binary Raster Image -----
    # Create a binary mask by mapping the domain to a grid.
    nx = resolution
    ny = int(Ly / Lx) * resolution
    x_lin = np.linspace(0, Lx, nx)
    y_lin = np.linspace(0, Ly, ny)
    xv, yv = np.meshgrid(x_lin, y_lin)
    mask = np.zeros((ny, nx), dtype=np.uint8)
    
    # For each circle, mark pixels as 1 if they lie inside the circle.
    for (cx, cy, r) in all_circles:
        circle_mask = (xv - cx)**2 + (yv - cy)**2 <= r**2
        mask[circle_mask] = 1
    
    binary_filename = "./outputs/RubyIceGrains_binary.png"
    # Save the binary image (using 'gray' so that 1 becomes white and 0 black)
    plt.imsave(binary_filename, mask, cmap='gray', format='png')
    print(f"Binary image saved as: {binary_filename}")

if __name__ == '__main__':
  # Domain dimensions (for example, 10 cm x 10 cm)
  Lx, Ly = 10e-2, 20e-2

  # Allowed regions for centers:
  # Bottom half: centers must have y in [0, Ly/2]
  # Top half: centers must have y in [Ly/2, Ly]
  bottom_region =       (0, Lx, 0, Ly/2)
  top_region =          (0, Lx, Ly/2, Ly)

  # Size ranges (for instance, in meters)
  top_size_range =      (0.45e-3, 0.55e-3)    # Top-half circles.
  bottom_size_range =   (1.3e-3, 2.7e-3) # Bottom-half circles.

  # Desired numbers and parameters.
  num_bottom =          600    # Generate bottom-half circles first.
  num_top =             6000      # Then generate top-half circles.
  max_attempts =        int(5e5)
  gap_top =             0.25e-3        # Minimum gap between circles.
  gap_bottom =          0.5e-3       # Minimum gap between circles.

  # Generate bottom-half circles.
  bottom_circles = generate_circle_packing(bottom_region, num_bottom, bottom_size_range, Lx, Ly, max_attempts, gap_bottom)
  
  # Generate top-half circles, ensuring no overlap with the bottom circles.
  top_circles = generate_circle_packing(top_region, num_top, top_size_range, Lx, Ly, max_attempts, gap_top, extra_circles=bottom_circles)

  # Save both the vectorized and binary images.
  save_vectorized_and_binary_images(top_circles, bottom_circles, Lx, Ly, resolution=1000)