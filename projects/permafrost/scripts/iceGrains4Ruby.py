#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

# Set a random seed for reproducibility.
np.random.seed(42)

def generate_circle_packing_kd(region, num_circles, size_range, Lx, Ly, 
  max_attempts=100000, gap=1e-6, extra_circles=None):
  """
  Generate a packing of circles (each represented as (x, y, r)) whose centers lie in
  the allowed 'region' using a cKDTree for fast neighbor queries.
  
  The entire circle must lie within the domain [0, Lx] x [0, Ly] and must not 
  overlap any previously accepted circles (or those in extra_circles) by at least
  the sum of the radii plus the given gap.
  
  Parameters:
    region        : (xmin, xmax, ymin, ymax) allowable center coordinates.
    num_circles   : desired number of circles to generate.
    size_range    : (r_min, r_max) tuple for circle radii.
    Lx, Ly        : overall domain dimensions.
    max_attempts  : maximum number of candidate attempts.
    gap           : enforced gap between circles.
    extra_circles : list of circles (tuples) that must not be overlapped.
    
  Returns:
    circles: list of tuples (x, y, r).
  """
  if extra_circles is None:
    extra_circles = []
  
  xmin, xmax, ymin, ymax = region
  circles = []       # Accepted circles (each as (x, y, r))
  centers = []       # List of accepted circle centers.
  radii = []         # List of accepted radii.

  # Build kd-tree data for extra circles if any.
  if extra_circles:
    extra_centers = np.array([[c[0], c[1]] for c in extra_circles])
    extra_radii = np.array([c[2] for c in extra_circles])
  else:
    extra_centers = np.empty((0, 2))
    extra_radii = np.empty((0,))

  attempts = 0
  while len(circles) < num_circles and attempts < max_attempts:
    attempts += 1
    r_new = np.random.uniform(*size_range)
    # Generate a candidate center uniformly within the allowed region,
    # ensuring that the circle lies completely within the domain.
    x = np.random.uniform(max(xmin, r_new), min(xmax, Lx - r_new))
    y = np.random.uniform(max(ymin, r_new), min(ymax, Ly - r_new))
    candidate = (x, y, r_new)
    
    # Verify candidate's circle is fully inside the domain.
    if (x - r_new < 0) or (x + r_new > Lx) or (y - r_new < 0) or (y + r_new > Ly):
      continue

    valid = True
    # Check against already accepted circles.
    if centers:
      # Rebuild kd-tree every time to include all accepted circles.
      kd_tree = cKDTree(np.array(centers))
      # Set search radius as candidate radius + maximum accepted radius + gap.
      search_radius = r_new + max(radii) + gap
      nearby_idx = kd_tree.query_ball_point([x, y], search_radius)
      for idx in nearby_idx:
        xi, yi = centers[idx]
        ri = radii[idx]
        if np.hypot(x - xi, y - yi) < (r_new + ri + gap):
          valid = False
          break
    # Check against extra circles if any.
    if valid and extra_centers.size:
      dists = np.hypot(extra_centers[:,0] - x, extra_centers[:,1] - y)
      if np.any(dists < (r_new + extra_radii + gap)):
        valid = False

    if valid:
      circles.append(candidate)
      centers.append([x, y])
      radii.append(r_new)
  
  if len(circles) < num_circles:
    print(f"Warning: Only generated {len(circles)} circles after {attempts} attempts.")
  else:
    print(f"Successfully generated {len(circles)} circles after {attempts} attempts.")
  return circles

def save_vectorized_and_binary_images(top_circles, bottom_circles, Lx, Ly, resolution=1000):
  """
  Save two outputs:
    1. A vectorized image (SVG) of the entire domain with filled white circles on a black background.
    2. A binary raster image (PNG) where pixels inside any circle are 1 (white) and outside are 0 (black).
  """
  all_circles = top_circles + bottom_circles
  
  # ----- Vectorized Image (SVG) -----
  fig, ax = plt.subplots(figsize=(6,6))
  fig.patch.set_facecolor('black')
  ax.set_facecolor('black')
  
  # Draw filled white circles for each generated circle.
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
  nx = resolution
  ny = int(Ly / Lx) * resolution
  x_lin = np.linspace(0, Lx, nx)
  y_lin = np.linspace(0, Ly, ny)
  xv, yv = np.meshgrid(x_lin, y_lin)
  mask = np.zeros((ny, nx), dtype=np.uint8)
  
  # Mark pixels inside each circle as 1.
  for (cx, cy, r) in all_circles:
    circle_mask = (xv - cx)**2 + (yv - cy)**2 <= r**2
    mask[circle_mask] = 1
  
  binary_filename = "./outputs/RubyIceGrains_binary.png"
  plt.imsave(binary_filename, mask, cmap='gray', format='png')
  print(f"Binary image saved as: {binary_filename}")

def save_all_circle_data(top_circles, bottom_circles, filename):
  """
  Save the data for all circles (positions and radii) to a CSV file.
  Each line contains: x, y, r, region
  """
  with open(filename, "w") as f:
    f.write("x,y,r,region\n")
    for circle in top_circles:
      f.write(f"{circle[0]},{circle[1]},{circle[2]},top\n")
    for circle in bottom_circles:
      f.write(f"{circle[0]},{circle[1]},{circle[2]},bottom\n")
  print(f"Circle data saved as: {filename}")

def save_pore_throat_data(all_circles, filename):
  """
  Compute the pore-throat for each circle and save the data to a CSV file.
  Here, the pore-throat is defined for each circle as the distance from its edge to the
  edge of its nearest neighbor (i.e. the Euclidean distance between centers minus the sum of
  the two radii). The output CSV has columns: circle_index, neighbor_index, pore_throat.
  """
  centers = np.array([[c[0], c[1]] for c in all_circles])
  radii = np.array([c[2] for c in all_circles])
  tree = cKDTree(centers)
  # Query the 2 nearest neighbors for each circle (first is self).
  dists, indices = tree.query(centers, k=2)
  with open(filename, "w") as f:
    f.write("circle_index,neighbor_index,pore_throat\n")
    for i in range(len(all_circles)):
      j = indices[i,1]  # nearest neighbor (other than itself)
      pore = dists[i,1] - (radii[i] + radii[j])
      f.write(f"{i},{j},{pore}\n")
  print(f"Pore-throat data saved as: {filename}")

if __name__ == '__main__':
  # Domain dimensions (for example, 10 cm x 20 cm)
  Lx, Ly = 10e-2, 20e-2

  # Allowed regions for centers:
  # Bottom half: centers have y in [0, Ly/2]
  # Top half: centers have y in [Ly/2, Ly]
  bottom_region = (0, Lx, 0, Ly/2)
  top_region = (0, Lx, Ly/2, Ly)

  # Size ranges (in meters)
  top_size_range = (0.4e-3, 0.6e-3)    # Top-half circles.
  bottom_size_range = (1.5e-3, 2.5e-3)  # Bottom-half circles.

  # Desired numbers and parameters.
  num_bottom = 900    # Generate bottom-half circles first.
  num_top = 9000      # Then generate top-half circles.
  max_attempts = int(1e6)
  gap_top = 0.05e-3        # Minimum gap between circles.
  gap_bottom = 0.15e-3     # Minimum gap between circles.

  print("Generating bottom-half circles...")
  bottom_circles = generate_circle_packing_kd(bottom_region, num_bottom, bottom_size_range, Lx, Ly, max_attempts, gap_bottom)
  
  print("Generating top-half circles (avoiding bottom circles)...")
  top_circles = generate_circle_packing_kd(top_region, num_top, top_size_range, Lx, Ly, max_attempts, gap_top, extra_circles=bottom_circles)
  
  # Save both the vectorized (SVG) and binary (PNG) images.
  save_vectorized_and_binary_images(top_circles, bottom_circles, Lx, Ly, resolution=1000)
  
  # Save circle data to CSV.
  save_all_circle_data(top_circles, bottom_circles, "./outputs/circle_data.csv")
  
  # Combine circles and compute & save pore-throat data.
  all_circles = top_circles + bottom_circles
  save_pore_throat_data(all_circles, "./outputs/pore_throat_data.csv")