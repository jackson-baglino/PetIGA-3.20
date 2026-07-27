import numpy as np
import matplotlib.pyplot as plt

def generate_fcc_with_all_extended_edges(radius, N_grains_x):
    """
    Generate 2D FCC packing in a square domain.
    Grains touch all four domain edges, with offset rows extending past left/right,
    and the top row extending past the top edge.

    Parameters:
        radius (float): Grain radius
        N_grains_x (int): Number of grains in aligned (even) rows

    Returns:
        positions (np.ndarray): (x, y) coordinates of grain centers
        domain_length (float): Side length of the square domain
    """
    dx = 2 * radius
    dy = np.sqrt(3) * radius

    domain_length = (N_grains_x - 1) * dx

    # Estimate number of rows needed so that last row touches or exceeds top boundary
    N_rows = 0
    y = 0
    while y + radius <= domain_length:
        N_rows += 1
        y = N_rows * dy

    positions = []
    for j in range(N_rows + 1):  # Add one extra row to ensure overshoot
        y = j * dy
        if j % 2 == 0:
            for i in range(N_grains_x):
                x = i * dx
                positions.append((x, y))
        else:
            for i in range(N_grains_x + 1):
                x = -radius + i * dx
                positions.append((x, y))

    return np.array(positions), domain_length

# Example usage
if __name__ == "__main__":
    radius = 100.0e-6
    N_grains_x = 4

    positions, domain_length = generate_fcc_with_all_extended_edges(radius, N_grains_x)

    print(f"Square domain side length: {domain_length:.2f}")
    print(f"Total grains placed: {len(positions)}")

    # Plotting
    fig, ax = plt.subplots()
    for x, y in positions:
        circle = plt.Circle((x, y), radius, edgecolor='black', facecolor='lightblue')
        ax.add_patch(circle)

    ax.set_aspect('equal')
    buffer = 2 * radius
    ax.set_xlim(-radius - buffer/2, domain_length + radius + buffer/2)
    ax.set_ylim(-radius, domain_length + radius + buffer/2)
    ax.set_title("2D FCC Packing Touching All Domain Edges")
    plt.xlabel("x")
    plt.ylabel("y")

    # Draw domain box
    domain = plt.Rectangle((0, 0), domain_length, domain_length,
                           linewidth=1.5, edgecolor='red', facecolor='none', linestyle='--')
    ax.add_patch(domain)

    plt.grid(True)
    plt.show()

        # Save to .dat file
    z_center = 2 * radius
    output_file = f"./inputs/grainReadFile-{len(positions)}FCC.dat"

    with open(output_file, 'w') as f:
        for x, y in positions:
            f.write(f"{x:.6f} {y:.6f} {z_center:.6f} {radius:.6f}\n")

    print(f"Grain data saved to '{output_file}' with format: x y z radius")
    print(f"Domain length: {domain_length:.6f} m")