import numpy as np

def generate_circle_phase_field(Lx, Ly, Nx, Ny, radius, filename="circle_phase_field.dat"):
    """
    Generates a phase field for a circle centered in a rectangular domain.
    
    Parameters:
    - Lx, Ly: Physical domain size in x and y directions (meters).
    - Nx, Ny: Number of grid points in x and y directions.
    - radius: Radius of the circle (meters).
    - epsilon: Interface width for smooth transition (meters).
    - filename: Output file name for saving the phase field.
    """
    
    # Grid spacing
    dx = Lx / Nx
    dy = Ly / Ny

    epsilon = min(dx, dy)
    
    # Create grid points (including boundaries, hence Nx+1 and Ny+1)
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    
    # Center of the domain
    xc, yc = Lx / 2, Ly / 2
    
    # Compute the phase field
    phase_field = np.zeros((Ny, Nx))
    for j in range(Ny):
        for i in range(Nx):
            d = np.sqrt((x[i] - xc)**2 + (y[j] - yc)**2) - radius
            phase_field[j, i] = 0.5 * (1 - np.tanh(2.0 * d / epsilon))  # Smooth transition using tanh
    
    # Save to .dat file
    with open(filename, "w") as f:
        for row in phase_field:
            for value in row:
                f.write(f"{value:.6f}\n")  # Write each value on a new line
    
    print(f"✅ Phase field saved to {filename}")

# Example Usage:
generate_circle_phase_field(
    Lx=0.5e-3,  # 0.5 mm domain
    Ly=0.5e-3,  # 0.5 mm domain
    Nx=32,     # 100 grid points
    Ny=32,     # 100 grid points
    radius=0.1e-3,  # 0.2 mm circle
    filename="/Users/jacksonbaglino/PetIGA-3.20/demo/input/Thermal_IO/circle_phase_field.dat"
)
