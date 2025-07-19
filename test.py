import numpy as np
import matplotlib.pyplot as plt

# Parameters
nx, ny = 50, 50  # Grid size
V_top, V_bottom = 5, -5  # Boundary conditions
V_left, V_right = 0, 0
max_iter = 10000
tolerance = 1e-5

# Initialize potential grid
V = np.zeros((nx, ny))

# Set boundary conditions
V[:, 0] = V_left
V[:, -1] = V_right
V[0, :] = V_top
V[-1, :] = V_bottom

# Iterative solver
for iteration in range(max_iter):
    V_old = V.copy()
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            V[i, j] = 0.25 * (V_old[i + 1, j] + V_old[i - 1, j] + V_old[i, j + 1] + V_old[i, j - 1])

    # Check for convergence
    if np.linalg.norm(V - V_old) < tolerance:
        print(f"Converged after {iteration} iterations")
        break

# Calculate checksum
checksum = np.sum(V)
print(f"Checksum of the potential grid: {checksum}")

# Plot the potential
plt.imshow(V, cmap='viridis', extent=[0, 1, 0, 1])
plt.colorbar(label='Potential (V)')
plt.title('2-D Laplace Equation Solution')
plt.xlabel('x')
plt.ylabel('y')
plt.show()