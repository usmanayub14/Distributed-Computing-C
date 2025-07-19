import numpy as np
import matplotlib.pyplot as plt

# Boundary values
top = 5.0
bottom = -5.0
left = 0.0
right = 0.0

# Grid size
rows, cols = 50, 50
grid = np.zeros((rows, cols))

# Apply boundary conditions
grid[0, :] = top
grid[-1, :] = bottom
grid[:, 0] = left
grid[:, -1] = right

# Solve using Jacobi iteration
def solve_laplace(grid, max_iterations=1000, tol=1e-5):
    for _ in range(max_iterations):
        old = grid.copy()
        grid[1:-1, 1:-1] = 0.25 * (
            old[0:-2, 1:-1] + old[2:, 1:-1] +
            old[1:-1, 0:-2] + old[1:-1, 2:]
        )
        if np.max(np.abs(grid - old)) < tol:
            break
    return grid

solved = solve_laplace(grid)

# Plot contour
x = np.linspace(0, 1, cols)
y = np.linspace(0, 1, rows)
X, Y = np.meshgrid(x, y)

plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, solved, 50, cmap='coolwarm')
plt.colorbar(contour, label='Voltage (V)')
plt.title('2D Laplace Equation Solution (Contour Plot)')
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.tight_layout()
plt.show()
