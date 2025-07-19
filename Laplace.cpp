#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <numeric>

using namespace std;

#define N 100               // Grid size
#define TOLERANCE 1e-6      // Convergence criteria
#define MAX_ITER 10000      // Maximum iterations

// Function to initialize the grid with boundary conditions
void initializeGrid(vector<vector<double>> &grid) {
    for (int i = 0; i < N; i++) {
        grid[0][i] = 5.0;    // Top boundary
        grid[N - 1][i] = -5.0; // Bottom boundary
        grid[i][0] = 0.0;     // Left boundary
        grid[i][N - 1] = 0.0; // Right boundary
    }
}

// Compute checksum for result verification (ignore boundary)
double computeChecksum(const vector<vector<double>> &grid) {
    double sum = 0.0;
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            sum += grid[i][j];
        }
    }
    return sum;
}

// Serial Laplace solver
double solveLaplaceSerial(vector<vector<double>> &grid) {
    vector<vector<double>> newGrid = grid;
    double error;
    int iter = 0;

    do {
        error = 0.0;
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                newGrid[i][j] = 0.25 * (grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1]);
                error = max(error, fabs(newGrid[i][j] - grid[i][j]));
            }
        }
        grid.swap(newGrid);
        iter++;
    } while (error > TOLERANCE && iter < MAX_ITER);

    return error;
}

// OpenMP Parallel Laplace solver
double solveLaplaceOpenMP(vector<vector<double>> &grid) {
    vector<vector<double>> newGrid = grid;
    double error;
    int iter = 0;

    do {
        error = 0.0;
        #pragma omp parallel for reduction(max:error)
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                newGrid[i][j] = 0.25 * (grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1]);
                error = max(error, fabs(newGrid[i][j] - grid[i][j]));
            }
        }
        grid.swap(newGrid);
        iter++;
    } while (error > TOLERANCE && iter < MAX_ITER);

    return error;
}

// MPI Distributed Laplace solver
double solveLaplaceMPI(vector<vector<double>> &grid, int rank, int size) {
    int start = (N / size) * rank;
    int end = (rank == size - 1) ? N : start + (N / size);

    vector<vector<double>> newGrid = grid;
    double error, globalError;
    int iter = 0;

    do {
        error = 0.0;
        for (int i = start + 1; i < end - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                newGrid[i][j] = 0.25 * (grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1]);
                error = max(error, fabs(newGrid[i][j] - grid[i][j]));
            }
        }

        MPI_Allreduce(&error, &globalError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        grid.swap(newGrid);
        iter++;
    } while (globalError > TOLERANCE && iter < MAX_ITER);

    return globalError;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <method>\n";
        cout << "Methods: serial | openmp | mpi\n";
        return 1;
    }

    string method = argv[1];

    if (method == "serial" || method == "openmp") {
        vector<vector<double>> grid(N, vector<double>(N, 0.0));
        initializeGrid(grid);

        double error;
        if (method == "serial") {
            cout << "Running Serial Solver...\n";
            error = solveLaplaceSerial(grid);
        } else {
            cout << "Running OpenMP Solver...\n";
            error = solveLaplaceOpenMP(grid);
        }

        cout << "Converged with error: " << error << endl;
        cout << "Checksum: " << computeChecksum(grid) << endl;
    }

    else if (method == "mpi") {
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        vector<vector<double>> grid(N, vector<double>(N, 0.0));
        if (rank == 0) initializeGrid(grid);
        MPI_Bcast(&grid[0][0], N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        double error = solveLaplaceMPI(grid, rank, size);

        if (rank == 0) {
            cout << "MPI Laplace Solver converged with error: " << error << endl;
            cout << "Checksum: " << computeChecksum(grid) << endl;
        }

        MPI_Finalize();
    } 

    else {
        cout << "Invalid method! Choose from: serial, openmp, mpi\n";
        return 1;
    }

    return 0;
}
