#include <iostream>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <mpi.h>

using namespace std;

const int N = 1000; // Large matrix size

//  Function to allocate a matrix dynamically
double** allocateMatrix(int size) {
    double** matrix = new double*[size];
    for (int i = 0; i < size; i++)
        matrix[i] = new double[size];
    return matrix;
}

//  Function to free allocated memory
void freeMatrix(double** matrix, int size) {
    for (int i = 0; i < size; i++)
        delete[] matrix[i];
    delete[] matrix;
}

//  Function to initialize matrix with random values
void initializeMatrix(double** matrix, int size) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix[i][j] = rand() % 10;
}

void printMatrix(double** matrix, int size) {
    int displaySize = (size > 5) ? 5 : size;  // Print only last 5x5 if matrix is large
    for (int i = size - displaySize; i < size; i++) { // Loop through all rows
        for (int j = size - displaySize; j < size; j++) { // Last 5 columns
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << "-----------------------------\n";
}

// 🚀 Function for Serial Matrix Multiplication
void serialMatrixMultiplication(double** A, double** B, double** C, int size) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            C[i][j] = 0;
            for (int k = 0; k < size; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
}

// 🚀 Function for Parallel Matrix Multiplication (OpenMP)
void parallelMatrixMultiplication(double** A, double** B, double** C, int size) {
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++) {
            C[i][j] = 0;
            for (int k = 0; k < size; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
}

// 🚀 Function for Distributed Matrix Multiplication (MPI)
void mpiMatrixMultiplication(double** A, double** B, double** C, int size, int rank, int numProcs) {
    int rowsPerProcess = size / numProcs;
    double** subA = allocateMatrix(rowsPerProcess);
    double** subC = allocateMatrix(rowsPerProcess);

    // Scatter matrix A among processes
    for (int i = 0; i < rowsPerProcess; i++)
        MPI_Scatter(A[i], size, MPI_DOUBLE, subA[i], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Broadcast matrix B to all processes
    for (int i = 0; i < size; i++)
        MPI_Bcast(B[i], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Compute partial matrix multiplication
    for (int i = 0; i < rowsPerProcess; i++)
        for (int j = 0; j < size; j++) {
            subC[i][j] = 0;
            for (int k = 0; k < size; k++)
                subC[i][j] += subA[i][k] * B[k][j];
        }

    // Gather results from all processes
    for (int i = 0; i < rowsPerProcess; i++)
        MPI_Gather(subC[i], size, MPI_DOUBLE, C[i], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    freeMatrix(subA, rowsPerProcess);
    freeMatrix(subC, rowsPerProcess);
}

// 🚀 Main Function
int main(int argc, char* argv[]) {
    srand(time(0));
    double** A = allocateMatrix(N);
    double** B = allocateMatrix(N);
    double** C1 = allocateMatrix(N);
    double** C2 = allocateMatrix(N);
    double** C3 = allocateMatrix(N);

    initializeMatrix(A, N);
    initializeMatrix(B, N);

    // Serial Execution
    cout << "Running Serial Matrix Multiplication...\n";
    serialMatrixMultiplication(A, B, C1, N);
    printMatrix(C1, N);

    // Parallel Execution
    cout << "Running Parallel Matrix Multiplication (OpenMP)...\n";
    parallelMatrixMultiplication(A, B, C2, N);
    printMatrix(C2, N);

    // MPI Execution
    cout << "Running Distributed Matrix Multiplication (MPI)...\n";
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpiMatrixMultiplication(A, B, C2, N, rank, numProcs);
    MPI_Finalize();
    printMatrix(C2, N);


    freeMatrix(A, N);
    freeMatrix(B, N);
    freeMatrix(C1, N);
    freeMatrix(C2, N);
    freeMatrix(C3, N);

    return 0;
}
