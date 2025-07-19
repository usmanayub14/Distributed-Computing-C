#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <mpi.h>

using namespace std;

#define ARRAY_SIZE 1000000  // Large array size

// Function for serial sum
long long serialSum(vector<int>& arr) {
    long long sum = 0;
    for (int num : arr) {
        sum += num;
    }
    return sum;
}

// Function for parallel sum using OpenMP
long long parallelSum(vector<int>& arr) {
    long long sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < arr.size(); i++) {
        sum += arr[i];
    }
    return sum;
}

// Function for distributed sum using MPI
long long distributedSum(vector<int>& arr) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunk_size = arr.size() / size;
    vector<int> local_chunk(chunk_size);

    // Scatter data to different processes
    MPI_Scatter(arr.data(), chunk_size, MPI_INT, local_chunk.data(), chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Compute local sum
    long long local_sum = 0;
    for (int num : local_chunk) {
        local_sum += num;
    }

    // Reduce to global sum
    long long global_sum = 0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    return global_sum;
}

int main(int argc, char* argv[]) {
    srand(time(0));
    vector<int> arr(ARRAY_SIZE);

    // Fill array with random numbers
    for (int& num : arr) {
        num = rand() % 100;
    }

    long long sum = 0;

    // Serial execution
    cout << "Running Serial Sum...\n";
    sum = serialSum(arr);
    cout << "Serial Sum: " << sum << endl;

    // Parallel execution with OpenMP
    cout << "Running Parallel Sum with OpenMP...\n";
    sum = parallelSum(arr);
    cout << "OpenMP Sum: " << sum << endl;

    // Distributed execution with MPI
    cout << "Running Distributed Sum with MPI...\n";
    MPI_Init(&argc, &argv);
    sum = distributedSum(arr);
    if (sum != 0) {
        cout << "MPI Sum: " << sum << endl;
    }
    MPI_Finalize();

    return 0;
}
