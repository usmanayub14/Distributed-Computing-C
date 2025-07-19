# Distributed Computing: MPI and Socket Programming in C/C++

This repository showcases my work in distributed computing, focusing on Message Passing Interface (MPI) for high-performance distributed applications and fundamental socket programming for network communication. These projects highlight techniques for data distribution, synchronization, and collective operations across multiple compute nodes.

## Key Learning Outcomes & Skills:

* **MPI (Message Passing Interface):**
    * `MPI_Init`, `MPI_Comm_rank`, `MPI_Comm_size`, `MPI_Finalize`.
    * Point-to-point communication (`MPI_Send`, `MPI_Recv`).
    * Collective operations (`MPI_Bcast`, `MPI_Scatterv`, `MPI_Gatherv`, `MPI_Reduce`).
    * `Master-Worker paradigm`.
    * `Halo Exchange` for domain decomposition in iterative solvers.
    * `Host file configuration` for multi-node execution.
* **Network Socket Programming:**
    * `TCP client/server model`.
    * `Asynchronous chat application (non-blocking I/O or multi-threading)`.
    * `File transfer over sockets`.
* **High-Performance Computing (HPC) Environment Interaction:** Compiling and running distributed applications on HPC clusters.
* **Scalability Analysis:** Evaluating how performance changes with increasing number of processes.
* **C/C++ Programming:** Robust and efficient C/C++ code for distributed environments.

## Projects:

### 1. Asynchronous Chat Application (TCP Sockets)
* **Description:** Implements a client/server chat application using TCP sockets, allowing asynchronous communication where both client and server can send and receive messages at any time.
* **Skills Showcased:** `TCP sockets (socket, bind, listen, accept, connect)`, `multi-threading for concurrent send/receive`, `non-blocking I/O (select/poll) or separate threads for send/receive operations`.
* **Source:** `PDP-Spring2025-Assignment-02-21022025.pdf` (Question 1)
* **Details:** `chat_app_server.cpp` and `chat_app_client.cpp`.
* **Results:** Screenshots demonstrating simultaneous message exchange.

### 2. Distributed Array Summation (Sockets)
* **Description:** Computes the sum of a large array using a distributed approach over a network via sockets, allowing processing on multiple machines.
* **Skills Showcased:** `Distributed computing using sockets`, `data partitioning and distribution across machines`, `aggregation of partial results`, `inter-machine communication for computational tasks`.
* **Source:** `PDP-Spring2025-Assignment-03-10032025.pdf` (Question 1)
* **Details:** `distributed_array_sum_server.cpp`, `distributed_array_sum_client.cpp`.
    * **Approach:** A server distributes portions of the array to multiple clients (machines). Clients compute partial sums and send them back to the server for final aggregation. Can optionally integrate OpenMP/Pthreads on each client to fully utilize its local cores.
* **Results:** Performance data comparing sequential, OpenMP, and distributed (socket-based) implementations.

### 3. Distributed Matrix Multiplication (Sockets)
* **Description:** Performs matrix multiplication by distributing computation tasks across multiple machines using TCP/IP sockets.
* **Skills Showcased:** `Distributed matrix multiplication over network`, `data transfer protocols for matrices`, `handling large data transfers over sockets`, `scalability of network-based distributed computing`.
* **Source:** `PDP-Spring2025-Assignment-03-10032025.pdf` (Question 2)
* **Details:** `distributed_matrix_mul_server.cpp`, `distributed_matrix_mul_client.cpp`.
    * **Approach:** Similar to array sum, the server orchestrates the distribution of matrix blocks and collection of results from client machines.
* **Results:** Performance comparison between sequential, OpenMP, and distributed (socket) versions.

### 4. MPI Hello World
* **Description:** A basic MPI program to verify the setup of the MPI environment and demonstrate fundamental communication concepts.
* **Skills Showcased:** `MPI initialization and finalization`, `determining process rank and size`, `retrieving processor name`.
* **Source:** `CS435 Lab 12 Manual - Distributed Computing with MPI.docx` (Task 1)
* **Details:** `mpi_hello_world.cpp`.

### 5. MPI Distributed Matrix Multiplication
* **Description:** Implements a distributed matrix multiplication algorithm where matrix computations are parallelized across multiple MPI processes.
* **Skills Showcased:** `MPI_Bcast`, `MPI_Scatterv`, `MPI_Gatherv`, `master-worker architecture`, `matrix partitioning (row-wise)`.
* **Source:** `CS435 Lab 12 Manual - Distributed Computing with MPI.docx` (Task 2)
* **Details:** `mpi_matrix_multiplication.cpp`.
* **Results:** Terminal output showing computation time and (for small matrices) the results, along with hostfile configuration used.

### 6. MPI Distributed 2D Convolution (Conceptual Outline for Implementation)
* **Description:** This project would extend the 2D convolution problem to a distributed setting using MPI, partitioning the input image/matrix across multiple processes.
* **Skills Showcased:** `MPI data partitioning strategies (e.g., block-wise)`, `ghost cell exchange/halo exchange` (if filter overlaps boundaries), `point-to-point communication (MPI_Send/Recv)` for boundary data, `collective operations for initialization/finalization`.
* **Source:** `CS435 Lab 12 Manual - Distributed Computing with MPI.docx` (Task 3, original task was sequential C++)
* **Details:** (Your implementation of `mpi_2d_convolution.cpp`)

### 7. MPI Numerical Integration
* **Description:** Calculates definite integrals of mathematical functions using a distributed numerical integration approach with MPI.
* **Skills Showcased:** `MPI_Reduce` for aggregating results, `load balancing` for distributing integration steps, `function pointers` for flexible integral calculation.
* **Source:** `CS435 Lab 12 Manual - Distributed Computing with MPI.docx` (Task 4)
* **Details:** `mpi_numerical_integration.cpp`.

### 8. MPI Distributed 2D Laplace Solver
* **Description:** Solves the 2D Laplace equation using the Jacobi iterative method, distributed across multiple MPI processes using a row-wise domain decomposition and halo exchange.
* **Skills Showcased:** `MPI domain decomposition`, `halo exchange (MPI_Sendrecv)`, `handling boundary conditions for distributed grids`, `iterative solvers in distributed environments`, `MPI_Barrier` for synchronization.
* **Source:** `Project Proposal.docx` (Part II) and `PDP-Spring2025-Assignment-03-10032025.pdf` (Question 3).
* **Details:** `laplace_mpi.cpp`.
    * **Boundary Conditions:** Explicitly mentions using `top=5V, bottom=-5V, left=right=0V` (from Assignment 3).
    * **Approach:** Each process manages a subset of rows, and `MPI_Sendrecv` handles the exchange of boundary data (ghost cells) between neighbors for accurate stencil computations.
    * **Correctness:** Verified by visual inspection for small grids and consistency with sequential/OpenMP versions.
* **Results:** Performance analysis included in the `HPC-Cluster-Implementations` repository, but mention it here as well.
