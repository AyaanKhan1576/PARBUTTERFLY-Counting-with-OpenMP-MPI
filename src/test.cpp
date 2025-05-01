// src/test.cpp
#include <iostream>
#include <omp.h>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #pragma omp parallel
    {
        printf("Hello from thread %d in process %d\n", omp_get_thread_num(), rank);
    }
    MPI_Finalize();
    return 0;
}
