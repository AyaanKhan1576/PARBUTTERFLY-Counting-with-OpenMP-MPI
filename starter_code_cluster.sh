#!/bin/bash

# Allow MPI to run as root (necessary inside Docker)
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Explicit prefix for OpenMPI (so remote nodes can find orted)
OMPI_PREFIX="--prefix /usr"

# Step 1: Compile preprocessing tools
g++ -std=c++17 -O3 preprocessing/normalize_bipartite.cpp -o normalize
g++ -std=c++17 -O3 preprocessing/preprocess_and_partition.cpp -lmetis -o preprocess

# Step 2: Normalize and partition the dataset
./normalize data/out.dblp-author processed
./preprocess processed/normalized_edges.txt 4 processed

# Step 3: Compile butterfly counting programs
g++ -std=c++17 -O3 src/butter_core.cpp -fopenmp -o butter_core
mpic++ -std=c++17 -O3 src/butter_mpi.cpp -fopenmp -o butter_mpi

# Step 4: Run single-node OpenMP tests
./butter_core --mode par --file processed/subgraph_0.txt --L 1953085 --rank degree --agg hash --threads 4
./butter_core --mode seq --file processed/subgraph_0.txt --L 1953085 --rank degree --agg hash

# Step 5: Run distributed MPI tests (sequential and parallel)
mpirun $OMPI_PREFIX -np 4 --hostfile machinefile ./butter_mpi --mode par --subgraph_dir ./processed --L 1953085 --rank degree --agg hash --threads 2
mpirun $OMPI_PREFIX -np 4 --hostfile machinefile ./butter_mpi --mode seq --subgraph_dir ./processed --L 1953085 --rank degree --agg hash

# Step 6: Benchmarking (OpenMP only)
export OMP_NUM_THREADS=1
for i in {1..20}
do
  { time ./butter_core; } >> dataOpenMPSeq.txt 2>&1
done

# Step 7: Benchmarking (MPI + OpenMP, sequential)
for i in {1..20}
do
  { time mpirun $OMPI_PREFIX -perhost 2 -np 4 --hostfile machinefile ./butter_mpi; } >> dataOpenMPMPISeq.txt 2>&1
done

# Step 8: Benchmarking (OpenMP parallel)
export OMP_NUM_THREADS=4
for i in {1..20}
do
  { time ./butter_core; } >> dataOpenMPPara.txt 2>&1
done

# Step 9: Benchmarking (MPI + OpenMP, parallel)
for i in {1..20}
do
  { time mpirun $OMPI_PREFIX -perhost 2 -np 4 --hostfile machinefile ./butter_mpi; } >> dataOpenMPMPIPara.txt 2>&1
done
