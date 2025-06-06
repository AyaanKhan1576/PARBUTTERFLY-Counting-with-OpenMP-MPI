# Parallel Butterfly Computation in Bipartite Graphs

## Contributors

* **Ayaan Khan (22i-0832)**
* **Ayaan Mughal (22i-0861)**
* **Salar Shoaib (20i-0830)**

---

## Overview

This project implements a **parallel algorithm** to count **butterfly motifs** in large bipartite graphs. Butterfly motifs are fundamental structures in bipartite networks and crucial for identifying dense local interactions, such as co-authorship networks and co-purchase patterns.

The algorithm is based on the **PARBUTTERFLY** framework, adapted to leverage a **hybrid parallelization approach** using **MPI** (inter-node) and **OpenMP** (intra-node) for high-performance computing. **METIS** is employed for balanced partitioning.

**Note**: This project is designed for execution on a **cluster environment** with support for **MPI and OpenMP**. Running it on a single machine without a proper cluster configuration may not yield the intended performance results.

---

## Key Features

* **Scalable butterfly counting** in large bipartite graphs
* Supports global butterfly counting with optional tip/wing decomposition
* **Hybrid parallelism** combining MPI and OpenMP
* Partition-aware computation using **METIS** to minimize cross-node communication
* Flexible ranking and aggregation strategies for performance tuning

---

## Getting Started

### Quickstart with starter\_code.sh

For ease of use, a single file—**starter\_code.sh**—is provided that automates the compilation, preprocessing, and execution steps. Simply place your raw dataset in the `data/` directory and run:

```bash
bash starter_code.sh
```

This script:

* Compiles the preprocessing and counting modules
* Normalizes and partitions the graph
* Runs the sequential, OpenMP, MPI, and hybrid butterfly counting routines
* Performs benchmark tests and logs the results

If **starter\_code.sh** does not fit your setup (e.g. you want more control or debugging), please follow the **Manual Setup** section below.

---

## Manual Setup

### Prerequisites

* **Cluster environment** with support for:

  * **C++17 compiler** (g++ ≥ 7 recommended)
  * **OpenMP**
  * **MPI** (e.g. OpenMPI)
  * **METIS** library

* **Data**: DBLP-author bipartite dataset (or other bipartite graphs)

### Installation

1. **Download & Prepare Data**

   ```bash
   wget http://konect.cc/files/download.tsv.dblp-author.tar.bz2
   tar -xvjf download.tsv.dblp-author.tar.bz2 -C data/
   ```

2. **Preprocessing**

   ```bash
   g++ -std=c++17 -O3 preprocessing/normalize_bipartite.cpp -o normalize
   g++ -std=c++17 -O3 preprocessing/preprocess_and_partition.cpp -lmetis -o preprocess

   ./normalize data/out.dblp-author processed
   ./preprocess processed/normalized_edges.txt 4 processed
   ```

3. **Compiling Butterfly Counting Engines**

   ```bash
   g++ -std=c++17 -O3 src/butter_core.cpp -fopenmp -o butter_core
   mpic++ -std=c++17 -O3 src/butter_mpi.cpp -fopenmp -o butter_mpi
   ```

---

## Running the Code

### Sequential Execution

```bash
./butter_core --mode seq --file processed/subgraph_0.txt --L 1953085 --rank degree --agg hash
```

### OpenMP (Intra-node Parallelism)

```bash
./butter_core --mode par --file processed/subgraph_0.txt --L 1953085 --rank degree --agg hash --threads 4
```

### MPI (Inter-node Parallelism)

```bash
mpirun -np 4 ./butter_mpi --mode seq --subgraph_dir ./processed --L 1953085 --rank degree --agg hash
mpirun -np 4 ./butter_mpi --mode par --subgraph_dir ./processed --L 1953085 --rank degree --agg hash --threads 2
```

### Hybrid MPI + OpenMP

```bash
mpirun -np 4 ./butter_mpi --mode par --subgraph_dir ./processed --L 1953085 --rank degree --agg hash --threads 2
```

---

## System Architecture

```plaintext
      ┌────────────┐
      │  METIS     │
      └────┬───────┘
           ↓
 ┌─────────────────────┐
 │  MPI: Partitioned   │
 │  Processes          │
 └────┬─────┬─────┬────┘
      ↓     ↓     ↓
┌───────┐ ┌───────┐ ┌───────┐
│OpenMP │ │OpenMP │ │OpenMP │
└───────┘ └───────┘ └───────┘
           ↓
     MPI Reduce
           ↓
   Global Butterfly Count
```

**Workflow Summary:**

1. **Normalize** the bipartite graph with 0-based IDs
2. **Partition** with METIS into balanced subgraphs
3. **Parallel Counting** with MPI + OpenMP
4. **Aggregate** results for global butterfly counts

---

## Performance Highlights

* **Sequential vs MPI**: Up to 4x speedup on large graphs with 4 MPI ranks
* **MPI vs Hybrid**: 15–30% faster with OpenMP threads per MPI process, especially on multi-core nodes
* **OpenMP Scalability**: Good scaling up to a point, with diminishing returns beyond optimal thread count

---

## Notes & Troubleshooting

* Designed for **cluster execution**; single-machine setups may not reflect the intended performance.
* Use `mpirun --oversubscribe` if you encounter slot allocation issues:

  ```bash
  mpirun --oversubscribe -np 4 ./butter_mpi ...
  ```
* Validate outputs:

  ```bash
  cat processed/normalized_edges.txt | wc -l
  sum(subgraph_*.txt lines) == normalized_edges.txt lines
  ```

---

## Acknowledgments

* Inspired by **“Parallel Algorithms for Butterfly Computations”** by Jessica Shi and Julian Shun
* **METIS** library for partitioning: [METIS Overview](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
* **DBLP-author dataset** from KONECT

---
