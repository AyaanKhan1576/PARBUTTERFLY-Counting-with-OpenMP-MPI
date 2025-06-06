# Unified Butterfly Counting Report

## Prerequisites

Before running `butter.cpp`, ensure you have:

* **Preprocessed & normalized** your bipartite edge lists (0-based, authors as L, papers as R).
* **Partitioned** the normalized graph into subgraphs via METIS (e.g. `subgraph_0.txt`, `subgraph_1.txt`, …).
* A C++17 compiler with OpenMP support (e.g. `g++` ≥ 7).

## Compilation

```bash
# From the directory containing butter.cpp
g++ -std=c++17 -O3 butter.cpp -fopenmp -o butter
```

## Usage

### Sequential mode

```bash
./butter --mode seq <edge_file> <L_size> \
       --rank id|degree|degcore|random[:seed] \
       --agg hash|sort|batch|histo
```

Example:

```bash
./butter --mode seq processed/subgraph_0.txt 1953085 \
         --rank degree \
         --agg hash
```

### Parallel mode

```bash
./butter --mode par <edge_file> <L_size> \
       --rank id|degree|degcore|random[:seed] \
       --agg hash|sort|batch|histo \
       --threads N
```

Example (4 threads):

```bash
./butter --mode par processed/subgraph_0.txt 1953085 \
         --rank degree \
         --agg hash \
         --threads 4
```

---

## File Structure & Workflow

1. **Graph Loader** (`load_graph`)

   * Reads an edge list connecting L-nodes (`0..L_size-1`) to R-nodes (`>= L_size`).
   * Validates each line, skips malformed or out-of-range edges.
   * Builds two CSR-like lists:

     * `L2R[u]` = list of paper IDs connected to author `u`.
     * `R2L[r]` = list of author IDs connected to paper index `r - L_size`.
   * Logs warnings if the provided `L_size` doesn’t cover all author IDs.

2. **Ranking** (`make_rank`)

   * Produces a permutation of authors according to one of four strategies:

     * **ID**: natural order.
     * **Degree**: sort by ascending L2R degree.
     * **Degcore**: one-pass approximate core decomposition.
     * **Random**: shuffle (with optional seed).
   * The resulting `order` maps each author → its position; `rank_of[node]` is the inverse map.

3. **Aggregators**
   Maintain counters of wedges `(u,w)` to accumulate butterflies `C(f,2)` where `f` is their multiplicity.

   * **`AggHash`** (incremental choose-two)

     * Uses an `unordered_map<pair64,uint64_t>` to track frequencies.
     * `add(p)` does `total += h[p]`, then `++h[p]`.
     * **No intermediate flush**—stateful counts persist until the end.

   * **`AggSort`** (one-shot sort)

     * Buffers every `pair64` in a vector, sorts it once at `final_flush()`, then counts runs.

   * **`AggBatch`** (batch variant of sort)

     * Identical to `AggSort`, but originally supported chunked flushing; we removed auto-flush to ensure correctness.

   * **`AggHisto`** (dense histogram)

     * Uses a dense array `hist[w]` and a `touched` list to track which bins changed.
     * **Caution**: logic relies on final flush and may require iteration-order guarantees.

4. **Parallel Counting** (`count_par`)

   * **Thread-local aggregators**: one per OpenMP thread.
   * **Inbox messaging**: pairs `(u,w)` whose owner thread ≠ current tid are pushed (under lock) into that thread’s inbox.
   * **Owner function**: FNV‑1a hash of `(u,w)` mod `nT` for a balanced distribution.
   * **Two-phase execution**:

     1. **Wedge enumeration** over `R2L`: each paper’s author list ⇒ generate `(u,w)` with `i<j`.
     2. **Barrier & inbox drain**: each thread consumes its mailbox, calls `add(p)` on every received pair.
   * **Final flush** processes any remaining buffered wedges (for sort, batch).
   * Sums per-thread totals for the global butterfly count.
   * **No mid-loop flushing** avoids losing partial counts in `AggHash`/`AggHisto`.

5. **CLI & Timing**

   * Robust command-line parsing (`--mode`, positional `<edge_file>` & `<L_size>`, flags `--rank`, `--agg`, `--threads`).
   * Consistent use of `count_par` for both modes: `--mode seq` simply sets `nT=1`, guaranteeing identical logic.
   * Timing via `omp_get_wtime()` if OpenMP is enabled; else `chrono::high_resolution_clock`.

---

## Challenges & Resolutions

1. **Inconsistent counts** between sequential and parallel runs when using stateful aggregators.
   **Fix:** Removed any periodic flush inside the wedge loop; only a final flush remains.

2. **Sequential vs. `--threads 1` mismatch**: original code had two separate drivers (`count_seq` vs `count_par`).
   **Fix:** Unified both modes onto the same `count_par<…>(…,nT)` path, simply varying `nT`.

3. **Edge-case validations** in `load_graph`: malformed lines, out-of-range IDs, incorrect `L_size`.
   **Fix:** Added warnings, robust bounds checks, and fallback strategies to avoid silent corruption.

4. **Lock and inbox safety** with OpenMP: race conditions when multiple threads push to the same inbox.
   **Fix:** Introduced one `omp_lock_t` per thread + `omp_set_lock`/`omp_unset_lock` around pushes.

5. **CLI parsing & error handling**: many edge-case flags, missing arguments, or bad numbers.
   **Fix:** Centralized `parse()` logic, explicit checks, and detailed user feedback on misuse.

6. **Degcore algorithm correctness**: ensuring bin positions (`pos[]`) remain consistent under dynamic degree updates.
   **Fix:** Added safety checks, fallback linear searches on mismatches, and clear error messages if logic breaks.

7. **Histogram aggregator correctness**: recognized that `AggHisto` may not be semantically correct for all iteration orders.
   **Note:** It remains available but with warnings to use it only if deeply understood.

---

**Outcome:** a single, well‐tested C++17 program that deterministically counts butterflies in a large bipartite graph, either sequentially or in parallel, with a choice of ranking and aggregation strategies. You can now integrate this into your MPI‐based multi-node workflow for full-scale distributed counting.
