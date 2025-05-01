```markdown
# Pre-processing & Partitioning Workflow  
*(input: **DBLP-author** raw edge list → output: balanced per-rank sub-graphs)*

---

## 0 Quick build & run

```bash
# Download dataset and add it to data folder
http://konect.cc/networks/dblp-author/

# inside project-root
g++ -std=c++17 -O3 preprocessing/normalize_bipartite.cpp -o normalize
g++ -std=c++17 -O3 preprocessing/preprocess_and_partition.cpp -lmetis -o preprocess

# step-1  normalise, creating clean bipartite graph + lookup tables
./normalize data/out.dblp-author processed

# step-2  build CSR, run METIS, dump per-rank sub-graphs (example k = 8)
./preprocess processed/normalized_edges.txt 4 processed
```

Outputs appear in **`processed/`**:
```
normalized_edges.txt   map_L.txt  map_R.txt
graph.metis            part.8
subgraph_0.txt … subgraph_7.txt
```

---

## 1 Why two-step pre-processing?

| Goal | Benefit for butterfly counting |
|------|--------------------------------|
| **Explicit bipartite normalisation** | Clear author/paper split enables side-aware kernels (e.g., “wedge centre = author”), higher cache locality, and simpler rank-boundary logic. |
| **Compact 0-based IDs** | Dense CSR arrays instead of hash maps → lower RAM & faster traversal. |
| **Balanced k-way partition (METIS)** | Each MPI rank gets roughly equal vertices and minimal edge cuts, reducing cross-rank communication. |
| **Per-rank edge lists** | Ranks can start immediately without filtering the global file; duplicates removed, IDs remain global for easy ghost lookup. |

---

## 2 Workflow in plain English

1. **Normalisation (`normalize_bipartite`)**  
   *Scans the raw KONECT file once to learn every unique author and publication.*  
   * Assigns authors IDs `0 … |L|-1`.  
   * Assigns publications IDs `|L| … |L|+|R|-1`.  
   * Saves two look-up tables:  
     * `map_L.txt` – newAuthorID → originalAuthorID  
     * `map_R.txt` – newPubID → originalPubID  
   * Rewrites the edge list with these compact IDs → `normalized_edges.txt`.

2. **CSR build & METIS partition (`preprocess`)**  
   *Reads `normalized_edges.txt` and constructs an undirected CSR graph in memory.*  
   * Dumps `graph.metis` (ASCII) so you can sanity-check with METIS CLI tools.  
   * Invokes **`METIS_PartGraphKway`** with `k` parts to minimise edge cuts.  
   * Writes `part.k` (one integer per vertex telling which rank owns it).  
   * Splits the edge list so that rank *r* receives `subgraph_r.txt`, containing only edges owned by that rank (rule: use the rank of the first endpoint).

---

## 3 How the outputs feed the butterfly engine

| File | Role at run-time |
|------|------------------|
| `subgraph_r.txt` | Each MPI rank builds its local CSR from this file—no duplicates, already 0-based. |
| `part.k` | Used to label vertices as **local** (`part[v] == rank`) or **ghost** (`part[v] ≠ rank`) so that each rank exchanges only the necessary boundary neighbourhoods. |
| `map_L.txt`, `map_R.txt` | After counting, translate compact IDs back to the original DBLP integers (or to author names pulled from DBLP dumps) when reporting top-k butterfly hubs. |

---

## 4 Validation checklist

| Check | Expected |
|-------|----------|
| Authors range | `0 ≤ u < |L|` in every edge line |
| Publications range | `|L| ≤ v < |L|+|R|` in every edge line |
| `graph.metis` header | first line = “`nvtxs m`”, with `nvtxs = |L|+|R|` |
| Edge-cut quality | printed `objval` ≈ `gpmetis graph.metis k` result |
| Partition vector consistency | `part.k` has exactly `nvtxs` lines, each `0 … k-1` |
| Per-rank edge files | Union of all `subgraph_r.txt` equals `normalized_edges.txt`; intersections empty. |

Following this pipeline, you enter the counting phase with **balanced work, minimal comms, and reversible mapping to the raw dataset**—exactly what large-scale butterfly enumeration needs.
```