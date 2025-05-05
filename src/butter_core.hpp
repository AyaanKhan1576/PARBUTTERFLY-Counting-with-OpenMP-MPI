#ifndef BUTTER_CORE_HPP
#define BUTTER_CORE_HPP

#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>
#include <atomic>
#include <numeric>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <random>
#include <unordered_map>
#include <cmath>
#include <filesystem> // Required for checking directory in MPI version

using vid    = uint32_t;
using pair64 = uint64_t;
static inline pair64 pack(vid a, vid b) { return (pair64)a << 32 | b; }

// ───────────────────── Graph loader ───────────────────────────
// Assumes node IDs in the file are GLOBAL IDs
struct Graph {
    vid L;
    std::vector<std::vector<vid>> L2R;
    std::vector<std::vector<vid>> R2L;
    vid L_param_original;
};

Graph load_graph(const std::string& path, vid L_param) {
    std::ifstream in(path);
    // Use std::filesystem to check existence clearly
    if (!std::filesystem::exists(path) || !std::filesystem::is_regular_file(path)) {
         throw std::runtime_error("Cannot open graph file (or not a file): " + path);
    }
     if (!in.is_open()) { // Double check after existence check
          throw std::runtime_error("Failed to open graph file stream: " + path);
     }

    std::vector<std::pair<vid, vid>> edges;
    size_t Rmax_idx = 0;
    vid Lmax_id = 0;
    bool first_r_node_found = false;

    for (std::string line; std::getline(in, line); ) {
        if (line.empty() || line[0] == '%' || line[0] == '#') continue;
        std::istringstream iss(line);
        vid u, v;
        if (!(iss >> u >> v)) {
             std::cerr << "Warning: Skipping malformed line in " << path << ": " << line << std::endl;
             continue;
        }
        // *** CRITICAL ASSUMPTION: u, v are GLOBAL IDs ***
        if (u >= L_param) {
             // This indicates an L node ID >= L_param, which breaks assumptions
             std::cerr << "Warning: Skipping edge (" << u << "," << v << ") in " << path
                       << " - L-node ID u (" << u << ") >= L_param (" << L_param << "). Check partitioning/L param." << std::endl;
             continue;
        }
        if (v < L_param) {
              // This indicates an R node ID < L_param
             std::cerr << "Warning: Skipping edge (" << u << "," << v << ") in " << path
                       << " - R-node ID v (" << v << ") < L_param (" << L_param << "). Check partitioning/L param." << std::endl;
             continue;
        }

        edges.emplace_back(u, v);
        size_t r_idx = v - L_param; // Index based on GLOBAL R node ID
        Rmax_idx = std::max(Rmax_idx, r_idx);
        Lmax_id = std::max(Lmax_id, u);
         if (!first_r_node_found) first_r_node_found = true;
    }
     in.close();

    // Effective sizes based on the MAX IDs seen GLOBALLY (important for consistent vector sizing if possible)
    // Note: For subgraph loading, these max values only reflect the local subgraph.
    // This might lead to smaller L/R2L vectors than needed if another subgraph has higher IDs.
    // A potential issue if ranking/counting relied on absolute max ID.
    vid effective_L_size = Lmax_id + 1;
     if (effective_L_size > L_param) {
        // This is still a potential issue, as L_param defines the L/R boundary
        std::cerr << "Warning in " << path << ": Max L vertex ID (" << Lmax_id << ") >= L_param (" << L_param
                  << "). L/R definition might be inconsistent." << std::endl;
        effective_L_size = Lmax_id + 1;
     } else {
         // If max L ID is < L_param, we still need L2R to be size L_param
         // if other subgraphs might reference nodes up to L_param-1.
         // Safest is usually to use L_param unless known otherwise. Let's stick to L_param.
          effective_L_size = L_param;
          // However, if Lmax_id+1 > L_param, we need the larger size. So:
          // effective_L_size = std::max(L_param, Lmax_id + 1); // Reconsider this based on global knowledge needs
          // Sticking to original logic for now:
          if (Lmax_id + 1 > L_param) effective_L_size = Lmax_id + 1; else effective_L_size = L_param;


     }
    size_t effective_R_vec_size = (first_r_node_found) ? (Rmax_idx + 1) : 0;

    Graph G;
    G.L = effective_L_size; // L size for vector allocation
    G.L_param_original = L_param; // The crucial boundary value
    G.L2R.assign(G.L, {});
    if (effective_R_vec_size > 0) G.R2L.assign(effective_R_vec_size, {});
    else G.R2L.clear();

    for (auto const& [u, v] : edges) {
        size_t r_idx = v - L_param;
        // Use allocated vector sizes for bounds checking
        if (u < G.L2R.size() && r_idx < G.R2L.size()) {
            G.L2R[u].push_back(v);
            G.R2L[r_idx].push_back(u);
        } else {
             std::cerr << "Internal Error in " << path << ": Indices (" << u << "," << v << ") -> (u=" << u << ", r_idx=" << r_idx
                       << ") out of bounds (L size=" << G.L2R.size() << ", R size=" << G.R2L.size() << ")." << std::endl;
        }
    }
    return G;
}


// ───────────────────── Ranking options ───────────────────────
// Operates on the Graph G passed to it (which is now a local subgraph in MPI case)
enum class RankType { ID, DEGREE, DEGCORE, RANDOM };
std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1) {
    vid n = G.L; // L size is based on loaded (sub)graph data and L_param
    if (n == 0) return {};

    std::vector<vid> ord(n);
    std::iota(ord.begin(), ord.end(), 0); // Fills with 0 to n-1

    // --- IMPORTANT: Ranking is now based on LOCAL subgraph structure ---
    // Degree/Degcore ranks might differ significantly from global ranks.

    if (typ == RankType::ID) return ord; // Still consistent

    if (typ == RankType::DEGREE) {
        std::sort(ord.begin(), ord.end(),
                  [&](vid a, vid b) {
                      // Degrees are based on LOCAL adjacencies in G.L2R
                      size_t deg_a = (a < G.L2R.size()) ? G.L2R[a].size() : 0;
                      size_t deg_b = (b < G.L2R.size()) ? G.L2R[b].size() : 0;
                      if (deg_a != deg_b) return deg_a < deg_b;
                      return a < b;
                  });
        return ord;
    }

    if (typ == RankType::RANDOM) {
        std::mt19937_64 rng(seed);
        std::shuffle(ord.begin(), ord.end(), rng); // Locally shuffled
        return ord;
    }

    if (typ == RankType::DEGCORE) {
        // Degcore calculation is based on LOCAL degrees and adjacencies
        std::vector<vid> deg(n);
        vid max_deg = 0;
        for (vid u = 0; u < n; ++u) {
            deg[u] = (u < G.L2R.size()) ? G.L2R[u].size() : 0;
            max_deg = std::max(max_deg, deg[u]);
        }
        // ... rest of Degcore logic is the same, but operates on local data ...
        std::vector<std::vector<vid>> bins(max_deg + 1);
        std::vector<vid> pos(n);
        for (vid u = 0; u < n; ++u) {
             if (deg[u] >= bins.size()) { throw std::runtime_error("Internal degcore degree error."); }
            pos[u] = bins[deg[u]].size(); bins[deg[u]].push_back(u);
        }
        std::vector<vid> core_order(n); vid current_rank = 0;
        for (vid d = 0; d <= max_deg; ++d) {
            size_t bin_ptr = 0;
            while(bin_ptr < bins[d].size()) {
                vid u = bins[d][bin_ptr++]; if (u >= n) continue;
                if (current_rank >= n) { throw std::runtime_error("Degcore rank overflow."); }
                core_order[current_rank++] = u;
                if (u >= G.L2R.size()) continue;
                for (vid r_orig : G.L2R[u]) {
                    vid L_offset = G.L_param_original; if (r_orig < L_offset) continue;
                    size_t r_idx = r_orig - L_offset; if (r_idx >= G.R2L.size()) continue;
                    for (vid w : G.R2L[r_idx]) {
                        if (w >= n) continue;
                        if (deg[w] > d) {
                            vid current_deg_w = deg[w]; if (current_deg_w >= bins.size() || current_deg_w == 0) continue;
                            size_t pos_w = pos[w];
                            if(pos_w >= bins[current_deg_w].size() || bins[current_deg_w][pos_w] != w) { /* find fallback */ }
                            vid last_node_in_bin = bins[current_deg_w].back(); bins[current_deg_w][pos_w] = last_node_in_bin;
                            pos[last_node_in_bin] = pos_w; bins[current_deg_w].pop_back(); deg[w]--;
                            pos[w] = bins[deg[w]].size(); bins[deg[w]].push_back(w);
                        }
                    }
                }
            }
        }
        if (current_rank != n) { std::cerr << "Warning: Local Degcore rank mismatch." << std::endl; }
        return core_order;
    }
    throw std::logic_error("Invalid RankType");
}


// ───────────────────── Aggregators ───────────────────────────
// (AggHash, AggSort, AggBatch, AggHisto structs remain exactly the same)
// ... [omitted for brevity - paste the aggregator structs from previous correct version here] ...
// HASH – incremental choose-two
struct AggHash {
    std::unordered_map<pair64, uint64_t> h; // Use uint64_t for counts
    uint64_t total = 0;
    AggHash()=default;
    explicit AggHash(size_t L_size){ /* Optional: h.reserve(estimate); */ }
    void add(pair64 p) {
        uint64_t& k = h[p]; // k is current count (0 if new)
        total += k;         // Add #butterflies involving previous k occurrences
        ++k;                // Increment count
    }
    void flush(){ /* DO NOTHING */ } void final_flush() { /* Optional clear */ }
    size_t mem() const { return h.bucket_count() * 1.5 * (sizeof(void*) + sizeof(pair64) + sizeof(uint64_t)) + h.size() * (sizeof(pair64) + sizeof(uint64_t)); }
};
// SORT – sort buffer at flush (only final flush is correct)
struct AggSort {
    std::vector<pair64> buf; uint64_t total = 0;
    AggSort()=default; explicit AggSort(size_t L_size){ /* reserve */ }
    void add(pair64 p){ buf.push_back(p); }
    void flush(){
        if (buf.empty()) return; std::sort(buf.begin(), buf.end());
        for (size_t i = 0; i < buf.size(); ) {
            size_t j = i; while (j < buf.size() && buf[j] == buf[i]) ++j;
            uint64_t f = j - i; if (f > 1) total += f * (f - 1) / 2; i = j;
        }
        buf.clear();
    }
    void final_flush() { flush(); }
    size_t mem() const { return buf.size() * sizeof(pair64) + sizeof(*this); }
};
// BATCH – Now behaves identically to AggSort
struct AggBatch : AggSort {
    AggBatch()=default; explicit AggBatch(size_t s):AggSort(s){ /* reserve */ }
    void add(pair64 p){ buf.push_back(p); /* NO internal flush */ }
};
// HISTO – incremental choose-two with dense array
struct AggHisto {
    std::vector<uint32_t> hist; std::vector<vid> touched; uint64_t total = 0; size_t L_size_internal = 0;
    explicit AggHisto(size_t L){
        if (L == 0) { L_size_internal = 0; return; } L_size_internal = L;
        try { hist.assign(L, 0); } catch (const std::bad_alloc& e) { throw std::runtime_error("Histo alloc failed"); }
    }
    void add(pair64 p){
        if (L_size_internal == 0) return; vid w = (vid)(p & 0xFFFFFFFFu);
        if (w >= L_size_internal) { throw std::runtime_error("Histo index OOB"); }
        if (hist[w] == 0) { touched.push_back(w); } total += hist[w]; hist[w]++;
    }
    void flush(){ for (vid w_idx : touched) { if (w_idx < L_size_internal) hist[w_idx] = 0; } touched.clear(); }
    void final_flush() { flush(); }
    size_t mem() const { return hist.size()*sizeof(uint32_t) + touched.size()*sizeof(vid) + sizeof(*this); }
};

// ──────────────── Parallel driver (SIMPLIFIED - no MPI range) ────────────────
template<class Agg>
uint64_t count_par(
    const Graph& G,                     // Now potentially a subgraph
    const std::vector<vid>& rank_of,    // Based on local subgraph ranking
    int nT                              // OpenMP threads per process
) {
    if (nT <= 0) nT = 1;
    if (G.L == 0 || G.R2L.empty()) return 0; // No work if no L nodes or no R nodes in this subgraph

    std::vector<Agg> local(nT, Agg(G.L));
    std::vector<std::vector<pair64>> inbox(nT);
    std::vector<omp_lock_t> locks(nT);
    for (int t = 0; t < nT; ++t) { omp_init_lock(&locks[t]); }

    auto owner = [&](vid u, vid w){ // Owner calculation based on (global) node IDs
        pair64 key = pack(u, w);
        uint64_t hash = 0xcbf29ce484222325ULL;
        const uint64_t prime = 0x100000001b3ULL;
        unsigned char* p = reinterpret_cast<unsigned char*>(&key);
        for(size_t i=0; i<sizeof(key); ++i) { hash ^= p[i]; hash *= prime; }
        return hash % nT;
    };

    // Peeling step: sort R-nodes by increasing degree
    std::vector<size_t> r_order(G.R2L.size());
    std::iota(r_order.begin(), r_order.end(), 0);
    std::sort(r_order.begin(), r_order.end(), [&](size_t a, size_t b) {
        return G.R2L[a].size() < G.R2L[b].size();
    });

    std::atomic<bool> critical_error_flag = false;

#pragma omp parallel num_threads(nT)
    {
        int tid = omp_get_thread_num();
        auto& A = local[tid];

        try {
            // *** Iterate over R-nodes in peeling order ***
#pragma omp for schedule(dynamic, 64) nowait
            for (size_t idx = 0; idx < r_order.size(); ++idx) {
                if (critical_error_flag.load()) continue;

                size_t rid = r_order[idx];
                const auto& aSet = G.R2L[rid];
                if (aSet.size() < 2) continue;

                for (size_t i = 0; i < aSet.size(); ++i) {
                    vid u0 = aSet[i];
                    // Use G.L for rank_of size check, as rank_of is sized based on the loaded G
                    if (u0 >= G.L) { // Check if u0 is within the L-partition boundary *and* vector bounds
                         if (u0 >= rank_of.size()) { // Check vector bounds explicitly
                            #pragma omp critical
                            { std::cerr << "Error: u0=" << u0 << " OOB rank_of (size=" << rank_of.size() << ", G.L=" << G.L << ")" << std::endl; }
                            critical_error_flag = true; continue;
                         }
                          // else: u0 >= G.L but < rank_of.size(). This implies G.L might be smaller than L_param? Log warning?
                          // This case shouldn't happen if G.L is set correctly based on L_param and max_L_id.
                    }


                    for (size_t j = i + 1; j < aSet.size(); ++j) {
                        vid w0 = aSet[j];
                         if (w0 >= G.L) { // Similar check for w0
                             if (w0 >= rank_of.size()) {
                                #pragma omp critical
                                { std::cerr << "Error: w0=" << w0 << " OOB rank_of (size=" << rank_of.size() << ", G.L=" << G.L << ")" << std::endl; }
                                critical_error_flag = true; continue;
                             }
                         }


                        vid u = u0, w = w0;
                        // Use rank_of for ordering. Requires u, w to be valid indices.
                        if (rank_of[u] > rank_of[w]) { std::swap(u, w); }

                        pair64 p = pack(u, w);
                        int o = owner(u, w);

                        if (o == tid) { A.add(p); }
                        else {
                            omp_set_lock(&locks[o]);
                            try { inbox[o].push_back(p); }
                            catch (const std::bad_alloc&) { critical_error_flag = true; /* handle OOM */ }
                            omp_unset_lock(&locks[o]);
                        }
                        if (critical_error_flag.load()) break;
                    }
                    if (critical_error_flag.load()) break;
                }
            } // end omp for

            if (!critical_error_flag.load()) {
#pragma omp barrier
                std::vector<pair64> local_inbox_copy;
                omp_set_lock(&locks[tid]); std::swap(local_inbox_copy, inbox[tid]); omp_unset_lock(&locks[tid]);
                for (const auto& p : local_inbox_copy) { if (critical_error_flag.load()) break; A.add(p); }
            }
            if (!critical_error_flag.load()) { A.final_flush(); }

        } catch (const std::exception& e) { critical_error_flag = true; /* handle error */ }
        catch (...) { critical_error_flag = true; /* handle unknown error */ }
    } // end omp parallel

    for (int t = 0; t < nT; ++t) { omp_destroy_lock(&locks[t]); }

    if (critical_error_flag.load()) { throw std::runtime_error("OMP computation error."); }

    uint64_t total_butterflies = 0;
    for (int t = 0; t < nT; ++t) { total_butterflies += local[t].total; }
    return total_butterflies;
}

// ───────────────────── CLI parsing ──────────────────────────
enum class ModeType { SEQ, PAR };
struct Params {
    ModeType mode = ModeType::PAR;
    std::string file;           // For standalone executable
    std::string subgraph_dir;   // For MPI executable
    vid L = 0;
    RankType rank = RankType::ID;
    uint64_t seed = 1;
    int agg = 0;
    int threads = 1;
};

// Common parsing logic, validation happens in specific main functions
Params parse(int argc, char** argv) {
    Params P;
    std::vector<std::string> args(argv + 1, argv + argc);
    bool L_set = false;

    for (size_t i = 0; i < args.size(); ++i) {
        std::string a = args[i];
        std::string next_arg = (i + 1 < args.size()) ? args[i+1] : "";

        try {
            if (a == "--mode" && !next_arg.empty()) {
                 if (next_arg == "seq") P.mode = ModeType::SEQ; else P.mode = ModeType::PAR; i++;
            } else if (a == "--file" && !next_arg.empty()) {
                 P.file = next_arg; i++; // Standalone will use this
            } else if (a == "--subgraph_dir" && !next_arg.empty()) {
                 P.subgraph_dir = next_arg; i++; // MPI will use this
            } else if (a == "--L" && !next_arg.empty()) {
                 P.L = std::stoul(next_arg); L_set = true; i++;
            } else if (a == "--rank" && !next_arg.empty()) {
                 std::string v = next_arg; /* Set P.rank */
                 if (v=="id") P.rank = RankType::ID; else if (v=="degree") P.rank = RankType::DEGREE;
                 else if (v=="degcore") P.rank = RankType::DEGCORE; else if (v.rfind("random",0)==0){/*set random + seed*/}
                 else { /* warn */ } i++;
            } else if (a == "--agg" && !next_arg.empty()) {
                 std::string v = next_arg; /* Set P.agg */
                 if (v=="hash") P.agg=0; else if (v=="sort") P.agg=1; else if (v=="batch") P.agg=2; else if (v=="histo") P.agg=3;
                 else { /* warn */ } i++;
            } else if (a == "--threads" && !next_arg.empty()) {
                 int nt = std::stoi(next_arg); P.threads = (nt >= 0) ? nt : 1; i++;
            }
            // Allow positional args for convenience? Requires careful order checking.
            // Example: Assume first non-flag is file/dir, second is L
             else if (P.file.empty() && P.subgraph_dir.empty() && a.rfind("--", 0) != 0 && !a.empty()) {
                 // Could be file OR dir - decide based on which executable is running
                 // Store in BOTH for now? Or check std::filesystem::is_directory? Risky here.
                 // Let's require flags --file or --subgraph_dir for clarity.
                  std::cerr << "Warning: Use --file or --subgraph_dir explicitly. Ignoring positional: " << a << std::endl;
             } else if (L_set == false && a.rfind("--", 0) != 0 && !a.empty()) {
                 // Assume second non-flag is L if file/dir was already set by flag
                 // P.L = std::stoul(a); L_set = true; // Disabled - require --L flag
                 std::cerr << "Warning: Use --L explicitly. Ignoring positional: " << a << std::endl;
             }
            else { /* Handle unknown flags/args */ }
        } catch (const std::exception& e) { /* Handle parsing errors */ }
    }

    // Final thread count adjustment (common logic)
    if (P.mode == ModeType::SEQ) P.threads = 1;
    else { // PAR
        if (P.threads <= 0) {
#ifdef _OPENMP
            P.threads = omp_get_max_threads(); if (P.threads <= 0) P.threads = 1;
#else
            P.threads = 1;
#endif
        }
    }
     if (P.threads < 1) P.threads = 1; // Ensure at least 1

    // L must always be set
    if (!L_set) throw std::runtime_error("Missing required argument: --L <L_size>");

    return P;
}

// ───────────────────── Timer helper ───────────────────────────
double get_time() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(now.time_since_epoch()).count();
#endif
}

#endif // BUTTER_CORE_HPP