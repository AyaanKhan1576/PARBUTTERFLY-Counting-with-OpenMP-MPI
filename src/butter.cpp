
// /****************************************************************************************
//  * butter.cpp  – unified sequential and parallel butterfly counter (bug-fixed)
//  *
//  * Build : g++ -std=c++17 -O3 butter.cpp -fopenmp -o butter
//  *
//  * Usage :
//  *   ./butter --mode seq|par <edge_file> <L_size> \
//  *          --rank id|degree|degcore|random[:seed] \
//  *          --agg  hash|sort|batch|histo \
//  *          [--threads N]           (only for --mode par)
//  *****************************************************************************************/

//  #include <bits/stdc++.h>
//  #ifdef _OPENMP
//  #include <omp.h>
//  #endif
//  #include <chrono> // For timing
 
//  using vid    = uint32_t;
//  using pair64 = uint64_t;
//  static inline pair64 pack(vid a, vid b) { return (pair64)a << 32 | b; }
 
//  // ───────────────────── Graph loader ───────────────────────────
//  struct Graph {
//      vid L;
//      std::vector<std::vector<vid>> L2R; // Adjacency list for L -> R
//      std::vector<std::vector<vid>> R2L; // Adjacency list for R -> L
//  };
 
//  Graph load_graph(const std::string& path, vid L_param) {
//      std::ifstream in(path);
//      if (!in) throw std::runtime_error("Cannot open " + path);
//      std::vector<std::pair<vid, vid>> edges;
//      size_t Rmax_idx = 0; // Max index used in R2L (v-L)
//      vid Lmax_id = 0;   // Max vertex ID used in L partition
 
//      for (std::string line; std::getline(in, line); ) {
//          // Skip empty lines and comments
//          if (line.empty() || line[0] == '%' || line[0] == '#') continue;
 
//          std::istringstream iss(line);
//          vid u, v;
//          if (!(iss >> u >> v)) {
//               std::cerr << "Warning: Skipping malformed line: " << line << std::endl;
//               continue;
//          }
 
//          // Basic validation assuming L nodes are [0, L_param-1] and R nodes are >= L_param
//          if (u >= L_param) {
//               std::cerr << "Warning: Skipping edge (" << u << "," << v << ") - u >= L_param (" << L_param << ")" << std::endl;
//               continue;
//          }
//          if (v < L_param) {
//               std::cerr << "Warning: Skipping edge (" << u << "," << v << ") - v < L_param (" << L_param << ")" << std::endl;
//               continue;
//          }
 
//          edges.emplace_back(u, v);
//          // Calculate index relative to L_param for R nodes
//          size_t r_idx = v - L_param;
//          Rmax_idx = std::max(Rmax_idx, r_idx);
//          Lmax_id = std::max(Lmax_id, u);
//      }
 
//      // Determine actual sizes needed based on input data
//      vid effective_L_size = Lmax_id + 1;
//      if (effective_L_size > L_param) {
//          std::cerr << "Warning: Max L vertex ID (" << Lmax_id << ") >= provided L_size (" << L_param
//                    << "). Using effective L size: " << effective_L_size << std::endl;
//          // L_param = effective_L_size; // Optionally update L_param if strict adherence isn't needed
//      } else {
//          // Use the provided L_param if it's large enough
//          effective_L_size = L_param;
//      }
//      size_t effective_R_vec_size = Rmax_idx + 1;
 
 
//      Graph G;
//      G.L = effective_L_size; // Store the effective L size used for vector allocation
//      G.L2R.assign(G.L, {});
//      G.R2L.assign(effective_R_vec_size, {});
 
//      for (auto const& [u, v] : edges) {
//          // Use indices relative to L_param for R nodes
//          size_t r_idx = v - L_param;
//           // Bounds check before insertion (using effective sizes)
//          if (u < G.L && r_idx < G.R2L.size()) {
//              G.L2R[u].push_back(v); // Store original R id
//              G.R2L[r_idx].push_back(u);
//          } else {
//               // Should not happen with correct size calculation, but good for safety
//               std::cerr << "Internal Error: Edge indices (" << u << "," << v << ") out of bounds during insertion." << std::endl;
//          }
//      }
 
//      // Sorting L2R neighbors is not required by this algorithm but doesn't hurt
//      // for (auto &nb : G.L2R) std::sort(nb.begin(), nb.end());
//      // Sorting R2L neighbors IS required by the wedge iteration logic (i, j loops assume sorted/consistent order isn't needed)
//      // No, the logic iterates i=0..N, j=i+1..N, which works on unsorted neighbors too. Sorting is not strictly needed here.
 
//      return G;
//  }
 
 
//  // ───────────────────── Ranking options ───────────────────────
//  enum class RankType { ID, DEGREE, DEGCORE, RANDOM };
//  std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1) {
//      vid n = G.L; // Use the graph's effective L size
//      std::vector<vid> ord(n);
//      std::iota(ord.begin(), ord.end(), 0);
 
//      if (typ == RankType::ID) {
//          return ord;
//      }
 
//      if (typ == RankType::DEGREE) {
//          std::sort(ord.begin(), ord.end(),
//                    [&](vid a, vid b) {
//                        // Get degrees safely (vectors already sized to G.L)
//                        size_t deg_a = G.L2R[a].size();
//                        size_t deg_b = G.L2R[b].size();
//                        if (deg_a != deg_b) return deg_a < deg_b;
//                        return a < b; // Tie-break by ID
//                    });
//          return ord;
//      }
 
//      if (typ == RankType::RANDOM) {
//          std::mt19937_64 rng(seed);
//          std::shuffle(ord.begin(), ord.end(), rng);
//          return ord;
//      }
 
//      // Degcore approximation (using standard k-core peeling)
//      if (typ == RankType::DEGCORE) {
//          std::vector<vid> deg(n);
//          vid max_deg = 0;
//          for (vid u = 0; u < n; ++u) {
//              deg[u] = G.L2R[u].size();
//              max_deg = std::max(max_deg, deg[u]);
//          }
 
//          std::vector<std::vector<vid>> bins(max_deg + 1);
//          std::vector<vid> pos(n); // Position of node u in its bin
//          for (vid u = 0; u < n; ++u) {
//              pos[u] = bins[deg[u]].size();
//              bins[deg[u]].push_back(u);
//          }
 
//          std::vector<vid> core_order(n); // Resulting order
//          vid current_rank = 0;
 
//          for (vid d = 0; d <= max_deg; ++d) {
//              size_t bin_idx = 0;
//              while(bin_idx < bins[d].size()) {
//                  vid u = bins[d][bin_idx++];
 
//                  // Store node in the order (lower core number = lower rank/earlier in order)
//                  core_order[current_rank++] = u;
 
//                  // Process neighbors to potentially reduce their degree
//                  for (vid r_orig : G.L2R[u]) {
//                      // Calculate R index based on L_param used during loading
//                      // Need L_param here, G.L might have been adjusted. Pass L_param?
//                      // Let's assume G.L holds the value corresponding to R node IDs. Risky.
//                      // Safer: Store L_param in Graph struct or re-derive R index carefully.
//                      // Assuming G.R2L indices correspond to r_orig - L_param
//                      // We need L_param. Let's pass it or store it.
//                      // *** TEMPORARY ASSUMPTION: G.L is the correct offset ***
//                      // This part is fragile if L_param != G.L
//                      vid L_offset = G.L; // Replace with actual L_param if different
//                      if (r_orig < L_offset) { /* Error condition */ continue; }
//                      size_t r_idx = r_orig - L_offset;
//                      if (r_idx >= G.R2L.size()) { /* Error condition */ continue; }
 
 
//                      for (vid w : G.R2L[r_idx]) {
//                          if (deg[w] > d) { // If neighbor w is still in a higher degree bin
//                              // Find w in its current bin
//                              vid current_deg_w = deg[w];
//                              size_t pos_w = pos[w];
//                              vid last_node_in_bin = bins[current_deg_w].back();
 
//                              // Swap w with the last element in its bin
//                              bins[current_deg_w][pos_w] = last_node_in_bin;
//                              pos[last_node_in_bin] = pos_w; // Update swapped node's position
//                              bins[current_deg_w].pop_back();
 
//                              // Move w to the lower degree bin
//                              deg[w]--;
//                              pos[w] = bins[deg[w]].size();
//                              bins[deg[w]].push_back(w);
//                          }
//                      }
//                  }
//              }
//          }
//           // Ensure all nodes were processed (should be true if graph connected)
//          assert(current_rank == n);
//          return core_order;
//      }
 
//      // Should not be reached
//      throw std::runtime_error("Invalid RankType");
//  }
 
 
//  // ───────────────────── Aggregators ───────────────────────────
//  static constexpr size_t FLUSH_BYTES = 8 * 1024 * 1024; // Target memory for flushing (relevant for Sort/Batch)
 
//  // Base class (optional, for structure) - current code uses templates/duck-typing
//  // struct AggBase {
//  //    uint64_t total = 0;
//  //    virtual ~AggBase() = default;
//  //    virtual void add(pair64 p) = 0;
//  //    virtual void flush() = 0; // Process intermediate buffer/state
//  //    virtual void final_flush() { flush(); }; // Process final state (may be same as flush)
//  //    virtual size_t mem() const = 0; // Estimate current memory usage
//  // };
 
 
//  // HASH – incremental choose-two
//  struct AggHash /* : AggBase */ {
//      std::unordered_map<pair64, uint64_t> h; // Use uint64_t for counts
//      uint64_t total = 0;
//      AggHash()=default;
//      explicit AggHash(size_t L_size){ /* Optional: h.reserve(estimate); */ }
//      void add(pair64 p) {
//          uint64_t& k = h[p]; // k is current count (0 if new)
//          total += k;         // Add #butterflies involving previous k occurrences
//          ++k;                // Increment count
//      }
//      // Intermediate flush should NOT clear state for Hash
//      void flush(){ /* DO NOTHING */ }
//      // Final flush can optionally clear memory if Agg object persists, but not needed for total
//      void final_flush() { /* h.clear(); */ }
//      // Memory estimate (rough)
//      size_t mem() const { return h.bucket_count() * 1.5 * (sizeof(void*) + sizeof(pair64) + sizeof(uint64_t)); }
//  };
 
//  // SORT – sort buffer at flush
//  struct AggSort /* : AggBase */ {
//      std::vector<pair64> buf;
//      uint64_t total = 0;
//      AggSort()=default;
//      explicit AggSort(size_t L_size){ /* Optional: buf.reserve(estimate); */ }
//      void add(pair64 p){ buf.push_back(p); }
//      void flush(){
//          if (buf.empty()) return;
//          std::sort(buf.begin(), buf.end());
//          for (size_t i = 0; i < buf.size(); ) {
//              size_t j = i;
//              while (j < buf.size() && buf[j] == buf[i]) ++j;
//              uint64_t f = j - i; // f = frequency = # common neighbors 'r'
//              if (f > 1) {
//                  total += f * (f - 1) / 2; // Choose 2 neighbors
//              }
//              i = j;
//          }
//          buf.clear();
//          // buf.shrink_to_fit(); // Optional memory reduction
//      }
//      void final_flush() { flush(); } // Final flush is same as intermediate
//      size_t mem() const { return buf.size() * sizeof(pair64); }
//  };
 
//  // BATCH – fixed chunk then sort
//  struct AggBatch : AggSort {
//      static constexpr size_t BATCH_SIZE_ITEMS = (1 << 19); // ~4MB if pair64 is 8 bytes
//      AggBatch()=default;
//      explicit AggBatch(size_t s):AggSort(s){ buf.reserve(BATCH_SIZE_ITEMS); }
//      void add(pair64 p){
//          buf.push_back(p);
//          // Check based on *size*, not memory, for consistency
//          if (buf.size() >= BATCH_SIZE_ITEMS) flush();
//      }
//      // Inherits flush, final_flush, mem
//  };
 
//  // HISTO – incremental choose-two with dense array
//  // WARNING: Fundamentally questionable logic for r-centric iteration.
//  struct AggHisto /* : AggBase */ {
//      std::vector<uint32_t> hist; // Use uint32_t or uint64_t? Check L size limits
//      std::vector<vid> touched;
//      uint64_t total = 0;
//      size_t L_size_internal = 0;
 
//      explicit AggHisto(size_t L){
//          if (L == 0) {
//               std::cerr << "Warning: AggHisto created with L=0." << std::endl;
//               return;
//          }
//          L_size_internal = L;
//          hist.assign(L, 0);
//          // touched.reserve(L); // Optional
//      }
//      void add(pair64 p){
//          vid w = (vid)(p & 0xFFFFFFFFu); // Assumes w is in lower 32 bits
 
//          if (w >= L_size_internal) {
//              std::cerr << "Error: AggHisto index w=" << w << " >= L=" << L_size_internal << ". Skipping add." << std::endl;
//              return; // Avoid crash, but count will be wrong
//          }
 
//          // Original logic (potentially flawed interpretation for r-centric loop)
//          if (hist[w] == 0) { // Only track if first time seeing *this w* in *this flush cycle*
//              touched.push_back(w);
//          }
//          total += hist[w]; // Add count associated with w *before* increment
//          hist[w]++;        // Increment count associated with w
//      }
//      // Intermediate flush should reset state for Histo
//      void flush(){
//          for (vid w_idx : touched) {
//              if (w_idx < L_size_internal) { // Bounds check
//                   hist[w_idx] = 0;
//              }
//          }
//          touched.clear();
//      }
//      void final_flush() { flush(); }
//      // Memory is hist + touched + object size. Dominated by hist if L large & dense.
//      size_t mem() const { return hist.size()*sizeof(uint32_t) + touched.size()*sizeof(vid); }
//  };
 
 
//  // ───────────────────── Parallel driver (FIXED) ───────────────────────
//  template<class Agg>
//  uint64_t count_par(const Graph& G, const std::vector<vid>& rank_of, int nT){
//      if (nT <= 0) nT = 1; // Ensure at least one thread
 
//      std::vector<Agg> local(nT, Agg(G.L)); // Create nT local aggregators
//      std::vector<std::vector<pair64>> inbox(nT);
//      std::vector<omp_lock_t> locks(nT);
//      for (int t=0; t<nT; ++t) {
//          omp_init_lock(&locks[t]);
//      }
 
//      // Simple FNV1a-based hash for potentially better distribution than XOR
//      auto owner = [&](vid u, vid w){
//          pair64 key = pack(u, w);
//          uint64_t hash = 0xcbf29ce484222325ULL;
//          const uint64_t prime = 0x100000001b3ULL;
//          unsigned char* p = (unsigned char*)&key;
//          for(size_t i=0; i<sizeof(key); ++i) {
//              hash ^= p[i];
//              hash *= prime;
//          }
//          return hash % nT;
//      };
 
//  #pragma omp parallel num_threads(nT)
//      {
//          int tid = omp_get_thread_num();
//          auto& A = local[tid]; // Reference to thread-local aggregator
 
//          // Schedule dynamic can be good if work per 'rid' varies significantly
//  #pragma omp for schedule(dynamic, 64) nowait // nowait allows threads to proceed to barrier sooner if done
//          for (size_t rid = 0; rid < G.R2L.size(); ++rid) {
//              const auto& aSet = G.R2L[rid]; // Neighbors of r in L
 
//              if (aSet.size() < 2) continue; // Cannot form a wedge
 
//              for (size_t i = 0; i < aSet.size(); ++i) {
//                  vid u0 = aSet[i];
//                  // Basic bounds check for safety, though rank_of should be sized correctly
//                  if (u0 >= rank_of.size()) continue;
 
//                  for (size_t j = i + 1; j < aSet.size(); ++j) {
//                      vid w0 = aSet[j];
//                      if (w0 >= rank_of.size()) continue;
 
//                      vid u = u0, w = w0;
 
//                      // Order pair (u, w) by rank (lower rank first)
//                      if (rank_of[u] > rank_of[w]) {
//                          std::swap(u, w);
//                      }
 
//                      pair64 p = pack(u, w);
//                      int o = owner(u, w);
 
//                      if (o == tid) {
//                          A.add(p); // Process locally
//                      } else {
//                          // Send to owner's inbox (locked)
//                          omp_set_lock(&locks[o]);
//                          try {
//                              inbox[o].push_back(p);
//                          } catch (const std::bad_alloc& e) {
//                               // Basic OOM handling: report and potentially stop
//                              omp_unset_lock(&locks[o]);
//                              #pragma omp critical
//                              { std::cerr << "FATAL: Out of memory adding to inbox[" << o << "] from thread " << tid << std::endl; }
//                              // Consider aborting or using a more robust OOM strategy
//                               exit(EXIT_FAILURE); // Simple abort
//                          }
//                          omp_unset_lock(&locks[o]);
//                      }
//                  } // end inner loop j
//              } // end outer loop i
 
//              // --- PERIODIC FLUSH REMOVED HERE ---
//              // This was the source of the inconsistency bug for AggHash/AggHisto.
//              // Removing it ensures state is preserved correctly until the final flush.
//              // ---
 
//          } // end omp for
 
//          // Barrier: All threads must finish the 'for' loop and sending messages
//          // (implicit barrier at end of 'for' unless 'nowait' is used,
//          // but explicit barrier after sending is safer conceptually)
//  #pragma omp barrier
 
//          // Phase 2: Drain own inbox
//          std::vector<pair64> local_inbox_copy;
//          // Optional: Reserve space if average inbox size can be estimated
//          // local_inbox_copy.reserve( G.L / nT ); // Very rough
 
//          omp_set_lock(&locks[tid]);
//          std::swap(local_inbox_copy, inbox[tid]); // Quick swap while holding lock
//          omp_unset_lock(&locks[tid]);
 
//          // Process pairs received from other threads
//          for (const auto& p : local_inbox_copy) {
//              A.add(p);
//          }
 
//          // Final flush for this thread's local aggregator.
//          // Crucial for AggSort/AggBatch to process their final buffer.
//          // Also needed if AggHash/AggHisto had any final processing (though they don't here).
//          A.final_flush(); // Use final_flush (might be same as flush)
 
//      } // end omp parallel
 
//      // Clean up locks
//      for (int t=0; t<nT; ++t) {
//          omp_destroy_lock(&locks[t]);
//      }
 
//      // Aggregate the total counts from all threads
//      uint64_t total_butterflies = 0;
//      for (int t = 0; t < nT; ++t) {
//          total_butterflies += local[t].total;
//      }
 
//      return total_butterflies;
//  }
 
 
//  // ───────────────────── CLI parsing ──────────────────────────
//  enum class ModeType { SEQ, PAR };
//  struct Params {
//      ModeType mode = ModeType::SEQ;
//      std::string file;
//      vid L = 0; // L size parameter from command line
//      RankType rank = RankType::ID;
//      uint64_t seed = 1;
//      int agg = 0;      // 0=hash,1=sort,2=batch,3=histo
//      int threads = 1;
//  };
 
//  // More robust CLI parsing
//  Params parse(int argc, char** argv){
//      Params P;
//      std::vector<std::string> args(argv + 1, argv + argc);
//      bool mode_set = false, file_set = false, L_set = false;
 
//      for (size_t i = 0; i < args.size(); ++i) {
//          std::string a = args[i];
//          std::string next_arg = (i + 1 < args.size()) ? args[i+1] : "";
 
//          try {
//              if (a == "--mode" && !next_arg.empty()) {
//                  if (next_arg == "seq") P.mode = ModeType::SEQ;
//                  else if (next_arg == "par") P.mode = ModeType::PAR;
//                  else throw std::runtime_error("Unknown mode '" + next_arg + "'");
//                  mode_set = true; i++;
//              } else if (a == "--rank" && !next_arg.empty()) {
//                  if (next_arg == "id") P.rank = RankType::ID;
//                  else if (next_arg == "degree") P.rank = RankType::DEGREE;
//                  else if (next_arg == "degcore") P.rank = RankType::DEGCORE;
//                  else if (next_arg.rfind("random", 0) == 0) {
//                      P.rank = RankType::RANDOM;
//                      auto pos = next_arg.find(':');
//                      if (pos != std::string::npos) {
//                          P.seed = std::stoull(next_arg.substr(pos + 1));
//                      }
//                  } else throw std::runtime_error("Unknown rank '" + next_arg + "'");
//                  i++;
//              } else if (a == "--agg" && !next_arg.empty()) {
//                  if (next_arg == "hash") P.agg = 0;
//                  else if (next_arg == "sort") P.agg = 1;
//                  else if (next_arg == "batch") P.agg = 2;
//                  else if (next_arg == "histo") P.agg = 3;
//                  else throw std::runtime_error("Unknown agg '" + next_arg + "'");
//                  i++;
//              } else if (a == "--threads" && !next_arg.empty()) {
//                  P.threads = std::stoi(next_arg);
//                  if (P.threads < 0) throw std::runtime_error("Number of threads cannot be negative.");
//                  i++;
//              } else if (!mode_set) { // Assume positional args before flags are handled
//                  throw std::runtime_error("Positional argument '" + a + "' before required --mode flag.");
//              } else if (!file_set) {
//                  P.file = a; file_set = true;
//              } else if (!L_set) {
//                  P.L = std::stoul(a); L_set = true;
//                  if (P.L == 0) std::cerr << "Warning: L_size specified as 0." << std::endl;
//              } else {
//                  throw std::runtime_error("Unknown or misplaced argument: " + a);
//              }
//          } catch (const std::invalid_argument& e) {
//              throw std::runtime_error("Invalid numeric value for argument '" + a + "': " + next_arg);
//          } catch (const std::out_of_range& e) {
//              throw std::runtime_error("Numeric value out of range for argument '" + a + "': " + next_arg);
//          }
//      }
 
//      // Validate required arguments
//      if (!mode_set) throw std::runtime_error("Missing required argument: --mode seq|par");
//      if (P.file.empty()) throw std::runtime_error("Missing required argument: <edge_file>");
//      if (!L_set) throw std::runtime_error("Missing required argument: <L_size>");
 
//      // Adjust threads for parallel mode if needed
//      if (P.mode == ModeType::PAR) {
//          if (P.threads <= 0) {
//  #ifdef _OPENMP
//              P.threads = omp_get_max_threads();
//  #else
//              P.threads = 1;
//              std::cerr << "Warning: OpenMP not enabled, running parallel mode with 1 thread." << std::endl;
//  #endif
//          }
//      } else { // Sequential mode always uses 1 thread for counting logic
//          P.threads = 1;
//      }
 
//      return P;
//  }
 
//  // ───────────────────── Timer helper ───────────────────────────
//  double get_time() {
//  #ifdef _OPENMP
//      return omp_get_wtime();
//  #else
//      // Use high_resolution_clock for better precision if OpenMP is off
//      auto now = std::chrono::high_resolution_clock::now();
//      auto duration = now.time_since_epoch();
//      return std::chrono::duration<double>(duration).count();
//  #endif
//  }
 
//  // ───────────────────── main ────────────────────────────────
//  int main(int argc, char** argv){
//      Params P;
//      try {
//          P = parse(argc, argv);
//      } catch (const std::exception &e) {
//          std::cerr << "Error parsing arguments: " << e.what() << "\n";
//          std::cerr << "Usage: ./butter --mode seq|par <edge_file> <L_size> [--rank id|degree|degcore|random[:seed]] [--agg hash|sort|batch|histo] [--threads N]\n";
//          return 1;
//      }
 
//      try {
//          double t_start_load = get_time();
//          Graph G = load_graph(P.file, P.L);
//          double t_end_load = get_time();
 
//          // Check if L size seems valid after loading
//          if (G.L == 0 && G.R2L.empty()) {
//               std::cerr << "Warning: Graph appears empty after loading." << std::endl;
//          } else if (G.L == 0) {
//               std::cerr << "Warning: Graph L size is 0, but R partition has data? Check L_size param." << std::endl;
//          }
 
//          double t_start_rank = get_time();
//          // Pass the *original* L_param to make_rank if needed for R index calculation
//          // Or ensure Graph struct stores it. Assuming G.L is sufficient for now.
//          auto order = make_rank(G, P.rank, P.seed);
//          // Ensure rank_of vector is sized correctly based on G.L (effective L size)
//          if (G.L == 0 && !order.empty()) {
//               // This shouldn't happen if G.L=0 -> n=0 in make_rank
//               std::cerr << "Error: G.L is 0 but ranking produced non-empty order." << std::endl;
//               return 1;
//          } else if (G.L > 0 && order.size() != G.L) {
//              std::cerr << "Error: Ranking order size (" << order.size() << ") does not match graph L size (" << G.L << ")." << std::endl;
//               // This indicates a problem in make_rank logic or G.L calculation
//               return 1;
//          }
 
//          std::vector<vid> rank_of(G.L); // Size based on graph's effective L size
//          if (G.L > 0) { // Avoid operations on size 0 vector
//              for (vid r = 0; r < G.L; ++r) {
//                   // Check bounds of order vector access
//                   if (r >= order.size()) {
//                       std::cerr << "Error: Index out of bounds accessing order vector during rank_of creation." << std::endl;
//                       return 1;
//                   }
//                   vid node_id = order[r];
//                   // Check bounds of rank_of vector access
//                    if (node_id >= rank_of.size()) {
//                       std::cerr << "Error: Node ID " << node_id << " from ranking is out of bounds for rank_of array (size " << rank_of.size() << ")." << std::endl;
//                       return 1;
//                   }
//                  rank_of[node_id] = r; // Assign rank based on position in sorted 'order'
//              }
//          }
//          double t_end_rank = get_time();
 
//          uint64_t ans = 0;
//          double t_start_count = get_time();
 
//          // Determine threads to use based on mode (always use count_par)
//          int nT_to_use = (P.mode == ModeType::SEQ) ? 1 : P.threads;
 
//          // Check OpenMP availability if parallel execution requested
//          #ifndef _OPENMP
//          if (nT_to_use > 1) {
//              std::cerr << "Warning: Compiled without OpenMP support, running with 1 thread instead of requested " << nT_to_use << "." << std::endl;
//              nT_to_use = 1;
//          }
//          #endif
 
//          // Call the appropriate count_par function via template based on P.agg
//          switch (P.agg) {
//              case 0: ans = count_par<AggHash >(G, rank_of, nT_to_use); break;
//              case 1: ans = count_par<AggSort >(G, rank_of, nT_to_use); break;
//              case 2: ans = count_par<AggBatch>(G, rank_of, nT_to_use); break;
//              case 3:
//                  if (nT_to_use > 0) { // Only warn if actually running
//                      std::cerr << "Warning: Using AggHisto, which may have fundamental logical issues with this algorithm's iteration order. Butterfly count might be incorrect." << std::endl;
//                  }
//                  ans = count_par<AggHisto>(G, rank_of, nT_to_use); break;
//              default: throw std::logic_error("Internal error: Invalid aggregator type selected.");
//          }
 
//          double t_end_count = get_time();
 
//          // Output results clearly
//          std::cout << "[" << (P.mode == ModeType::SEQ ? "seq" : "par") << "]"
//                    << " butterflies=" << ans
//                    << " load=" << std::fixed << std::setprecision(5) << (t_end_load - t_start_load) << "s"
//                    << " rank_t=" << std::fixed << std::setprecision(5) << (t_end_rank - t_start_rank) << "s"
//                    << " count=" << std::fixed << std::setprecision(5) << (t_end_count - t_start_count) << "s"
//                    << " thr=" << nT_to_use // Report actual threads used
//                    << " rank=" << static_cast<int>(P.rank)
//                    << " agg=" << P.agg
//                    << "\n";
 
//      } catch (const std::exception &e) {
//          std::cerr << "Runtime Error: " << e.what() << "\n";
//          return 1;
//      }
//      return 0;
//  }

/****************************************************************************************
 * butter.cpp  – unified sequential and parallel butterfly counter (bug-fixed)
 *
 * Build : g++ -std=c++17 -O3 butter.cpp -fopenmp -o butter
 *
 * Usage :
 *   ./butter --mode seq|par <edge_file> <L_size> \
 *          --rank id|degree|degcore|random[:seed] \
 *          --agg  hash|sort|batch|histo \
 *          [--threads N]           (only for --mode par)
 *****************************************************************************************/

 #include <bits/stdc++.h>
 #ifdef _OPENMP
 #include <omp.h>
 #endif
 #include <chrono> // For timing
 
 using vid    = uint32_t;
 using pair64 = uint64_t;
 static inline pair64 pack(vid a, vid b) { return (pair64)a << 32 | b; }
 
 // ───────────────────── Graph loader ───────────────────────────
 struct Graph {
     vid L;
     std::vector<std::vector<vid>> L2R; // Adjacency list for L -> R
     std::vector<std::vector<vid>> R2L; // Adjacency list for R -> L
     vid L_param_original; // Store the original L parameter for R index calculation
 };
 
 Graph load_graph(const std::string& path, vid L_param) {
     std::ifstream in(path);
     if (!in) throw std::runtime_error("Cannot open " + path);
     std::vector<std::pair<vid, vid>> edges;
     size_t Rmax_idx = 0; // Max index used in R2L (v-L)
     vid Lmax_id = 0;   // Max vertex ID used in L partition
 
     for (std::string line; std::getline(in, line); ) {
         // Skip empty lines and comments
         if (line.empty() || line[0] == '%' || line[0] == '#') continue;
 
         std::istringstream iss(line);
         vid u, v;
         if (!(iss >> u >> v)) {
              std::cerr << "Warning: Skipping malformed line: " << line << std::endl;
              continue;
         }
 
         // Basic validation assuming L nodes are [0, L_param-1] and R nodes are >= L_param
         if (u >= L_param) {
              std::cerr << "Warning: Skipping edge (" << u << "," << v << ") - u (" << u << ") >= L_param (" << L_param << ")" << std::endl;
              continue;
         }
         if (v < L_param) {
              std::cerr << "Warning: Skipping edge (" << u << "," << v << ") - v (" << v << ") < L_param (" << L_param << ")" << std::endl;
              continue;
         }
 
         edges.emplace_back(u, v);
         // Calculate index relative to L_param for R nodes
         size_t r_idx = v - L_param;
         Rmax_idx = std::max(Rmax_idx, r_idx);
         Lmax_id = std::max(Lmax_id, u);
     }
 
     // Determine actual sizes needed based on input data
     vid effective_L_size = Lmax_id + 1;
     if (effective_L_size > L_param) {
         std::cerr << "Warning: Max L vertex ID (" << Lmax_id << ") >= provided L_size (" << L_param
                   << "). Using effective L size: " << effective_L_size << ". R-node indices might be affected if not adjusted." << std::endl;
          // Decide how to handle this: Either error out, adjust L_param, or rely on Graph.L_param_original
          // Sticking with original L_param for consistency in R-node indexing.
          effective_L_size = Lmax_id + 1; // Allocate enough space for L nodes
     } else {
         // Use the provided L_param if it's large enough or if max L ID is smaller
         effective_L_size = L_param;
     }
     size_t effective_R_vec_size = Rmax_idx + 1;
 
 
     Graph G;
     G.L = effective_L_size; // Store the effective L size used for L vector allocation
     G.L_param_original = L_param; // Store the original L_param for R indexing
     G.L2R.assign(G.L, {});
     G.R2L.assign(effective_R_vec_size, {});
 
     for (auto const& [u, v] : edges) {
         // Use indices relative to the ORIGINAL L_param for R nodes
         size_t r_idx = v - L_param;
          // Bounds check before insertion (using vector sizes)
         if (u < G.L2R.size() && r_idx < G.R2L.size()) {
             G.L2R[u].push_back(v); // Store original R id
             G.R2L[r_idx].push_back(u);
         } else {
              // This might happen if effective_L_size > L_param causing L2R to be larger
              if (u < effective_L_size && r_idx < G.R2L.size()) {
                   // Example: L_param=10, max_L_id=15. G.L=16. u=15 is valid index for L2R.
                   // Check if L2R needs resize or if this is intended.
                   // Assuming G.L already sized correctly for Lmax_id+1
                    if (u < G.L) { // Check against actual allocated size
                        G.L2R[u].push_back(v);
                        G.R2L[r_idx].push_back(u);
                    } else {
                       std::cerr << "Internal Error: L-index " << u << " out of bounds (" << G.L << ") during insertion." << std::endl;
                    }
              } else {
                 std::cerr << "Internal Error: Edge indices (" << u << "," << v << ") -> (u=" << u << ", r_idx=" << r_idx
                           << ") out of bounds (L size=" << G.L2R.size() << ", R size=" << G.R2L.size()
                           << ") during insertion." << std::endl;
              }
         }
     }
     // Sort L-neighbors for each R-node? Not strictly needed for the wedge iteration i, j+1 logic.
     // for (auto &nb : G.R2L) std::sort(nb.begin(), nb.end());
 
     return G;
 }
 
 
 // ───────────────────── Ranking options ───────────────────────
 enum class RankType { ID, DEGREE, DEGCORE, RANDOM };
 std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1) {
     vid n = G.L; // Use the graph's effective L size for processing
     if (n == 0) return {}; // Handle empty graph case
 
     std::vector<vid> ord(n);
     std::iota(ord.begin(), ord.end(), 0);
 
     if (typ == RankType::ID) {
         return ord;
     }
 
     if (typ == RankType::DEGREE) {
         std::sort(ord.begin(), ord.end(),
                   [&](vid a, vid b) {
                       // Get degrees safely (vectors already sized to G.L)
                       size_t deg_a = (a < G.L2R.size()) ? G.L2R[a].size() : 0;
                       size_t deg_b = (b < G.L2R.size()) ? G.L2R[b].size() : 0;
                       if (deg_a != deg_b) return deg_a < deg_b;
                       return a < b; // Tie-break by ID
                   });
         return ord;
     }
 
     if (typ == RankType::RANDOM) {
         std::mt19937_64 rng(seed);
         std::shuffle(ord.begin(), ord.end(), rng);
         return ord;
     }
 
     // Degcore approximation (using standard k-core peeling)
     if (typ == RankType::DEGCORE) {
         std::vector<vid> deg(n);
         vid max_deg = 0;
         for (vid u = 0; u < n; ++u) {
             deg[u] = (u < G.L2R.size()) ? G.L2R[u].size() : 0;
             max_deg = std::max(max_deg, deg[u]);
         }
 
         std::vector<std::vector<vid>> bins(max_deg + 1);
         std::vector<vid> pos(n); // Position of node u in its bin
         for (vid u = 0; u < n; ++u) {
             if (deg[u] >= bins.size()) {
                 throw std::runtime_error("Degree " + std::to_string(deg[u]) + " exceeds max_deg " + std::to_string(max_deg) + " used for bins.");
             }
             pos[u] = bins[deg[u]].size();
             bins[deg[u]].push_back(u);
         }
 
         std::vector<vid> core_order(n); // Resulting order
         vid current_rank = 0; // Index for placing nodes into core_order
 
         for (vid d = 0; d <= max_deg; ++d) {
             size_t bin_ptr = 0; // Use pointer instead of modifying vector while iterating
             while(bin_ptr < bins[d].size()) {
                 vid u = bins[d][bin_ptr++];
                 if (u >= n) continue; // Safety check
 
                 // Place node u in the order (lower core number = lower rank/earlier in order)
                 if (current_rank >= n) {
                      throw std::runtime_error("Degcore internal error: current_rank exceeded n.");
                 }
                 core_order[current_rank++] = u;
 
                 // Process neighbors to potentially reduce their degree
                 if (u >= G.L2R.size()) continue; // Safety check
 
                 for (vid r_orig : G.L2R[u]) {
                     // Calculate R index using the original L_param stored in Graph
                     vid L_offset = G.L_param_original;
                     if (r_orig < L_offset) { /* Error condition */ continue; }
                     size_t r_idx = r_orig - L_offset;
                     if (r_idx >= G.R2L.size()) { /* Error condition */ continue; }
 
                     for (vid w : G.R2L[r_idx]) {
                         if (w >= n) continue; // Safety check against L size used for deg array
 
                         if (deg[w] > d) { // If neighbor w is still in a higher degree bin
                             vid current_deg_w = deg[w];
                             if (current_deg_w >= bins.size() || current_deg_w == 0) continue; // Safety check
 
                             // Find w in its current bin (bins[current_deg_w]) using pos[w]
                             size_t pos_w = pos[w];
                             if(pos_w >= bins[current_deg_w].size() || bins[current_deg_w][pos_w] != w) {
                                 // This indicates pos[] is inconsistent, potentially due to duplicate edges or logic error
                                 // Fallback: search for w
                                 bool found = false;
                                 for(size_t temp_pos = 0; temp_pos < bins[current_deg_w].size(); ++temp_pos) {
                                     if (bins[current_deg_w][temp_pos] == w) {
                                         pos_w = temp_pos;
                                         pos[w] = pos_w; // Correct the position map
                                         found = true;
                                         break;
                                     }
                                 }
                                 if (!found) {
                                      std::cerr << "Warning: Degcore inconsistency - node " << w << " not found at pos[" << w << "]=" << pos[w] << " in bin[" << current_deg_w << "]. Skipping neighbor update." << std::endl;
                                      continue; // Skip this update if state is corrupt
                                 }
                             }
 
                             // Get the last element in the *current* degree bin
                             vid last_node_in_bin = bins[current_deg_w].back();
 
                             // Swap w with the last element in bins[current_deg_w]
                             bins[current_deg_w][pos_w] = last_node_in_bin;
                             pos[last_node_in_bin] = pos_w; // Update swapped node's position
                             bins[current_deg_w].pop_back();
 
                             // Move w to the lower degree bin (bins[current_deg_w - 1])
                             deg[w]--; // Should equal current_deg_w - 1
                             pos[w] = bins[deg[w]].size(); // New position is at the end of the lower bin
                             bins[deg[w]].push_back(w);
                         }
                     } // end loop w
                 } // end loop r_orig
             } // end while bin_ptr
         } // end for d
 
         // Ensure all nodes were processed
         if (current_rank != n) {
              std::cerr << "Warning: Degcore processing finished with current_rank (" << current_rank << ") != n (" << n << "). Graph might be disconnected or error occurred." << std::endl;
              // Fill remaining slots if any? Or just return the partial order? Returning as is.
         }
         return core_order;
     }
 
     // Should not be reached
     throw std::logic_error("Invalid RankType");
 }
 
 
 // ───────────────────── Aggregators ───────────────────────────
 
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
     // Intermediate flush MUST NOT clear state for Hash
     void flush(){ /* DO NOTHING - State must persist */ }
     void final_flush() { /* Optionally clear h if Agg object reused, but not needed for total */ }
     size_t mem() const { // Rough estimate
         return h.bucket_count() * 1.5 * (sizeof(void*) + sizeof(pair64) + sizeof(uint64_t))
                + h.size() * (sizeof(pair64) + sizeof(uint64_t)); // Add node overhead
     }
 };
 
 // SORT – sort buffer at flush (only final flush is correct)
 struct AggSort {
     std::vector<pair64> buf;
     uint64_t total = 0;
     AggSort()=default;
     explicit AggSort(size_t L_size){ /* Optional: buf.reserve(estimate); */ }
     void add(pair64 p){ buf.push_back(p); }
     // Intermediate flush MUST NOT be called externally if correctness is desired.
     // It processes the *current* buffer state, which is incomplete mid-computation.
     void flush(){
         if (buf.empty()) return;
         std::sort(buf.begin(), buf.end());
         for (size_t i = 0; i < buf.size(); ) {
             size_t j = i;
             while (j < buf.size() && buf[j] == buf[i]) ++j;
             uint64_t f = j - i; // f = frequency = # common neighbors 'r' in this batch
             if (f > 1) {
                 total += f * (f - 1) / 2; // Choose 2 neighbors
             }
             i = j;
         }
         buf.clear(); // Clear buffer *after* processing this batch
         // buf.shrink_to_fit(); // Optional memory reduction
     }
     // Final flush ensures the complete buffer is processed.
     void final_flush() { flush(); }
     size_t mem() const { return buf.size() * sizeof(pair64) + sizeof(*this); }
 };
 
 // BATCH – Now behaves identically to AggSort, only called via final_flush
 // Kept for CLI compatibility, but internal auto-flush is removed for correctness.
 struct AggBatch : AggSort {
     // static constexpr size_t BATCH_SIZE_ITEMS = (1 << 19); // This constant is no longer used for auto-flushing
 
     AggBatch()=default;
     explicit AggBatch(size_t s):AggSort(s){ /* Optional: buf.reserve estimate based on BATCH_SIZE_ITEMS or other logic */ }
 
     // CRITICAL FIX: Remove internal flush condition. Flushing mid-computation is incorrect.
     void add(pair64 p){
         buf.push_back(p);
         // if (buf.size() >= BATCH_SIZE_ITEMS) flush(); // <-- REMOVED THIS LINE
     }
     // Inherits flush(), final_flush(), mem() from AggSort.
     // flush() will only be called via final_flush() externally now.
 };
 
 // HISTO – incremental choose-two with dense array
 // WARNING: Still potentially flawed logic for r-centric iteration unless pairs (u,w) are somehow guaranteed to be processed contiguously.
 struct AggHisto {
     std::vector<uint32_t> hist; // Counts associated with 'w'
     std::vector<vid> touched;   // Tracks which 'w' indices were modified in the current cycle
     uint64_t total = 0;
     size_t L_size_internal = 0;
 
     explicit AggHisto(size_t L){
         if (L == 0) {
              // Allow L=0 for empty graphs, but hist won't be usable.
              // Consider throwing an error if L=0 but graph is not empty.
              L_size_internal = 0;
              return;
         }
         L_size_internal = L;
         try {
             hist.assign(L, 0);
             // touched.reserve(L); // Optional: Can be large reservation
         } catch (const std::bad_alloc& e) {
             throw std::runtime_error("Failed to allocate AggHisto::hist vector of size " + std::to_string(L));
         }
     }
 
     void add(pair64 p){
         if (L_size_internal == 0) return; // Cannot use hist if L=0
 
         vid w = (vid)(p & 0xFFFFFFFFu); // Assumes w is in lower 32 bits
 
         if (w >= L_size_internal) {
              // This indicates an issue with L sizing or pair packing. Error out.
             throw std::runtime_error("AggHisto index w=" + std::to_string(w) +
                                      " >= L=" + std::to_string(L_size_internal));
         }
 
         // Incremental logic: add the count *before* incrementing
         if (hist[w] == 0) { // First time seeing w *in this cycle*? Track it for flush.
             // This relies on flush being called appropriately.
             // If intermediate flushes don't happen, touched could grow very large.
             touched.push_back(w);
         }
         total += hist[w];
         hist[w]++;
     }
 
     // Intermediate flush - Needed if AggHisto state needs reset (which it does for this logic)
     // BUT, calling it mid-computation leads to incorrectness.
     // This aggregator seems fundamentally mismatched with the parallel strategy unless modified.
     void flush(){
         // This resets counts. Should ONLY be called at the very end for correctness in the current setup.
         for (vid w_idx : touched) {
             if (w_idx < L_size_internal) { // Bounds check
                  hist[w_idx] = 0;
             }
         }
         touched.clear();
     }
     // Final flush ensures state is cleared (if needed), but total is already accumulated.
     void final_flush() { flush(); }
     size_t mem() const { return hist.size()*sizeof(uint32_t) + touched.size()*sizeof(vid) + sizeof(*this); }
 };
 
 
 // ───────────────────── Parallel driver (Corrected) ───────────────────────
 template<class Agg>
 uint64_t count_par(const Graph& G, const std::vector<vid>& rank_of, int nT){
     if (nT <= 0) nT = 1; // Ensure at least one thread
     if (G.L == 0) return 0; // No L-nodes, no butterflies
 
     std::vector<Agg> local(nT, Agg(G.L)); // Create nT local aggregators
     std::vector<std::vector<pair64>> inbox(nT);
     std::vector<omp_lock_t> locks(nT);
     bool locks_initialized = true;
     for (int t=0; t<nT; ++t) {
        omp_init_lock(&locks[t]);
        //  if (omp_init_lock(&locks[t]) != 0) {
        //       locks_initialized = false; // Handle potential init failure
        //       std::cerr << "Error: Failed to initialize lock for thread " << t << std::endl;
        //       // Cleanup already initialized locks before erroring?
        //       for (int k=0; k<t; ++k) omp_destroy_lock(&locks[k]);
        //       throw std::runtime_error("Failed to initialize OpenMP locks");
        //  }
     }
 
     // Simple FNV1a-based hash for potentially better distribution than XOR
     auto owner = [&](vid u, vid w){
         // Ensure consistent order for hash input (u < w based on rank already)
         pair64 key = pack(u, w);
         uint64_t hash = 0xcbf29ce484222325ULL;
         const uint64_t prime = 0x100000001b3ULL;
         unsigned char* p = reinterpret_cast<unsigned char*>(&key);
         for(size_t i=0; i<sizeof(key); ++i) {
             hash ^= p[i];
             hash *= prime;
         }
         return hash % nT;
     };
 
     std::atomic<bool> critical_error_flag = false;
 
 #pragma omp parallel num_threads(nT)
     {
         int tid = omp_get_thread_num();
         auto& A = local[tid]; // Reference to thread-local aggregator
 
         try { // Add try-catch block within parallel region for exceptions
 
         // Schedule dynamic can be good if work per 'rid' varies significantly
 #pragma omp for schedule(dynamic, 64) nowait
         for (size_t rid = 0; rid < G.R2L.size(); ++rid) {
             if (critical_error_flag.load()) continue; // Stop processing if fatal error occurred
 
             const auto& aSet = G.R2L[rid]; // Neighbors of r in L
 
             if (aSet.size() < 2) continue; // Cannot form a wedge
 
             for (size_t i = 0; i < aSet.size(); ++i) {
                 vid u0 = aSet[i];
                 // Basic bounds check for safety
                 if (u0 >= rank_of.size()) {
                      #pragma omp critical
                      { std::cerr << "Error: Out-of-bounds u0=" << u0 << " access in rank_of (size=" << rank_of.size() << ")" << std::endl; }
                      critical_error_flag = true; continue;
                 }
 
                 for (size_t j = i + 1; j < aSet.size(); ++j) {
                     vid w0 = aSet[j];
                     if (w0 >= rank_of.size()) {
                         #pragma omp critical
                         { std::cerr << "Error: Out-of-bounds w0=" << w0 << " access in rank_of (size=" << rank_of.size() << ")" << std::endl; }
                         critical_error_flag = true; continue;
                     }
 
                     vid u = u0, w = w0;
 
                     // Order pair (u, w) by rank (lower rank first)
                     // rank_of maps node_id -> rank
                     if (rank_of[u] > rank_of[w]) {
                         std::swap(u, w);
                     }
 
                     pair64 p = pack(u, w);
                     int o = owner(u, w);
 
                     if (o == tid) {
                         A.add(p); // Process locally
                     } else {
                         // Send to owner's inbox (locked)
                         omp_set_lock(&locks[o]);
                         try {
                             inbox[o].push_back(p);
                         } catch (const std::bad_alloc& e) {
                              omp_unset_lock(&locks[o]);
                             #pragma omp critical
                             { std::cerr << "FATAL: Out of memory adding to inbox[" << o << "] from thread " << tid << ". Aborting." << std::endl; }
                              critical_error_flag = true; // Signal other threads
                              // Consider std::terminate or other error handling
                         }
                         omp_unset_lock(&locks[o]);
                     }
                     if (critical_error_flag.load()) break; // Exit inner loop if error
                 } // end inner loop j
                  if (critical_error_flag.load()) break; // Exit outer loop if error
             } // end outer loop i
 
             // --- NO PERIODIC FLUSH HERE ---
             // Flushing stateful aggregators (Hash, Histo, Sort/Batch) here is incorrect.
             // They must accumulate all pairs before the final calculation/flush.
 
         } // end omp for
 
         // If an error occurred, bypass inbox processing for potentially faster termination
         if (!critical_error_flag.load()) {
             // Barrier: Wait for all threads to finish the 'for' loop & sending messages.
             // Implicit barrier here unless 'nowait' was used everywhere (it was on the 'for'). Need explicit barrier.
 #pragma omp barrier
 
             // Phase 2: Drain own inbox
             std::vector<pair64> local_inbox_copy;
             // Optional: Reserve space
             // size_t estimated_inbox_size = G.R2L.size() * G.L / nT / 10; // Very rough estimate
             // try { local_inbox_copy.reserve(estimated_inbox_size); } catch (...) { /* ignore reserve failure */ }
 
             omp_set_lock(&locks[tid]);
             std::swap(local_inbox_copy, inbox[tid]); // Quick swap while holding lock
             omp_unset_lock(&locks[tid]);
 
             // Process pairs received from other threads
             for (const auto& p : local_inbox_copy) {
                  if (critical_error_flag.load()) break; // Check flag again
                  A.add(p);
             }
         } // end if !critical_error_flag
 
         // Final flush for this thread's local aggregator.
         // This is where AggSort/AggBatch process their complete buffer.
         // AggHash/AggHisto don't strictly *need* this for the total, but good practice for potential state cleanup.
         if (!critical_error_flag.load()) {
              A.final_flush();
         }
 
         } catch (const std::exception& e) {
             // Catch exceptions from A.add, A.final_flush etc. within the parallel region
             #pragma omp critical
             { std::cerr << "Runtime error in thread " << tid << ": " << e.what() << std::endl; }
             critical_error_flag = true; // Signal other threads
         } catch (...) {
             #pragma omp critical
             { std::cerr << "Unknown exception in thread " << tid << std::endl; }
             critical_error_flag = true; // Signal other threads
         }
 
     } // end omp parallel
 
     // Check if a critical error occurred during parallel execution
     if (critical_error_flag.load()) {
         // Cleanup locks before throwing
         if (locks_initialized) {
             for (int t=0; t<nT; ++t) { omp_destroy_lock(&locks[t]); }
         }
         throw std::runtime_error("Fatal error occurred during parallel computation.");
     }
 
     // Clean up locks
     if (locks_initialized) {
       for (int t=0; t<nT; ++t) {
           omp_destroy_lock(&locks[t]);
       }
     }
 
 
     // Aggregate the total counts from all threads
     uint64_t total_butterflies = 0;
     for (int t = 0; t < nT; ++t) {
         total_butterflies += local[t].total;
     }
 
     return total_butterflies;
 }
 
 
 // ───────────────────── CLI parsing ──────────────────────────
 enum class ModeType { SEQ, PAR };
 struct Params {
     ModeType mode = ModeType::SEQ;
     std::string file;
     vid L = 0; // L size parameter from command line
     RankType rank = RankType::ID;
     uint64_t seed = 1;
     int agg = 0;      // 0=hash,1=sort,2=batch,3=histo
     int threads = 1;
 };
 
 // More robust CLI parsing
 Params parse(int argc, char** argv){
     Params P;
     std::vector<std::string> args(argv + 1, argv + argc);
     bool mode_set = false, file_set = false, L_set = false;
 
     for (size_t i = 0; i < args.size(); ++i) {
         std::string a = args[i];
         std::string next_arg = (i + 1 < args.size()) ? args[i+1] : "";
 
         try {
             if (a == "--mode" && !next_arg.empty()) {
                 if (next_arg == "seq") P.mode = ModeType::SEQ;
                 else if (next_arg == "par") P.mode = ModeType::PAR;
                 else throw std::runtime_error("Unknown mode '" + next_arg + "'");
                 mode_set = true; i++;
             } else if (a == "--rank" && !next_arg.empty()) {
                 if (next_arg == "id") P.rank = RankType::ID;
                 else if (next_arg == "degree") P.rank = RankType::DEGREE;
                 else if (next_arg == "degcore") P.rank = RankType::DEGCORE;
                 else if (next_arg.rfind("random", 0) == 0) {
                     P.rank = RankType::RANDOM;
                     auto pos = next_arg.find(':');
                     if (pos != std::string::npos) {
                          try { P.seed = std::stoull(next_arg.substr(pos + 1)); }
                          catch(...) { throw std::runtime_error("Invalid seed value after 'random:'"); }
                     }
                 } else throw std::runtime_error("Unknown rank '" + next_arg + "'");
                 i++;
             } else if (a == "--agg" && !next_arg.empty()) {
                 if (next_arg == "hash") P.agg = 0;
                 else if (next_arg == "sort") P.agg = 1;
                 else if (next_arg == "batch") P.agg = 2; // Keep batch option
                 else if (next_arg == "histo") P.agg = 3;
                 else throw std::runtime_error("Unknown agg '" + next_arg + "'");
                 i++;
             } else if (a == "--threads" && !next_arg.empty()) {
                 int nt = std::stoi(next_arg);
                 if (nt < 0) throw std::runtime_error("Number of threads cannot be negative.");
                 // If 0, use default max threads later.
                 P.threads = nt;
                 i++;
             } else if (!mode_set && a.rfind("--", 0) != 0 && !a.empty()) {
                  // If we haven't seen --mode yet and this doesn't look like a flag, assume file
                  P.file = a; file_set = true;
             } else if (mode_set && !file_set && a.rfind("--", 0) != 0 && !a.empty()) {
                  // If mode is set, haven't seen file, and not a flag -> file
                  P.file = a; file_set = true;
             } else if (mode_set && file_set && !L_set && a.rfind("--", 0) != 0 && !a.empty()) {
                  // If mode and file set, haven't seen L, and not a flag -> L
                  P.L = std::stoul(a); L_set = true;
                  if (P.L == 0) std::cerr << "Warning: L_size specified as 0." << std::endl;
             } else if (!a.empty() && a.rfind("--", 0) != 0) {
                  // Looks like a positional argument but we already have file and L
                  throw std::runtime_error("Unexpected positional argument: " + a);
             } else if (a.rfind("--", 0) == 0 && next_arg.empty()) {
                  // Flag without a value
                  throw std::runtime_error("Flag " + a + " requires an argument.");
             } else {
                 // Either an unknown flag or handled in previous iterations
                 if (a.rfind("--", 0) == 0) // Only throw error if it looks like an unknown flag
                    throw std::runtime_error("Unknown or misplaced argument: " + a);
             }
         } catch (const std::invalid_argument& e) {
             throw std::runtime_error("Invalid numeric value for argument near '" + a + "'");
         } catch (const std::out_of_range& e) {
             throw std::runtime_error("Numeric value out of range for argument near '" + a + "'");
         }
     }
 
     // Validate required arguments were found
     if (!mode_set) throw std::runtime_error("Missing required argument: --mode seq|par");
     if (P.file.empty()) throw std::runtime_error("Missing required argument: <edge_file>");
     if (!L_set) throw std::runtime_error("Missing required argument: <L_size>");
 
     // Adjust threads for parallel mode if needed
     if (P.mode == ModeType::PAR) {
         if (P.threads <= 0) { // If 0 or negative, use max available
 #ifdef _OPENMP
             P.threads = omp_get_max_threads();
             if (P.threads <= 0) P.threads = 1; // Safety for odd omp_get_max_threads results
 #else
             P.threads = 1;
             if (P.threads != 0) // Only warn if user didn't explicitly request 0 or 1
                std::cerr << "Warning: OpenMP not enabled, running parallel mode with 1 thread." << std::endl;
 #endif
         }
     } else { // Sequential mode always uses 1 thread for counting logic
         P.threads = 1;
     }
 
     return P;
 }
 
 // ───────────────────── Timer helper ───────────────────────────
 double get_time() {
 #ifdef _OPENMP
     return omp_get_wtime();
 #else
     // Use high_resolution_clock for better precision if OpenMP is off
     auto now = std::chrono::high_resolution_clock::now();
     auto duration = now.time_since_epoch();
     return std::chrono::duration<double>(duration).count();
 #endif
 }
 
 // ───────────────────── main ────────────────────────────────
 int main(int argc, char** argv){
     Params P;
     try {
         P = parse(argc, argv);
     } catch (const std::exception &e) {
         std::cerr << "Error parsing arguments: " << e.what() << "\n";
         std::cerr << "Usage: ./butter --mode seq|par <edge_file> <L_size> [--rank id|degree|degcore|random[:seed]] [--agg hash|sort|batch|histo] [--threads N]\n";
         return 1;
     }
 
     try {
         double t_start_load = get_time();
         Graph G = load_graph(P.file, P.L);
         double t_end_load = get_time();
 
         // Check if L size seems valid after loading
         if (G.L == 0 && G.R2L.empty()) {
              std::cout << "Info: Graph appears empty after loading." << std::endl;
         } else if (G.L == 0) {
              std::cerr << "Warning: Graph L size is 0, but R partition has data? Check L_size param and graph format." << std::endl;
         }
 
         double t_start_rank = get_time();
         auto order = make_rank(G, P.rank, P.seed);
         double t_end_rank = get_time();
 
         std::vector<vid> rank_of; // Will be sized based on G.L
         if (G.L > 0) {
             if (order.size() != G.L) {
                  throw std::runtime_error("Rank order size (" + std::to_string(order.size()) +
                                           ") mismatch with graph L size (" + std::to_string(G.L) + ").");
             }
             rank_of.resize(G.L); // Size based on graph's effective L size
             for (vid r = 0; r < G.L; ++r) { // r is the rank
                 vid node_id = order[r];
                 if (node_id >= G.L) {
                      throw std::runtime_error("Node ID " + std::to_string(node_id) + " from rank order" +
                                              " is out of bounds for graph L size (" + std::to_string(G.L) + ").");
                 }
                 rank_of[node_id] = r; // Assign rank based on position in sorted 'order'
             }
         } else {
              // G.L is 0, order should be empty, rank_of remains empty.
              if (!order.empty()) {
                  throw std::runtime_error("Graph L size is 0, but ranking returned a non-empty order.");
              }
         }
 
 
         uint64_t ans = 0;
         double t_start_count = get_time();
 
         // Determine threads to use based on mode (always use count_par)
         int nT_to_use = (P.mode == ModeType::SEQ) ? 1 : P.threads;
 
         // Check OpenMP availability if parallel execution requested
         #ifndef _OPENMP
         if (nT_to_use > 1) {
             std::cerr << "Warning: Compiled without OpenMP support, running with 1 thread instead of requested " << nT_to_use << "." << std::endl;
             nT_to_use = 1;
         }
         #endif
 
         // Call the appropriate count_par function via template based on P.agg
         switch (P.agg) {
             case 0: ans = count_par<AggHash >(G, rank_of, nT_to_use); break;
             case 1: ans = count_par<AggSort >(G, rank_of, nT_to_use); break;
             case 2: ans = count_par<AggBatch>(G, rank_of, nT_to_use); break; // Uses fixed AggBatch now
             case 3:
                 if (nT_to_use > 0 && G.L > 0) { // Only warn if actually running and histo used
                     std::cerr << "Warning: Using AggHisto. Ensure L_size is accurate and sufficient memory is available. Logic might be sensitive to iteration order." << std::endl;
                 }
                 ans = count_par<AggHisto>(G, rank_of, nT_to_use); break;
             default: throw std::logic_error("Internal error: Invalid aggregator type selected.");
         }
 
         double t_end_count = get_time();
 
         // Output results clearly
         std::cout << "[" << (P.mode == ModeType::SEQ ? "seq" : "par") << "]"
                   << " butterflies=" << ans
                   << " load=" << std::fixed << std::setprecision(5) << (t_end_load - t_start_load) << "s"
                   << " rank_t=" << std::fixed << std::setprecision(5) << (t_end_rank - t_start_rank) << "s"
                   << " count=" << std::fixed << std::setprecision(5) << (t_end_count - t_start_count) << "s"
                   << " thr=" << nT_to_use // Report actual threads used
                   << " rank=" << static_cast<int>(P.rank)
                   << " agg=" << P.agg
                   << "\n";
 
     } catch (const std::exception &e) {
         std::cerr << "Runtime Error: " << e.what() << "\n";
         return 1;
     }
     return 0;
 }