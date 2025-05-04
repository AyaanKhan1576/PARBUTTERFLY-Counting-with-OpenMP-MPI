// /****************************************************************************************
//  * butter.cpp  – unified sequential and parallel butterfly counter
//  *
//  * Build : g++ -std=c++17 -O3 butter.cpp -fopenmp -o butter
//  *
//  * Usage :
//  *   ./butter --mode seq|par <edge_file> <L_size>
//  *          --rank id|degree|degcore|random[:seed]
//  *          --agg  hash|sort|batch|histo
//  *          [--threads N]           (only for --mode par)
//  *****************************************************************************************/

//  #include <bits/stdc++.h>
//  #ifdef _OPENMP
//  #include <omp.h>
//  #endif
 
//  using vid    = uint32_t;
//  using pair64 = uint64_t;
//  static inline pair64 pack(vid a, vid b) { return (pair64)a << 32 | b; }
 
//  // ───────────────────── Graph loader ───────────────────────────
//  // We read the normalized bipartite subgraph (L-size = number of L-nodes).
//  // Then we build L2R (not used here) and R2L (hyperedges: each R is a hyperedge
//  // connecting its neighbors in L).
//  struct Graph {
//      vid L;
//      std::vector<std::vector<vid>> L2R, R2L;
//  };
//  Graph load_graph(const std::string& path, vid L) {
//      std::ifstream in(path);
//      if (!in) throw std::runtime_error("Cannot open " + path);
//      std::vector<std::pair<vid, vid>> edges;
//      size_t Rmax = 0;
//      for (std::string line; std::getline(in, line); ) {
//          if (line.empty() || line[0] == '%') continue;
//          vid u, v;
//          std::istringstream iss(line);
//          iss >> u >> v;
//          edges.emplace_back(u, v);
//          Rmax = std::max(Rmax, size_t(v - L));
//      }
//      Graph G;
//      G.L = L;
//      G.L2R.assign(L, {});
//      G.R2L.assign(Rmax + 1, {});
//      for (auto &e : edges) {
//          G.L2R[e.first].push_back(e.second);
//          G.R2L[e.second - L].push_back(e.first);
//      }
//      for (auto &nb : G.L2R) std::sort(nb.begin(), nb.end());
//      return G;
//  }
 
//  // ───────────────────── Ranking options ───────────────────────
//  enum class RankType { ID, DEGREE, DEGCORE, RANDOM };
//  std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1) {
//      vid n = G.L;
//      std::vector<vid> ord(n);
//      std::iota(ord.begin(), ord.end(), 0);
//      if (typ == RankType::ID) return ord;
//      if (typ == RankType::DEGREE) {
//          std::sort(ord.begin(), ord.end(),
//                    [&](vid a, vid b){ return G.L2R[a].size() < G.L2R[b].size(); });
//          return ord;
//      }
//      if (typ == RankType::RANDOM) {
//          std::mt19937_64 rng(seed);
//          std::shuffle(ord.begin(), ord.end(), rng);
//          return ord;
//      }
//      // one‐pass approximate core order (degcore)
//      std::vector<vid> deg(n);
//      for (vid u = 0; u < n; ++u) deg[u] = G.L2R[u].size();
//      vid maxd = *std::max_element(deg.begin(), deg.end());
//      std::vector<std::vector<vid>> bins(maxd + 1);
//      for (vid u = 0; u < n; ++u) bins[deg[u]].push_back(u);
//      vid idx = 0;
//      for (vid d = 0; d <= maxd; ++d) {
//          auto &bin = bins[d];
//          while (!bin.empty()) {
//              vid u = bin.back(); bin.pop_back();
//              ord[idx++] = u;
//              for (vid r : G.L2R[u]) {
//                  for (vid w : G.R2L[r - G.L]) {
//                      if (deg[w] > d) {
//                          --deg[w];
//                          bins[deg[w]].push_back(w);
//                      }
//                  }
//              }
//          }
//      }
//      return ord;
//  }
 
//  // ───────────────────── Aggregators ───────────────────────────
//  static constexpr size_t FLUSH_BYTES = 8 * 1024 * 1024;
 
//  // HASH – incremental choose-two
//  struct AggHash {
//      std::unordered_map<pair64,int> h;
//      uint64_t total = 0;
//      AggHash()=default; AggHash(size_t){}  // dummy ctor
//      void add(pair64 p) {
//          int &k = h[p];
//          total += k;
//          ++k;
//      }
//      void flush(){ h.clear(); }
//      size_t mem() const { return h.size() * sizeof(*h.begin()); }
//  };
 
//  // SORT – sort buffer at flush
//  struct AggSort {
//      std::vector<pair64> buf;
//      uint64_t total = 0;
//      AggSort()=default; AggSort(size_t){}
//      void add(pair64 p){ buf.push_back(p); }
//      void flush(){
//          std::sort(buf.begin(), buf.end());
//          for (size_t i = 0; i < buf.size(); ) {
//              size_t j = i;
//              while (j < buf.size() && buf[j] == buf[i]) ++j;
//              uint64_t f = j - i;
//              total += f * (f-1) / 2;
//              i = j;
//          }
//          buf.clear();
//      }
//      size_t mem() const { return buf.size() * sizeof(pair64); }
//  };
 
//  // BATCH – fixed chunk then sort
//  struct AggBatch : AggSort {
//      static constexpr size_t BATCH = 1 << 19;
//      AggBatch()=default; AggBatch(size_t s):AggSort(s){}
//      void add(pair64 p){
//          buf.push_back(p);
//          if (buf.size() == BATCH) flush();
//      }
//  };
 
//  // HISTO – incremental choose-two with dense array
//  struct AggHisto {
//      std::vector<int> hist;
//      std::vector<vid> touched;
//      uint64_t total = 0;
//      explicit AggHisto(size_t L){ hist.assign(L,0); }
//      void add(pair64 p){
//          vid w = (vid)(p & 0xFFFFFFFFu);
//          total += hist[w];
//          if (hist[w]++ == 0) touched.push_back(w);
//      }
//      void flush(){
//          for (vid w : touched) hist[w] = 0;
//          touched.clear();
//      }
//      size_t mem() const { return touched.size()*sizeof(vid); }
//  };
 
//  // ───────────────────── Sequential driver ─────────────────────
//  template<class Agg>
//  uint64_t count_seq(const Graph& G, const std::vector<vid>& rank){
//      Agg A(G.L);
//      for (size_t rid = 0; rid < G.R2L.size(); ++rid) {
//          const auto &aSet = G.R2L[rid];
//          for (size_t i = 0; i < aSet.size(); ++i) {
//              vid u = aSet[i];
//              for (size_t j = i+1; j < aSet.size(); ++j) {
//                  vid w = aSet[j];
//                  if (rank[u] > rank[w]) std::swap(u,w);
//                  A.add(pack(u,w));
//              }
//          }
//          if (A.mem() > FLUSH_BYTES) A.flush();
//      }
//      A.flush();
//      return A.total;
//  }
 
//  // ───────────────────── Parallel driver ───────────────────────
//  // “Owner‐hashing” across threads
//  template<class Agg>
//  uint64_t count_par(const Graph& G, const std::vector<vid>& rank, int nT){
//      std::vector<Agg> local(nT, Agg(G.L));
//      std::vector<std::vector<pair64>> inbox(nT);
//      std::vector<omp_lock_t> locks(nT);
//      for (int t=0; t<nT; ++t) omp_init_lock(&locks[t]);
 
//      auto owner = [&](vid u, vid w){ return (u ^ w) % nT; };
 
//  #pragma omp parallel num_threads(nT)
//      {
//          int tid = omp_get_thread_num();
//          auto &A = local[tid];
//  #pragma omp for schedule(dynamic,64) nowait
//          for (size_t rid=0; rid<G.R2L.size(); ++rid){
//              const auto &aSet = G.R2L[rid];
//              for (size_t i=0; i<aSet.size(); ++i){
//                  vid u = aSet[i];
//                  for (size_t j=i+1; j<aSet.size(); ++j){
//                      vid w = aSet[j];
//                      if (rank[u] > rank[w]) std::swap(u,w);
//                      int o = owner(u,w);
//                      if (o==tid) {
//                          A.add(pack(u,w));
//                      } else {
//                          omp_set_lock(&locks[o]);
//                          inbox[o].push_back(pack(u,w));
//                          omp_unset_lock(&locks[o]);
//                      }
//                  }
//              }
//              if (A.mem() > FLUSH_BYTES) A.flush();
//          }
//  #pragma omp barrier
//          // drain my inbox
//          std::vector<pair64> tmp;
//          omp_set_lock(&locks[tid]);
//          tmp.swap(inbox[tid]);
//          omp_unset_lock(&locks[tid]);
//          for (auto &p : tmp) A.add(p);
//          A.flush();
//      }
//      uint64_t sum = 0;
//      for (auto &A : local) sum += A.total;
//      return sum;
//  }
 
//  // ───────────────────── CLI parsing ──────────────────────────
//  enum class ModeType { SEQ, PAR };
//  struct Params {
//      ModeType mode = ModeType::SEQ;
//      std::string file;
//      vid L = 0;
//      RankType rank = RankType::ID;
//      uint64_t seed = 1;
//      int agg = 0;      // 0=hash,1=sort,2=batch,3=histo
//      int threads = 1;
//  };
 
//  Params parse(int argc, char** argv){
//      if (argc < 5 || std::string(argv[1])!="--mode")
//          throw std::runtime_error("Usage: butter --mode seq|par <edge_file> <L_size> [--rank …] [--agg …] [--threads N]");
//      Params P;
//      std::string m = argv[2];
//      if (m=="seq") P.mode = ModeType::SEQ;
//      else if (m=="par") P.mode = ModeType::PAR;
//      else throw std::runtime_error("Unknown mode "+m);
 
//      P.file = argv[3];
//      P.L    = std::stoul(argv[4]);
 
//      for (int i=5; i<argc; ++i){
//          std::string a = argv[i];
//          if (a=="--rank"){
//              std::string v = argv[++i];
//              if (v=="id")        P.rank = RankType::ID;
//              else if (v=="degree")  P.rank = RankType::DEGREE;
//              else if (v=="degcore") P.rank = RankType::DEGCORE;
//              else if (v.rfind("random",0)==0){
//                  P.rank = RankType::RANDOM;
//                  auto pos = v.find(':');
//                  if (pos!=std::string::npos) P.seed = std::stoull(v.substr(pos+1));
//              } else throw std::runtime_error("Unknown rank "+v);
//          }
//          else if (a=="--agg"){
//              std::string v = argv[++i];
//              if (v=="hash")  P.agg=0;
//              else if (v=="sort")  P.agg=1;
//              else if (v=="batch") P.agg=2;
//              else if (v=="histo") P.agg=3;
//              else throw std::runtime_error("Unknown agg "+v);
//          }
//          else if (a=="--threads"){
//              P.threads = std::stoi(argv[++i]);
//          }
//          else throw std::runtime_error("Unknown flag "+a);
//      }
//      if (P.mode==ModeType::PAR && P.threads<1)
//  #ifdef _OPENMP
//          P.threads = omp_get_max_threads();
//  #else
//          P.threads = 1;
//  #endif
//      return P;
//  }
 
//  // ───────────────────── main ────────────────────────────────
//  int main(int argc, char** argv){
//      try {
//          auto P = parse(argc,argv);
 
//          double t0, t1, t2;
//  #ifdef _OPENMP
//          t0 = omp_get_wtime();
//  #else
//          t0 = clock()/double(CLOCKS_PER_SEC);
//  #endif
//          Graph G = load_graph(P.file, P.L);
//  #ifdef _OPENMP
//          t1 = omp_get_wtime();
//  #else
//          t1 = clock()/double(CLOCKS_PER_SEC);
//  #endif
 
//          auto order = make_rank(G, P.rank, P.seed);
//          std::vector<vid> rank_of(G.L);
//          for (vid i=0; i<G.L; ++i) rank_of[order[i]] = i;
 
//          uint64_t ans = 0;
//          if (P.mode == ModeType::SEQ) {
//              switch(P.agg){
//                case 0: ans = count_seq<AggHash >(G, rank_of); break;
//                case 1: ans = count_seq<AggSort >(G, rank_of); break;
//                case 2: ans = count_seq<AggBatch>(G, rank_of); break;
//                case 3: ans = count_seq<AggHisto>(G, rank_of); break;
//              }
//  #ifdef _OPENMP
//              t2 = omp_get_wtime();
//  #else
//              t2 = clock()/double(CLOCKS_PER_SEC);
//  #endif
//              std::cout << "[seq] butterflies="<<ans
//                        <<" load="<<(t1-t0)<<"s"
//                        <<" count="<<(t2-t1)<<"s"
//                        <<" rank="<<(int)P.rank
//                        <<" agg="<<P.agg
//                        <<"\n";
//          }
//          else {
//  #ifdef _OPENMP
//              t2 = omp_get_wtime(); // reuse t2 as end of load for par
//  #else
//              t2 = clock()/double(CLOCKS_PER_SEC);
//  #endif
//              switch(P.agg){
//                case 0: ans = count_par<AggHash >(G, rank_of, P.threads); break;
//                case 1: ans = count_par<AggSort >(G, rank_of, P.threads); break;
//                case 2: ans = count_par<AggBatch>(G, rank_of, P.threads); break;
//                case 3: ans = count_par<AggHisto>(G, rank_of, P.threads); break;
//              }
//  #ifdef _OPENMP
//              double t3 = omp_get_wtime();
//              std::cout << "[par] butterflies="<<ans
//                        <<" load="<<(t2-t0)<<"s"
//                        <<" count="<<(t3-t2)<<"s"
//                        <<" thr="<<P.threads
//                        <<" rank="<<(int)P.rank
//                        <<" agg="<<P.agg
//                        <<"\n";
//  #else
//              double t3 = clock()/double(CLOCKS_PER_SEC);
//              std::cout << "[par] butterflies="<<ans
//                        <<" load="<<(t2-t0)<<"s"
//                        <<" count="<<(t3-t2)<<"s"
//                        <<" thr="<<P.threads
//                        <<" rank="<<(int)P.rank
//                        <<" agg="<<P.agg
//                        <<"\n";
//  #endif
//          }
//      }
//      catch(const std::exception &e){
//          std::cerr<<"Error: "<<e.what()<<"\n";
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
 *   ./butter --mode seq|par <edge_file> <L_size>
 *          --rank id|degree|degcore|random[:seed]
 *          --agg  hash|sort|batch|histo
 *          [--threads N]           (only for --mode par)
 *****************************************************************************************/

 #include <bits/stdc++.h>
 #ifdef _OPENMP
 #include <omp.h>
 #endif
 
 using vid    = uint32_t;
 using pair64 = uint64_t;
 static inline pair64 pack(vid a, vid b) { return (pair64)a << 32 | b; }
 
 // ───────────────────── Graph loader ───────────────────────────
 struct Graph {
     vid L;
     std::vector<std::vector<vid>> L2R, R2L;
 };
 Graph load_graph(const std::string& path, vid L) {
     std::ifstream in(path);
     if (!in) throw std::runtime_error("Cannot open " + path);
     std::vector<std::pair<vid, vid>> edges;
     size_t Rmax = 0;
     for (std::string line; std::getline(in, line); ) {
         if (line.empty() || line[0] == '%') continue;
         vid u, v;
         std::istringstream iss(line);
         iss >> u >> v;
         edges.emplace_back(u, v);
         Rmax = std::max(Rmax, size_t(v - L));
     }
     Graph G;
     G.L = L;
     G.L2R.assign(L, {});
     G.R2L.assign(Rmax + 1, {});
     for (auto &e : edges) {
         G.L2R[e.first].push_back(e.second);
         G.R2L[e.second - L].push_back(e.first);
     }
     for (auto &nb : G.L2R) std::sort(nb.begin(), nb.end());
     return G;
 }
 
 // ───────────────────── Ranking options ───────────────────────
 enum class RankType { ID, DEGREE, DEGCORE, RANDOM };
 std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1) {
     vid n = G.L;
     std::vector<vid> ord(n);
     std::iota(ord.begin(), ord.end(), 0);
     if (typ == RankType::ID) return ord;
     if (typ == RankType::DEGREE) {
         std::sort(ord.begin(), ord.end(),
                   [&](vid a, vid b){ return G.L2R[a].size() < G.L2R[b].size(); });
         return ord;
     }
     if (typ == RankType::RANDOM) {
         std::mt19937_64 rng(seed);
         std::shuffle(ord.begin(), ord.end(), rng);
         return ord;
     }
     // one‐pass approximate core order (degcore)
     std::vector<vid> deg(n);
     for (vid u = 0; u < n; ++u) deg[u] = G.L2R[u].size();
     vid maxd = *std::max_element(deg.begin(), deg.end());
     std::vector<std::vector<vid>> bins(maxd + 1);
     for (vid u = 0; u < n; ++u) bins[deg[u]].push_back(u);
     vid idx = 0;
     for (vid d = 0; d <= maxd; ++d) {
         auto &bin = bins[d];
         while (!bin.empty()) {
             vid u = bin.back(); bin.pop_back();
             ord[idx++] = u;
             for (vid r : G.L2R[u]) {
                 for (vid w : G.R2L[r - G.L]) {
                     if (deg[w] > d) {
                         --deg[w];
                         bins[deg[w]].push_back(w);
                     }
                 }
             }
         }
     }
     return ord;
 }
 
 // ───────────────────── Aggregators ───────────────────────────
 static constexpr size_t FLUSH_BYTES = 8 * 1024 * 1024;
 
 // HASH – incremental choose-two
 struct AggHash {
     std::unordered_map<pair64,int> h;
     uint64_t total = 0;
     AggHash()=default; AggHash(size_t){}  // dummy ctor
     void add(pair64 p) {
         int &k = h[p];
         total += k;
         ++k;
     }
     void flush(){ h.clear(); }
     size_t mem() const { return h.size() * sizeof(*h.begin()); }
 };
 
 // SORT – sort buffer at flush
 struct AggSort {
     std::vector<pair64> buf;
     uint64_t total = 0;
     AggSort()=default; AggSort(size_t){}
     void add(pair64 p){ buf.push_back(p); }
     void flush(){
         std::sort(buf.begin(), buf.end());
         for (size_t i = 0; i < buf.size(); ) {
             size_t j = i;
             while (j < buf.size() && buf[j] == buf[i]) ++j;
             uint64_t f = j - i;
             total += f * (f-1) / 2;
             i = j;
         }
         buf.clear();
     }
     size_t mem() const { return buf.size() * sizeof(pair64); }
 };
 
 // BATCH – fixed chunk then sort
 struct AggBatch : AggSort {
     static constexpr size_t BATCH = 1 << 19;
     AggBatch()=default; AggBatch(size_t s):AggSort(s){}
     void add(pair64 p){
         buf.push_back(p);
         if (buf.size() == BATCH) flush();
     }
 };
 
 // HISTO – incremental choose-two with dense array
 struct AggHisto {
     std::vector<int> hist;
     std::vector<vid> touched;
     uint64_t total = 0;
     explicit AggHisto(size_t L){ hist.assign(L,0); }
     void add(pair64 p){
         vid w = (vid)(p & 0xFFFFFFFFu);
         total += hist[w];
         if (hist[w]++ == 0) touched.push_back(w);
     }
     void flush(){
         for (vid w : touched) hist[w] = 0;
         touched.clear();
     }
     size_t mem() const { return touched.size()*sizeof(vid); }
 };
 
 // ───────────────────── Sequential driver ─────────────────────
 template<class Agg>
 uint64_t count_seq(const Graph& G, const std::vector<vid>& rank){
     Agg A(G.L);
     for (size_t rid = 0; rid < G.R2L.size(); ++rid) {
         const auto &aSet = G.R2L[rid];
         for (size_t i = 0; i < aSet.size(); ++i) {
             vid u0 = aSet[i];
             for (size_t j = i+1; j < aSet.size(); ++j) {
                 vid w0 = aSet[j];
                 vid u = u0, w = w0;              // copy before swap
                 if (rank[u] > rank[w]) std::swap(u,w);
                 A.add(pack(u,w));
             }
         }
         if (A.mem() > FLUSH_BYTES) A.flush();
     }
     A.flush();
     return A.total;
 }
 
 // ───────────────────── Parallel driver ───────────────────────
 template<class Agg>
 uint64_t count_par(const Graph& G, const std::vector<vid>& rank, int nT){
     std::vector<Agg> local(nT, Agg(G.L));
     std::vector<std::vector<pair64>> inbox(nT);
     std::vector<omp_lock_t> locks(nT);
     for (int t=0; t<nT; ++t) omp_init_lock(&locks[t]);
     auto owner = [&](vid u, vid w){ return (u ^ w) % nT; };
 
 #pragma omp parallel num_threads(nT)
     {
         int tid = omp_get_thread_num();
         auto &A = local[tid];
 #pragma omp for schedule(dynamic,64) nowait
         for (size_t rid=0; rid<G.R2L.size(); ++rid){
             const auto &aSet = G.R2L[rid];
             for (size_t i=0; i<aSet.size(); ++i){
                 vid u0 = aSet[i];
                 for (size_t j=i+1; j<aSet.size(); ++j){
                     vid w0 = aSet[j];
                     vid u = u0, w = w0;          // copy before swap
                     if (rank[u] > rank[w]) std::swap(u,w);
                     int o = owner(u,w);
                     if (o==tid) {
                         A.add(pack(u,w));
                     } else {
                         omp_set_lock(&locks[o]);
                         inbox[o].push_back(pack(u,w));
                         omp_unset_lock(&locks[o]);
                     }
                 }
             }
             if (A.mem() > FLUSH_BYTES) A.flush();
         }
 #pragma omp barrier
         // drain my inbox
         std::vector<pair64> tmp;
         omp_set_lock(&locks[tid]);
         tmp.swap(inbox[tid]);
         omp_unset_lock(&locks[tid]);
         for (auto &p : tmp) A.add(p);
         A.flush();
     }
     uint64_t sum = 0;
     for (auto &A : local) sum += A.total;
     return sum;
 }
 
 // ───────────────────── CLI parsing ──────────────────────────
 enum class ModeType { SEQ, PAR };
 struct Params {
     ModeType mode = ModeType::SEQ;
     std::string file;
     vid L = 0;
     RankType rank = RankType::ID;
     uint64_t seed = 1;
     int agg = 0;      // 0=hash,1=sort,2=batch,3=histo
     int threads = 1;
 };
 
 Params parse(int argc, char** argv){
     if (argc < 5 || std::string(argv[1])!="--mode")
         throw std::runtime_error("Usage: butter --mode seq|par <edge_file> <L_size> [--rank …] [--agg …] [--threads N]");
     Params P;
     std::string m = argv[2];
     if (m=="seq") P.mode = ModeType::SEQ;
     else if (m=="par") P.mode = ModeType::PAR;
     else throw std::runtime_error("Unknown mode "+m);
 
     P.file = argv[3];
     P.L    = std::stoul(argv[4]);
 
     for (int i=5; i<argc; ++i){
         std::string a = argv[i];
         if (a=="--rank"){
             std::string v = argv[++i];
             if (v=="id")        P.rank = RankType::ID;
             else if (v=="degree")  P.rank = RankType::DEGREE;
             else if (v=="degcore") P.rank = RankType::DEGCORE;
             else if (v.rfind("random",0)==0){
                 P.rank = RankType::RANDOM;
                 auto pos = v.find(':');
                 if (pos!=std::string::npos) P.seed = std::stoull(v.substr(pos+1));
             } else throw std::runtime_error("Unknown rank "+v);
         }
         else if (a=="--agg"){
             std::string v = argv[++i];
             if (v=="hash")  P.agg=0;
             else if (v=="sort")  P.agg=1;
             else if (v=="batch") P.agg=2;
             else if (v=="histo") P.agg=3;
             else throw std::runtime_error("Unknown agg "+v);
         }
         else if (a=="--threads"){
             P.threads = std::stoi(argv[++i]);
         }
         else throw std::runtime_error("Unknown flag "+a);
     }
     if (P.mode==ModeType::PAR && P.threads<1)
 #ifdef _OPENMP
         P.threads = omp_get_max_threads();
 #else
         P.threads = 1;
 #endif
     return P;
 }
 
 // ───────────────────── main ────────────────────────────────
 int main(int argc, char** argv){
     try {
         auto P = parse(argc,argv);
 
 #ifdef _OPENMP
         double t0 = omp_get_wtime();
 #else
         double t0 = clock()/double(CLOCKS_PER_SEC);
 #endif
         Graph G = load_graph(P.file, P.L);
 #ifdef _OPENMP
         double t1 = omp_get_wtime();
 #else
         double t1 = clock()/double(CLOCKS_PER_SEC);
 #endif
 
         auto order = make_rank(G, P.rank, P.seed);
         std::vector<vid> rank_of(G.L);
         for (vid i=0; i<G.L; ++i) rank_of[order[i]] = i;
 
         uint64_t ans = 0;
         if (P.mode == ModeType::SEQ) {
             switch(P.agg){
               case 0: ans = count_seq<AggHash >(G, rank_of); break;
               case 1: ans = count_seq<AggSort >(G, rank_of); break;
               case 2: ans = count_seq<AggBatch>(G, rank_of); break;
               case 3: ans = count_seq<AggHisto>(G, rank_of); break;
             }
 #ifdef _OPENMP
             double t2 = omp_get_wtime();
 #else
             double t2 = clock()/double(CLOCKS_PER_SEC);
 #endif
             std::cout << "[seq] butterflies="<<ans
                       <<" load="<<(t1-t0)<<"s"
                       <<" count="<<(t2-t1)<<"s"
                       <<" rank="<<(int)P.rank
                       <<" agg="<<P.agg
                       <<"\n";
         }
         else {
 #ifdef _OPENMP
             double load_end = omp_get_wtime();
 #else
             double load_end = clock()/double(CLOCKS_PER_SEC);
 #endif
             switch(P.agg){
               case 0: ans = count_par<AggHash >(G, rank_of, P.threads); break;
               case 1: ans = count_par<AggSort >(G, rank_of, P.threads); break;
               case 2: ans = count_par<AggBatch>(G, rank_of, P.threads); break;
               case 3: ans = count_par<AggHisto>(G, rank_of, P.threads); break;
             }
 #ifdef _OPENMP
             double t3 = omp_get_wtime();
             std::cout << "[par] butterflies="<<ans
                       <<" load="<<(load_end-t0)<<"s"
                       <<" count="<<(t3-load_end)<<"s"
                       <<" thr="<<P.threads
                       <<" rank="<<(int)P.rank
                       <<" agg="<<P.agg
                       <<"\n";
 #else
             double t3 = clock()/double(CLOCKS_PER_SEC);
             std::cout << "[par] butterflies="<<ans
                       <<" load="<<(load_end-t0)<<"s"
                       <<" count="<<(t3-load_end)<<"s"
                       <<" thr="<<P.threads
                       <<" rank="<<(int)P.rank
                       <<" agg="<<P.agg
                       <<"\n";
 #endif
         }
     }
     catch(const std::exception &e){
         std::cerr<<"Error: "<<e.what()<<"\n";
         return 1;
     }
     return 0;
 }
 