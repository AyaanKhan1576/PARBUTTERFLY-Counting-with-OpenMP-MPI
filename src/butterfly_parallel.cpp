// /****************************************************************************************
//  * butterfly_parallel.cpp  –  deterministic ParButterfly counter with owner hashing
//  *
//  * Build : g++ -std=c++17 -O3 -fopenmp butterfly_parallel.cpp -o butter_par
//  * Usage : ./butter_par edge_file L_size
//  *         --rank  id | degree | degcore | random[:seed]
//  *         --agg   hash | sort | batch | histo
//  *         [--threads N]
//  *
//  * Output : [par] butterflies=<cnt> load=<s> count=<s> thr=<n> rank=<id> agg=<id>
//  *****************************************************************************************/
// #include <bits/stdc++.h>
// #include <omp.h>

// /* ---------------------------------------------------------------------------- types */
// using vid    = uint32_t;
// using pair64 = uint64_t;
// static inline pair64 pack(vid a, vid b) { return (uint64_t)a << 32 | b; }

// /* ---------------------------------------------------------------------------- graph */
// struct Graph {
//     vid L;
//     std::vector<std::vector<vid>> L2R, R2L;
// };

// Graph load_graph(const std::string& path, vid L)
// {
//     std::ifstream in(path);
//     if (!in) throw std::runtime_error("Cannot open " + path);

//     std::vector<std::pair<vid,vid>> E;
//     size_t Rmax = 0;
//     for (std::string line; std::getline(in, line);) {
//         if (line.empty() || line[0]=='%') continue;
//         vid u,v; std::istringstream(line) >> u >> v;
//         E.emplace_back(u,v);
//         Rmax = std::max<size_t>(Rmax, v - L);
//     }
//     Graph G; G.L = L;
//     G.L2R.assign(L, {});
//     G.R2L.assign(Rmax+1, {});
//     for (auto [u,v] : E) { G.L2R[u].push_back(v); G.R2L[v-L].push_back(u); }
//     for (auto& nb : G.L2R) std::sort(nb.begin(), nb.end());
//     return G;
// }

// /* -------------------------------------------------------------------- ranking */
// enum class RankType { ID, DEGREE, DEGCORE, RANDOM };

// std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1)
// {
//     vid n = G.L;
//     std::vector<vid> order(n); std::iota(order.begin(), order.end(), 0);

//     if (typ == RankType::ID) return order;

//     if (typ == RankType::DEGREE) {
//         std::sort(order.begin(), order.end(),
//                   [&](vid a, vid b){ return G.L2R[a].size() < G.L2R[b].size(); });
//         return order;
//     }
//     if (typ == RankType::RANDOM) {
//         std::mt19937_64 rng(seed);
//         std::shuffle(order.begin(), order.end(), rng);
//         return order;
//     }

//     /* DEGCORE – one-pass approximate core ordering */
//     std::vector<vid> deg(n); for (vid u=0; u<n; ++u) deg[u]=G.L2R[u].size();
//     vid maxd=*std::max_element(deg.begin(),deg.end());
//     std::vector<std::vector<vid>> bin(maxd+1);
//     for (vid u=0; u<n; ++u) bin[deg[u]].push_back(u);

//     vid out=0;
//     for (vid d=0; d<=maxd; ++d)
//         while (!bin[d].empty()) {
//             vid u = bin[d].back(); bin[d].pop_back();
//             order[out++] = u;
//             for (vid r: G.L2R[u])
//                 for (vid w: G.R2L[r-G.L]) if (deg[w] > d) {
//                     --deg[w]; bin[deg[w]].push_back(w);
//                 }
//         }
//     return order;
// }

// /* ---------------------------------------------------------------- aggregators */
// static constexpr size_t FLUSH_BYTES = 8*1024*1024;

// /* --- hash */
// struct AggHash {
//     std::unordered_map<pair64,int> h; uint64_t local = 0;
//     AggHash() = default; AggHash(size_t){}          /* dummy ctor */
//     void  add(pair64 p){ ++h[p]; }
//     void  flush(){ for(auto& kv:h) if(kv.second>1) local += (uint64_t)kv.second*(kv.second-1)/2; h.clear(); }
//     size_t mem() const { return h.size()*sizeof(std::unordered_map<pair64,int>::value_type); }
// };

// /* --- sort */
// struct AggSort {
//     std::vector<pair64> buf; uint64_t local = 0;
//     AggSort() = default; AggSort(size_t){}          /* dummy ctor */
//     void  add(pair64 p){ buf.push_back(p); }
//     void  flush(){
//         std::sort(buf.begin(), buf.end());
//         for(size_t i=0;i<buf.size();){
//             size_t j=i; while(j<buf.size() && buf[j]==buf[i]) ++j;
//             size_t f=j-i; if(f>1) local += (uint64_t)f*(f-1)/2;
//             i=j;
//         }
//         buf.clear();
//     }
//     size_t mem() const { return buf.size()*sizeof(pair64); }
// };

// /* --- batch */
// struct AggBatch {
//     static constexpr size_t BATCH = 1<<19;
//     std::vector<pair64> buf; uint64_t local=0;
//     AggBatch() = default; AggBatch(size_t){}        /* dummy ctor */
//     void add(pair64 p){ buf.push_back(p); if(buf.size()==BATCH) flush(); }
//     void flush(){
//         if(buf.empty()) return;
//         std::sort(buf.begin(), buf.end());
//         for(size_t i=0;i<buf.size();){
//             size_t j=i; while(j<buf.size()&&buf[j]==buf[i]) ++j;
//             size_t f=j-i; if(f>1) local += (uint64_t)f*(f-1)/2;
//             i=j;
//         }
//         buf.clear();
//     }
//     size_t mem() const { return buf.size()*sizeof(pair64); }
// };

// /* --- histo */
// struct AggHisto {
//     std::vector<int> hist; std::vector<vid> touched; uint64_t local=0;
//     explicit AggHisto(size_t L=0){ hist.assign(L,0); }
//     void add(pair64 p){
//         vid w = p & 0xFFFFFFFFu;
//         if(hist[w]++ == 0) touched.push_back(w);
//     }
//     void flush(){
//         for(vid w:touched){ int c=hist[w]; if(c>1) local+=(uint64_t)c*(c-1)/2; hist[w]=0; }
//         touched.clear();
//     }
//     size_t mem() const { return touched.size()*sizeof(vid); }
// };

// /* -------------------------------------------------- owner-hashing driver */
// template<class Agg>
// uint64_t owner_hash(const Graph& G,
//                     const std::vector<vid>& rank_of,
//                     int nT)
// {
//     /* per-thread state */
//     std::vector<Agg>              local(nT, Agg(G.L));
//     std::vector<std::vector<pair64>> inbox(nT);
//     std::vector<omp_lock_t>       lock(nT);
//     for(int t=0;t<nT;++t) omp_init_lock(&lock[t]);

//     auto owner = [=](vid u, vid w){ return (u ^ w) % nT; };

// #pragma omp parallel num_threads(nT)
//     {
//         int tid = omp_get_thread_num();
//         auto& L  = local[tid];

// #pragma omp for schedule(dynamic,64) nowait
//         for(size_t rid=0; rid<G.R2L.size(); ++rid){
//             const auto& Aset = G.R2L[rid];
//             for(size_t i=0;i<Aset.size();++i){
//                 vid u=Aset[i];
//                 for(size_t j=i+1;j<Aset.size();++j){
//                     vid w=Aset[j];
//                     if(rank_of[u]>rank_of[w]) std::swap(u,w);
//                     int ow = owner(u,w);
//                     if(ow==tid) L.add(pack(u,w));
//                     else{
//                         omp_set_lock(&lock[ow]);
//                         inbox[ow].push_back(pack(u,w));
//                         omp_unset_lock(&lock[ow]);
//                     }
//                 }
//             }
//             if(L.mem()>FLUSH_BYTES) L.flush();
//         }

// #pragma omp barrier
//         /* drain my inbox */
//         {
//             std::vector<pair64> tmp;
//             omp_set_lock(&lock[tid]);  tmp.swap(inbox[tid]);  omp_unset_lock(&lock[tid]);
//             for(pair64 p:tmp) L.add(p);
//             L.flush();
//         }
//     }
//     uint64_t tot=0; for(auto& L:local) tot += L.local;
//     return tot;
// }

// /* -------------------------------------------------- CLI */
// struct Params{
//     std::string file; vid L;
//     RankType rank=RankType::ID; uint64_t seed=1;
//     enum {HASH,SORT,BATCH,HISTO} agg=HASH;
//     int threads=omp_get_max_threads();
// };
// Params parse(int argc,char**argv){
//     if(argc<3) throw std::runtime_error("syntax");
//     Params P; P.file=argv[1]; P.L=std::stoul(argv[2]);
//     for(int i=3;i<argc;++i){
//         std::string a=argv[i];
//         if(a=="--threads") { P.threads=std::stoi(argv[++i]); }
//         else if(a=="--rank"){
//             std::string v=argv[++i];
//             if(v=="id") P.rank=RankType::ID;
//             else if(v=="degree") P.rank=RankType::DEGREE;
//             else if(v=="degcore")P.rank=RankType::DEGCORE;
//             else if(v.rfind("random",0)==0){
//                 P.rank=RankType::RANDOM;
//                 auto c=v.find(':'); if(c!=std::string::npos) P.seed=std::stoull(v.substr(c+1));
//             }else throw std::runtime_error("bad rank");
//         }else if(a=="--agg"){
//             std::string v=argv[++i];
//             if(v=="hash") P.agg=Params::HASH;
//             else if(v=="sort") P.agg=Params::SORT;
//             else if(v=="batch")P.agg=Params::BATCH;
//             else if(v=="histo")P.agg=Params::HISTO;
//             else throw std::runtime_error("bad agg");
//         }else throw std::runtime_error("unknown flag "+a);
//     }
//     return P;
// }

// /* -------------------------------------------------- main */
// int main(int argc,char**argv)
// {
//     try{
//         Params P=parse(argc,argv);

//         double t0=omp_get_wtime();
//         Graph G = load_graph(P.file,P.L);
//         double t1=omp_get_wtime();

//         auto order = make_rank(G,P.rank,P.seed);
//         std::vector<vid> rank_of(G.L);
//         for(vid i=0;i<G.L;++i) rank_of[order[i]]=i;

//         uint64_t ans=0;
//         switch(P.agg){
//             case Params::HASH : ans=owner_hash<AggHash >(G,rank_of,P.threads); break;
//             case Params::SORT : ans=owner_hash<AggSort >(G,rank_of,P.threads); break;
//             case Params::BATCH: ans=owner_hash<AggBatch>(G,rank_of,P.threads); break;
//             case Params::HISTO: ans=owner_hash<AggHisto>(G,rank_of,P.threads); break;
//         }
//         double t2=omp_get_wtime();

//         std::cout<<"[par] butterflies="<<ans
//                  <<" load="<<t1-t0<<"s"
//                  <<" count="<<t2-t1<<"s"
//                  <<" thr="<<P.threads
//                  <<" rank="<<(int)P.rank
//                  <<" agg="<<P.agg
//                  <<'\n';
//     }catch(const std::exception& e){
//         std::cerr<<"error: "<<e.what()
//                  <<"\nUsage: butter_par edge_file L_size "
//                  "--rank id|degree|degcore|random[:seed] "
//                  "--agg hash|sort|batch|histo "
//                  "[--threads N]\n";
//         return 1;
//     }
// }


/****************************************************************************************
 * butterfly_parallel.cpp  –  deterministic ParButterfly counter (owner hashing)
 * Build : g++ -std=c++17 -O3 -fopenmp butterfly_parallel.cpp -o butter_par
 * Usage : ./butter_par edge_file L_size
 *         --rank  id | degree | degcore | random[:seed]
 *         --agg   hash | sort | batch | histo
 *         [--threads N]
 *****************************************************************************************/
#include <bits/stdc++.h>
#include <omp.h>

/* ─────────────────────────────────── typedefs & helpers ─────────────────────────── */
using vid   = uint32_t;
using pair64= uint64_t;
static inline pair64 pack(vid a, vid b){ return (uint64_t)a<<32 | b; }

/* ─────────────────────────────────── Graph loader ───────────────────────────────── */
struct Graph{
    vid L;
    std::vector<std::vector<vid>> L2R, R2L;
};
Graph load_graph(const std::string& path, vid L)
{
    std::ifstream in(path); if(!in) throw std::runtime_error("open "+path);
    std::vector<std::pair<vid,vid>> E; size_t Rmax=0;
    for(std::string line; std::getline(in,line);){
        if(line.empty()||line[0]=='%') continue;
        vid u,v; std::istringstream(line)>>u>>v;
        E.emplace_back(u,v); Rmax=std::max<size_t>(Rmax,v-L);
    }
    Graph G; G.L=L; G.L2R.assign(L,{}); G.R2L.assign(Rmax+1,{});
    for(auto [u,v]:E){ G.L2R[u].push_back(v); G.R2L[v-L].push_back(u); }
    for(auto& nb:G.L2R) std::sort(nb.begin(),nb.end());
    return G;
}

/* ─────────────────────────────────── Ranking ────────────────────────────────────── */
enum class RankType{ ID, DEGREE, DEGCORE, RANDOM };
std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1)
{
    vid n=G.L; std::vector<vid> ord(n); std::iota(ord.begin(),ord.end(),0);
    if(typ==RankType::ID) return ord;
    if(typ==RankType::DEGREE){
        std::sort(ord.begin(),ord.end(),
                  [&](vid a,vid b){return G.L2R[a].size()<G.L2R[b].size();});
        return ord;
    }
    if(typ==RankType::RANDOM){
        std::mt19937_64 rng(seed); std::shuffle(ord.begin(),ord.end(),rng); return ord;
    }
    /* one-pass degcore */
    std::vector<vid> deg(n); for(vid u=0;u<n;++u) deg[u]=G.L2R[u].size();
    vid maxd=*std::max_element(deg.begin(),deg.end());
    std::vector<std::vector<vid>> bin(maxd+1);
    for(vid u=0;u<n;++u) bin[deg[u]].push_back(u);
    vid out=0;
    for(vid d=0; d<=maxd; ++d)
        while(!bin[d].empty()){
            vid u=bin[d].back(); bin[d].pop_back();
            ord[out++]=u;
            for(vid r:G.L2R[u])
                for(vid w:G.R2L[r-G.L]) if(deg[w]>d){
                    --deg[w]; bin[deg[w]].push_back(w);
                }
        }
    return ord;
}

/* ─────────────────────────────── Aggregation back-ends ──────────────────────────── */
static constexpr size_t FLUSH_BYTES = 8*1024*1024;

/* HASH – incremental choose-two */
struct AggHash{
    std::unordered_map<pair64,int> h; uint64_t local=0;
    AggHash()=default; AggHash(size_t){}                /* dummy ctor */
    void add(pair64 p){
        int& k=h[p]; local+=k; ++k;                     // incremental Δ
    }
    void flush(){ h.clear(); }
    size_t mem() const { return h.size()*sizeof(std::unordered_map<pair64,int>::value_type); }
};

/* SORT (vector) */
struct AggSort{
    std::vector<pair64> buf; uint64_t local=0;
    AggSort()=default; AggSort(size_t){}                /* dummy ctor */
    void add(pair64 p){ buf.push_back(p); }
    void flush(){
        std::sort(buf.begin(),buf.end());
        for(size_t i=0;i<buf.size();){
            size_t j=i; while(j<buf.size()&&buf[j]==buf[i]) ++j;
            size_t f=j-i; local += (uint64_t)f*(f-1)/2;
            i=j;
        }
        buf.clear();
    }
    size_t mem() const { return buf.size()*sizeof(pair64); }
};

/* BATCH (fixed chunk) */
struct AggBatch{
    static constexpr size_t BATCH=1<<19;
    std::vector<pair64> buf; uint64_t local=0;
    AggBatch()=default; AggBatch(size_t){}              /* dummy ctor */
    void add(pair64 p){ buf.push_back(p); if(buf.size()==BATCH) flush(); }
    void flush(){
        if(buf.empty()) return;
        std::sort(buf.begin(),buf.end());
        for(size_t i=0;i<buf.size();){
            size_t j=i; while(j<buf.size()&&buf[j]==buf[i]) ++j;
            size_t f=j-i; local += (uint64_t)f*(f-1)/2; i=j;
        }
        buf.clear();
    }
    size_t mem() const { return buf.size()*sizeof(pair64); }
};

/* HISTO */
struct AggHisto{
    std::vector<int> hist; std::vector<vid> touched; uint64_t local=0;
    explicit AggHisto(size_t L=0){ hist.assign(L,0); }
    void add(pair64 p){
        vid w=p & 0xFFFFFFFFu; local+=hist[w];
        if(hist[w]++==0) touched.push_back(w);
    }
    void flush(){ for(vid w:touched) hist[w]=0; touched.clear(); }
    size_t mem() const { return touched.size()*sizeof(vid); }
};

/* ───────────────────────── Owner-hashing driver (deterministic) ─────────────────── */
template<class Agg>
uint64_t owner_hash(const Graph& G,
                    const std::vector<vid>& rank_of,
                    int nT)
{
    std::vector<Agg>                  local(nT, Agg(G.L));
    std::vector<std::vector<pair64>>  inbox(nT);
    std::vector<omp_lock_t>           lock(nT);
    for(int t=0;t<nT;++t) omp_init_lock(&lock[t]);

    auto owner=[=](vid u,vid w){ return (u^w)%nT; };

#pragma omp parallel num_threads(nT)
    {
        int tid=omp_get_thread_num();
        auto& L=local[tid];

#pragma omp for schedule(dynamic,64) nowait
        for(size_t rid=0; rid<G.R2L.size(); ++rid){
            const auto& Aset=G.R2L[rid];
            for(size_t i=0;i<Aset.size();++i){
                vid u=Aset[i];
                for(size_t j=i+1;j<Aset.size();++j){
                    vid w=Aset[j];
                    if(rank_of[u]>rank_of[w]) std::swap(u,w);
                    int ow=owner(u,w);
                    if(ow==tid) L.add(pack(u,w));
                    else{
                        omp_set_lock(&lock[ow]);
                        inbox[ow].push_back(pack(u,w));
                        omp_unset_lock(&lock[ow]);
                    }
                }
            }
            if(L.mem()>FLUSH_BYTES) L.flush();
        }

#pragma omp barrier
        /* drain inbox */
        std::vector<pair64> tmp;
        omp_set_lock(&lock[tid]); tmp.swap(inbox[tid]); omp_unset_lock(&lock[tid]);
        for(pair64 p:tmp) L.add(p); L.flush();
    }
    uint64_t tot=0; for(auto& L:local) tot+=L.local; return tot;
}

/* ───────────────────────────── CLI & main ───────────────────────────────────────── */
struct Params{
    std::string file; vid L;
    RankType rank=RankType::ID; uint64_t seed=1;
    enum {HASH,SORT,BATCH,HISTO} agg=HASH;
    int threads=omp_get_max_threads();
};
Params parse(int argc,char**argv){
    if(argc<3) throw std::runtime_error("syntax");
    Params P; P.file=argv[1]; P.L=std::stoul(argv[2]);
    for(int i=3;i<argc;++i){
        std::string a=argv[i];
        if(a=="--threads") P.threads=std::stoi(argv[++i]);
        else if(a=="--rank"){
            std::string v=argv[++i];
            if(v=="id") P.rank=RankType::ID;
            else if(v=="degree") P.rank=RankType::DEGREE;
            else if(v=="degcore")P.rank=RankType::DEGCORE;
            else if(v.rfind("random",0)==0){
                P.rank=RankType::RANDOM;
                auto c=v.find(':'); if(c!=std::string::npos) P.seed=std::stoull(v.substr(c+1));
            }else throw std::runtime_error("bad rank");
        }else if(a=="--agg"){
            std::string v=argv[++i];
            if(v=="hash") P.agg=Params::HASH;
            else if(v=="sort") P.agg=Params::SORT;
            else if(v=="batch")P.agg=Params::BATCH;
            else if(v=="histo")P.agg=Params::HISTO;
            else throw std::runtime_error("bad agg");
        }else throw std::runtime_error("unknown "+a);
    }
    return P;
}

template<class A> uint64_t run(const Graph& G, const std::vector<vid>& rank, int nT){
    return owner_hash<A>(G,rank,nT);
}

int main(int argc,char**argv)
{
    try{
        Params P=parse(argc,argv);

        double t0=omp_get_wtime();
        Graph G=load_graph(P.file,P.L);
        double t1=omp_get_wtime();

        auto order=make_rank(G,P.rank,P.seed);
        std::vector<vid> rank_of(G.L); for(vid i=0;i<G.L;++i) rank_of[order[i]]=i;

        uint64_t ans=0;
        switch(P.agg){
            case Params::HASH : ans=run<AggHash >(G,rank_of,P.threads); break;
            case Params::SORT : ans=run<AggSort >(G,rank_of,P.threads); break;
            case Params::BATCH: ans=run<AggBatch>(G,rank_of,P.threads); break;
            case Params::HISTO: ans=run<AggHisto>(G,rank_of,P.threads); break;
        }
        double t2=omp_get_wtime();

        std::cout<<"[par] butterflies="<<ans
                 <<" load="<<t1-t0<<"s"
                 <<" count="<<t2-t1<<"s"
                 <<" thr="<<P.threads
                 <<" rank="<<(int)P.rank
                 <<" agg="<<P.agg
                 <<'\n';
    }catch(const std::exception& e){
        std::cerr<<"error: "<<e.what()
                 <<"\nUsage: butter_par edge_file L_size "
                 "--rank id|degree|degcore|random[:seed] "
                 "--agg hash|sort|batch|histo "
                 "[--threads N]\n";
        return 1;
    }
}
