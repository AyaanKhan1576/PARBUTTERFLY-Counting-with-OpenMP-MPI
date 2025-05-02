
// /****************************************************************************************
//  * butterfly_sequential.cpp  –  single-thread, deterministic ParButterfly counter
//  *
//  * Build : g++ -std=c++17 -O3 butterfly_sequential.cpp -o butter_seq
//  * Usage : ./butter_seq edge_file L_size
//  *         --rank  id | degree | degcore | random[:seed]
//  *         --agg   hash | sort | batch | histo
//  *****************************************************************************************/
// #include <bits/stdc++.h>

// /* ───────────────────── typedefs & helpers ───────────────────── */
// using vid   = uint32_t;
// using pair64= uint64_t;
// static inline pair64 pack(vid a, vid b){ return (uint64_t)a<<32 | b; }

// /* ───────────────────── Graph loader ─────────────────────────── */
// struct Graph{
//     vid L;
//     std::vector<std::vector<vid>> L2R, R2L;
// };
// Graph load_graph(const std::string& path, vid L)
// {
//     std::ifstream in(path); if(!in) throw std::runtime_error("open "+path);
//     std::vector<std::pair<vid,vid>> E; size_t Rmax=0;
//     for(std::string line; std::getline(in,line);){
//         if(line.empty()||line[0]=='%') continue;
//         vid u,v; std::istringstream(line)>>u>>v;
//         E.emplace_back(u,v); Rmax=std::max<size_t>(Rmax,v-L);
//     }
//     Graph G; G.L=L;
//     G.L2R.assign(L,{}); G.R2L.assign(Rmax+1,{});
//     for(auto [u,v]:E){ G.L2R[u].push_back(v); G.R2L[v-L].push_back(u); }
//     for(auto& nb:G.L2R) std::sort(nb.begin(),nb.end());
//     return G;
// }

// /* ───────────────────── Ranking options ─────────────────────── */
// enum class RankType{ ID, DEGREE, DEGCORE, RANDOM };

// std::vector<vid> make_rank(const Graph& G, RankType typ, uint64_t seed=1)
// {
//     vid n=G.L; std::vector<vid> ord(n); std::iota(ord.begin(),ord.end(),0);
//     if(typ==RankType::ID) return ord;

//     if(typ==RankType::DEGREE){
//         std::sort(ord.begin(),ord.end(),
//                   [&](vid a,vid b){ return G.L2R[a].size()<G.L2R[b].size();});
//         return ord;
//     }
//     if(typ==RankType::RANDOM){
//         std::mt19937_64 rng(seed);
//         std::shuffle(ord.begin(),ord.end(),rng);
//         return ord;
//     }
//     /* one-pass approximate core order */
//     std::vector<vid> deg(n); for(vid u=0;u<n;++u) deg[u]=G.L2R[u].size();
//     vid maxd=*std::max_element(deg.begin(),deg.end());
//     std::vector<std::vector<vid>> bins(maxd+1);
//     for(vid u=0;u<n;++u) bins[deg[u]].push_back(u);
//     vid out=0;
//     for(vid d=0; d<=maxd; ++d)
//         while(!bins[d].empty()){
//             vid u=bins[d].back(); bins[d].pop_back();
//             ord[out++]=u;
//             for(vid r:G.L2R[u])
//                 for(vid w:G.R2L[r-G.L]) if(deg[w]>d){
//                     --deg[w]; bins[deg[w]].push_back(w);
//                 }
//         }
//     return ord;
// }

// /* ───────────────────── Aggregators ─────────────────────────── */
// static constexpr size_t FLUSH_BYTES = 8*1024*1024;

// /* HASH – incremental Δ choose-two */
// struct AggHash{
//     std::unordered_map<pair64,int> h; uint64_t local=0;
//     AggHash()=default; AggHash(size_t){}                /* dummy ctor */
//     void add(pair64 p){ int& k=h[p]; local+=k; ++k; }
//     void flush(){ h.clear(); }
//     size_t mem()const{return h.size()*sizeof(std::unordered_map<pair64,int>::value_type);}
// };

// /* SORT */
// struct AggSort{
//     std::vector<pair64> buf; uint64_t local=0;
//     AggSort()=default; AggSort(size_t){}                /* dummy ctor */
//     void add(pair64 p){ buf.push_back(p); }
//     void flush(){
//         std::sort(buf.begin(),buf.end());
//         for(size_t i=0;i<buf.size();){
//             size_t j=i; while(j<buf.size()&&buf[j]==buf[i]) ++j;
//             size_t f=j-i; local+=(uint64_t)f*(f-1)/2; i=j;
//         }
//         buf.clear();
//     }
//     size_t mem() const { return buf.size()*sizeof(pair64); }
// };

// /* BATCH – fixed chunk */
// struct AggBatch: AggSort{
//     static constexpr size_t B = 1<<19;
//     AggBatch()=default; AggBatch(size_t s):AggSort(s){}
//     void add(pair64 p){ buf.push_back(p); if(buf.size()==B) flush(); }
// };

// /* HISTO – incremental Δ choose-two */
// struct AggHisto{
//     std::vector<int> hist; std::vector<vid> touched; uint64_t local=0;
//     explicit AggHisto(size_t L){ hist.assign(L,0); }
//     void add(pair64 p){
//         vid w = p & 0xFFFFFFFFu;
//         local += hist[w];
//         if(hist[w]++==0) touched.push_back(w);
//     }
//     void flush(){ for(vid w:touched) hist[w]=0; touched.clear(); }
//     size_t mem()const{return touched.size()*sizeof(vid);}
// };

// /* ───────────────────── Sequential driver ───────────────────── */
// template<class Agg>
// uint64_t count_seq(const Graph& G,
//                    const std::vector<vid>& rank)
// {
//     Agg A(G.L);
//     for(size_t rid=0; rid<G.R2L.size(); ++rid){
//         const auto& aSet = G.R2L[rid];
//         for(size_t i=0;i<aSet.size();++i){
//             vid u=aSet[i];
//             for(size_t j=i+1;j<aSet.size();++j){
//                 vid w=aSet[j];
//                 if(rank[u] > rank[w]) std::swap(u,w);
//                 A.add(pack(u,w));
//             }
//         }
//         if(A.mem()>FLUSH_BYTES) A.flush();
//     }
//     A.flush();
//     return A.local;
// }

// /* ───────────────────── CLI / main ─────────────────────────── */
// struct Params{
//     std::string file; vid L;
//     RankType rank=RankType::ID; uint64_t seed=1;
//     enum {HASH,SORT,BATCH,HISTO} agg=HASH;
// };

// Params parse(int argc,char**argv){
//     if(argc<3) throw std::runtime_error("syntax");
//     Params P; P.file=argv[1]; P.L=std::stoul(argv[2]);
//     for(int i=3;i<argc;++i){
//         std::string a=argv[i];
//         if(a=="--rank"){
//             std::string v=argv[++i];
//             if(v=="degree") P.rank=RankType::DEGREE;
//             else if(v=="degcore")P.rank=RankType::DEGCORE;
//             else if(v.rfind("random",0)==0){
//                 P.rank=RankType::RANDOM;
//                 auto c=v.find(':'); if(c!=std::string::npos) P.seed=std::stoull(v.substr(c+1));
//             }
//         }else if(a=="--agg"){
//             std::string v=argv[++i];
//             if(v=="sort") P.agg=Params::SORT;
//             else if(v=="batch")P.agg=Params::BATCH;
//             else if(v=="histo")P.agg=Params::HISTO;
//         }else throw std::runtime_error("unknown flag "+a);
//     }
//     return P;
// }

// int main(int argc,char**argv)
// {
//     try{
//         Params P=parse(argc,argv);

//         double t0=clock()/double(CLOCKS_PER_SEC);
//         Graph  G = load_graph(P.file,P.L);
//         double t1=clock()/double(CLOCKS_PER_SEC);

//         auto order=make_rank(G,P.rank,P.seed);
//         std::vector<vid> rank_of(G.L); for(vid i=0;i<G.L;++i) rank_of[order[i]]=i;

//         uint64_t ans=0;
//         if(P.agg==Params::HASH ) ans=count_seq<AggHash >(G,rank_of);
//         else if(P.agg==Params::SORT ) ans=count_seq<AggSort >(G,rank_of);
//         else if(P.agg==Params::BATCH) ans=count_seq<AggBatch>(G,rank_of);
//         else                           ans=count_seq<AggHisto>(G,rank_of);

//         double t2=clock()/double(CLOCKS_PER_SEC);

//         std::cout<<"[seq] butterflies="<<ans
//                  <<" load="<<t1-t0<<"s"
//                  <<" count="<<t2-t1<<"s"
//                  <<" rank="<<(int)P.rank
//                  <<" agg="<<P.agg
//                  <<'\n';
//     }catch(const std::exception& e){
//         std::cerr<<"error: "<<e.what()
//                  <<"\nUsage: butter_seq edge_file L_size "
//                  "--rank id|degree|degcore|random[:seed] "
//                  "--agg hash|sort|batch|histo\n";
//         return 1;
//     }
// }


/****************************************************************************************
 * butterfly_sequential.cpp  –  simple, fully-deterministic butterfly counter
 *
 * NO incremental tricks, NO flushing: we store the whole (u,w) frequency map
 * and do one C(freq,2) scan at the end.  This is the ground-truth reference
 * for validating the OpenMP/MPI versions.
 *
 * Build : g++ -std=c++17 -O3 butterfly_sequential.cpp -o butter_seq
 * Usage : ./butter_seq edge_file L_size
 *         --rank  id | degree | degcore | random[:seed]
 *****************************************************************************************/
#include <bits/stdc++.h>

using vid   = uint32_t;
using pair64= uint64_t;
static inline pair64 pack(vid a, vid b){ return (uint64_t)a<<32 | b; }

/* ─────────── Graph ─────────── */
struct Graph{
    vid L;
    std::vector<std::vector<vid>> L2R, R2L;
};
Graph load_graph(const std::string& p, vid L)
{
    std::ifstream in(p); if(!in) throw std::runtime_error("open "+p);
    std::vector<std::pair<vid,vid>> E; size_t Rmax=0;
    for(std::string l; std::getline(in,l);){
        if(l.empty()||l[0]=='%') continue;
        vid u,v; std::istringstream(l)>>u>>v;
        E.emplace_back(u,v); Rmax=std::max<size_t>(Rmax,v-L);
    }
    Graph G; G.L=L; G.L2R.assign(L,{}); G.R2L.assign(Rmax+1,{});
    for(auto [u,v]:E){ G.L2R[u].push_back(v); G.R2L[v-L].push_back(u); }
    for(auto& nb:G.L2R) std::sort(nb.begin(),nb.end());
    return G;
}

/* ─────────── Ranking (ID, Degree, DegCore, Random) ─────────── */
enum class RankType{ID,DEGREE,DEGCORE,RANDOM};
std::vector<vid> make_rank(const Graph& G, RankType t, uint64_t seed=1)
{
    vid n=G.L; std::vector<vid> o(n); std::iota(o.begin(),o.end(),0);
    if(t==RankType::ID) return o;
    if(t==RankType::DEGREE){
        std::sort(o.begin(),o.end(),
                  [&](vid a,vid b){return G.L2R[a].size()<G.L2R[b].size();});
        return o;
    }
    if(t==RankType::RANDOM){
        std::mt19937_64 rng(seed);
        std::shuffle(o.begin(),o.end(),rng);
        return o;
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
            o[out++]=u;
            for(vid r:G.L2R[u])
                for(vid w:G.R2L[r-G.L]) if(deg[w]>d){
                    --deg[w]; bin[deg[w]].push_back(w);
                }
        }
    return o;
}

/* ─────────── Counting (single pass, full map) ─────────── */
uint64_t butterflies_seq(const Graph& G,
                         const std::vector<vid>& rank_of)
{
    std::unordered_map<pair64,uint32_t> freq;
    freq.reserve(1<<22);                         // pre-size ~4M entries

    for(size_t rid=0; rid<G.R2L.size(); ++rid){
        const auto& Aset = G.R2L[rid];
        for(size_t i=0;i<Aset.size();++i){
            vid u=Aset[i];
            for(size_t j=i+1;j<Aset.size();++j){
                vid w=Aset[j];
                if(rank_of[u]>rank_of[w]) std::swap(u,w);
                ++freq[pack(u,w)];
            }
        }
    }
    uint64_t total=0;
    for(auto& kv:freq){
        uint64_t f=kv.second;
        total += f*(f-1)/2;
    }
    return total;
}

/* ─────────── CLI + main ─────────── */
int main(int argc,char**argv)
{
    if(argc<3){
        std::cerr<<"Usage: butter_seq edge_file L_size "
                 "--rank id|degree|degcore|random[:seed]\n";
        return 1;
    }
    std::string file = argv[1];  vid L = std::stoul(argv[2]);

    RankType rT = RankType::ID;  uint64_t seed=1;
    for(int i=3;i<argc;++i){
        std::string a=argv[i];
        if(a=="--rank"){
            std::string v=argv[++i];
            if(v=="degree")  rT=RankType::DEGREE;
            else if(v=="degcore") rT=RankType::DEGCORE;
            else if(v.rfind("random",0)==0){
                rT=RankType::RANDOM;
                auto c=v.find(':'); if(c!=std::string::npos) seed=std::stoull(v.substr(c+1));
            }
        }
    }

    double t0=clock()/double(CLOCKS_PER_SEC);
    Graph  G = load_graph(file,L);
    double t1=clock()/double(CLOCKS_PER_SEC);

    auto order = make_rank(G,rT,seed);
    std::vector<vid> rank_of(G.L); for(vid i=0;i<G.L;++i) rank_of[order[i]]=i;

    uint64_t ans = butterflies_seq(G,rank_of);
    double t2=clock()/double(CLOCKS_PER_SEC);

    std::cout<<"[seq] butterflies="<<ans
             <<" load="<<t1-t0<<"s"
             <<" count="<<t2-t1<<"s"
             <<" rank="<<(int)rT
             <<'\n';
    return 0;
}
