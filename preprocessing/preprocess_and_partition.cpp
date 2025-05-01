
//  /*
//  * preprocess_and_partition.cpp
//  * ----------------------------
//  * Build : g++ -std=c++17 -O3 preprocess_and_partition.cpp -lmetis -o preprocess
//  * Usage : ./preprocess <edge_file> <k_parts> <output_dir>
//  */

// #include <bits/stdc++.h>
// #include <metis.h>
// #include <filesystem>

// namespace fs = std::filesystem;
// using idx = idx_t;          // METIS integer
// using vid_t = int;          // matches normalize_bipartite

// int main(int argc,char**argv)
// {
//     if (argc!=4){
//         std::cerr<<"Usage: "<<argv[0]<<" <edge_file> <k_parts> <out_dir>\n";
//         return 1;
//     }
//     const std::string edgeFile = argv[1];
//     const idx kParts = std::stoi(argv[2]);
//     const fs::path outDir = argv[3];
//     if (!fs::exists(outDir)) fs::create_directories(outDir);

//     /* -------- 1. read edges -------- */
//     std::ifstream in(edgeFile);
//     if(!in){ std::cerr<<"cannot open "<<edgeFile<<'\n'; return 1; }

//     std::vector<std::pair<vid_t,vid_t>> edges;
//     vid_t maxId = 0; std::string line;
//     while(std::getline(in,line)){
//         if(line.empty()||line[0]=='%') continue;
//         std::istringstream iss(line);
//         vid_t u,v; iss>>u>>v;
//         if(u==v) continue;
//         edges.emplace_back(u,v);
//         maxId = std::max(maxId, std::max(u,v));
//     }
//     const idx nv = maxId+1;

//     /* CSR build */
//     std::vector<std::vector<idx>> adj(nv);
//     for(auto [u,v]:edges){ adj[u].push_back(v); adj[v].push_back(u); }
//     for(auto& v:adj){ std::sort(v.begin(),v.end()); v.erase(std::unique(v.begin(),v.end()),v.end()); }

//     std::vector<idx> xadj(nv+1);
//     std::vector<idx> adjncy; adjncy.reserve(edges.size()*2);
//     for(idx i=0;i<nv;++i){ xadj[i]=adjncy.size(); adjncy.insert(adjncy.end(),adj[i].begin(),adj[i].end()); }
//     xadj[nv]=adjncy.size(); adj.clear(); edges.clear();

//     /* graph.metis */
//     {
//         std::ofstream g(outDir/"graph.metis");
//         g<<nv<<' '<<adjncy.size()/2<<'\n';
//         for(idx i=0;i<nv;++i){
//             for(idx j=xadj[i];j<xadj[i+1];++j) g<<adjncy[j]+1<<' ';
//             g<<'\n';
//         }
//     }

//     /* METIS */
//     std::vector<idx> part(nv,0);
//     idx nv_copy=nv,ncon=1,objval;
//     idx opts[METIS_NOPTIONS]; METIS_SetDefaultOptions(opts); opts[METIS_OPTION_NUMBERING]=0;
//     int rc = METIS_PartGraphKway(&nv_copy,&ncon,xadj.data(),adjncy.data(),
//                                  nullptr,nullptr,nullptr,
//                                  const_cast<idx*>(&kParts),
//                                  nullptr,nullptr,opts,&objval,part.data());
//     if(rc!=METIS_OK){ std::cerr<<"METIS error "<<rc<<'\n'; return 1; }

//     /* part.<k> */
//     {
//         std::ofstream p(outDir/("part."+std::to_string(kParts)));
//         for(idx v:part) p<<v<<'\n';
//     }

//     /* per-rank edges */
//     std::vector<std::ofstream> f(kParts);
//     for(idx r=0;r<kParts;++r) f[r].open(outDir/("subgraph_"+std::to_string(r)+".txt"));
//     for(idx u=0;u<nv;++u){
//         idx r=part[u];
//         for(idx j=xadj[u];j<xadj[u+1];++j){ idx v=adjncy[j]; if(u<v) f[r]<<u<<' '<<v<<'\n'; }
//     }

//     std::cout<<"[preprocess] nv="<<nv<<", cut="<<objval<<", outputs→"<<outDir<<'\n';
//     return 0;
// }


/*
 * preprocess_and_partition.cpp
 * ----------------------------
 * Build : g++ -std=c++17 -O3 preprocess_and_partition.cpp -lmetis -o preprocess
 * Usage : ./preprocess <edge_file> <k_parts> <output_dir>
 */

 #include <bits/stdc++.h>
 #include <metis.h>
 #include <filesystem>
 
 namespace fs = std::filesystem;
 using idx   = idx_t;   // METIS integer
 using vid_t = int;     // compact vertex IDs produced by normalize_bipartite
 
 int main(int argc, char** argv)
 {
     if (argc != 4) {
         std::cerr << "Usage: " << argv[0]
                   << " <edge_file> <k_parts> <out_dir>\n";
         return 1;
     }
     const std::string edgeFile = argv[1];
     const idx   kParts = std::stoi(argv[2]);
     const fs::path outDir = argv[3];
     if (!fs::exists(outDir)) fs::create_directories(outDir);
 
     /* ------------------------------------------------------------------ */
     std::cout << "[1] Reading edge list …\n";
 
     std::ifstream in(edgeFile);
     if (!in) { std::cerr << "Cannot open " << edgeFile << '\n'; return 1; }
 
     std::vector<std::pair<vid_t,vid_t>> edges;
     vid_t maxId = 0;  std::string line;
     while (std::getline(in, line)) {
         if (line.empty() || line[0] == '%') continue;
         std::istringstream iss(line);
         vid_t u, v;  iss >> u >> v;
         if (u == v) continue;
         edges.emplace_back(u, v);
         maxId = std::max(maxId, std::max(u, v));
     }
     const idx nv = maxId + 1;
     std::cout << "    parsed " << edges.size() << " edges, nv = "
               << nv << '\n';
 
     /* ------------------------------------------------------------------ */
     std::cout << "[2] Building CSR …\n";
 
     std::vector<std::vector<idx>> adj(nv);
     for (auto [u, v] : edges) {
         adj[u].push_back(v);
         adj[v].push_back(u);
     }
     for (auto &nbrs : adj) {
         std::sort(nbrs.begin(), nbrs.end());
         nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
     }
 
     std::vector<idx> xadj(nv + 1);
     std::vector<idx> adjncy;  adjncy.reserve(edges.size() * 2ULL);
     for (idx i = 0; i < nv; ++i) {
         xadj[i] = adjncy.size();
         adjncy.insert(adjncy.end(), adj[i].begin(), adj[i].end());
     }
     xadj[nv] = adjncy.size();
     adj.clear();  edges.clear();
     std::cout << "    CSR arrays sizes: xadj = " << xadj.size()
               << ", adjncy = " << adjncy.size() << '\n';
 
     /* ------------------------------------------------------------------ */
     std::cout << "[3] Writing graph.metis …\n";
 
     {
         std::ofstream g(outDir / "graph.metis");
         g << nv << ' ' << adjncy.size() / 2 << '\n';
         for (idx i = 0; i < nv; ++i) {
             for (idx j = xadj[i]; j < xadj[i + 1]; ++j) g << adjncy[j] + 1 << ' ';
             g << '\n';
         }
     }
 
     /* ------------------------------------------------------------------ */
     std::cout << "[4] Calling METIS (k = " << kParts << ") …\n";
 
     std::vector<idx> part(nv, 0);
     idx nv_copy = nv, ncon = 1, objval;
     idx opts[METIS_NOPTIONS];  METIS_SetDefaultOptions(opts);
     opts[METIS_OPTION_NUMBERING] = 0;
 
     int rc = METIS_PartGraphKway(&nv_copy, &ncon,
                                  xadj.data(), adjncy.data(),
                                  nullptr, nullptr, nullptr,
                                  const_cast<idx*>(&kParts),
                                  nullptr, nullptr,
                                  opts, &objval, part.data());
 
     if (rc != METIS_OK) { std::cerr << "METIS error " << rc << '\n'; return 1; }
     std::cout << "    done   (edge-cut = " << objval << ")\n";
 
     /* ------------------------------------------------------------------ */
     std::cout << "[5] Saving part." << kParts << " …\n";
 
     {
         std::ofstream p(outDir / ("part." + std::to_string(kParts)));
         for (idx v : part) p << v << '\n';
     }
 
     /* ------------------------------------------------------------------ */
     std::cout << "[6] Emitting subgraph files …\n";
 
     std::vector<std::ofstream> fout(kParts);
     for (idx r = 0; r < kParts; ++r)
         fout[r].open(outDir / ("subgraph_" + std::to_string(r) + ".txt"));
 
     for (idx u = 0; u < nv; ++u) {
         idx r = part[u];
         for (idx j = xadj[u]; j < xadj[u + 1]; ++j) {
             idx v = adjncy[j];
             if (u < v) fout[r] << u << ' ' << v << '\n';
         }
     }
     std::cout << "    wrote subgraph_0.txt … subgraph_" << kParts - 1
               << ".txt\n";
 
     /* ------------------------------------------------------------------ */
     std::cout << "[✓] Pre-processing complete.  Outputs in " << outDir << '\n';
     return 0;
 }
 