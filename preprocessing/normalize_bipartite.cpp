/*
 * normalize_bipartite.cpp
 * -----------------------
 * Converts KONECT dblp-author edge list into:
 *   • normalized_edges.txt   (0-based; authors first, papers shifted by |L|)
 *   • map_L.txt              (authorNewID  authorOrigID)
 *   • map_R.txt              (paperNewID   paperOrigID)
 *
 * Build : g++ -std=c++17 -O3 normalize_bipartite.cpp -o normalize
 * Usage : ./normalize <edge_file> <output_dir>
 */

 #include <bits/stdc++.h>
 #include <filesystem>
 
 namespace fs = std::filesystem;
 using  orig_t = int;   // original IDs in raw file (1-based)
 using  vid_t  = int;   // compact 0-based IDs we assign  (renamed!)
 
 /* ---------- Pass-1 : collect unique authors / papers ------------------ */
 struct Stats {
     std::unordered_map<orig_t,vid_t> mapL, mapR;
     size_t edges = 0;
 };
 
 static Stats first_pass(const std::string& file)
 {
     Stats st;
     std::ifstream in(file);
     if (!in) throw std::runtime_error("cannot open " + file);
 
     vid_t nextA = 0, nextP = 0;
     std::string line;
     while (std::getline(in,line)) {
         if (line.empty() || line[0]=='%') continue;
         std::istringstream iss(line);
         orig_t a,p; iss>>a>>p;
         if (st.mapL.find(a)==st.mapL.end()) st.mapL[a]=nextA++;
         if (st.mapR.find(p)==st.mapR.end()) st.mapR[p]=nextP++;
         ++st.edges;
     }
     return st;
 }
 
 /* ------------------------------ main ---------------------------------- */
 int main(int argc,char**argv)
 {
     if (argc!=3) {
         std::cerr<<"Usage: "<<argv[0]<<" <edge_file> <output_dir>\n";
         return 1;
     }
     const std::string inFile = argv[1];
     const fs::path outDir   = argv[2];
     if (!fs::exists(outDir)) fs::create_directories(outDir);
 
     Stats s = first_pass(inFile);
     const vid_t L = s.mapL.size();          // number of authors (shift for R)
 
     /* write maps */
     {
         std::ofstream lm(outDir/"map_L.txt");
         for (auto& [orig,newID] : s.mapL) lm<<newID<<' '<<orig<<'\n';
     }
     {
         std::ofstream rm(outDir/"map_R.txt");
         for (auto& [orig,newID] : s.mapR) rm<<(newID+L)<<' '<<orig<<'\n';
     }
 
     /* second pass: emit normalized edges */
     std::ifstream in(inFile);
     std::ofstream norm(outDir/"normalized_edges.txt");
     std::string line;
     while (std::getline(in,line)) {
         if (line.empty()||line[0]=='%') continue;
         std::istringstream iss(line);
         orig_t a,p; iss>>a>>p;
         vid_t ua = s.mapL[a];
         vid_t vp = s.mapR[p] + L;
         norm<<ua<<' '<<vp<<'\n';
     }
 
     std::cout<<"[normalize] Authors(L): "<<L
              <<"\n[normalize] Papers(R) : "<<s.mapR.size()
              <<"\n[normalize] Edges    : "<<s.edges
              <<"\n[normalize] Output   : "<<outDir<<'\n';
     return 0;
 }
 