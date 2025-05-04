// butter.cpp: Standalone Sequential/OpenMP version

#include "butter_core.hpp" // Include all core logic

#include <iostream>
#include <iomanip>
#include <stdexcept>

int main(int argc, char** argv){
    int mpi_rank = 0; // Standalone is effectively rank 0
    Params P;
    try {
        P = parse(argc, argv);
        // --- Validation for Standalone ---
        if (P.file.empty()) {
            throw std::runtime_error("Missing required argument: --file <edge_file>");
        }
        if (!P.subgraph_dir.empty()) {
             std::cerr << "[" << mpi_rank << "] Warning: --subgraph_dir provided but ignored in standalone mode." << std::endl;
        }
        // --- End Validation ---

    } catch (const std::exception &e) {
        std::cerr << "[" << mpi_rank << "] Error parsing arguments: " << e.what() << "\n";
        std::cerr << "[" << mpi_rank << "] Usage: ./butter --file <edge_file> --L <L_size> [--mode seq|par] [--rank ...] [--agg ...] [--threads N]\n";
        return 1;
    }

    try {
        double t_start_load = get_time();
        // Load the single graph file specified
        Graph G = load_graph(P.file, P.L);
        double t_end_load = get_time();

        if (G.L == 0 && G.R2L.empty()) { /* Info msg */ }
        else if (G.L == 0) { /* Warning msg */ }

        double t_start_rank = get_time();
        auto order = make_rank(G, P.rank, P.seed); // Rank the loaded graph
        double t_end_rank = get_time();

        std::vector<vid> rank_of; // Compute rank_of array
        if (G.L > 0) {
            if (order.size() != G.L) throw std::runtime_error("Rank size mismatch.");
            rank_of.resize(G.L);
            for (vid r = 0; r < G.L; ++r) {
                vid node_id = order[r];
                if (node_id >= G.L) throw std::runtime_error("Rank node ID OOB.");
                rank_of[node_id] = r;
            }
        } else if (!order.empty()) throw std::runtime_error("L=0 but rank non-empty.");

        uint64_t ans = 0;
        double t_start_count = get_time();
        int nT_to_use = P.threads;

        // ** Call SIMPLIFIED count_par (no range needed) **
        switch (P.agg) {
            case 0: ans = count_par<AggHash >(G, rank_of, nT_to_use); break;
            case 1: ans = count_par<AggSort >(G, rank_of, nT_to_use); break;
            case 2: ans = count_par<AggBatch>(G, rank_of, nT_to_use); break;
            case 3:
                 if (nT_to_use > 0 && G.L > 0) { /* Histo warning */ }
                 ans = count_par<AggHisto>(G, rank_of, nT_to_use); break;
            default: throw std::logic_error("Invalid aggregator type.");
        }
        double t_end_count = get_time();

        // Output results
        std::cout << "[" << (P.mode == ModeType::SEQ ? "seq" : "par") << "]"
                  << " butterflies=" << ans
                  // ... rest of output fields ...
                  << " load=" << std::fixed << std::setprecision(5) << (t_end_load - t_start_load) << "s"
                  << " rank_t=" << std::fixed << std::setprecision(5) << (t_end_rank - t_start_rank) << "s"
                  << " count=" << std::fixed << std::setprecision(5) << (t_end_count - t_start_count) << "s"
                  << " thr=" << nT_to_use
                  << " rank=" << static_cast<int>(P.rank)
                  << " agg=" << P.agg
                  << "\n";

    } catch (const std::exception &e) {
        std::cerr << "[" << mpi_rank << "] Runtime Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}