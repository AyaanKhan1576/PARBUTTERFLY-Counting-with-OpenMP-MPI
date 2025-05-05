// butter_mpi.cpp: MPI wrapper for pre-partitioned subgraphs with individual rank output

#include <mpi.h>
#include "butter_core.hpp" // Include all core logic

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <filesystem> // Required
#include <sstream>    // Required for formatting local output string

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    double mpi_start_time = MPI_Wtime();

    Params P;
    try {
        P = parse(argc, argv);
        if (P.subgraph_dir.empty()) throw std::runtime_error("--subgraph_dir required");
        if (!std::filesystem::exists(P.subgraph_dir) || !std::filesystem::is_directory(P.subgraph_dir)) {
             throw std::runtime_error("Subgraph directory not found: " + P.subgraph_dir);
        }
    } catch (const std::exception &e) {
        if (mpi_rank == 0) { /* Print usage */ }
        MPI_Abort(MPI_COMM_WORLD, 1); return 1;
    }

    Graph G;
    std::vector<vid> rank_of;
    double t_load_rank_start = 0, t_load_rank_end = 0;
    double t_start_count = 0, t_end_count = 0;
    uint64_t local_ans = 0;
    std::string subgraph_filename;

    try {
        t_load_rank_start = MPI_Wtime();
        subgraph_filename = P.subgraph_dir + "/subgraph_" + std::to_string(mpi_rank) + ".txt";
        if (mpi_rank == 0 && world_size > 1) { /* Print info messages */ } // Only rank 0 prints assumptions once
         if (!std::filesystem::exists(subgraph_filename)) {
            std::cerr << "[MPI Rank " << mpi_rank << "] Error: Subgraph file not found: " << subgraph_filename << std::endl;
            // Graceful exit? Or abort? Let's abort for simplicity.
            MPI_Abort(MPI_COMM_WORLD, 2); return 2;
         }
        G = load_graph(subgraph_filename, P.L);
        auto order = make_rank(G, P.rank, P.seed);
        if (G.L > 0) { rank_of.resize(G.L); if(order.size()!=G.L) throw std::runtime_error("rank size"); for(vid r=0;r<G.L;++r){ vid nid=order[r]; if(nid>=G.L) throw std::runtime_error("rank id oob"); rank_of[nid]=r; } }
        else if (!order.empty()) throw std::runtime_error("L=0 rank non-empty");
        t_load_rank_end = MPI_Wtime();

        t_start_count = MPI_Wtime();
        int nT_to_use = P.threads;
        switch (P.agg) { /* Call count_par<...>(G, rank_of, nT_to_use) */
            case 0: local_ans = count_par<AggHash >(G, rank_of, nT_to_use); break;
            case 1: local_ans = count_par<AggSort >(G, rank_of, nT_to_use); break;
            case 2: local_ans = count_par<AggBatch>(G, rank_of, nT_to_use); break;
            case 3: local_ans = count_par<AggHisto>(G, rank_of, nT_to_use); break;
            default: throw std::logic_error("Invalid agg type.");
        }
        t_end_count = MPI_Wtime();

        // --- Calculate Local Timings ---
        double local_load_rank_time = t_load_rank_end - t_load_rank_start;
        double local_count_time = t_end_count - t_start_count;

        // --- Print Local Results (Each Rank) ---
        std::stringstream ss;
        ss << "[MPI Rank " << std::setw(3) << mpi_rank << "/" << world_size << "]"
           << " file=" << subgraph_filename.substr(subgraph_filename.find_last_of('/') + 1) // Just filename
           << " butterflies=" << local_ans
           << " load+rank=" << std::fixed << std::setprecision(5) << local_load_rank_time << "s"
           << " count=" << std::fixed << std::setprecision(5) << local_count_time << "s"
           << " thr=" << nT_to_use
           << "\n";
        std::cout << ss.str(); // Print the formatted string

        // --- Barrier before Reduction/Final Output ---
        // Ensures all ranks finish printing their local info before rank 0 prints summary
        std::cout.flush(); // Ensure buffer is flushed before barrier
        MPI_Barrier(MPI_COMM_WORLD);

        // --- Phase 3: Reduce Results ---
        uint64_t global_ans = 0;
        MPI_Reduce(&local_ans, &global_ans, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

        double max_count_time = 0;
        MPI_Reduce(&local_count_time, &max_count_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        double max_load_rank_time = 0;
        MPI_Reduce(&local_load_rank_time, &max_load_rank_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


        // --- Phase 4: Output from Rank 0 ---
        if (mpi_rank == 0) {
            double mpi_end_time = MPI_Wtime();
            double total_wall_time = mpi_end_time - mpi_start_time;

            std::cout << "----------------------------------------------------------------------\n"; // Separator
            std::cout << "[mpi-SUMMARY]"
                      << " butterflies=" << global_ans
                      << " load+rank(max)=" << std::fixed << std::setprecision(5) << max_load_rank_time << "s"
                      << " count(max)=" << std::fixed << std::setprecision(5) << max_count_time << "s"
                      << " total_wall=" << std::fixed << std::setprecision(5) << total_wall_time << "s"
                      << " nodes=" << world_size
                      << " thr_per_node=" << P.threads // Use P.threads as it reflects config
                      << " rank=" << static_cast<int>(P.rank)
                      << " agg=" << P.agg
                      << "\n";
             std::cout << "[mpi-SUMMARY] Info: Processed subgraphs from directory: " << P.subgraph_dir << std::endl;
        }

    } catch (const std::exception &e) { /* Handle errors, MPI_Abort */
        std::cerr << "[MPI Rank " << mpi_rank << "] Runtime Error processing " << subgraph_filename << ": " << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1); return 1;
    } catch (...) { /* Handle unknown errors, MPI_Abort */
         std::cerr << "[MPI Rank " << mpi_rank << "] Unknown Runtime Error processing " << subgraph_filename << "." << "\n";
         MPI_Abort(MPI_COMM_WORLD, 1); return 1;
    }

    MPI_Finalize();
    return 0;
}