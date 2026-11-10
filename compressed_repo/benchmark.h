// Utilities for creating main files that read graphs.
#pragma once

#include "assert.h"
#include "graph_io.h"

#define run_app(G, APP, mutates, rounds)    \
  double total_time = 0.0;                  \
  for (size_t r = 0; r < rounds; r++) {     \
    if (mutates) {                          \
      auto G_copy = G;                      \
      total_time += APP(G_copy, P);         \
    } else {                                \
      total_time += APP(G, P);              \
    }                                       \
  }                                         \
  auto time_per_iter = total_time / rounds; \
  std::cout << "# time per iter: " << time_per_iter << "\n";

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * either symmetric or asymmetric graph inputs */
#define generate_main(APP, mutates)                                            \
  int main(int argc, char* argv[]) {                                           \
    gbbs::commandLine P(argc, argv, " [-s] <inFile>");                         \
    char* iFile = P.getArgument(0);                                            \
    bool symmetric = P.getOptionValue("-s");                                   \
    bool compressed = P.getOptionValue("-c");                                  \
    bool binary = P.getOptionValue("-b");                                      \
    bool mmap = P.getOptionValue("-m");                                        \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                        \
    if (compressed) {                                                          \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(  \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>( \
            iFile, mmap);                                                      \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    } else {                                                                   \
      if (symmetric) {                                                         \
        auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(iFile, mmap,   \
                                                                binary);       \
        run_app(G, APP, mutates, rounds)                                       \
      } else {                                                                 \
        auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(iFile, mmap,  \
                                                                 binary);      \
        run_app(G, APP, mutates, rounds)                                       \
      }                                                                        \
    }                                                                          \
  }
