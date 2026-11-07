# GBBS core (`core-slim` branch)

This branch cuts [GBBS](https://github.com/ParAlg/gbbs) down to its core primitives. The goal is to integrate it with EM-graph and PLAID. Following `note_from_prof_dhulipala.txt`, we keep only the neighbor-access API and the Ligra-style traversal layer. Anything else should come back only after careful thought. Upstream GBBS is on `master` (`git diff master` shows every cut).

## What's here

| Layer | Files | API |
| --- | --- | --- |
| Parallel runtime | `gbbs/bridge.h`, `macros.h`, `flags.h` | parlaylib wrappers, `uintE`/`uintT`, `gbbs::empty`, edgeMap flags |
| Neighbors | `gbbs/vertex.h`, `compressed_vertex.h`, `encodings/`, `intersect.h` | `G.get_vertex(v).out_neighbors()` → `get_degree()`, `get_neighbor(i)`, `map(f)`, `reduce(m, monoid)`, `count(f)`, `intersect(...)` |
| Graphs | `gbbs/graph.h`, `edge_array.h` | `symmetric_graph` / `asymmetric_graph` (and `_ptr_` variants), `from_edges`, `mapEdges`, `reduceEdges`, `edges()` |
| Frontiers | `gbbs/vertex_subset.h` | `vertexSubset`, `vertexSubsetData<D>` (sparse ↔ dense) |
| Traversal | `gbbs/edge_map_data.h`, `edge_map_blocked.h`, `edge_map_utils.h`, `interface.h` | `edgeMap(G, vs, F, threshold, flags)`, `edgeMapData<D>`, `nghMapData` |
| I/O | `gbbs/graph_io.h/.cc`, `io.h/.cc`, `benchmark.h` | `read_unweighted_symmetric_graph(...)` etc. (adjacency, binary, and compressed bytePDA), `generate_*_main` macros |
| Helpers | `gbbs/helpers/` | `parse_command_line`, `assert`, `directed_edge`/`undirected_edge` (tests) |

The core library depends only on parlaylib. googletest is used only by the tests.

Benchmarks kept as smoke tests:
- `benchmarks/BFS/NonDeterministicBFS`: sparse and dense edgeMap.
- `benchmarks/PageRank`: power iteration over dense edgeMap. Personalized sources are supported.

Utilities kept: `//utils:compressor` (adjacency → bytePDA/binary), `//utils:snap_converter`, and `//utils:converter`.

## What was removed

Every other benchmark. Julienne bucketing (`bucket.h`, `julienne.h`). `edge_map_reduce.h` and the histogram, sequential-hash, and counting-sort helpers. Graph mutation and contraction (`graph_mutation.h`, `contract.h`, sparse and resizable tables). `speculative_for`. Atomic counters. abseil-based progress reporting and status macros. The PageRank delta and edgeMapReduce variants. The Make build, `scripts/`, the graph generators, `random_reorder`, `to_edge_list`, and the abseil, fuzztest, and google_benchmark dependencies.

## Next step

Define a minimal `get_neighbors()` (`fetch_edges`) interface that both the in-memory graphs here and EM-graph's SSD-backed `SemiEMGraph` can implement. Then rebuild edgeMap on top of it, so that switching between DRAM and SSD costs nothing compared with the specialized code.

## Building

```sh
bazel build //...
bazel test //...
bazel run //benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -src 0 $PWD/inputs/rMatGraph_J_5_100
bazel run //benchmarks/PageRank:PageRank_main -- -s $PWD/inputs/rMatGraph_J_5_100
```

Common flags: `-s` means the graph is symmetric, `-c` means compressed (bytePDA), `-b` means binary CSR, `-m` means mmap the file, and `-rounds N` sets the number of rounds. Build options are in `.bazelrc` (`--config=serial|openmp|cilk|asan`).

To make a compressed input:

```sh
bazel run //utils:compressor -- -s -o $PWD/inputs/rMatGraph_J_5_100.bytepda $PWD/inputs/rMatGraph_J_5_100
```

## Input format

Use the PBBS/Ligra adjacency format, which is plain text:

```
AdjacencyGraph
<n>
<m>
<o0> ... <o(n-1)>      # offsets into the edge array, one per vertex
<e0> ... <e(m-1)>      # edge targets, grouped by source
```

Weighted graphs start with `WeightedAdjacencyGraph` and put the m weights after the targets. To convert a SNAP edge list, run `bazel run //utils:snap_converter -- -s -i <in.txt> -o <out.adj>`.

Citation for the original suite: Dhulipala, Blelloch, Shun. *Theoretically Efficient Parallel Graph Algorithms Can Be Fast and Scalable.* SPAA 2018.
