// Runs the unchanged GBBS BFS and PageRank on a graph stored on SSD, and
// checks the SSD backend against the same graph loaded into DRAM.
//
//   ssd_bench bfs      <prefix> [-src N] [-rounds R] [-io_threads T]
//   ssd_bench pagerank <prefix> [-iters K] [-rounds R] [-io_threads T]
//   ssd_bench check    <prefix> <graph.adj> [-src N] [-iters K]
//
// <prefix> is an EM-graph shard prefix (`ls /mnt/ssd0/*_metadata0`; the name
// before `_metadata0`), written by `graph_main shard`. Only symmetric,
// unweighted shards are supported. The block format is read from the shard's
// metadata, as graph_main does.
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "gbbs/edge_map_data.h"
#include "gbbs/fetch_edges.h"
#include "gbbs/graph_io.h"
#include "gbbs/ssd_fetch_edges.h"
#include "gbbs/ssd_graph.h"
#include "gbbs/vertex_subset.h"

#include "gbbs/ssd/apps/BFS.h"
#include "gbbs/ssd/apps/PageRank.h"

#include "graph/block_format/metadata_header.h"
#include "graph/block_format/uncompressed_format.h"
#include "graph/block_format/varint_delta_format.h"
#include "graph/block_format/varint_delta_interval_format.h"

namespace {

using gbbs::uintE;

// ---------------------------------------------------------------- arguments

struct Args {
  std::string command, prefix, adj;
  uintE src = 0;
  size_t rounds = 3;
  size_t iters = 10;
  size_t io_threads = 0;  // 0 = per-command default (12 BFS, 4 PageRank)
};

Args ParseArgs(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: ssd_bench bfs|pagerank|check <prefix> [graph.adj] "
                 "[-src N] [-rounds R] [-iters K] [-io_threads T]\n";
    std::exit(1);
  }
  Args a;
  a.command = argv[1];
  a.prefix = argv[2];
  int i = 3;
  if (a.command == "check") {
    if (argc < 4) {
      std::cerr << "check needs the .adj file of the same graph\n";
      std::exit(1);
    }
    a.adj = argv[3];
    i = 4;
  }
  for (; i + 1 < argc; i += 2) {
    std::string flag = argv[i];
    size_t value = std::stoull(argv[i + 1]);
    if (flag == "-src") a.src = static_cast<uintE>(value);
    else if (flag == "-rounds") a.rounds = value;
    else if (flag == "-iters") a.iters = value;
    else if (flag == "-io_threads") a.io_threads = value;
    else {
      std::cerr << "unknown flag " << flag << "\n";
      std::exit(1);
    }
  }
  return a;
}

// ---------------------------------------------------------------- benchmarks

template <class G>
void PrintHeader(const char* app, const Args& a, const G& g) {
  std::cout << "### Application: " << app << " (SSD backend)\n"
            << "### Graph: " << a.prefix << "\n"
            << "### Threads: " << gbbs::num_workers()
            << "  (reader threads for dense rounds: " << g.dense_io_threads
            << ")\n"
            << "### n: " << g.n << "  m: " << g.m << "\n"
            << "### Degree pass at load: " << g.degree_pass_seconds << " s\n";
}

template <class Format>
void RunBfs(const Args& a) {
  gbbs::ssd_graph<Format> g(a.prefix);
  g.dense_io_threads = a.io_threads ? a.io_threads : 12;
  PrintHeader("BFS", a, g);
  for (size_t r = 0; r < a.rounds; r++) {
    gbbs::timer t;
    auto parents = gbbs::BFS(g, a.src);
    std::cout << "### Running Time: " << t.next_time() << " s\n";
  }
}

template <class Format>
void RunPageRank(const Args& a) {
  gbbs::ssd_graph<Format> g(a.prefix);
  g.dense_io_threads = a.io_threads ? a.io_threads : 4;
  PrintHeader("PageRank", a, g);
  std::cout << "### Iterations per run: " << a.iters << " (eps = 0)\n";
  for (size_t r = 0; r < a.rounds; r++) {
    gbbs::timer t;
    auto pr = gbbs::PageRank_edgeMap(g, /*eps=*/0.0, /*sources=*/{},
                                     /*damping_factor=*/0.85, a.iters);
    double s = t.next_time();
    std::cout << "### Running Time: " << s << " s  (" << s / a.iters
              << " s per iteration)\n";
  }
}

// ---------------------------------------------------------------- check

bool g_ok = true;
void Report(const std::string& what, bool ok) {
  std::cout << (ok ? "  OK    " : "  FAIL  ") << what << "\n";
  g_ok = g_ok && ok;
}

// For every vertex: how many neighbors fetch_edges delivered, their sum, and
// their sum of squares. Equal digests on SSD and DRAM mean (with very high
// probability) that the same neighbor lists were delivered.
struct Digest {
  gbbs::sequence<uint64_t> count, sum, sum_sq;
};

template <class G, class VS, class Keep = gbbs::keep_all>
Digest NeighborDigest(G& g, VS& frontier, Keep keep = {}) {
  Digest d{gbbs::sequence<uint64_t>(g.n, 0), gbbs::sequence<uint64_t>(g.n, 0),
           gbbs::sequence<uint64_t>(g.n, 0)};
  gbbs::fetch_edges(g, frontier, [&](uintE v, const auto& piece) {
    uint64_t c = 0, s = 0, q = 0;
    for (const auto& [u, w] : piece) {
      c++;
      s += u;
      q += static_cast<uint64_t>(u) * u;
    }
    __atomic_fetch_add(&d.count[v], c, __ATOMIC_RELAXED);
    __atomic_fetch_add(&d.sum[v], s, __ATOMIC_RELAXED);
    __atomic_fetch_add(&d.sum_sq[v], q, __ATOMIC_RELAXED);
  }, keep);
  return d;
}

bool SameDigest(const Digest& x, const Digest& y) {
  size_t bad = parlay::count_if(parlay::iota<size_t>(x.count.size()),
                                [&](size_t v) {
                                  return x.count[v] != y.count[v] ||
                                         x.sum[v] != y.sum[v] ||
                                         x.sum_sq[v] != y.sum_sq[v];
                                });
  if (bad) std::cout << "        " << bad << " vertices differ\n";
  return bad == 0;
}

// Level of every vertex from `src`, by a plain level-by-level BFS over the
// DRAM graph (the reference answer).
template <class G>
std::vector<uint32_t> ReferenceLevels(G& g, uintE src) {
  constexpr uint32_t kInf = std::numeric_limits<uint32_t>::max();
  std::vector<uint32_t> level(g.n, kInf);
  std::vector<uintE> current = {src}, next;
  level[src] = 0;
  for (uint32_t l = 0; !current.empty(); l++) {
    next.clear();
    gbbs::vertexSubset frontier(
        g.n, gbbs::sequence<uintE>(current.begin(), current.end()));
    // split=false: each vertex's list is read on one thread, so the
    // non-atomic writes below are safe per vertex; vertices still run in
    // parallel, so claim neighbors with a compare-and-swap.
    std::vector<std::vector<uintE>> found(gbbs::num_workers());
    gbbs::fetch_edges(g, frontier, [&](uintE, const auto& piece) {
      auto& out = found[gbbs::worker_id()];
      for (const auto& [u, w] : piece) {
        if (level[u] == kInf &&
            gbbs::atomic_compare_and_swap(&level[u], kInf, l + 1)) {
          out.push_back(u);
        }
      }
    }, gbbs::keep_all{}, gbbs::fetch_opts{false, false});
    for (auto& f : found) next.insert(next.end(), f.begin(), f.end());
    std::swap(current, next);
  }
  return level;
}

template <class Format>
void RunCheck(const Args& a) {
  std::cout << "Loading " << a.adj << " into DRAM\n";
  auto dram = gbbs::gbbs_io::read_unweighted_symmetric_graph(
      a.adj.c_str(), /*mmap=*/false, /*binary=*/false);
  std::cout << "Opening SSD shards " << a.prefix << "\n";
  gbbs::ssd_graph<Format> ssd(a.prefix);
  const size_t n = ssd.n;

  std::cout << "Graph size and degrees\n";
  Report("n matches", dram.n == ssd.n);
  Report("m matches", dram.m == ssd.m);
  if (dram.n != ssd.n) std::exit(1);
  Report("every degree matches",
         parlay::count_if(parlay::iota<size_t>(n), [&](size_t v) {
           return dram.get_vertex(v).out_degree() !=
                  ssd.get_vertex(v).out_degree();
         }) == 0);

  std::cout << "fetch_edges delivers the same neighbors\n";
  {
    gbbs::all_vertices everyone{n};
    Report("dense frontier, every vertex",
           SameDigest(NeighborDigest(dram, everyone),
                      NeighborDigest(ssd, everyone)));

    auto keep = [](uintE v) { return v % 3 != 0; };
    Report("dense frontier, keep(v) skips every third vertex",
           SameDigest(NeighborDigest(dram, everyone, keep),
                      NeighborDigest(ssd, everyone, keep)));

    auto every100 = parlay::filter(parlay::iota<uintE>(n),
                                   [](uintE v) { return v % 100 == 7; });
    gbbs::vertexSubset list_a(n, gbbs::sequence<uintE>(every100));
    gbbs::vertexSubset list_b(n, gbbs::sequence<uintE>(every100));
    Report("list frontier, 1% of vertices",
           SameDigest(NeighborDigest(dram, list_a),
                      NeighborDigest(ssd, list_b)));

    // The highest-degree vertices: their lists span many blocks on SSD, so
    // their pieces arrive concurrently.
    auto ids = parlay::tabulate(n, [](size_t i) { return static_cast<uintE>(i); });
    auto by_degree = parlay::sort(ids, [&](uintE x, uintE y) {
      return ssd.get_vertex(x).out_degree() > ssd.get_vertex(y).out_degree();
    });
    auto top = gbbs::sequence<uintE>(by_degree.begin(),
                                     by_degree.begin() + std::min<size_t>(n, 64));
    std::sort(top.begin(), top.end());
    size_t multi_block = parlay::count_if(top, [&](uintE v) {
      return ssd.sem().BlocksForVertex(v) > 1;
    });
    std::cout << "        (" << multi_block
              << " of the 64 top-degree vertices span several blocks)\n";
    gbbs::vertexSubset hubs_a(n, gbbs::sequence<uintE>(top));
    gbbs::vertexSubset hubs_b(n, gbbs::sequence<uintE>(top));
    Report("list frontier, 64 highest-degree vertices",
           SameDigest(NeighborDigest(dram, hubs_a),
                      NeighborDigest(ssd, hubs_b)));
  }

  std::cout << "BFS from " << a.src << " on SSD\n";
  {
    auto level = ReferenceLevels(dram, a.src);
    auto parents = gbbs::BFS(ssd, a.src);
    constexpr uint32_t kInf = std::numeric_limits<uint32_t>::max();
    size_t bad_reach = 0, bad_parent = 0;
    for (size_t v = 0; v < n; v++) {
      bool reached = parents[v] != UINT_E_MAX;
      if (reached != (level[v] != kInf)) bad_reach++;
      else if (reached && v != a.src &&
               level[parents[v]] + 1 != level[v]) bad_parent++;
    }
    if (bad_reach) std::cout << "        " << bad_reach << " vertices\n";
    Report("same reachable set as the DRAM reference", bad_reach == 0);
    if (bad_parent) std::cout << "        " << bad_parent << " vertices\n";
    Report("every parent is exactly one level closer to the source",
           bad_parent == 0);
  }

  std::cout << "PageRank, " << a.iters << " iterations\n";
  {
    auto p_dram = gbbs::PageRank_edgeMap(dram, 0.0, {}, 0.85, a.iters);
    auto p_ssd = gbbs::PageRank_edgeMap(ssd, 0.0, {}, 0.85, a.iters);
    double max_diff = parlay::reduce(
        parlay::delayed_seq<double>(n, [&](size_t v) {
          return std::fabs(p_dram[v] - p_ssd[v]);
        }),
        parlay::maximum<double>());
    std::cout << "        largest difference: " << max_diff << "\n";
    Report("SSD and DRAM PageRank agree to 1e-9", max_diff < 1e-9);
  }

  std::cout << (g_ok ? "ALL OK" : "FAILED") << "\n";
  if (!g_ok) std::exit(1);
}

// ---------------------------------------------------------------- dispatch

template <class Format>
void Run(const Args& a) {
  if (a.command == "bfs") RunBfs<Format>(a);
  else if (a.command == "pagerank") RunPageRank<Format>(a);
  else if (a.command == "check") RunCheck<Format>(a);
  else {
    std::cerr << "unknown command " << a.command << "\n";
    std::exit(1);
  }
}

}  // namespace

int main(int argc, char** argv) {
  Args a = ParseArgs(argc, argv);
  using V = DefaultVertex;
  using Edge = Graph<V>::Edge;
  switch (ReadFormatIdFromMetadata(a.prefix)) {
    case BlockFormatId::kUncompressed:
      Run<UncompressedFormat<V, Edge>>(a);
      break;
    case BlockFormatId::kVarintDelta:
      Run<VarintDeltaFormat<V, Edge>>(a);
      break;
    case BlockFormatId::kVarintDeltaInterval:
      Run<VarintDeltaIntervalFormat<V, Edge>>(a);
      break;
    default:
      std::cerr << "unsupported block format\n";
      return 1;
  }
  return 0;
}
