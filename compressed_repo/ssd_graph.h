// ssd_graph: a graph whose edges live on SSD, in EM-graph's 4 KiB-block
// shards (EM-graph/semi_em_graph.h). It is the graph type the GBBS algorithms
// see; the reads themselves happen in ssd_fetch_edges.h.
//
// Only what edgeMap and the algorithms ask of a graph is provided: n, m,
// weight_type, and get_vertex(v).out_degree() / in_degree(). Neighbor lists
// are reached only through fetch_edges.
//
// First version: symmetric, unweighted graphs. In-edges and out-edges are the
// same list, so a pull (edgeMapDense) reads the same shards as a push.
#pragma once

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "bridge.h"
#include "macros.h"

#include "graph/semi_em_graph.h"

namespace gbbs {

template <class Format>
struct ssd_graph {
  using vertex_id = DefaultVertex;  // uint32_t, the same type as uintE
  using sem_graph = SemiEMGraph<vertex_id, Empty, unsigned int, Format>;
  using format = Format;
  using weight_type = gbbs::empty;
  static_assert(std::is_same_v<vertex_id, uintE>,
                "the SSD shards and GBBS must use the same vertex id type");

  // What get_vertex returns: just the degree. There is no neighbor access
  // here on purpose; use fetch_edges.
  struct vertex_view {
    uintE deg;
    uintE out_degree() const { return deg; }
    uintE in_degree() const { return deg; }  // symmetric
  };

  size_t n = 0;
  size_t m = 0;
  double degree_pass_seconds = 0;  // time spent filling `degrees` at load
  // Reader threads for dense (streaming) reads; the remaining workers decode
  // and run the callbacks. EM-graph uses 12 for BFS and 4 for PageRank.
  size_t dense_io_threads = sem_graph::kDenseIoThreads;

  // Opens the shards at `prefix` (metadata on /mnt/ssd0, edge files across
  // the SSDs) and computes every vertex's degree with one streaming pass,
  // because the shards record which blocks a vertex owns, not its degree.
  explicit ssd_graph(const std::string& prefix)
      : sem_(std::make_unique<sem_graph>(prefix)) {
    sem_->OpenFiles(prefix);
    n = sem_->n;
    m = sem_->m;

    timer t;  // starts on construction
    degrees_ = sequence<uintE>(n, 0);
    uintE* deg = degrees_.data();
    sem_->StreamActiveSlots(nullptr, [&](vertex_id u, const vertex_id*,
                                         size_t cnt) {
      // A hub's blocks arrive on different threads, so add atomically.
      if (cnt) __atomic_fetch_add(&deg[u], static_cast<uintE>(cnt),
                                  __ATOMIC_RELAXED);
    });
    degree_pass_seconds = t.next_time();

    size_t total = parlay::reduce(parlay::delayed_seq<size_t>(
        n, [&](size_t i) { return static_cast<size_t>(deg[i]); }));
    if (total != m) {
      std::cerr << "ssd_graph: degrees sum to " << total << " but the shard "
                << "metadata says m = " << m << std::endl;
      abort();
    }
  }

  vertex_view get_vertex(uintE v) const { return {degrees_[v]}; }

  // The underlying EM-graph object, for the SSD fetch_edges backend. The
  // algorithms hold the graph by const reference, but reading from SSD
  // changes the reader's internal state, hence the unique_ptr.
  sem_graph& sem() const { return *sem_; }

 private:
  std::unique_ptr<sem_graph> sem_;
  sequence<uintE> degrees_;
};

}  // namespace gbbs
