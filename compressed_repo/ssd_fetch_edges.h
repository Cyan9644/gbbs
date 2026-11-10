// fetch_edges for graphs on SSD (ssd_graph.h). Same call, same contract as the
// DRAM version in fetch_edges.h; only the way neighbors are read differs.
//
//   frontier is a list  -> read only the 4 KiB blocks that hold those
//                          vertices (sorted by block, many reads in flight).
//   frontier is dense   -> stream the edge files from start to end, skipping
//                          every region that holds no vertex we still need.
//
// Either way the callback gets one piece per (vertex, block): the part of v's
// neighbor list stored in that block. A small vertex fits in one block, so it
// gets one piece. A high-degree vertex spans several blocks, and those blocks
// can be decoded on different threads at once, so its pieces are marked
// `concurrent` (edgeMap then uses updateAtomic for it).
//
// Overload resolution picks these functions over the generic DRAM template
// because ssd_graph<F> is more specific than any Graph, so edge_map_data.h
// and the algorithms need no change.
//
// Ignored options: in_edges (symmetric graphs only, for now) and split (a
// large list is already split, one piece per block).
#pragma once

#include <algorithm>
#include <atomic>
#include <tuple>
#include <type_traits>
#include <vector>

#include "bridge.h"
#include "fetch_edges.h"
#include "macros.h"
#include "ssd_graph.h"

namespace gbbs {

// A piece of v's neighbor list, decoded into a per-thread buffer. Valid only
// until the callback returns. Iterates as (neighbor, weight) tuples like the
// DRAM pieces; the weight is gbbs::empty (unweighted shards).
struct ssd_piece {
  using edge = std::tuple<uintE, gbbs::empty>;
  const uintE* first;
  size_t len;
  bool concurrent;  // other blocks of this vertex may be in a callback now

  struct iterator {
    const uintE* p;
    edge operator*() const { return edge(*p, gbbs::empty()); }
    iterator& operator++() {
      ++p;
      return *this;
    }
    bool operator!=(const iterator& o) const { return p != o.p; }
  };
  iterator begin() const { return {first}; }
  iterator end() const { return {first + len}; }
  size_t size() const { return len; }
};

namespace ssd_fetch_internal {

// Edges a list-frontier block can hold before the per-thread buffer must
// grow. The interval format can expand one block into many more edges; its
// decoders grow the buffer themselves.
constexpr size_t kInitialScratch = 1 << 16;

template <class V>
inline std::vector<V>& thread_scratch() {
  static thread_local std::vector<V> scratch(kInitialScratch);
  return scratch;
}

// The frontier's vertex list, without a copy when it is already a plain list
// of ids (vertexSubset); built from vtx(i) otherwise (vertexSubsetData).
template <class VS>
decltype(auto) vertex_list(VS& frontier) {
  if constexpr (std::is_same_v<std::decay_t<decltype(frontier.s)>,
                               sequence<uintE>>) {
    return (frontier.s);
  } else {
    return sequence<uintE>::from_function(
        frontier.size(), [&](size_t i) { return frontier.vtx(i); });
  }
}

// Dense frontier: stream the files, skipping regions with no active vertex.
template <class F, class VS, class Fn, class Keep>
void fetch_dense(const ssd_graph<F>& G, VS& frontier, Fn& fn, Keep& keep) {
  using V = typename ssd_graph<F>::vertex_id;
  auto& sem = G.sem();
  const size_t n = G.n;

  // skip bit = "v needs no reading": not in the frontier, or keep(v) is
  // already false. One task builds one 64-bit word, so no atomics.
  FrontierBitset skip(n);
  std::atomic<bool> any_skipped{false};
  parallel_for(0, skip.num_words, [&](size_t w) {
    const size_t base = w * 64, end = std::min(base + 64, n);
    uint64_t bits = 0;
    for (size_t v = base; v < end; v++) {
      if (!(frontier.isIn(v) && keep(v))) bits |= 1ULL << (v - base);
    }
    if (bits != 0) any_skipped.store(true, std::memory_order_relaxed);
    skip.data[w] = bits;
  });
  // Nothing to skip (e.g. PageRank): a plain full scan, which also avoids
  // building the per-region bitmap.
  const FrontierBitset* skip_ptr = any_skipped.load() ? &skip : nullptr;

  sem.StreamActiveBlocks(
      skip_ptr,
      [&](unsigned char* block, std::vector<V>& scratch) {
        // keep(v) is asked again here: it can turn false during the round
        // (e.g. BFS claimed v from an earlier block).
        F::decode_to_buffer(
            [&](V v) { return !(skip_ptr && skip.test(v)) && keep(v); },
            block, scratch, [&](V v, const V* edges, size_t len) {
              fn(v, ssd_piece{edges, len, sem.BlocksForVertex(v) > 1});
            });
      },
      G.dense_io_threads);
}

// List frontier: read only the blocks holding the frontier's vertices.
template <class F, class VS, class Fn, class Keep>
void fetch_sparse(const ssd_graph<F>& G, VS& frontier, Fn& fn, Keep& keep) {
  using V = typename ssd_graph<F>::vertex_id;
  using Edge = typename ssd_graph<F>::sem_graph::Edge;
  using BlockId = unsigned int;
  auto& sem = G.sem();

  const auto& vertices = vertex_list(frontier);
  // (block, the frontier vertices in it), sorted by block.
  auto blocks = sem.GetBlockVertices(vertices);
  const size_t num_blocks = blocks.size();
  if (num_blocks == 0) return;

  auto emit = [&](V v, const V* edges, size_t len)
                  __attribute__((always_inline)) {
    if (len == 0 || !keep(v)) return;
    fn(v, ssd_piece{edges, len, sem.BlocksForVertex(v) > 1});
  };

  // Decode one block that just arrived, for the frontier vertices in it.
  auto on_block = [&](BlockId i, void* data) {
    const auto& wanted = blocks[i].second;
    auto* buf = static_cast<unsigned char*>(data);
    auto& scratch = thread_scratch<V>();
    // Does the format decode a whole slot at a time for a vertex list?
    // (The interval format does; the others give one edge at a time.)
    constexpr bool kSlotDecoder = requires(const sequence<V>& vs,
                                           unsigned char* b,
                                           std::vector<V>& s) {
      F::decode_for_to_buffer(vs, b, s, [](V, const V*, size_t) {});
    };
    if constexpr (kSlotDecoder) {
      F::decode_for_to_buffer(wanted, buf, scratch, emit);
    } else {
      // The format hands over one edge at a time, grouped by vertex;
      // collect each vertex's edges into the buffer and emit them together.
      size_t len = 0;
      V cur = 0;
      F::decode_for(wanted, buf, [&](V v, Edge e) {
        if (len > 0 && v != cur) {
          emit(cur, scratch.data(), len);
          len = 0;
        }
        cur = v;
        if (len == scratch.size()) scratch.resize(2 * len);
        scratch[len++] = e.vertex;
      });
      if (len > 0) emit(cur, scratch.data(), len);
    }
  };

  // Each task submits a contiguous range of blocks to its thread's io_uring
  // ring and decodes them as they complete (in any order).
  const size_t num_chunks = std::min(num_blocks, num_workers() * 8);
  const size_t chunk_size = (num_blocks + num_chunks - 1) / num_chunks;
  parallel_for(0, num_chunks, [&](size_t c) {
    const size_t start = c * chunk_size;
    const size_t end = std::min(start + chunk_size, num_blocks);
    if (start >= end) return;
    auto& io = GetThreadIOModule<BlockId>();
    // One captured pointer, so the std::function the ring stores per
    // request needs no heap allocation.
    auto callback = [&on_block](BlockId i, void* data) { on_block(i, data); };
    for (size_t i = start; i < end; i++) {
      auto [fd_index, physical] =
          sem.GetBlockMapper().logical_to_physical(blocks[i].first);
      io.submitRequest(static_cast<BlockId>(i), sem.GetFd(fd_index),
                       physical * GRAPH_BLOCK_SIZE, callback);
    }
    io.Drain();
  }, 1);
}

}  // namespace ssd_fetch_internal

template <class F, class VS, class Fn, class Keep = keep_all>
inline void fetch_edges(const ssd_graph<F>& G, VS& frontier, Fn fn,
                        Keep keep = {}, fetch_opts opts = {}) {
  (void)opts;
  // A frontier with no list form (all_vertices, used by pulls) is always
  // read densely.
  if constexpr (requires { frontier.s; }) {
    if (!frontier.dense()) {
      ssd_fetch_internal::fetch_sparse(G, frontier, fn, keep);
      return;
    }
  }
  ssd_fetch_internal::fetch_dense(G, frontier, fn, keep);
}

// Non-const graphs too: without this, the generic DRAM template (Graph&)
// would be an exact match and win over the const& overload above.
template <class F, class VS, class Fn, class Keep = keep_all>
inline void fetch_edges(ssd_graph<F>& G, VS& frontier, Fn fn, Keep keep = {},
                        fetch_opts opts = {}) {
  fetch_edges(static_cast<const ssd_graph<F>&>(G), frontier, std::move(fn),
              std::move(keep), opts);
}

}  // namespace gbbs
