// fetch_edges: the single neighbor-access primitive of the GBBS rewrite.
// See gbbs/neighbor_fetch_design.tex.
//
//   fetch_edges(G, frontier, fn, keep, opts)
//
// For every vertex v in `frontier` with keep(v), calls fn(v, piece) for the
// pieces of v's neighbor list. A piece is a sorted, contiguous run of
// (neighbor, weight) pairs, valid only until fn returns. The backend decides
// how to read: this file is the DRAM backend, which loops in parallel; an SSD
// backend reads blocks (list frontier) or streams files (dense frontier).
//
// Contract for callers:
//   * No assumption about the order in which vertices or pieces arrive.
//   * piece.concurrent is true when other pieces of the same v may be running
//     at the same moment, so updates to v-owned state must be atomic. It is
//     false when v's pieces arrive one after another on one thread (always
//     the case when opts.split is false).
//   * keep(v) is checked before each piece; once it is false, the rest of v
//     is skipped (this is how a pull stops early).
//   * A piece is a range of (neighbor, weight) tuples. Iterate it at most
//     once: a compressed piece decodes while you iterate, so breaking out of
//     the loop also stops the decoding.
#pragma once

#include <algorithm>
#include <tuple>

#include "bridge.h"
#include "compressed_vertex.h"
#include "macros.h"
#include "vertex.h"

namespace gbbs {

// A run of v's neighbors. Points into graph memory or a scratch buffer.
template <class W>
struct edge_piece {
  using edge = std::tuple<uintE, W>;
  const edge* first;
  size_t len;
  bool concurrent;  // other pieces of this vertex may run at the same time
  const edge* begin() const { return first; }
  const edge* end() const { return first + len; }
  size_t size() const { return len; }
  const edge& operator[](size_t i) const { return first[i]; }
};

// A piece that decodes as it is iterated. `Cursor` yields the piece's edges
// via first() and next() (see bytepd_amortized::block_cursor). Single pass:
// begin() hands the decoder state to the iterator by value, so it can live in
// registers for the whole loop.
template <class W, class Cursor>
struct decoding_piece {
  using edge = std::tuple<uintE, W>;
  Cursor cursor;
  size_t len;
  bool concurrent;

  struct sentinel {};
  struct iterator {
    Cursor cursor;
    size_t i;
    size_t len;
    edge e;
    const edge& operator*() const { return e; }
    iterator& operator++() {
      if (++i < len) e = cursor.next();
      return *this;
    }
    bool operator!=(sentinel) const { return i != len; }
  };
  iterator begin() const {
    iterator it{cursor, 0, len, edge()};
    if (len > 0) it.e = it.cursor.first();
    return it;
  }
  sentinel end() const { return {}; }
  size_t size() const { return len; }
};

struct fetch_opts {
  bool in_edges = false;  // read in-neighbors instead of out-neighbors
  bool split = true;      // allow a large list to be split across threads
};

struct keep_all {
  bool operator()(uintE) const { return true; }
};

// The set {0, ..., n-1}. Used by pulls, where `keep` does the filtering.
struct all_vertices {
  size_t n;
  bool dense() const { return true; }
  bool isIn(uintE) const { return true; }
  size_t size() const { return n; }
  uintE vtx(size_t i) const { return i; }
};

namespace fetch_internal {

// Edges per piece when an uncompressed list is split. (Compressed lists are
// split at their encoding blocks.)
constexpr size_t kPieceSize = 1024;
// Vertices per task when scanning a dense frontier.
constexpr size_t kDenseGranularity = 2048;

// The split paths below are kept out of line (noinline): they are rare, and
// inlining a nested parallel_for into the per-vertex loop slows the common
// path down.

// Uncompressed, large list: slices of kPieceSize run in parallel.
template <class W, class Fn, class Keep>
__attribute__((noinline)) void split_pieces(uintE v, const std::tuple<uintE, W>* e,
                                            size_t d, Fn& fn, Keep& keep) {
  parallel_for(0, parlay::num_blocks(d, kPieceSize), [&](size_t b) {
    if (!keep(v)) return;
    size_t s = b * kPieceSize, t = std::min(s + kPieceSize, d);
    fn(v, edge_piece<W>{e + s, t - s, true});
  }, 1);
}

// Uncompressed: pieces are slices of the edge array, no copying.
template <class W, class Fn, class Keep>
__attribute__((always_inline)) inline void for_each_piece(
    uncompressed_neighbors<W> nghs, Fn& fn, Keep& keep, bool split) {
  const uintE v = nghs.id;
  const size_t d = nghs.degree;
  if (!split || d <= kPieceSize) [[likely]] {
    fn(v, edge_piece<W>{nghs.neighbors, d, false});
  } else {
    split_pieces<W>(v, nghs.neighbors, d, fn, keep);
  }
}

// Compressed. Each encoding block is one piece that decodes while the caller
// iterates it, so decoding is fused with the caller's loop and a `break`
// stops it. Blocks of one list run one after another, or in parallel with
// `split`. A codec without blocks (no block_cursor) is one whole-list piece.
template <class C, class W>
concept has_block_cursor = requires(uchar* p, uintE x) {
  C::template block_cursor<W>(p, x, x, x);
};

// Adapts a codec's whole-list iter<W> (cur/next) to the cursor interface.
template <class It>
struct list_cursor {
  It it;
  auto first() { return it.cur(); }
  auto next() { return it.next(); }
};

template <class W, class C, class Fn>
__attribute__((always_inline)) inline void block_piece(uintE v, uchar* nghs,
                                                      uintE d, uintE b,
                                                      bool concurrent, Fn& fn) {
  auto cursor = C::template block_cursor<W>(nghs, v, d, b);
  const size_t len = cursor.size();
  if (len > 0) fn(v, decoding_piece<W, decltype(cursor)>{cursor, len, concurrent});
}

// Lists with more than one encoding block (degree above PARALLEL_DEGREE).
template <class W, class C, class Fn, class Keep>
__attribute__((noinline)) void multi_block_pieces(uintE v, uchar* nghs,
                                                  uintE d, size_t num_blocks,
                                                  bool split, Fn& fn,
                                                  Keep& keep) {
  if (split) {
    parallel_for(0, num_blocks, [&](size_t b) {
      if (keep(v)) block_piece<W, C>(v, nghs, d, b, /*concurrent=*/true, fn);
    }, 1);
  } else {
    for (size_t b = 0; b < num_blocks && keep(v); b++) {
      block_piece<W, C>(v, nghs, d, b, /*concurrent=*/false, fn);
    }
  }
}

template <class W, class C, class Fn, class Keep>
__attribute__((always_inline)) inline void for_each_piece(
    compressed_neighbors<W, C> nghs, Fn& fn, Keep& keep, bool split) {
  const uintE v = nghs.id;
  const uintE d = nghs.degree;
  if (d == 0) return;
  if constexpr (has_block_cursor<C, W>) {
    const size_t num_blocks = C::get_num_blocks(nghs.neighbors, d);
    if (num_blocks == 1) [[likely]] {
      block_piece<W, C>(v, nghs.neighbors, d, 0, /*concurrent=*/false, fn);
    } else {
      multi_block_pieces<W, C>(v, nghs.neighbors, d, num_blocks, split, fn,
                               keep);
    }
  } else {
    using It = decltype(C::template iter<W>(nghs.neighbors, d, v));
    using Cursor = list_cursor<It>;
    fn(v, decoding_piece<W, Cursor>{
              Cursor{C::template iter<W>(nghs.neighbors, d, v)}, d, false});
  }
}

}  // namespace fetch_internal

template <class Graph, class VS, class Fn, class Keep = keep_all>
inline void fetch_edges(Graph& G, VS& frontier, Fn fn, Keep keep = {},
                        fetch_opts opts = {}) {
  // The per-vertex path (visit -> for_each_piece -> fn) is forced inline;
  // edgeMap's callbacks are marked the same way.
  auto visit = [&](uintE v) __attribute__((always_inline)) {
    if (!keep(v)) return;
    auto vertex = G.get_vertex(v);
    fetch_internal::for_each_piece(
        opts.in_edges ? vertex.in_neighbors() : vertex.out_neighbors(), fn,
        keep, opts.split);
  };
  if (frontier.dense()) {
    parallel_for(0, G.n, [&](size_t v) {
      if (frontier.isIn(v)) visit(v);
    }, fetch_internal::kDenseGranularity);
  } else {
    parallel_for(0, frontier.size(), [&](size_t i) {
      visit(frontier.vtx(i));
    }, 1);
  }
}

}  // namespace gbbs
