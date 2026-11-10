// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// edgeMap, written once on top of fetch_edges (fetch_edges.h). Nothing here
// touches a vertex's neighbor list directly, so swapping the backend under
// fetch_edges changes no code in this file or in the algorithms.
//
// An edgeMap struct F provides:
//   cond(d)                  should d still be considered?
//   update(s, d, w)          non-atomic update, s -> d
//   updateAtomic(s, d, w)    atomic update, s -> d
// update/updateAtomic return bool (or std::optional<Data>): "add d to the
// output frontier".
#pragma once

#include <optional>
#include <tuple>
#include <type_traits>
#include <vector>

#include "bridge.h"
#include "fetch_edges.h"
#include "flags.h"
#include "vertex_subset.h"

namespace gbbs {

namespace edge_map_internal {

// Was update/updateAtomic's result "add to output"? (bool or optional<Data>)
template <class M>
inline bool wants_output(const M& m) {
  if constexpr (std::is_same_v<M, bool>) {
    return m;
  } else {
    return m.has_value();
  }
}

// Entry for the dense output array (bool, or (bool, Data)).
template <class Data, class M>
inline auto dense_entry(const M& m) {
  if constexpr (std::is_same_v<Data, gbbs::empty>) {
    return true;
  } else {
    return std::make_tuple(true, *m);
  }
}

// Entry for the sparse output list (uintE, or (uintE, Data)).
template <class Data, class M>
inline auto sparse_entry(uintE d, const M& m) {
  if constexpr (std::is_same_v<Data, gbbs::empty>) {
    return d;
  } else {
    return std::make_tuple(d, *m);
  }
}

}  // namespace edge_map_internal

// Pull. Every v with cond(v) scans its in-neighbors (out-neighbors if
// `in_edges` is set) for members of vs, stopping once cond(v) turns false.
// update is used when v's pieces arrive one at a time, and updateAtomic when
// the piece says other pieces of v may be running concurrently. Large lists
// are split across threads only with `dense_parallel`.
template <class Data, class Graph, class VS, class F>
inline vertexSubsetData<Data> edgeMapDense(Graph& G, VS& vs, F& f,
                                           const flags fl) {
  using D = typename vertexSubsetData<Data>::D;
  const size_t n = G.n;
  const bool output = should_output(fl);
  sequence<D> next = output ? sequence<D>(n, D()) : sequence<D>();
  fetch_opts opts{/*in_edges=*/!(fl & in_edges),
                  /*split=*/static_cast<bool>(fl & dense_parallel)};

  auto scan = [&](uintE v, const auto& nghs, auto atomic)
                  __attribute__((always_inline)) {
    for (const auto& [u, w] : nghs) {
      if (vs.isIn(u)) {
        auto m = [&] {
          if constexpr (decltype(atomic)::value) {
            return f.updateAtomic(u, v, w);
          } else {
            return f.update(u, v, w);
          }
        }();
        if (output && edge_map_internal::wants_output(m)) {
          next[v] = edge_map_internal::dense_entry<Data>(m);
        }
      }
      if (!f.cond(v)) break;
    }
  };
  all_vertices everyone{n};
  fetch_edges(G, everyone, [&](uintE v, const auto& nghs)
                                __attribute__((always_inline)) {
    if (nghs.concurrent) {
      scan(v, nghs, std::true_type{});
    } else {
      scan(v, nghs, std::false_type{});
    }
  }, [&](uintE v) __attribute__((always_inline)) { return f.cond(v); },
  opts);
  return output ? vertexSubsetData<Data>(n, std::move(next))
                : vertexSubsetData<Data>(n);
}

// Push from a dense frontier: every member s calls updateAtomic(s, d) on its
// out-neighbors d with cond(d). Output is written into a dense array.
template <class Data, class Graph, class VS, class F>
inline vertexSubsetData<Data> edgeMapDenseForward(Graph& G, VS& vs, F& f,
                                                  const flags fl) {
  using D = typename vertexSubsetData<Data>::D;
  const size_t n = G.n;
  const bool output = should_output(fl);
  sequence<D> next = output ? sequence<D>(n, D()) : sequence<D>();
  fetch_opts opts{/*in_edges=*/static_cast<bool>(fl & in_edges),
                  /*split=*/true};

  fetch_edges(G, vs, [&](uintE s, const auto& nghs) {
    for (const auto& [d, w] : nghs) {
      if (f.cond(d)) {
        auto m = f.updateAtomic(s, d, w);
        if (output && edge_map_internal::wants_output(m)) {
          next[d] = edge_map_internal::dense_entry<Data>(m);
        }
      }
    }
  }, keep_all{}, opts);
  return output ? vertexSubsetData<Data>(n, std::move(next))
                : vertexSubsetData<Data>(n);
}

// Push from a sparse frontier. Each worker appends the vertices it adds to its
// own buffer, and the buffers are concatenated at the end.
template <class Data, class Graph, class VS, class F>
inline vertexSubsetData<Data> edgeMapSparse(Graph& G, VS& vs, F& f,
                                            const flags fl) {
  using S = typename vertexSubsetData<Data>::S;
  const size_t n = G.n;
  fetch_opts opts{/*in_edges=*/static_cast<bool>(fl & in_edges),
                  /*split=*/true};

  if (!should_output(fl)) {
    fetch_edges(G, vs, [&](uintE s, const auto& nghs) {
      for (const auto& [d, w] : nghs) {
        if (f.cond(d)) f.updateAtomic(s, d, w);
      }
    }, keep_all{}, opts);
    return vertexSubsetData<Data>(n);
  }

  struct alignas(64) out_buffer {
    std::vector<S> items;
  };
  std::vector<out_buffer> buffers(num_workers());
  fetch_edges(G, vs, [&](uintE s, const auto& nghs) {
    auto& out = buffers[worker_id()].items;
    for (const auto& [d, w] : nghs) {
      if (f.cond(d)) {
        auto m = f.updateAtomic(s, d, w);
        if (edge_map_internal::wants_output(m)) {
          out.push_back(edge_map_internal::sparse_entry<Data>(d, m));
        }
      }
    }
  }, keep_all{}, opts);

  auto offsets = sequence<size_t>::from_function(
      buffers.size(), [&](size_t i) { return buffers[i].items.size(); });
  size_t total = parlay::scan_inplace(offsets);
  if (total == 0) return vertexSubsetData<Data>(n);
  auto out = sequence<S>::uninitialized(total);
  parallel_for(0, buffers.size(), [&](size_t i) {
    std::copy(buffers[i].items.begin(), buffers[i].items.end(),
              out.begin() + offsets[i]);
  }, 1);
  return vertexSubsetData<Data>(n, std::move(out));
}

// Decides on sparse or dense based on the number of nonzeros in the active
// vertices.
template <
    class Data /* data associated with vertices in the output vertex_subset */,
    class Graph /* graph type */, class VS /* vertex_subset type */,
    class F /* edgeMap struct */>
inline vertexSubsetData<Data> edgeMapData(Graph& GA, VS& vs, F f,
                                          intT threshold = -1,
                                          const flags& fl = 0) {
  size_t numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
  size_t dense_threshold = threshold;
  if (threshold == -1) dense_threshold = numEdges / 20;
  if (vs.size() == 0) return vertexSubsetData<Data>(numVertices);

  auto run_dense = [&] {
    vs.toDense();
    return (fl & dense_forward) ? edgeMapDenseForward<Data>(GA, vs, f, fl)
                                : edgeMapDense<Data>(GA, vs, f, fl);
  };

  if ((fl & dense_only) || (vs.isDense && vs.size() > numVertices / 10)) {
    return run_dense();
  }

  size_t out_degrees = 0;
  if (vs.out_degrees_set()) {
    out_degrees = vs.get_out_degrees();
  } else {
    vs.toSparse();
    auto degree_f = [&](size_t i) {
      return (fl & in_edges) ? GA.get_vertex(vs.vtx(i)).in_degree()
                             : GA.get_vertex(vs.vtx(i)).out_degree();
    };
    auto degree_im = parlay::delayed_seq<size_t>(vs.size(), degree_f);
    out_degrees = parlay::reduce(degree_im);
    vs.set_out_degrees(out_degrees);
  }

  if (out_degrees == 0) return vertexSubsetData<Data>(numVertices);
  if (m + out_degrees > dense_threshold && !(fl & no_dense)) {
    return run_dense();
  }
  return edgeMapSparse<Data>(GA, vs, f, fl);
}

// Regular edgeMap, where no extra data is stored per vertex.
template <class Graph /* graph type */, class VS /* vertex_subset type */,
          class F /* edgeMap struct */>
inline vertexSubset edgeMap(Graph& GA, VS& vs, F f, intT threshold = -1,
                            const flags& fl = 0) {
  return edgeMapData<gbbs::empty>(GA, vs, f, threshold, fl);
}

// cond function that always returns true
inline bool cond_true(intT d) { return 1; }

}  // namespace gbbs
