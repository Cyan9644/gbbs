#pragma once

#include "edge_array.h"
#include "edge_map_data.h"
#include "flags.h"
#include "graph.h"
#include "macros.h"
#include "vertex_subset.h"

namespace gbbs {

// =================== VertexSubset Operators: Mapping =====================
// Applies the the supplied (update, cond) operators on edges incident to the
// input vertex_subset, aggregating results at the neighbors of this vset.
// This is the specialized version of the general nghMap function where the
// output is a plain vertex_subset.
// F is a struct providing the following functions:
//   update : (uintE * uintE * W) -> bool
//   updateAtomic : (uintE * uintE * W) -> bool
//   cond : uintE -> bool
template <class Graph, class VS, class F>
inline vertexSubsetData<gbbs::empty> neighbor_map(Graph& G, VS& vs, F f,
                                                  intT threshold = -1,
                                                  flags fl = 0) {
  return edgeMapData<gbbs::empty>(G, vs, f, threshold, fl);
}

// Generalized version of edgeMap. Takes an input vertex_subset and
// aggregates results at the neighbors of the vset. See above for the
// requirements on F.
template <class Graph,  // underlying graph class
          class Data,   // data associated with vertices in the output
          class VS,     // input vertex_subset type
          class F>      // map function type
inline vertexSubsetData<Data>
nghMapData(Graph& G, VS& vs, F f, intT threshold = -1, flags fl = 0) {
  return edgeMapData<Data>(G, vs, f, threshold, fl);
}

}  // namespace gbbs
