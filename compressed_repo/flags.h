#pragma once

#include <cstdint>

namespace gbbs {

typedef uint32_t flags;
const flags no_output = 1;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags no_dense = 64;
const flags in_edges = 128;        // map over in edges instead of out edges
const flags dense_only = 1024;
inline bool should_output(const flags& fl) { return !(fl & no_output); }

}  // namespace gbbs
