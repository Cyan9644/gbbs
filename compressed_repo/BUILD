licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

# main header for GBBS
cc_library(
    name = "gbbs",
    hdrs = ["gbbs.h"],
    deps = [
        ":benchmark",
        ":bridge",
        ":graph",
        ":interface",
        ":macros",
        ":vertex_subset",
        "//gbbs/helpers:parse_command_line",
    ],
)

# main header for Julienne

cc_library(
    name = "bridge",
    hdrs = ["bridge.h"],
    linkopts = ["-pthread"],
    deps = [
        "@parlaylib//parlay:delayed_sequence",
        "@parlaylib//parlay:io",
        "@parlaylib//parlay:monoid",
        "@parlaylib//parlay:parallel",
        "@parlaylib//parlay:primitives",
        "@parlaylib//parlay:random",
        "@parlaylib//parlay:range",
        "@parlaylib//parlay:sequence",
        "@parlaylib//parlay:slice",
        "@parlaylib//parlay:utilities",
        "@parlaylib//parlay/internal:binary_search",
        "@parlaylib//parlay/internal:get_time",
    ],
)

cc_library(
    name = "compressed_vertex",
    hdrs = ["compressed_vertex.h"],
    deps = [
        ":bridge",
        ":macros",
        "//gbbs/encodings:decoders",
    ],
)

cc_library(
    name = "fetch_edges",
    hdrs = ["fetch_edges.h"],
    deps = [
        ":bridge",
        ":compressed_vertex",
        ":macros",
        ":vertex",
    ],
)

cc_library(
    name = "edge_array",
    hdrs = ["edge_array.h"],
    deps = [
        ":bridge",
        ":macros",
    ],
)

cc_library(
    name = "flags",
    hdrs = ["flags.h"],
)

cc_library(
    name = "edge_map_data",
    hdrs = ["edge_map_data.h"],
    deps = [
        ":bridge",
        ":fetch_edges",
        ":flags",
        ":vertex_subset",
    ],
)

cc_library(
    name = "graph",
    hdrs = ["graph.h"],
    deps = [
        ":bridge",
        ":compressed_vertex",
        ":edge_array",
        ":flags",
        ":macros",
        ":vertex",
    ],
)

cc_library(
    name = "graph_generator",
    hdrs = ["graph_generator.h"],
    deps = [
        ":bridge",
        ":graph",
        ":macros",
        ":vertex",
    ],
)

cc_library(
    name = "interface",
    hdrs = ["interface.h"],
    deps = [
        ":edge_array",
        ":edge_map_data",
        ":flags",
        ":graph",
        ":macros",
        ":vertex_subset",
    ],
)

cc_library(
    name = "graph_io",
    srcs = ["graph_io.cc"],
    hdrs = ["graph_io.h"],
    deps = [
        ":graph",
        ":io",
        ":macros",
        ":vertex",
        "@parlaylib//parlay:io",
    ],
)

cc_library(
    name = "io",
    srcs = ["io.cc"],
    hdrs = ["io.h"],
    deps = [
        ":macros",
    ],
)

cc_library(
    name = "benchmark",
    hdrs = ["benchmark.h"],
    deps = [
        ":graph_io",
    ],
)

cc_library(
    name = "macros",
    hdrs = ["macros.h"],
    deps = [
        ":bridge",
    ],
)

cc_library(
    name = "intersect",
    hdrs = ["intersect.h"],
    deps = [
        ":macros",
    ],
)

cc_library(
    name = "vertex",
    hdrs = ["vertex.h"],
    deps = [
        ":bridge",
        ":intersect",
        ":macros",
    ],
)

cc_library(
    name = "vertex_subset",
    hdrs = ["vertex_subset.h"],
    deps = [
        ":bridge",
        ":flags",
        ":macros",
    ],
)
