# Benchmarks and checks for the SSD backend. Builds only inside PLAID; run
# ssd/install_into_plaid.sh first (it also copies the BFS and PageRank
# headers into ssd/apps/). Everything is tagged manual so the gbbs repo's
# `bazel build //...` skips it.

package(
    default_visibility = ["//visibility:public"],
)

# The unchanged GBBS BFS and PageRank headers, copied in by the install script.
cc_library(
    name = "apps",
    hdrs = [
        "apps/BFS.h",
        "apps/PageRank.h",
    ],
    tags = ["manual"],
    deps = [
        "//gbbs",
        "//gbbs:bridge",
        "//gbbs:edge_map_data",
        "//gbbs:flags",
        "//gbbs:macros",
        "//gbbs:vertex_subset",
        "//gbbs/helpers:assert",
    ],
)

cc_binary(
    name = "ssd_bench",
    srcs = ["ssd_bench.cc"],
    tags = ["manual"],
    deps = [
        ":apps",
        "//gbbs:edge_map_data",
        "//gbbs:fetch_edges",
        "//gbbs:graph_io",
        "//gbbs:ssd_graph",
        "//gbbs:vertex_subset",
        "//graph:block_format",
        "//graph:semi_em_graph",
    ],
)
