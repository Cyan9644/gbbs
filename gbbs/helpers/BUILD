licenses(["notice"])


package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "assert",
    hdrs = ["assert.h"],
    deps = [],
)

cc_library(
    name = "directed_edge",
    srcs = ["directed_edge.cc"],
    hdrs = ["directed_edge.h"],
    deps = [
        "//gbbs:macros",
    ],
)

cc_library(
    name = "parse_command_line",
    hdrs = ["parse_command_line.h"],
)

cc_library(
    name = "undirected_edge",
    srcs = ["undirected_edge.cc"],
    hdrs = ["undirected_edge.h"],
    deps = [
        "//gbbs:macros",
    ],
)
