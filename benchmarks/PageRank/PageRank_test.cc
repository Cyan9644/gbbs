// TODO: Add tests for `PageRank_edgeMap` with `has_in_edges = false`.

#include "benchmarks/PageRank/PageRank.h"

#include <cmath>
#include <cstddef>
#include <unordered_set>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "gbbs/bridge.h"
#include "gbbs/graph.h"
#include "gbbs/helpers/undirected_edge.h"
#include "gbbs/macros.h"
#include "gbbs/unit_tests/graph_test_utils.h"
#include "gbbs/vertex.h"

using ::testing::DoubleNear;
using ::testing::ElementsAre;
using ::testing::ElementsAreArray;
using ::testing::IsEmpty;
using ::testing::Pointwise;

namespace gbbs {

template <class Graph>
sequence<double> compute_pagerank(Graph& G, double eps = 0.000001,
                                  const sequence<uintE>& sources = {},
                                  double damping_factor = 0.85,
                                  size_t max_iters = 100) {
  return PageRank_edgeMap(G, eps, sources, damping_factor, max_iters,
                          /*has_in_edges=*/true);
}

TEST(PageRankTest, EmptyGraph) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_THAT(compute_pagerank(graph), IsEmpty());
}

TEST(PageRankTest, EdgelessGraph) {
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_THAT(compute_pagerank(graph),
              ElementsAre(1.0 / 3, 1.0 / 3, 1.0 / 3));
}

TEST(PageRankTest, ZeroIterations) {
  // Graph diagram:
  //     0 - 1 - 2
  //     3 (isolated)
  //
  constexpr uintE kNumVertices{4};
  const std::unordered_set<UndirectedEdge> kEdges{{0, 1}, {1, 2}};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_THAT(compute_pagerank(graph, /*eps=*/0.1, /*sources=*/{},
                                     /*damping_factor=*/0.85, /*max_iters=*/0),
              ElementsAre(0.25, 0.25, 0.25, 0.25));
}

TEST(PageRankTest, Cycle) {
  // Graph diagram:
  //     0 - 1 - 2 - 3 - 4 (loops back to 0)
  //
  constexpr uintE kNumVertices{5};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    EXPECT_THAT(compute_pagerank(graph),
                ElementsAre(0.2, 0.2, 0.2, 0.2, 0.2));
  }
}

TEST(PageRankTest, Path) {
  // Graph diagram:
  //     0 - 1 - 2
  //
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1},
      {1, 2},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_THAT(
      compute_pagerank(graph),
      Pointwise(DoubleNear(1e-4),
                             {0.2567570878, 0.4864858243, 0.2567570878}));
}

TEST(PageRankTest, RespectsEps) {
  // Graph diagram:
  //     0 - 1 - 2
  //
  constexpr uintE kNumVertices{3};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1},
      {1, 2},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    // Between iterations 4 and 5 and between iterations 5 and 6, the L1
    // distance between the PageRank values are approximately 0.295804 and
    // 0.251433 respectively. The next lines verify that assumption.
    const sequence<double> result_with_4_iterations =
        compute_pagerank(graph, /*eps=*/1e-6,
                               /*sources=*/{},
                               /*damping_factor=*/0.85, /*max_iters=*/4);
    const sequence<double> result_with_5_iterations =
        compute_pagerank(graph, /*eps=*/1e-6,
                               /*sources=*/{},
                               /*damping_factor=*/0.85, /*max_iters=*/5);
    const sequence<double> result_with_6_iterations =
        compute_pagerank(graph, /*eps=*/1e-6,
                               /*sources=*/{},
                               /*damping_factor=*/0.85, /*max_iters=*/6);

    double l1_distance_iteration_4_to_5 = 0.0;
    for (int i = 0; i < 3; ++i) {
      l1_distance_iteration_4_to_5 +=
          fabs(result_with_4_iterations[i] - result_with_5_iterations[i]);
    }
    ASSERT_GT(l1_distance_iteration_4_to_5, 0.28);
    ASSERT_LT(l1_distance_iteration_4_to_5, 0.3);

    double l1_distance_iteration_5_to_6 = 0.0;
    for (int i = 0; i < 3; ++i) {
      l1_distance_iteration_5_to_6 +=
          fabs(result_with_5_iterations[i] - result_with_6_iterations[i]);
    }
    ASSERT_GT(l1_distance_iteration_5_to_6, 0.24);
    ASSERT_LT(l1_distance_iteration_5_to_6, 0.26);

    // Run the algorithm with `eps` small enough to reach iteration 5, but not
    // iteration 6.
    EXPECT_THAT(
        compute_pagerank(graph,
                               /*eps=*/0.31 / 3, /*sources=*/{},
                               /*damping_factor=*/0.85,
                               /*max_iters=*/9),
        ElementsAreArray(result_with_5_iterations));

    // Run the algorithm with `eps` small enough to reach iteration
    // 6, but not iteration 7.
    EXPECT_THAT(
        compute_pagerank(graph,
                               /*eps=*/0.27 / 3, /*sources=*/{},
                               /*damping_factor=*/0.85,
                               /*max_iters=*/9),
        ElementsAreArray(result_with_6_iterations));
  }
}

TEST(PageRankTest, BasicUndirected) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6
  constexpr uintE kNumVertices{7};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<double> expected{0.1428571429, 0.1428571429, 0.0802049539,
                                  0.2074472351, 0.1389813363, 0.2074472351,
                                  0.0802049539};

  EXPECT_THAT(compute_pagerank(graph),
              Pointwise(DoubleNear(1e-4), expected));
}

TEST(PageRankTest, TwoSources) {
  // Graph diagram:
  //     0 - 1    2 - 3 - 4
  //                    \ |
  //                      5 -- 6
  constexpr uintE kNumVertices{7};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {2, 3}, {3, 4}, {3, 5}, {4, 5}, {5, 6},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  // The results below were computed by using NetworkX's source-based PageRank
  // as a reference implementation.
  const sequence<double> expected{0.2702716450, 0.2297283550, 0.0384308511,
                                  0.1356382979, 0.1518617021, 0.1356382979,
                                  0.0384308511};
  EXPECT_THAT(compute_pagerank(graph,
                                     /*eps=*/0.000001, {0, 4}),
              Pointwise(DoubleNear(1e-4), expected));
}

TEST(PageRankTest, AllButOneSource) {
  // Graph diagram:
  //     0 - 1 - 2 - 3 - 4 (loops back to 0)
  //
  constexpr uintE kNumVertices{5};
  const std::unordered_set<UndirectedEdge> kEdges{
      {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0},
  };
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  const sequence<double> expected{0.2137115266, 0.2137115266, 0.2009045292,
                                  0.1707678884, 0.2009045292};
  EXPECT_THAT(compute_pagerank(graph,
                                     /*eps=*/0.000001, {0, 1, 2, 4}),
              Pointwise(DoubleNear(1e-4), expected));
}

TEST(PageRankTest, TwoNodesNoEdges) {
  constexpr uintE kNumVertices{2};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    // No mass is migrated.
    const sequence<double> expected{0.5, 0.5};
    EXPECT_THAT(compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 1}),
                Pointwise(DoubleNear(1e-4), expected));
  }

  {
    // No mass is migrated.
    const sequence<double> expected{1.0, 0};
    EXPECT_THAT(compute_pagerank(graph,
                                       /*eps=*/0.000001, {0}),
                Pointwise(DoubleNear(1e-4), expected));
  }
}

TEST(PageRankTest, TwoNodesOneEdge) {
  // Graph diagram:
  //     0 - 1
  constexpr uintE kNumVertices{2};
  const std::unordered_set<UndirectedEdge> kEdges{{0, 1}};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    // Both have equal mass; again no mass gets migrated.
    const sequence<double> expected{0.5, 0.5};
    EXPECT_THAT(compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 1}),
                Pointwise(DoubleNear(1e-4), expected));
  }

  {
    // The steady state for this setting is when P[0] = 20/37 and P[1] = 17/37.
    // This can be verified by plugging in the values above into the PPR
    // equations
    // P[0] = 0.85*P[1] + 0.15
    // P[1] = 0.85*P[0]
    const sequence<double> expected{0.5405409317, 0.4594590684};
    EXPECT_THAT(compute_pagerank(graph,
                                       /*eps=*/0.000001, {0}),
                Pointwise(DoubleNear(1e-4), expected));
  }
}

TEST(PageRankTest, FourNodePathOneSource) {
  // Graph diagram:
  //     0 - 1 - 2 - 3
  constexpr uintE kNumVertices{4};
  const std::unordered_set<UndirectedEdge> kEdges{{0, 1}, {1, 2}, {2, 3}};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  {
    // The values below can be verified as being stationary for the following
    // equations (note that with only one source, the addedConstant is only
    // present for node 0)
    // P[0] = 0.85*P[1]/2 + 0.15/4
    // P[1] = 0.85*(P[0] + P[2]/2) + 0.15/4
    // P[2] = 0.85*(P[1]/2 + P[3]) + 0.15/4
    // P[3] = 0.85*P[2]/2 + 0.15/4
    const sequence<double> expected{0.17543856058862931, 0.32456143941137061,
                                    0.32456143941137061, 0.17543856058862931};
    EXPECT_THAT(compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 1, 2, 3}),
                Pointwise(DoubleNear(1e-4), expected));
  }

  {
    // The values below can be verified as being stationary for the following
    // equations (note that with only one source, the addedConstant is only
    // present for node 0)
    // P[0] = 0.85*P[1]/2 + 0.15
    // P[1] = 0.85*(P[0] + P[2]/2)
    // P[2] = 0.85*(P[1]/2 + P[3])
    // P[3] = 0.85*P[2]/2
    const sequence<double> expected{0.3022241274, 0.3581756964, 0.2383155317,
                                    0.1012846445};
    EXPECT_THAT(compute_pagerank(graph,
                                       /*eps=*/0.000001, {0}),
                Pointwise(DoubleNear(1e-4), expected));
  }

  {
    // The values below can be verified as being stationary for the following
    // equations (now that we have two sources, the added constant is split
    // between the two sources).
    // P[0] = 0.85*P[1]/2 + 0.15/2
    // P[1] = 0.85*(P[0] + P[2]/2)
    // P[2] = 0.85*(P[1]/2 + P[3])
    // P[3] = 0.85*P[2]/2 + 0.15/2
    const sequence<double> expected{0.2017542424, 0.2982457576, 0.2982457576,
                                    0.2017542424};
    EXPECT_THAT(compute_pagerank(graph,
                                       /*eps=*/0.000001, {0, 3}),
                Pointwise(DoubleNear(1e-4), expected));
  }
}

TEST(PageRankTest, WeightedDigraph) {
  // Graph diagram (each number wrapped in parentheses denotes the weight of its
  // nearest edge):
  //
  //            3
  //           ^ |
  //       (3) | | (4)
  //           | v
  //   0 ---->  1  <---- 2
  //      (1)   |    (2)
  //            | (5)
  //            v
  //            4
  //

  using WeightedDigraph = asymmetric_graph<asymmetric_vertex, float>;
  sequence<WeightedDigraph::edge> edges(
      {{0, 1, 1.0}, {1, 3, 3.0}, {1, 4, 5.0}, {2, 1, 2.0}, {3, 1, 4.0}});
  WeightedDigraph graph = WeightedDigraph::from_edges(edges, /*n=*/5);

  // The values below can be verified as being stationary for the following
  // equations:
  //
  // P[0] = (0.15 + 0.85 * P[4]) / 5
  // P[1] = (0.15 + 0.85 * P[4]) / 5
  //        + 0.85 * (P[0] * 1/1 + P[2] * 2/2 + P[3] * 4/4)
  // P[2] = (0.15 + 0.85 * P[4]) / 5
  // P[3] = (0.15 + 0.85 * P[4]) / 5 + 0.85 * (P[1] * 3/8)
  // P[4] = (0.15 + 0.85 * P[4]) / 5 + 0.85 * (P[1] * 5/8)
  EXPECT_THAT(compute_pagerank(graph),
              Pointwise(
                  DoubleNear(1e-4),
                  {2333.0 / 30348.0, 11360.0 / 30348.0, 2333.0 / 30348.0,
                   5954.0 / 30348.0, 8368.0 / 30348.0}));
}

// Tests the case where a node (in this case, node 2) does have some outgoing
// edges, but they all have weight zero. An implementation that doesn't handle
// that case carefully could end up producing NaNs due to computing the ratio
// 0/0.
TEST(PageRankTest, NodeWithAllOutgoingEdgesHavingZeroWeight) {
  // Graph diagram (each number wrapped in parentheses denotes the weight of its
  // nearest edge):
  //
  //        (1)
  //     0 -----> 3
  //     |        ^
  // (2) |        | (0)
  //     v  (1)   |
  //     1 -----> 2
  //       <-----
  //         (0)
  //

  using WeightedDigraph = asymmetric_graph<asymmetric_vertex, float>;
  sequence<WeightedDigraph::edge> edges(
      {{0, 1, 2.0}, {0, 3, 1.0}, {1, 2, 1.0}, {2, 1, 0.0}, {2, 3, 0.0}});
  WeightedDigraph graph = WeightedDigraph::from_edges(edges, /*n=*/4);

  // The values below can be verified as being stationary for the following
  // equations:
  //
  // P[0] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4
  // P[1] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4 + 0.85 * (P[0] * 2/3)
  // P[2] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4 + 0.85 * (P[1] * 1/1)
  // P[3] = (0.15 + 0.85 * P[2] + 0.85 * P[3]) / 4 + 0.85 * (P[0] * 1/3)
  EXPECT_THAT(compute_pagerank(graph),
              Pointwise(DoubleNear(1e-4),
                                     {600.0 / 3709.0, 940.0 / 3709.0,
                                      1399.0 / 3709.0, 770.0 / 3709.0}));
}

TEST(PageRankDeathTest, EpsNegative) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_DEATH(auto result = compute_pagerank(graph, /*eps=*/-1.0),
               "Failed assertion `eps >= 0.0`");
}

TEST(PageRankDeathTest, DampingFactorNegative) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_DEATH(
      auto result =
          compute_pagerank(graph,
                                 /*eps=*/0.0, {}, /*damping_factor=*/-1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TEST(PageRankDeathTest, DampingFactorEqualToOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_DEATH(
      auto result =
          compute_pagerank(graph,
                                 /*eps=*/0.0, {}, /*damping_factor=*/1.0),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

TEST(PageRankDeathTest, DampingFactorGreaterThanOne) {
  constexpr uintE kNumVertices{0};
  const std::unordered_set<UndirectedEdge> kEdges{};
  auto graph{graph_test::MakeUnweightedSymmetricGraph(kNumVertices, kEdges)};

  EXPECT_DEATH(
      auto result =
          compute_pagerank(graph,
                                 /*eps=*/0.0, {}, /*damping_factor=*/1.1),
      "Failed assertion `0.0 <= damping_factor && damping_factor < 1.0`");
}

}  // namespace gbbs