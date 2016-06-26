#include "igraph.h"
#include "ivector.h"

#include <gtest/gtest.h>
#include <iostream>
#include <vector>

namespace {
using cmplx::common::IGraph;
using cmplx::common::IVector;

class IGraphTest : public testing::Test {};

TEST_F(IGraphTest, IGraphLatticeSanity) {
  std::vector<int> dims({3, 3});
  IGraph* g = IGraph::UndirectedLattice(dims);
  EXPECT_EQ(g->vertices(), 9);
  EXPECT_EQ(g->diameter(), 4);
  EXPECT_EQ(g->adj_list(4).size(), 4);
  delete g;
}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
}  // namespace
