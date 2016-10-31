#include "simulator.h"

#include "../common/ivector.h"
#include "../common/igraph.h"
#include "../common/bit_array.h"
#include "../common/realization.h"

#include <iostream>
#include <gtest/gtest.h>
#include <vector>

namespace {
using cmplx::common::IGraph;
using cmplx::common::IVector;
using cmplx::common::BitArray;
using cmplx::common::Realization;
using cmplx::simul::Simulator;

int Tmax = 100;
class SimulatorTest : public testing::Test {

 protected:
  Realization dummyRealization(double p, double q, int n) {
    BitArray infected(n);
    BitArray susceptible(n);
    for (int i = 0; i < n; ++i) {
      if (i == 4) {
        infected.set(i, 1);
      }
      susceptible.set(i, 1);
    }
    return Realization(p, q, Tmax, susceptible, infected, BitArray::zeros(n));
  }
};

TEST_F(SimulatorTest, SIROnLattice) {
  std::vector<int> dims({10, 10});
  IGraph* g = IGraph::UndirectedLattice(dims);
  EXPECT_EQ(g->vertices(), 100);
  EXPECT_EQ(g->diameter(), 18);

  Realization exp_sp = dummyRealization(0, 0, g->vertices());
  for (double p = 0; p < 1; p += 0.1)
    for (double q = 0; q < 1; q += 0.1) {
      for (int t = 0; t < Tmax; ++t) {
        Realization sp = dummyRealization(p, q, g->vertices());
        if (!p) {
          EXPECT_EQ(exp_sp.infected(), sp.infected());
          EXPECT_EQ(exp_sp.recovered(), sp.recovered());
        }

        if (!q) {
          EXPECT_EQ(exp_sp.recovered(), sp.recovered());
        }

        if (p == 1 && t >= g->diameter()) {
          BitArray zeros(g->vertices());
          EXPECT_EQ(zeros, sp.susceptible());
          if (q == 1) {
            EXPECT_EQ(zeros, sp.infected());
          }
        }
      }
    }

  delete g;
}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
}  // namespace
