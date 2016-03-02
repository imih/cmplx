#include "simulator.h"

#include "../common/ivector.h"
#include "../common/igraph.h"
#include "../common/bit_array.h"
#include "../common/sir_params.h"

#include <iostream>
#include <gtest/gtest.h>
#include <vector>

namespace {
using cmplx::common::IGraph;
using cmplx::common::IVector;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::simul::Simulator;

int Tmax = 100;
class SimulatorTest : public testing::Test {

protected:
  SirParams dummySirParams(double p, double q, int n) {
    BitArray infected(n);
    BitArray susceptible(n);
    for (int i = 0; i < n; ++i) {
      if (i == 4) {
        infected.set(i, 1);
      } else {
        susceptible.set(i, 1);
      }
    }
    return SirParams(p, q, Tmax, infected, susceptible);
  }
};

TEST_F(SimulatorTest, SIROnLattice) {
  std::vector<int> dims({10, 10});
  IGraph g = IGraph::UndirectedLattice(dims);
  EXPECT_EQ(g.vertices(), 100);
  EXPECT_EQ(g.diameter(), 18);

  SirParams exp_sp = dummySirParams(0, 0, g.vertices());
  for (double p = 0; p < 1; p += 0.1)
    for (double q = 0; q < 1; q += 0.1) {
      for (int t = 0; t < Tmax; ++t) {
        SirParams sp = dummySirParams(p, q, g.vertices());
        /*
        std::cout << p << " " << q << std::endl;
        std::cout << "S " << sp.susceptible() << std::endl;
        std::cout << "I " << sp.infected() << std::endl;
        std::cout << "R " << sp.recovered() << std::endl << std::endl;
        */
        if (!p) {
          EXPECT_EQ(exp_sp.infected(), sp.infected());
          EXPECT_EQ(exp_sp.recovered(), sp.recovered());
        }

        if (!q) {
          EXPECT_EQ(exp_sp.recovered(), sp.recovered());
        }

        if (p == 1 && t >= g.diameter()) {
          BitArray zeros(g.vertices());
          EXPECT_EQ(zeros, sp.susceptible());
          if (q == 1) {
            EXPECT_EQ(zeros, sp.infected());
          }
        }
      }
    }
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
} // namespace