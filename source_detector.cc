#include "source_detector.h"

#include "common/bit_array.h"
#include "simul/simulator.h"

using cmplx::common::IGraph;
using cmplx::common::SirParams;
using cmplx::common::BitArray;
using cmplx::simul::Simulator;

using std::vector;

namespace cmplx {
vector<double> SourceDetector::directMonteCarloDetection(IGraph &g,
                                                         SirParams &sir_params,
                                                         int no_simulations) {
  SirParams expected = sir_params;
  BitArray exp_active_nodes = expected.infected() | expected.recovered();
  int maxT = sir_params.maxT();
  int vertices = g.vertices();
  std::vector<double> outcomes_prob;
  for (int v = 0; v < vertices; ++v) {
    int outcomes = 0;
    for (int sim_id = 0; sim_id < no_simulations; ++sim_id) {
      SirParams params0 = SourceDetector::paramsForSingleSource(v, sir_params);
      for (int t = 0; t < maxT; ++t) {
        Simulator::NaiveSIROneStep(g, params0);
        // TODO pruning
      }
      if (expected.infected() == params0.infected() &&
          expected.susceptible() == params0.susceptible()) {
        outcomes++;
      }
    }
    outcomes_prob.push_back((double)outcomes / no_simulations);
  }
  return outcomes_prob;
}

SirParams
SourceDetector::paramsForSingleSource(int vertex,
                                      common::SirParams &ending_params) {
  int vertices = ending_params.infected().bits_num();
  BitArray infected(vertices);
  BitArray susceptible(vertices);
  for (int i = 0; i < vertices; ++i) {
    if (i == vertex) {
      infected.set(i, true);
      susceptible.set(i, false);
    } else if (ending_params.recovered().bit(i) ||
               ending_params.infected().bit(i) ||
               ending_params.susceptible().bit(i)) {
      susceptible.set(i, true);
    }
  }
  return SirParams(ending_params.p(), ending_params.q(), ending_params.maxT(),
                   infected, susceptible);
}

} // namespace cmplx
