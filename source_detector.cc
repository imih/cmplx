#include "source_detector.h"

#include <mpi.h>
#include <cstring>

#include "common/bit_array.h"
#include "simul/simulator.h"

using cmplx::common::IGraph;
using cmplx::common::SirParams;
using cmplx::common::BitArray;
using cmplx::simul::Simulator;

using std::vector;

namespace cmplx {
vector<double> SourceDetector::directMonteCarloDetection(
    const IGraph &g, const SirParams &sir_params, int no_simulations) {
  assert(g.vertices() == sir_params.population_size());
  std::vector<double> outcomes_prob;
  int population_size = g.vertices();
  for (int v = 0; v < population_size; ++v) {
    // P(source = v | snapshot)
    int outcomes = 0;
    for (int sim_id = 0; sim_id < no_simulations; ++sim_id) {
      outcomes += SSSirSimulation(v, g, sir_params);
    }
    outcomes_prob.push_back((double)outcomes / no_simulations);
  }
  return outcomes_prob;
}

int SourceDetector::SSSirSimulation(int source_id, const IGraph &g,
                                    const SirParams &sir_params) {
  SirParams snapshot = sir_params;
  BitArray exp_active_nodes = snapshot.infected() | snapshot.recovered();
  BitArray banned_nodes = ~exp_active_nodes;
  int maxT = sir_params.maxT();
  BitArray zeros(sir_params.population_size());

  SirParams params0 =
      SourceDetector::paramsForSingleSource(source_id, snapshot);
  for (int t = 0; t < maxT; ++t) {
    Simulator::NaiveSIROneStep(g, params0);
    if (!((banned_nodes & params0.infected()) == zeros)) {
      break;
    }
  }
  return (snapshot.infected() == params0.infected() &&
          snapshot.recovered() == params0.recovered())
             ? 1
             : 0;
}

SirParams
SourceDetector::paramsForSingleSource(int source_vertex,
                                      common::SirParams &ending_params) {
  int population_size = ending_params.population_size();
  BitArray infected(population_size);
  // If the vertex was susceptible at some point.
  BitArray susceptible = ending_params.recovered() | ending_params.infected() |
                         ending_params.susceptible();
  infected.set(source_vertex, true);
  susceptible.set(source_vertex, false);
  return SirParams(ending_params.p(), ending_params.q(), ending_params.maxT(),
                   infected, susceptible);
}

} // namespace cmplx
