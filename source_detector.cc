#include "source_detector.h"

#include <cstring>
#include <mpi.h>

#include "common/bit_array.h"
#include "common/realization.h"
#include "simul/simulator.h"

using cmplx::common::IGraph;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::common::Random;
using cmplx::common::BitArray;
using cmplx::simul::Simulator;

#include <set>

using std::vector;

namespace cmplx {
vector<double> SourceDetector::directMonteCarloDetection(
    const IGraph &g, const Realization &realization, int no_simulations,
    const Random &random) {
  assert(g.vertices() == realization.population_size());
  std::vector<double> outcomes_prob;
  int population_size = g.vertices();
  double sum = 0;
  for (int v = 0; v < population_size; ++v) {
    int outcomes = 0;
    if (!realization.realization().bit(v)) {
      outcomes = 0;
    } else {
      // P(source = v | snapshot)
      for (int sim_id = 0; sim_id < no_simulations; ++sim_id) {
        outcomes += SingleSourceSirSimulation(v, g, realization, random);
      }
    }
    sum += outcomes;
    outcomes_prob.push_back((double)outcomes);
  }
  for(int v = 0; v < population_size; ++v) {
    outcomes_prob[v] /= sum;
  }
  return outcomes_prob;
}

int SourceDetector::SingleSourceSirSimulation(int source_id, const IGraph &g,
                                    const Realization &realization,
                                    const Random &random) {
  int maxT = realization.maxT();
  SirParams params0 =
      SourceDetector::paramsForSingleSource(source_id, realization);
  for (int t = 0; t < maxT; ++t) {
    Simulator::NaiveSIROneStep(g, params0, random);
    if ((realization.realization() | params0.infected()).bitCount() !=
        realization.realization().bitCount()) {
      return 0;
    }
  }
  return realization.realization().bitCount() ==
         (params0.infected() | params0.recovered()).bitCount();
}

SirParams
SourceDetector::paramsForSingleSource(int source_vertex,
                                      const Realization &realization) {
  int population_size = realization.population_size();
  // If the vertex was susceptible at some point.
  BitArray infected = BitArray::zeros(population_size);
  BitArray susceptible = realization.realization();
  infected.set(source_vertex, true);
  susceptible.set(source_vertex, false);
  return SirParams(realization.p(), realization.q(), realization.maxT(),
                   infected, susceptible);
}

Realization
SourceDetector::createBenchmarkLatticeSnapshot(const IGraph &graph) {
  int vertices = graph.vertices();
  BitArray inf = BitArray::zeros(vertices);
  inf.set(vertices / 2, true);
  BitArray succ = BitArray::ones(vertices);
  succ.set(vertices / 2, false);
  SirParams snapshot(0.4 /*p*/, 0 /*q*/, 5 /* T*/, inf, succ);
  // snapshot.printForLattice((int)sqrt(vertices));
  Random random(time(NULL));
  Simulator::NaiveSIR(graph, snapshot, random);
  return Realization(snapshot);
}

} // namespace cmplx
