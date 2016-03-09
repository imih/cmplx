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
    const IGraph &g, const Realization &realization, int no_simulations) {
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
        outcomes += DMCSingleSourceSirSimulation(v, g, realization);
      }
    }
    sum += outcomes;
    outcomes_prob.push_back((double)outcomes);
  }
  for (int v = 0; v < population_size; ++v) {
    outcomes_prob[v] /= sum;
  }
  return outcomes_prob;
}

int SourceDetector::DMCSingleSourceSirSimulation(
    int source_id, const IGraph &g, const Realization &realization) {
  int maxT = realization.maxT();
  SirParams params0 = paramsForSingleSource(source_id, realization);
  for (int t = 0; t < maxT; ++t) {
    Simulator::NaiveSIROneStep(g, params0);
    if ((realization.realization() | params0.infected()).bitCount() !=
        realization.realization().bitCount()) {
      return 0;
    }
  }
  return realization.realization().bitCount() ==
         (params0.infected() | params0.recovered()).bitCount();
}

vector<double>
SourceDetector::softMarginDetection(const IGraph &g,
                                    const Realization &realization,
                                    int no_simulations, double a) {
  vector<double> P;
  int population_size = g.vertices();
  for (int v = 0; v < population_size; ++v) {
    vector<double> fi;
    for (int s = 0; s < no_simulations; ++s) {
      fi.push_back(SMSingleSourceSirSimulation(v, g, realization));
    }
    P.push_back(likelihood(fi, a));
  }
  return P;
}

double SourceDetector::SMSingleSourceSirSimulation(
    int source_id, const common::IGraph &g,
    const common::Realization &realization) {
  SirParams params0 = paramsForSingleSource(source_id, realization);
  Simulator::NaiveSIR(g, params0);
  BitArray observed = params0.infected() | params0.recovered();
  return JaccardSimilarity(realization.realization(), observed);
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

double SourceDetector::likelihood(vector<double> fi, double a) {
  int n = (int)fi.size();
  double P = 0;
  for (int i = 0; i < n; ++i) {
    P += w_(fi[i], a);
  }
  return P / n;
}

} // namespace cmplx
