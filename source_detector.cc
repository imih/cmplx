#include "source_detector.h"

#include <cstring>
#include <mpi.h>

#include "seq_mc_stat.h"
#include "common/bit_array.h"
#include "common/realization.h"

using cmplx::common::IGraph;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::common::BitArray;
using cmplx::simul::Simulator;

#include <set>

using std::vector;

namespace cmplx {
vector<double> SourceDetector::directMonteCarloDetection(
    const Realization &realization, int no_simulations) {
  std::vector<double> outcomes_prob;
  int population_size = realization.population_size();
  double sum = 0;
  for (int v = 0; v < population_size; ++v) {
    int outcomes = 0;
    if (!realization.realization().bit(v)) {
      outcomes = 0;
    } else {
      // P(source = v | snapshot)
      for (int sim_id = 0; sim_id < no_simulations; ++sim_id) {
        outcomes += DMCSingleSourceSirSimulation(v, realization);
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
    int source_id, const Realization &realization) {
  SirParams params0 = paramsForSingleSource(source_id, realization);
  bool prunned =
      simulator_.NaiveSIR(params0, true, (realization.realization()));
  if (prunned) return 0;
  return realization.realization().bitCount() ==
         (params0.infected() | params0.recovered()).bitCount();
}

namespace {
vector<double> normalize(vector<double> P) {
  double sum = 0;
  for (double p : P) sum += p;
  for (int i = 0; i < (int)P.size(); ++i) P[i] /= sum;
  return P;
}
}  // namespace

vector<double> SourceDetector::softMarginDetection(
    const Realization &realization, int no_simulations, double a) {
  vector<double> P;
  int population_size = realization.population_size();
  for (int v = 0; v < population_size; ++v) {
    vector<double> fi;
    for (int s = 0; s < no_simulations; ++s) {
      fi.push_back(SMSingleSourceSirSimulation(v, realization));
    }
    P.push_back(likelihood(fi, a));
  }
  return normalize(P);
}

double SourceDetector::SMSingleSourceSirSimulation(
    int source_id, const common::Realization &realization) {
  SirParams params0 = paramsForSingleSource(source_id, realization);
  bool prunned = simulator_.NaiveSIR(params0);
  BitArray observed = params0.infected() | params0.recovered();
  return JaccardSimilarity(realization.realization(), observed);
}

SirParams SourceDetector::paramsForSingleSource(
    int source_vertex, const Realization &realization) {
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

SeqSample SourceDetector::forwardSeqSample(const SeqSample &seqSample) {
  SeqSample novi = seqSample;
  double u = simulator_.NaiveSIROneStep(novi.sir_params());
  novi.addW(u);
  return novi;
}

double SourceDetector::sequentialMCPosterior(
    int v, const common::Realization &realization) {
  std::vector<SeqSample> stat;
  stat.assign(10000, SeqSample(v, realization));
  int maxT = realization.maxT();
  std::string file_name = "seq_vulgaris_weights" + std::to_string(v);
  FILE *f = fopen(file_name.c_str(), "w");
  for (int t = 1; t <= maxT; ++t) {
    for (const SeqSample &seq : stat) {
      assert((int)seq.w().size() > 0);
      fprintf(f, "%.10lf ", seq.w().back());
    }
    fprintf(f, "\n");
    for (int i = 0; i < (int)stat.size(); ++i) {
      stat[i] = forwardSeqSample(stat[i]);
    }
  }

  fclose(f);

  double p = 0;
  for (const SeqSample &seq : stat) {
    if (seq.match(realization)) {
      // p += seq.w().back();
      p++;
    }
  }
  return p / (int)stat.size();
}

std::vector<double> SourceDetector::sequentialMCDetection(
    const common::Realization &realization) {
  int nodes = realization.population_size();
  std::vector<double> P;
  double sum = 0;
  for (int v = 0; v < nodes; ++v) {
    if (realization.realization().bit(v) == false)
      P.push_back(0);
    else
      P.push_back(sequentialMCPosterior(v, realization));
    sum += P.back();
  }
  for (double &P_v : P) {
    P_v /= sum;
  }
  return P;
}

}  // namespace cmplx
