#include "source_detector.h"

#include <cstring>
#include <mpi.h>
#include <set>

#include "common/bit_array.h"
#include "common/realization.h"
#include "common/ivector.h"

using cmplx::common::IGraph;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::common::BitArray;
using cmplx::simul::Simulator;

#include <set>

using std::vector;

namespace cmplx {
SirParams SourceDetector::paramsForSingleSource(
    int source_vertex, const Realization& realization) {
  int population_size = realization.population_size();
  // If the vertex was susceptible at some point.
  BitArray infected = BitArray::zeros(population_size);
  BitArray susceptible = realization.realization();
  infected.set(source_vertex, true);
  susceptible.set(source_vertex, false);
  return SirParams(realization.p(), realization.q(), realization.maxT(),
                   infected, susceptible);
}

vector<double> DirectMonteCarloDetector::directMonteCarloDetection(
    const Realization& realization, int no_simulations, ModelType model_type) {
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
        outcomes += DMCSingleSourceSimulation(v, realization, model_type);
      }
      sum += outcomes;
      outcomes_prob.push_back((double)outcomes);
    }
  }
  for (int v = 0; v < population_size; ++v) {
    outcomes_prob[v] /= sum;
  }
  return outcomes_prob;
}

int DirectMonteCarloDetector::DMCSingleSourceSimulation(
    int source_id, const Realization& realization, ModelType model_type) {
  SirParams params0 = paramsForSingleSource(source_id, realization);
  bool prunned = false;
  if (model_type == ModelType::SIR) {
    prunned = simulator_.NaiveSIR(params0, true, (realization.realization()));
  } else if (model_type == ModelType::ISS) {
    prunned = simulator_.NaiveISS(params0, true, (realization.realization()));
  }
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

vector<double> SoftMarginDetector::softMarginDetection(
    const Realization& realization, int no_simulations, double a,
    ModelType model_type) {
  vector<double> P;
  int population_size = realization.population_size();
  for (int v = 0; v < population_size; ++v) {
    vector<double> fi;
    for (int s = 0; s < no_simulations; ++s) {
      fi.push_back(SMSingleSourceSimulation(v, realization, model_type));
    }
    P.push_back(likelihood(fi, a));
  }
  return normalize(P);
}

double SoftMarginDetector::SMSingleSourceSimulation(
    int source_id, const common::Realization& realization,
    ModelType model_type) {
  SirParams params0 = paramsForSingleSource(source_id, realization);
  bool prunned = false;
  if (model_type == ModelType::SIR) {
    prunned = simulator_.NaiveSIR(params0);
  } else if (model_type == ModelType::ISS) {
    prunned = simulator_.NaiveISS(params0);
  }
  BitArray observed = params0.infected() | params0.recovered();
  return JaccardSimilarity(realization.realization(), observed);
}

double SoftMarginDetector::likelihood(vector<double> fi, double a) {
  int n = (int)fi.size();
  double P = 0;
  for (int i = 0; i < n; ++i) {
    P += w_(fi[i], a);
  }
  return P / n;
}

std::vector<double> SequentialMCDetector::seqMonteCarloDetectionSIR(
    const common::Realization& realization, int sample_size) {
  vector<double> P;
  int population_size = realization.population_size();
  double sum = 0;
  for (int v = 0; v < population_size; ++v) {
    if (realization.realization().bit(v)) {
      P.push_back(seqPosterior(v, sample_size, realization));
    } else
      P.push_back(0);
    sum += P.back();
  }

  for (int v = 0; v < population_size; ++v) P[v] /= sum;
  return P;
}

double SequentialMCDetector::seqPosterior(
    int v, int sample_size, const common::Realization& target_realization) {
  std::vector<SeqSample> prev_samples;
  prev_samples.clear();
  std::vector<SeqSample> samples;
  samples.assign(sample_size, SeqSample(v, target_realization));

  std::vector<int> target_infected_idx_ =
      target_realization.realization().positions();
  double p = target_realization.p();
  double q = target_realization.q();

  for (int t = 0; t < target_realization.maxT(); ++t) {
    prev_samples = samples;
    samples.clear();

    for (const SeqSample& sample : prev_samples) {
      BitArray prev_inf = sample.infected();
      BitArray prev_rec = sample.recovered();
      SeqSample newSample = sample;
      NewSample ns = drawSample(p, q, target_infected_idx_, prev_inf, prev_rec);
      newSample.update(ns.new_inf, ns.new_rec, ns.new_g, ns.new_pi);
      samples.push_back(newSample);
    }
    prev_samples.clear();
    printvc2(samples);
  }

  double pos_P = 0;
  for (const SeqSample& sample : samples) {
    if (sample.match(target_realization)) {
      pos_P += sample.w();
    }
  }
  printf("Post: %.10lf\n", pos_P / (int)samples.size());
  return pos_P / (int)samples.size();
}

SequentialMCDetector::NewSample SequentialMCDetector::drawSample(
    double p, double q, const std::vector<int>& target_infected_idx,
    const BitArray& prev_inf, const BitArray& prev_rec) {
  SequentialMCDetector::NewSample sample;
  std::set<int> reachable = buildReachable(prev_inf);
  BitArray next_inf = prev_inf;
  BitArray next_rec = prev_rec;
  double G = 1;
  for (int t : target_infected_idx) {
    if (prev_inf.bit(t)) {
      if (simulator_.eventDraw(q)) {
        G *= q;
        // try I -> R
        next_rec.set(t, true);
        next_inf.set(t, false);
      } else {
        G *= (1 - q);
      }
    } else if (reachable.count(t) && !prev_inf.bit(t) && !prev_rec.bit(t)) {
      G *= 0.5;
      if (simulator_.eventDraw(0.5)) {
        // S -> I
        next_inf.set(t, true);
      }
    }
  }
  sample.new_inf = next_inf;
  sample.new_rec = next_rec;
  sample.new_g = G;

  for (int p : prev_inf.positions()) {
    reachable.insert(p);
  }

  {
    double P = 1;
    for (int b : reachable) {
      if (prev_inf.bit(b) && !next_inf.bit(b)) {
        // printf("b: %d ", b);
        // puts("I -> R");
        P *= q;
      } else if (prev_inf.bit(b) && next_inf.bit(b)) {
        // printf("b: %d ", b);
        // puts("I -> I");
        P *= (1 - q);
      } else if (!prev_inf.bit(b) && !prev_rec.bit(b)) {
        // S ->
        const common::IGraph& graph = simulator_.graph();
        const common::IVector<int>& adj_list = graph.adj_list(b);

        int deg = 0;
        for (int i = 0; i < (int)adj_list.size(); ++i) {
          if (prev_inf.bit(adj_list[i])) deg++;
        }
        // printf("b: %d deg: %d ", b, deg);

        if (!next_inf.bit(b)) {
          // S -> S
          P *= pow(1 - p, deg);
        }
        if (next_inf.bit(b)) {
          // S -> I
          if (!deg) P = 0;
          P *= (1 - pow(1 - p, deg));
        }
      }
    }
    sample.new_pi = P;
  }

  return sample;
}

std::set<int> SequentialMCDetector::buildReachable(const BitArray& infected) {
  std::set<int> s;
  for (int pos : infected.positions()) {
    const IGraph& graph = simulator_.graph();
    const common::IVector<int>& adj_list = graph.adj_list(pos);
    for (int i = 0; i < (int)adj_list.size(); ++i) {
      s.insert(adj_list[i]);
    }
  }
  return s;
}

void SequentialMCDetector::printvc2(const std::vector<SeqSample>& samples) {
  printf("vc2 %.10lf\n", vc2(samples));
  printf("ESS %.10lf\n", ESS(samples));
}

double SequentialMCDetector::vc2(const std::vector<cmplx::SeqSample>& samples) {
  double avg_w = 0;
  double vc2 = 0;
  for (const SeqSample& sample : samples) {
    avg_w += sample.w();
  }
  for (const SeqSample& sample : samples) {
    vc2 += (sample.w() - avg_w) * (sample.w() - avg_w);
  }
  vc2 /= ((int)samples.size() - 1);
  vc2 /= (avg_w * avg_w);
  return vc2;
}

double SequentialMCDetector::ESS(const std::vector<cmplx::SeqSample>& samples) {
  double p = 0;
  double q = 0;
  for (const SeqSample& sample : samples) {
    p += sample.w();
    q += sample.w() * sample.w();
  }
  return p * p / q;
}

}  // namespace cmplx
