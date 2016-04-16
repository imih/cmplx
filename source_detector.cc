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
vector<double> SourceDetector::directMonteCarloDetection(
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

int SourceDetector::DMCSingleSourceSimulation(int source_id,
                                              const Realization& realization,
                                              ModelType model_type) {
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

vector<double> SourceDetector::softMarginDetection(
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

double SourceDetector::SMSingleSourceSimulation(
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

std::vector<double> SourceDetector::seqMonteCarloDetectionSIR(
    const common::Realization& realization) {
  vector<double> P;
  int population_size = realization.population_size();
  double sum = 0;
  for (int v = 0; v < population_size; ++v) {
    if (realization.realization().bit(v)) {
      P.push_back(seqPosterior(v, realization));
    } else
      P.push_back(0);
    sum += P.back();
  }

  for (int v = 0; v < population_size; ++v) P[v] /= sum;
  return P;
}

double SourceDetector::seqPosterior(
    int v, const common::Realization& target_realization) {
  std::vector<SeqSample> prev_samples;
  prev_samples.clear();
  std::vector<SeqSample> samples;
  samples.assign(100000, SeqSample(v, target_realization));

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
      std::pair<BitArray, BitArray> ir_pair =
          draw_g_1(p, q, target_infected_idx_, prev_inf, prev_rec);
      std::cout << ir_pair.first.to_string() << std::endl
                << ir_pair.second.to_string() << std::endl;
      double newG = g_1_cond(p, q, target_infected_idx_, ir_pair.first,
                             ir_pair.second, prev_inf, prev_rec);
      double newPi =
          Pi_1_cond(p, q, ir_pair.first, ir_pair.second, prev_inf, prev_rec);
      printf("\nG: %.10lf Pi: %.10lf\n", newG, newPi);
      newSample.update(ir_pair.first, ir_pair.second, newG, newPi);
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

std::set<int> SourceDetector::buildReachable(const BitArray& infected) {
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

std::pair<BitArray, BitArray> SourceDetector::draw_g_1(
    double p, double q, const vector<int>& target_infected_idx,
    const BitArray& cur_inf, const BitArray& cur_rec) {
  std::set<int> reachable = buildReachable(cur_inf);
  BitArray next_inf = cur_inf;
  BitArray next_rec = cur_rec;
  for (int t : target_infected_idx) {
    if (cur_inf.bit(t) && simulator_.eventDraw(q)) {
      // try I -> R
      next_rec.set(t, true);
      next_inf.set(t, false);
    } else if (reachable.count(t) && !cur_inf.bit(t) && !cur_rec.bit(t) &&
               simulator_.eventDraw(0.5)) {
      // S -> I
      next_inf.set(t, true);
    }
  }
  return std::make_pair(next_inf, next_rec);
}

double SourceDetector::g_1_cond(double p, double q,
                                const std::vector<int> target_infected_idx,
                                const common::BitArray& new_i,
                                const common::BitArray& new_rec,
                                const common::BitArray& old_i,
                                const common::BitArray& old_rec) {
  std::set<int> reachable = buildReachable(old_i);
  double P = 1;
  for (int b : target_infected_idx) {
    if (old_i.bit(b) && !new_i.bit(b)) {
      // printf("b: %d ", b);
      // puts("I -> R");
      P *= q;
    }
    if (old_i.bit(b) && new_i.bit(b)) {
      // printf("b: %d ", b);
      // puts("I -> I");
      P *= (1 - q);
    }
    if (reachable.count(b)) {
      if (!old_i.bit(b) && !new_i.bit(b) && !old_rec.bit(b)) {
        // printf("b: %d ", b);
        // puts("S -> S");
        P *= (1 - 0.5);
      }
      if (!old_i.bit(b) && new_i.bit(b) && !old_rec.bit(b)) {
        // printf("b: %d ", b);
        // puts("S -> I");
        P *= 0.5;
      }
    }
  }
  printf("\n");
  return P;
}

double SourceDetector::Pi_1_cond(double p, double q,
                                 const common::BitArray& new_i,
                                 const common::BitArray& new_rec,
                                 const common::BitArray& old_i,
                                 const common::BitArray& old_rec) {
  int n = new_i.bits_num();
  double P = 1;
  for (int b = 0; b < n; ++b) {
    if (old_i.bit(b) && !new_i.bit(b)) {
      //printf("b: %d ", b);
      //puts("I -> R");
      P *= q;
    } else if (old_i.bit(b) && new_i.bit(b)) {
      //printf("b: %d ", b);
      //puts("I -> I");
      P *= (1 - q);
    } else if (!old_i.bit(b) && !old_rec.bit(b)) {
      // S ->
      const common::IGraph& graph = simulator_.graph();
      const common::IVector<int>& adj_list = graph.adj_list(b);

      int deg = 0;
      for (int i = 0; i < (int)adj_list.size(); ++i) {
        if (old_i.bit(adj_list[i])) deg++;
      }
      //printf("b: %d deg: %d ", b, deg);

      if (!new_i.bit(b)) {
        // S -> S
        P *= pow(1 - p, deg);
      }
      if (new_i.bit(b)) {
        // S -> I
        if (!deg) P = 0;
        P *= (1 - pow(1 - p, deg));
      }
    }
  }

  return P;
}

void SourceDetector::printvc2(const std::vector<SeqSample>& samples) {
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
  printf("vc2 %.10lf\n", vc2);
}

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

double SourceDetector::likelihood(vector<double> fi, double a) {
  int n = (int)fi.size();
  double P = 0;
  for (int i = 0; i < n; ++i) {
    P += w_(fi[i], a);
  }
  return P / n;
}

}  // namespace cmplx
