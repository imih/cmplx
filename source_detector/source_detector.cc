#include "source_detector.h"

#include <set>

#include <cstring>
#include <map>
#include <thread>

#include "../common/bit_array.h"
#include "../common/ivector.h"
#include "../common/ivector.cc"

using cmplx::common::IGraph;
using cmplx::common::RealizationRead;
using cmplx::common::Realization;
using cmplx::common::BitArray;
using cmplx::simul::Simulator;

using std::vector;

namespace cmplx {
Realization SourceDetector::paramsForSingleSource(
    int source_vertex, const common::RealizationRead& realization) {
  int population_size = realization.population_size();
  // If the vertex was susceptible at some point.
  BitArray infected = BitArray::zeros(population_size);
  infected.set(source_vertex, true);
  BitArray susceptible = BitArray::ones(population_size);
  BitArray recovered = BitArray::zeros(population_size);
  return Realization(realization.p(), realization.q(), realization.maxT(),
                     susceptible, infected, recovered);
}

vector<double> DirectMonteCarloDetector::directMonteCarloDetection(
    const RealizationRead& realization, int no_simulations,
    ModelType model_type) {
  std::vector<double> outcomes_prob;
  outcomes_prob.clear();
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
    int source_id, const common::RealizationRead& realization,
    ModelType model_type) {
  Realization params0 = paramsForSingleSource(source_id, realization);
  bool prunned = false;
  switch (model_type) {
    case ModelType::SIR:
      prunned = simulator_->NaiveSIR(params0, true, realization.realization());
      break;
    case ModelType::ISS:
      prunned = simulator_->NaiveISS(params0, true, realization.realization());
      break;
  }
  if (prunned) return 0;
  return realization.bitCount() ==
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
    const RealizationRead& realization, int no_simulations, double a,
    ModelType model_type) {
  vector<double> P;
  P.clear();
  int population_size = realization.population_size();
  for (int v = 0; v < population_size; ++v) {
    vector<double> fi;
    fi.clear();
    for (int s = 0; s < no_simulations; ++s) {
      fi.push_back(SMSingleSourceSimulation(v, realization, model_type));
    }
    P.push_back(likelihood(fi, a));
  }
  return normalize(P);
}

double SoftMarginDetector::SMSingleSourceSimulation(
    int source_id, const common::RealizationRead& realization,
    ModelType model_type) {
  Realization params0 = paramsForSingleSource(source_id, realization);
  bool prunned = false;
  if (model_type == ModelType::SIR) {
    prunned = simulator_->NaiveSIR(params0);
  } else if (model_type == ModelType::ISS) {
    prunned = simulator_->NaiveISS(params0);
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
  return P;
}

std::vector<double> SequentialMCDetector::seqMonteCarloDetectionSIR(
    const common::RealizationRead& realization, int sample_size,
    ResamplingType resampling_type) {
  vector<double> P;
  P.clear();
  int population_size = realization.population_size();
  double sum = 0;
  for (int v = 0; v < population_size; ++v) {
    if (realization.realization().bit(v)) {
      P.push_back(seqPosterior(v, sample_size, realization, resampling_type,
                               true /* maximize hits*/));
    } else
      P.push_back(0);
    sum += P.back();
  }

  for (int v = 0; v < population_size; ++v) P[v] /= sum;
  return P;
}

double SequentialMCDetector::seqPosterior(
    int v, int sample_size, const common::RealizationRead& target_realization,
    cmplx::ResamplingType resampling_type, bool maximize_hits) {
  std::vector<SeqSample> samples(
      sample_size, SeqSample(v, target_realization.population_size()));

  std::vector<int> target_infected_idx_ =
      target_realization.realization().positions();
  double p = target_realization.p();
  double q = target_realization.q();

  for (int t = 0; t < target_realization.maxT(); ++t) {
    samples = resampling(t, resampling_type, samples);

    for (SeqSample& sample : samples) {
      BitArray prev_inf = sample.infected();
      BitArray prev_rec = sample.recovered();
      NewSample ns =
          drawSample(t, target_realization.maxT(), p, q, target_infected_idx_,
                     prev_inf, prev_rec, maximize_hits);
      sample.update(ns.new_inf, ns.new_rec, ns.new_g, ns.new_pi);
    }
  }
  return posteriorFromSamples(samples, target_realization.realization());
}

double SequentialMCDetector::posteriorFromSamples(
    const std::vector<SeqSample>& samples,
    const common::BitArray& realization) {
  double pos_P = 0;
  double sum = 0;
  for (const SeqSample& sample : samples) {
    if (sample.match(realization)) {
      pos_P += sample.w();
    }
    sum += sample.w();
  }
  fprintf(stderr, "\nPost: %.10lf\n", pos_P);
  return pos_P;
}

std::set<int> SequentialMCDetector::buildReachable(const BitArray& infected) {
  std::set<int> s;
  s.clear();
  std::vector<int> infected_positions = infected.positions();
  for (int pos : infected_positions) {
    const IGraph* graph = simulator_->graph();
    const common::IVector<int>& adj_list = graph->adj_list(pos);
    for (int i = 0; i < (int)adj_list.size(); ++i) {
      s.insert(adj_list[i]);
    }
  }
  return s;
}

SequentialMCDetector::NewSample SequentialMCDetector::drawSample(
    int T, int tMAX, double p, double q,
    const std::vector<int>& target_infected_idx, const BitArray& prev_inf,
    const BitArray& prev_rec, bool maximize_hits) {
  SequentialMCDetector::NewSample sample;
  std::set<int> reachable = buildReachable(prev_inf);
  BitArray next_inf = prev_inf;
  BitArray next_rec = prev_rec;
  const common::IGraph* graph = simulator_->graph();
  double G = 1;
  for (int t : target_infected_idx) {
    if (prev_inf.bit(t)) {
      double q2 = q;
      if (simulator_->eventDraw(q2)) {
        G *= q2;
        // try I -> R
        next_rec.set(t, true);
        next_inf.set(t, false);
      } else {
        G *= (1 - q2);
      }
    } else if (reachable.count(t) && !prev_inf.bit(t) && !prev_rec.bit(t)) {
      const common::IVector<int>& adj_list = graph->adj_list(t);
      int D = 0;
      for (int i = 0; i < (int)adj_list.size(); ++i) {
        if (prev_inf.bit(adj_list[i])) D++;
      }

      double p2 = 1 - pow((1 - p), D);
      if (maximize_hits && T == tMAX - 1) p2 = 1;
      if (simulator_->eventDraw(p2)) {
        // S -> I
        next_inf.set(t, true);
        G *= p2;
      } else {
        G *= (1 - p2);
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
        P *= q;
      } else if (prev_inf.bit(b) && next_inf.bit(b)) {
        P *= (1 - q);
      } else if (!prev_inf.bit(b) && !prev_rec.bit(b)) {
        // S ->
        const common::IVector<int>& adj_list = graph->adj_list(b);
        int D = 0;
        for (int i = 0; i < (int)adj_list.size(); ++i) {
          if (prev_inf.bit(adj_list[i])) D++;
        }

        if (!next_inf.bit(b)) {
          // S -> S
          P *= pow(1 - p, D);
        }
        if (next_inf.bit(b)) {
          // S -> I
          if (!D) P = 0;
          P *= (1 - pow(1 - p, D));
        }
      }
    }
    sample.new_pi = P;
  }

  return sample;
}

void SequentialMCDetector::printvc2(const std::vector<SeqSample>& samples) {
  printf("vc2 %.10lf\n", vc2(samples));
}

std::vector<SeqSample> SequentialMCDetector::resampling(
    int t, const ResamplingType& resampling_type,
    const vector<SeqSample>& samplesOrg) {
  int sample_size = (int)samplesOrg.size();
  if ((vc2(samplesOrg) <= (1LL << t)) ||
      resampling_type == ResamplingType::NONE) {
    return samplesOrg;
  }

  puts("resampling...");
  std::vector<SeqSample> prev_samples = samplesOrg;
  std::vector<SeqSample> samples;
  samples.clear();
  switch (resampling_type) {
    case(ResamplingType::NONE) : {
      return samplesOrg;
      break;
    }
    case(ResamplingType::SIMPLE_RANDOM_SAMPLING) : {
      double sum = 0;
      vector<double> P;
      P.clear();
      P.push_back(0);
      for (const SeqSample& sample : prev_samples) {
        sum += sample.w();
        P.push_back(P.back() + sample.w());
      }

      for (int i = 0; i < sample_size; ++i) {
        double p = simulator_->P() * sum;
        int lo = 0;
        int hi = sample_size - 1;
        while (lo < hi) {
          int mid = (lo + hi + 1) / 2;
          if (P[i] <= p) {
            lo = mid;
          } else
            hi = mid - 1;
        }
        samples.push_back(prev_samples[lo]);
        samples.back().setW(sum / sample_size);
      }
      break;
    }
    case(ResamplingType::RESIDUAL_SAMPLING) : {
      vector<double> Wnorm;
      Wnorm.clear();
      double sum = 0;
      for (const SeqSample& sample : prev_samples) {
        Wnorm.push_back(sample.w());
        sum += sample.w();
      }
      for (int i = 0; i < sample_size; ++i) {
        Wnorm[i] /= sum;
        int k = (int)(Wnorm[i] * sample_size);
        Wnorm[i] = Wnorm[i] * sample_size - k;
        for (int j = 0; j < k; ++j) samples.push_back(prev_samples[i]);
      }
      int m_r = sample_size - (int)samples.size();
      vector<double> P;
      P.clear();
      P.push_back(0);
      for (double wn : Wnorm) {
        P.push_back(P.back() + wn);
      }
      double maxP = P.back();
      for (int i = 0; i < m_r; ++i) {
        double p = simulator_->P() * maxP;
        int lo = 0;
        int hi = sample_size - 1;
        while (lo < hi) {
          int mid = (lo + hi + 1) / 2;
          if (P[i] <= p) {
            lo = mid;
          } else
            hi = mid - 1;
        }
        samples.push_back(prev_samples[lo]);
      }
      for (SeqSample& sample : samples) {
        sample.setW(sum / (double)samples.size());
      }
      break;
    }
  }

  return samples;
}

double SequentialMCDetector::vc2(const std::vector<cmplx::SeqSample>& samples) {
  double avg_w = 0;
  double vc2 = 0;
  for (const SeqSample& sample : samples) {
    avg_w += sample.w();
  }
  avg_w /= (double)samples.size();

  for (const SeqSample& sample : samples) {
    vc2 += (sample.w() - avg_w) * (sample.w() - avg_w);
  }
  vc2 /= (((double)samples.size() - 1) * (avg_w * avg_w));
  return vc2;
}

double SequentialMCDetector::ESS(const std::vector<cmplx::SeqSample>& samples) {
  return (int)samples.size() / (1 + vc2(samples));
}

std::vector<double> SequentialSoftMCDetector::seqMonteCarloDetectionSIR(
    const common::Realization& realization, int sample_size,
    ResamplingType resampling_type) {
  vector<double> P;
  P.clear();
  int population_size = realization.population_size();
  double sum = 0;
  for (int v = 0; v < population_size; ++v) {
    if (realization.realization().bit(v)) {
      P.push_back(
          seqPosterior(v, sample_size, realization, resampling_type, false));
    } else
      P.push_back(0);
    sum += P.back();
  }

  for (int v = 0; v < population_size; ++v) P[v] /= sum;
  return P;
}

double SequentialSoftMCDetector::posteriorFromSamples(
    const std::vector<SeqSample>& samples,
    const common::BitArray& target_realization) {
  double a = pow(2, -5);
  double pos_P = 0;
  for (const SeqSample& sample : samples) {
    pos_P += sample.w() *
             w_(JaccardSimilarity(target_realization, sample.realization()), a);
  }
  fprintf(stderr, "\nPost: %.10lf\n", pos_P);
  return pos_P;
}

}  // namespace cmplx
