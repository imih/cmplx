#include "seq_mc_stat.h"

#include <algorithm>

using std::vector;
using cmplx::common::BitArray;

namespace cmplx {

SeqSample::SeqSample(int v, const common::Realization& realization)
    : infected_(BitArray::zeros(realization.population_size())),
      recovered_(BitArray::zeros(realization.population_size())) {
  infected_.set(v, true);
  t_ = 0;
  w_ = 1;
  g_ = 1;
  pi_ = 1;
}

SeqSample::SeqSample(const SeqSample& seqSample)
    : infected_(seqSample.infected()),
      recovered_(seqSample.recovered()),
      t_(seqSample.t()),
      w_(seqSample.w()),
      pi_(seqSample.pi()),
      g_(seqSample.g()) {}

void SeqSample::update(const BitArray& new_infected,
                       const BitArray& new_recovered, double newG,
                       double newPi) {
  t_++;
  infected_ = new_infected;
  recovered_ = new_recovered;
  pi_ = newPi;
  g_ = newG;
  w_ *= pi_ / g_;
}

bool SeqSample::match(const common::Realization& realization) const {
  int realizationBC = realization.realization().bitCount();
  if ((this->realization() | realization.realization()).bitCount() !=
      realizationBC) {
    return false;
  }
  int newBC = this->realization().bitCount();
  return newBC == realizationBC;
}

SeqMCStat::SeqMCStat(int v, const common::Realization& target_realization,
                     const common::IGraph& graph)
    : target_realization_(target_realization), v_(v), simulator_(graph) {
  target_infected_idx_ = target_realization_.realization().positions();
  samples_.assign(10000, SeqSample(v, target_realization));
  prev_samples_.clear();
}

void SeqMCStat::SISStep() {
  prev_samples_ = samples_;
  samples_.clear();

  for (const SeqSample& sample : prev_samples_) {
    BitArray prev_inf = sample.infected();
    BitArray prev_rec = sample.recovered();
    SeqSample newSample = sample;
    std::pair<BitArray, BitArray> ir_pair = draw_g_1(prev_inf, prev_rec);

    double newG = g_1_cond(ir_pair.first, ir_pair.second, prev_inf, prev_rec);
    double newPi = Pi_1_cond(ir_pair.first, ir_pair.second, prev_inf, prev_rec);
    newSample.update(ir_pair.first, ir_pair.second, newG, newPi);
    samples_.push_back(newSample);
  }
  prev_samples_.clear();
  printvc2();
}

std::pair<BitArray, BitArray> SeqMCStat::draw_g_1(const BitArray& cur_inf,
                                                  const BitArray& cur_rec) {
  double p = target_realization_.p();
  double q = target_realization_.q();
  BitArray next_inf = cur_inf;
  BitArray next_rec = cur_rec;
  for (int i = 0; i < (int)target_infected_idx_.size(); ++i) {
    if (cur_inf.bit(i) && simulator_.eventDraw(q)) {
      // try I -> R
      next_rec.set(i, true);
      next_inf.set(i, false);
    } else if (!cur_inf.bit(i) && !cur_rec.bit(i) && simulator_.eventDraw(p)) {
      // S -> I
      next_inf.set(target_infected_idx_[i], true);
    }
  }
  return std::make_pair(next_inf, next_rec);
}

double SeqMCStat::g_1_cond(const common::BitArray& new_i,
                           const common::BitArray& new_rec,
                           const common::BitArray& old_i,
                           const common::BitArray& old_rec) {
  int n = new_i.bits_num();
  double p = target_realization_.p();
  double q = target_realization_.q();
  double P = 1;
  for (int b = 0; b < n; ++b) {
    if (!old_i.bit(b) && !new_i.bit(b) && !old_rec.bit(b)) {
      // S -> S
      P *= (1 - p);
    }
    if (old_i.bit(b) && !new_i.bit(b)) {
      // I -> R
      P *= q;
    }
    if (old_i.bit(b) && new_i.bit(b)) {
      // I -> I
      P *= (1 - q);
    }
    if (!old_i.bit(b) && new_i.bit(b) && !old_rec.bit(b)) {
      // S -> I
      P *= p;
    }
  }
  return P;
}

double SeqMCStat::Pi_1_cond(const common::BitArray& new_i,
                            const common::BitArray& new_rec,
                            const common::BitArray& old_i,
                            const common::BitArray& old_rec) {
  int n = new_i.bits_num();
  double p = target_realization_.p();
  double q = target_realization_.q();
  double P = 1;
  for (int b = 0; b < n; ++b) {
    if (old_i.bit(b) && !new_i.bit(b)) {
      // I -> R
      P *= q;
    } else if (old_i.bit(b) && new_i.bit(b)) {
      // I -> I
      P *= (1 - q);
    } else if (!old_i.bit(b) && !old_rec.bit(b)) {
      // S ->
      const common::IGraph& graph = simulator_.graph();
      const common::IVector<int>& adj_list = graph.adj_list(b);

      int deg = 0;
      for (int i = 0; i < (int)adj_list.size(); ++i) {
        if (old_i.bit(adj_list[i])) deg++;
      }

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

void SeqMCStat::printvc2() {
  double avg_w = 0;
  double vc2 = 0;
  for (const SeqSample& sample : samples_) {
    avg_w += sample.w();
  }
  for (const SeqSample& sample : samples_) {
    vc2 += (sample.w() - avg_w) * (sample.w() - avg_w);
  }
  vc2 /= ((int)samples_.size() - 1);
  vc2 /= (avg_w * avg_w);
  printf("%.10lf\n", vc2);
}

}  // namespace cmplx
