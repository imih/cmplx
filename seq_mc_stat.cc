#include "seq_mc_stat.h"

#include <algorithm>

using std::vector;
using cmplx::common::BitArray;

namespace cmplx {

SeqSample::SeqSample(int v, const common::Realization& realization)
    : realization_(BitArray::zeros(realization.population_size())) {
  realization_.set(v, true);
  t_ = 0;
  w_ = 1;
  g_ = 1;
  pi_ = 1;
}

SeqSample::SeqSample(const SeqSample& seqSample)
    : realization_(seqSample.realization()),
      t_(seqSample.t()),
      w_(seqSample.w()),
      pi_(seqSample.pi()),
      g_(seqSample.g()) {}

void SeqSample::update(const common::BitArray new_realization, double newG,
                       double newPi) {
  t_++;
  realization_ = new_realization;
  pi_ = newPi;
  g_ = newG;
  w_ *= pi_ / g_;
}

bool SeqSample::match(const common::Realization& realization) const {
  int realizationBC = realization.realization().bitCount();
  if ((realization_ | realization.realization()).bitCount() != realizationBC) {
    return false;
  }
  int newBC = realization_.bitCount();
  return newBC == realizationBC;
}

SeqMCStat::SeqMCStat(int v, const common::Realization& target_realization,
                     const common::IGraph& graph)
    : target_realization_(target_realization), v_(v), simulator_(graph) {
  target_infected_idx_ = target_realization_.realization().positions();
  samples_.assign(1000, SeqSample(v, target_realization));
  prev_samples_.clear();
}

void SeqMCStat::SISStep() {
  prev_samples_ = samples_;
  samples_.clear();

  for (const SeqSample& sample : prev_samples_) {
    common::BitArray prevR = sample.realization();
    SeqSample newSample = sample;
    common::BitArray newR = draw_g_1(prevR);
    double newG = g_1_cond(newR, prevR);
    double newPi = Pi_1_cond(newR, prevR);
    newSample.update(newR, newG, newPi);
    samples_.push_back(newSample);
  }
}

common::BitArray SeqMCStat::draw_g_1(const BitArray& seq_sample) {
  double p = target_realization_.p();
  BitArray next_sample = seq_sample;
  for (int i = 0; i < (int)target_infected_idx_.size(); ++i) {
    if (seq_sample.bit(i)) continue;
    if (simulator_.eventDraw(p)) {
      next_sample.set(target_infected_idx_[i], true);
    }
  }
  return next_sample;
}

double SeqMCStat::g_1_cond(const common::BitArray& new_r,
                           const common::BitArray& old_r) {
  int diff_all_events = (int)target_infected_idx_.size() - old_r.bitCount();
  int diff_pos = new_r.bitCount() - old_r.bitCount();
  return pow(target_realization_.p(), diff_pos) *
         pow(1 - target_realization_.p(), diff_all_events - diff_pos);
}

double Pi_1_cond(const common::BitArray& new_r,
                   const common::BitArray& old_r) {
  //TODO
  return 0;
}

void printvc2() {
  //TODO
}


}  // namespace cmplx
