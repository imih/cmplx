#include "seq_sample.h"

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
  // w_ *= newPi / (pi_ * newG);
  w_ *= newPi / newG;
  pi_ = newPi;
  g_ = newG;
  // printf("%.10lf %.10lf %.10lf\n", pi_, g_, w_);
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

}  // namespace cmplx
