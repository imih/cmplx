#include "seq_sample.h"

#include <algorithm>

using std::vector;
using cmplx::common::BitArray;

namespace cmplx {

SeqSample::SeqSample(int v, int population_size)
    : infected_(BitArray::zeros(population_size)),
      recovered_(BitArray::zeros(population_size)) {
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

bool SeqSample::match(const common::BitArray& realization) const {
  int realizationBC = realization.bitCount();
  common::BitArray thisRealization = this->realization();
  if ((thisRealization | realization).bitCount() != realizationBC) {
    return false;
  }
  int newBC = thisRealization.bitCount();
  return newBC == realizationBC;
}

void SeqSample::update(const BitArray& new_infected,
                       const BitArray& new_recovered, double newG,
                       double newPi) {
  t_++;
  infected_ = new_infected;
  recovered_ = new_recovered;
  w_ *= newPi / newG;
  pi_ = newPi;
  g_ = newG;
}

}  // namespace cmplx
