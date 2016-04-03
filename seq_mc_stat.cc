#include "seq_mc_stat.h"

using std::vector;
using cmplx::common::BitArray;

namespace cmplx {

SeqSample::SeqSample(int v, const common::Realization& realization) {
  common::BitArray infected = BitArray::zeros(realization.population_size());
  infected.set(v, true);
  sir_params_ = common::SirParams(
      realization.p(), realization.q(), realization.maxT(), infected,
      BitArray::ones(realization.population_size()));
  t_ = 0;
  w_ = vector<double>();
  w_.push_back(1);
}

SeqSample::SeqSample(const SeqSample& seqSample)
    : sir_params_(seqSample.sir_params_const()),
      t_(seqSample.t()),
      w_(seqSample.w()) {}

}  // namespace cmplx
