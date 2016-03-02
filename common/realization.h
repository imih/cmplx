#ifndef CMPLX_COMMON_REALIZATION_H
#define CMPLX_COMMON_REALIZATION_H

#include "bit_array.h"
#include "sir_params.h"

namespace cmplx {
namespace common {
class Realization {
public:
  Realization(const SirParams &snapshot)
      : Realization(snapshot.p(), snapshot.q(), snapshot.maxT(),
                    snapshot.infected() | snapshot.recovered()) {}

  Realization(double p, double q, int maxT, const BitArray &realization)
      : p_(p), q_(q), T_(maxT), realization_(realization) {}

  ~Realization() = default;

  int maxT() const { return T_; }

  double p() const { return p_; }
  double q() const { return q_; }

  int population_size() const { return realization_.bits_num(); }

  BitArray realization() const { return realization_; }

  void print() const { std::cout << "R: " << realization_ << std::endl; }

private:
  const double p_;
  const double q_;
  const int T_;
  const BitArray realization_;
};
} // namespace common
} // namespace cmplx

#endif // CMPLX_COMMON_SIR_PARAMS_H
