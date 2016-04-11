#ifndef CMPLX_COMMON_SIR_PARAMS_H
#define CMPLX_COMMON_SIR_PARAMS_H

#include "bit_array.h"
#include "idqueue.h"

namespace cmplx {
namespace common {
class SirParams {
 public:
  SirParams(double p, double q, int T, const BitArray &infected,
            const BitArray &susceptible);

  SirParams() : SirParams(0, 0, 0, BitArray::zeros(0), BitArray::zeros(0)) {}

  SirParams(const SirParams &other)
      : SirParams(other.p(), other.q(), other.maxT(), other.infected(),
                  other.susceptible()) {
    *recovered_ = other.recovered();
  }

  ~SirParams() {
    delete infected_;
    delete susceptible_;
    delete recovered_;
  }

  int maxT() const { return T_; }

  double p() const { return p_; }
  double q() const { return q_; }

  int population_size() const { return susceptible_->bits_num(); }

  const BitArray &infected() const { return *infected_; }
  const BitArray &susceptible() const { return *susceptible_; }
  const BitArray &recovered() const { return *recovered_; }

  void set_infected(const BitArray &infected) {
    *infected_ = BitArray(infected);
  }

  void set_susceptible(const BitArray &susceptible) {
    *susceptible_ = BitArray(susceptible);
  }

  void set_recovered(const BitArray &recovered) {
    *recovered_ = BitArray(recovered);
  }

  void print() const;
  void printForLattice(int n) const;

 private:
  double p_;
  double q_;
  int T_;

  BitArray *infected_;
  BitArray *susceptible_;
  BitArray *recovered_;
};
}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_SIR_PARAMS_H
