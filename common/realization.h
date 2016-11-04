#ifndef CMPLX_COMMON_REALIZATION_H
#define CMPLX_COMMON_REALIZATION_H

#include "bit_array.h"

namespace cmplx {
namespace common {

class Realization {
 public:
  Realization(double p, double q, int maxT, const BitArray& susceptible,
              const BitArray& infected, const BitArray& recovered);

  Realization(double p, double q, int maxT, const BitArray& susceptible)
      : Realization(p, q, maxT, susceptible,
                    BitArray::zeros(susceptible.bits_num()),
                    BitArray::zeros(susceptible.bits_num())) {}

  Realization()
      : Realization(0, 0, 0, BitArray::zeros(0), BitArray::zeros(0),
                    BitArray::zeros(0)) {}

  Realization(const Realization& other)
      : Realization(other.p(), other.q(), other.maxT(), other.susceptible(),
                    other.infected(), other.recovered()) {}

  ~Realization() = default;

  int maxT() const { return T_; }

  double p() const { return p_; }
  double q() const { return q_; }

  int population_size() const { return realization_.bits_num(); }

  const BitArray& susceptible() const { return susceptible_; }
  const BitArray& infected() const { return infected_; }
  const BitArray& recovered() const { return recovered_; }
  const BitArray& realization() const { return realization_; }

  void set_infected(const BitArray& infected);
  void set_recovered(const BitArray& recovered);

  void update(const BitArray& infected, const BitArray& recovered);

  void print() const { std::cout << "R: " << realization_ << std::endl; }
  void printForLattice(int n) const;

 private:
  const double p_;
  const double q_;
  const int T_;

  const BitArray susceptible_;
  BitArray infected_;
  BitArray recovered_;
  BitArray realization_;

  void sir_sanity_check();
  void calcRealization();
};

class RealizationRead {
 public:
  RealizationRead(const BitArray& realization, double p, double q, int maxT)
      : realization_(realization), p_(p), q_(q), maxT_(maxT) {}

  RealizationRead(const Realization& realization)
      : realization_(realization.realization()),
        p_(realization.p()),
        q_(realization.q()),
        maxT_(realization.maxT()) {}

  const BitArray& realization() const { return realization_; }
  double p() const { return p_; }
  double q() const { return q_; }
  double maxT() const { return maxT_; }
  int population_size() const { return realization_.bits_num(); }
  std::vector<int> positions() const { return realization_.positions(); }
  int bitCount() const { return realization_.bitCount(); }

  void setRealization(const BitArray& realization) {
    realization_ = realization;
  }

 private:
  BitArray realization_;
  double p_;
  double q_;
  int maxT_;
};

}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_SIR_PARAMS_H
