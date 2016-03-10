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

  SirParams(const SirParams &other)
      : infected_q_(other.infected()),
        time_steps_(other.time_steps()), p_(other.p()), q_(other.q()),
        T_(other.maxT()), infected_(other.infected()),
        susceptible_(other.susceptible()), recovered_(other.recovered()) {}

  ~SirParams() = default;

  IDqueue &infected_q() { return infected_q_; }

  int maxT() const { return T_; }
  void incrTime() { time_steps_++; }
  int time_steps() const { return time_steps_; }

  double p() const { return p_; }
  double q() const { return q_; }

  int population_size() const { return susceptible_.bits_num(); }

  BitArray infected() const { return infected_; }
  BitArray susceptible() const { return susceptible_; }
  BitArray recovered() const { return recovered_; }

  void set_infected(const BitArray &infected) { infected_ = infected; }
  void set_susceptible(const BitArray &susceptible) {
    susceptible_ = susceptible;
  }

  void set_recovered(const BitArray &recovered) { recovered_ = recovered; }

  void print() const;
  void printForLattice(int n) const;

private:
  IDqueue infected_q_;

  int time_steps_;
  double p_;
  double q_;
  int T_;

  BitArray infected_;
  BitArray susceptible_;
  BitArray recovered_;
};
} // namespace common
} // namespace cmplx

#endif // CMPLX_COMMON_SIR_PARAMS_H
