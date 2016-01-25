#ifndef CMPLX_COMMON_SIR_PARAMS_H
#define CMPLX_COMMON_SIR_PARAMS_H

#include "bit_array.h"
#include "idqueue.h"
#include "random.h"

namespace cmplx {
namespace common {
class SirParams {
public:
  SirParams(double p, double q, int T, const BitArray &infected,
            const BitArray &susceptible);

  IDqueue &infected_q() { return infected_q_; }

  bool drawP() { return random_.eventDraw(p_); }
  bool drawQ() { return random_.eventDraw(q_); }

  int maxT() const { return T_; }
  void incrTime() { time_steps_++; }
  int time_steps() const { return time_steps_; }

  double p() const { return p_; }
  double q() const { return q_; }

  int population_size() const { return susceptible_.bits_num(); }

  BitArray &infected() const { return const_cast<BitArray &>(infected_); }
  BitArray &susceptible() const { return const_cast<BitArray &>(susceptible_); }
  BitArray &recovered() const { return const_cast<BitArray &>(recovered_); }

private:
  Random random_;

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
