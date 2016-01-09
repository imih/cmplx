#ifndef CMPLX_SIMUL_SIMULATOR_H
#define CMPLX_SIMUL_SIMULATOR_H

#include "../common/bit_array.h"
#include "../common/idqueue.h"
#include "../common/igraph.h"
#include "../common/ivector.h"
#include "../common/random.h"

namespace cmplx {
namespace simul {
class SirParams {
public:
  SirParams(double p, double q, int T, const common::BitArray infected,
            const common::BitArray susceptible)
      : random_(), infected_q_(infected), time_steps_(0), p_(p), q_(q), T_(T),
        infected_(infected), susceptible_(susceptible),
        recovered_(infected.bits_num()) {}

  common::IDqueue &infected_q() { return infected_q_; }

  bool drawP() { return random_.eventDraw(p_); }
  bool drawQ() { return random_.eventDraw(q_); }

  int maxT() { return T_; }
  void incrTime() { time_steps_++; }
  int time_steps() { return time_steps_; }

  // double p() { return p_; }
  // double q() { return q_; }

  common::BitArray &infected() { return infected_; }
  common::BitArray &susceptible() { return susceptible_; }
  common::BitArray &recovered() { return recovered_; }

private:
  common::Random random_;

  common::IDqueue infected_q_;

  int time_steps_;
  double p_;
  double q_;
  int T_;

  common::BitArray infected_;
  common::BitArray susceptible_;
  common::BitArray recovered_;
};

class Simulator {
public:
  static void NaiveSIROneStep(common::IGraph& graph, SirParams& sir_params);
  static void NaiveSIR(common::IGraph &graph, SirParams &sir_params);
};
} // namespace simul
} // namespace cmplx

#endif // CMPLX_SIMUL_SIMULATOR_H
