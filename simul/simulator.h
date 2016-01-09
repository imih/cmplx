#ifndef CMPLX_SIMUL_SIMULATOR_H
#define CMPLX_SIMUL_SIMULATOR_H

#include "../common/bit_array.h"
#include "../common/igraph.h"
#include "../common/ivector.h"

namespace cmplx {
namespace simul {
class SirParams {
public:
  SirParams(double p, double q, const common::BitArray infected,
            const common::BitArray susceptible)
      : p_(p), q_(q), infected_(infected), susceptible_(susceptible) {}

  double p() { return p_; }
  double q() { return q_; }

  const common::BitArray &infected() { return infected_; }
  const common::BitArray &susceptible() { return susceptible_; }

private:
  double p_;
  double q_;
  common::BitArray infected_;
  common::BitArray susceptible_;
};

class SimulationStats {
public:
  SimulationStats(SirParams &sir_params)
      : infected_(sir_params.infected()),
        susceptible_(sir_params.susceptible()),
        recovered_(sir_params.infected().bits_num()) {}

  common::BitArray &infected() { return infected_; }
  common::BitArray &susceptible() { return susceptible_; }
  common::BitArray &recovered() { return recovered_; }

private:
  common::BitArray infected_;
  common::BitArray susceptible_;
  common::BitArray recovered_;
};

class Simulator {
public:
  static SimulationStats NaiveSIR(common::IGraph &graph, SirParams &sir_params);
};
} // namespace simul
} // namespace cmplx

#endif // CMPLX_SIMUL_SIMULATOR_H
