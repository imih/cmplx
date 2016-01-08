#ifndef CMPLX_SIMULATOR_H
#define CMPLX_SIMULATOR_H

#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/ivector.h"

namespace cmplx {
class SirParams {
 public:
  SirParams(double p, double q) : p_(p), q_(q) {}
  double p() { return p_; }
  double q() { return q_; }

  const common::BitArray& infected() { return infected_; }
  const common::BitArray& susceptible_() { return susceptible_;}

 private:
  double p_;
  double q_;
  common::BitArray infected_;
  common::BitArray susceptible_;
};

class SimulationStats {
 public:
  SimulationStats() {}

  common::BitArray& infected() { return infected_; }
  common::BitArray& susceptible() { return susceptible_;}
  common::BitArray& recovered() {  return recovered_; }
  // TODO
 private:
  common::BitArray infected_;
  common::BitArray susceptible_;
  common::BitArray recovered_;
};

class Simulator {
 public:
  static SimulationStats NaiveSIRStep(const common::IGraph& graph,
                                      const SirParams& sir_params);

  /*
  static SimulationStats simulateSIRNaive(
      const common::IGraph& graph, const common::IVector<int>& infected_nodes,
      const SirParams& sir_params);
      */
};
}  // namespace cmplx

#endif  // CMPLX_SIMULATOR_H
