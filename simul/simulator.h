#ifndef CMPLX_SIMUL_SIMULATOR_H
#define CMPLX_SIMUL_SIMULATOR_H

#include "../common/bit_array.h"
#include "../common/igraph.h"
#include "../common/realization.h"
#include <sys/time.h>
#include <random>
#include <unistd.h>

#include <string>

namespace cmplx {
namespace simul {

class Simulator {
 public:
  Simulator(const common::IGraph *graph);

  bool NaiveSIR(common::Realization &sir_params, bool prunning = false,
                const common::BitArray &allowed_nodes =
                    common::BitArray::zeros(1));

  bool NaiveISS(common::Realization &sir_params, bool prunning = false,
                const common::BitArray &allowed_nodes =
                    common::BitArray::zeros(1));

  bool eventDraw(double probability) {
    return prob_distribution_(generator_) <= probability;
  }

  double P() { return prob_distribution_(generator_); }

  const common::IGraph *graph() { return graph_; }

 private:
  // NOT owned!
  const common::IGraph *graph_;

  std::mt19937_64 generator_;
  std::uniform_real_distribution<double> prob_distribution_;
};
}  // namespace simul
}  // namespace cmplx

#endif  // CMPLX_SIMUL_SIMULATOR_H
