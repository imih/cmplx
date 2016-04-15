#ifndef CMPLX_SIMUL_SIMULATOR_H
#define CMPLX_SIMUL_SIMULATOR_H

#include "../common/bit_array.h"
#include "../common/idqueue.h"
#include "../common/igraph.h"
#include "../common/ivector.h"
#include "../common/sir_params.h"
#include <sys/time.h>
#include <random>
#include <unistd.h>

#include <string>

namespace cmplx {
namespace simul {

class Simulator {
 public:
  Simulator(const common::IGraph &graph)
      : graph_(graph), prob_distribution_(0, 1) {
    struct timeval tv;
    gettimeofday(&tv, 0);
    srand(tv.tv_usec * getpid());
    generator_.seed(tv.tv_usec * getpid());
  }

  bool NaiveSIR(common::SirParams &sir_params, bool prunning = false,
                const common::BitArray &allowed_nodes =
                    common::BitArray::zeros(1));

  // Returns  probability of drawn sample.
  double NaiveSIROneStep(common::SirParams &sir_params);
  
  bool NaiveISS(common::SirParams &sir_params, bool prunning = false,
      const common::BitArray &allowed_nodes = common::BitArray::zeros(1));

  bool eventDraw(double probability) {
    return prob_distribution_(generator_) <= probability;
  }

  const common::IGraph graph() { return graph_; }

 private:
  const common::IGraph &graph_;
  std::mt19937_64 generator_;
  std::uniform_real_distribution<double> prob_distribution_;

  // static void calcCummulativeInfecting(int n, const std::string&
  // file_name);
};
}  // namespace simul
}  // namespace cmplx

#endif  // CMPLX_SIMUL_SIMULATOR_H
