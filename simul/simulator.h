#ifndef CMPLX_SIMUL_SIMULATOR_H
#define CMPLX_SIMUL_SIMULATOR_H

#include "../common/bit_array.h"
#include "../common/idqueue.h"
#include "../common/igraph.h"
#include "../common/ivector.h"
#include "../common/random.h"
#include "../common/sir_params.h"

#include <string>

namespace cmplx {
namespace simul {

class Simulator {
public:
  Simulator(const common::IGraph &graph)
      : graph_(graph), random_(1LL * time(NULL) * getpid()) {}

  void NaiveSIR(common::SirParams &sir_params);

  void NaiveSIROneStep(common::SirParams &sir_params);

private:
  bool draw(double p) { return random_.eventDraw(p); }

  const common::IGraph &graph_;
  common::Random random_;

  // static void calcCummulativeInfecting(int n, const std::string& file_name);
};
} // namespace simul
} // namespace cmplx

#endif // CMPLX_SIMUL_SIMULATOR_H
