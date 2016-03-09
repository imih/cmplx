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
  static void NaiveSIR(const common::IGraph &graph,
                       common::SirParams &sir_params);

  static void NaiveSIROneStep(const common::IGraph &graph,
                              common::SirParams &sir_params);

  // static void calcCummulativeInfecting(int n, const std::string& file_name);
};
} // namespace simul
} // namespace cmplx

#endif // CMPLX_SIMUL_SIMULATOR_H
