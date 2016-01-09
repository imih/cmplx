#ifndef CMPLX_SIMUL_SIMULATOR_H
#define CMPLX_SIMUL_SIMULATOR_H

#include "../common/bit_array.h"
#include "../common/idqueue.h"
#include "../common/igraph.h"
#include "../common/ivector.h"
#include "../common/random.h"
#include "../common/sir_params.h"

namespace cmplx {
namespace simul {

class Simulator {
public:
  static void NaiveSIROneStep(common::IGraph& graph, common::SirParams& sir_params);
  static void NaiveSIR(common::IGraph &graph, common::SirParams &sir_params);
};
} // namespace simul
} // namespace cmplx

#endif // CMPLX_SIMUL_SIMULATOR_H
