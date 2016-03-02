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

  /* Calculates cummulative distribution of the probabilty the infected node
   * infects k neighbours
   * out of total n susceptible neigbours in the limit of discrete time under
   * SIR model(p, q).
   * Save in the file "file_name"  values of CDF for k = [1...n], p =
   * [0..0.1..1], q = [0..0.1..1].
   * TODO use precision arithmetic library.
   */

  static void NaiveSIROneStep(const common::IGraph &graph,
                              common::SirParams &sir_params);

  // static void calcCummulativeInfecting(int n, const std::string& file_name);
};
} // namespace simul
} // namespace cmplx

#endif // CMPLX_SIMUL_SIMULATOR_H
