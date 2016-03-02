#ifndef DIRECT_MC_PARAMS_H
#define DIRECT_MC_PARAMS_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/bit_array.h"

namespace cmplx {
class DirectMCParams {
public:
  static DirectMCParams SupFig2Params();

  const common::IGraph &graph() { return graph_; }
  const common::Realization &realization() { return realization_; }
  int simulations() { return simulations_; }

private:
  DirectMCParams(const common::IGraph &graph, const common::Realization &r,
                 int simulations)
      : graph_(graph), realization_(r), simulations_(simulations) {}

  common::IGraph graph_;
  common::Realization realization_;
  int simulations_;
};

} // namespace cmplx

#endif // DIRECT_MC_PARAMS_H
