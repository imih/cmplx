#ifndef SOURCE_DETECTION_PARAMS_H
#define SOURCE_DETECTION_PARAMS_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/bit_array.h"

namespace cmplx {
class SourceDetectionParams {
public:
  // DirectMC
  static SourceDetectionParams SupFig2Params();
  static SourceDetectionParams LatticeCenter();
  static SourceDetectionParams BenchmarkParams(int realization_no);

  const common::IGraph &graph() { return graph_; }
  const common::Realization &realization() { return realization_; }
  int simulations() { return simulations_; }
  double a() { return a_; }

private:
  SourceDetectionParams(const common::IGraph &graph,
                        const common::Realization &r, int simulations,
                        double a = 0.05)
      : graph_(graph), realization_(r), a_(a), simulations_(simulations) {}

  common::IGraph graph_;
  common::Realization realization_;
  double a_;
  int simulations_;
};

} // namespace cmplx

#endif // SOURCE_DETECTION_PARAMS_H
