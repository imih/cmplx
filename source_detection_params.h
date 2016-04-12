#ifndef SOURCE_DETECTION_PARAMS_H
#define SOURCE_DETECTION_PARAMS_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/bit_array.h"
#include <string>

namespace cmplx {
class SourceDetectionParams {
 public:
  // DirectMC
  static SourceDetectionParams SupFig2Params();
  static SourceDetectionParams LatticeCenter();
  static SourceDetectionParams BenchmarkParams(int realization_no);
  static SourceDetectionParams ParamsFromGrid(double p, double q, int n);

  const common::IGraph &graph() const { return graph_; }
  const common::Realization &realization() const { return realization_; }
  int simulations() const { return simulations_; }
  double a() const { return a_; }

  void setSimulations(int simulations) { simulations_ = simulations; }
  void setA(double a) { a_ = a; }
  void setRealization(const common::BitArray &r) {
    realization_.setRealization(r);
  }

  std::string summary() const;

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

}  // namespace cmplx

#endif  // SOURCE_DETECTION_PARAMS_H
