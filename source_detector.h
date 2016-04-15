#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "simul/simulator.h"

#include <vector>

namespace cmplx {
enum ModelType {
  SIR,
  ISS
};

class SourceDetector {
 public:
  SourceDetector(const common::IGraph &g) : simulator_(g) {}

  // Return distribution of nodes being the source of SIR epidemic simulation
  // based on epidemic snapshot defined by sir_params.
  std::vector<double> directMonteCarloDetection(
      const common::Realization &realization, int no_simulations,
      ModelType model_type = ModelType::SIR);

  int DMCSingleSourceSimulation(int source_id,
                                const common::Realization &realization,
                                ModelType model_type = ModelType::SIR);

  // Return starting parameters for the epidemic that starts with a single
  // source defined by source_vertex and is capable of producing a snapshot
  // defined by
  // ending_params.
  common::SirParams paramsForSingleSource(
      int source_vertex, const common::Realization &realization);

  std::vector<double> softMarginDetection(
      const common::Realization &realization, int no_simulations, double a,
      ModelType model_type = ModelType::SIR);

  double SMSingleSourceSimulation(int source_id,
                                  const common::Realization &realization,
                                  ModelType model_type = ModelType::SIR);

  double likelihood(std::vector<double> fi, double a);

 private:
  double w_(double x, double a) {
    return exp(-1.0 * (x - 1) * (x - 1) / (a * a));
  }

  simul::Simulator simulator_;
};
}  // namespace cmplx
#endif  // CMPLX_SOURCE_DETECTOR_H
