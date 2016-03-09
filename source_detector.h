#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "common/random.h"

#include <vector>

namespace cmplx {
class SourceDetector {
public:
  // Return distribution of nodes being the source of SIR epidemic simulation
  // based on epidemic snapshot defined by sir_params.
  std::vector<double>
  directMonteCarloDetection(const common::IGraph &g,
                            const common::Realization &realization,
                            int no_simulations, const common::Random &random);

  int DMCSingleSourceSirSimulation(int source_id, const common::IGraph &g,
                                   const common::Realization &realization,
                                   const common::Random &random);

  // Return starting parameters for the epidemic that starts with a single
  // source defined by source_vertex and is capable of producing a snapshot
  // defined by
  // ending_params.
  common::SirParams
  paramsForSingleSource(int source_vertex,
                        const common::Realization &realization);

  std::vector<double> softMarginDetection(
      const common::IGraph &g, const common::Realization &realization,
      int no_simulations, double a, const common::Random &random);

  double SMSingleSourceSirSimulation(int source_id, const common::IGraph &g,
                                     const common::Realization &realization,
                                     const common::Random &random);

private:
  double w_(double x, double a) {
    return exp(-1.0 * (x - 1) * (x - 1) / (a * a));
  }

  double calcP(std::vector<double> fi, double a);
};
} // namespace cmplx
#endif // CMPLX_SOURCE_DETECTOR_H
