#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"

#include <vector>

namespace cmplx {
class SourceDetector {
public:
  // Return distribution of nodes being the source of SIR epidemic simulation
  // based on epidemic snapshot defined by sir_params.
  std::vector<double>
  directMonteCarloDetection(const common::IGraph &g,
                            const common::Realization &realization,
                            int no_simulations, bool debug_print = false);

  int SSSirSimulation(int source_id, const common::IGraph &g,
                      const common::Realization &realization);

  // Return starting parameters for the epidemic that starts with a single
  // source define by source_vertex and is capable of producing a snapshot
  // defined by
  // ending_params.
  static common::SirParams
  paramsForSingleSource(int source_vertex,
                        const common::Realization &realization);

  static common::Realization
  createBenchmarkLatticeSnapshot(const common::IGraph &graph);
};
} // namespace cmplx
#endif // CMPLX_SOURCE_DETECTOR_H
