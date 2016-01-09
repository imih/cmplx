#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "common/igraph.h"
#include "common/sir_params.h"

#include <vector>

namespace cmplx {
class SourceDetector {
public:
  // Return distribution of nodes being the source of sir epidemic simulation
  // with ending parameters sir_params.
  static std::vector<double> directMonteCarloDetection(common::IGraph &g,
                                        common::SirParams &sir_params,
                                        int no_simulations);

private:
  static common::SirParams
  paramsForSingleSource(int vertex, common::SirParams &ending_params);
};
} // namespace cmplx
#endif // CMPLX_SOURCE_DETECTOR_H
