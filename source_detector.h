#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "common/igraph.h"
#include "simul/simulator.h"

namespace cmplx {
class SourceDetector {
public:
  static void directMonteCarloDetection(common::IGraph &g,
                                        simul::SirParams &sir_params);
};
} // namespace cmplx
#endif // CMPLX_SOURCE_DETECTOR_H
