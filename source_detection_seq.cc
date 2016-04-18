#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "source_detection_params.h"

#include <iostream>
#include <vector>

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::SourceDetectionParams;

int main() {
  // SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  auto params = SourceDetectionParams::ParamsFromGrid(0.5, 0.5, 5);

  cmplx::SequentialMCDetector sd(params->graph().get());

  std::vector<double> probs =
      sd.seqMonteCarloDetectionSIR(params->realization(), 160000);

  for (double p : probs) {
    printf("%.10lf\n", p);
  }
  return 0;
}
