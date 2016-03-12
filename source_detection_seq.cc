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
using cmplx::common::Random;
using cmplx::SourceDetectionParams;

int main() {
  SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  // DirectMCParams params = DirectMCParams::BenchmarkParams(1);
  //SourceDetectionParams params = SourceDetectionParams::LatticeCenter();

  SourceDetector sd(params.graph());
  /*
  std::vector<double> probs =
      sd.directMonteCarloDetection(params.realization(), params.simulations());
      */
  std::vector<double> probs = sd.softMarginDetection(
      params.realization(), params.simulations(), params.a());

  for (double p : probs) {
    printf("%.10lf\n", p);
  }
  return 0;
}
