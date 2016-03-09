#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "simul/simulator.h"
#include "common/random.h"
#include "source_detection_params.h"

#include <iostream>
#include <vector>

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::common::Random;
using cmplx::simul::Simulator;
using cmplx::SourceDetectionParams;

int main() {
  // DirectMCParams params = DirectMCParams::SupFig2Params();
  // DirectMCParams params = DirectMCParams::BenchmarkParams(1);
  SourceDetectionParams params = SourceDetectionParams::LatticeCenter();

  SourceDetector sd;
  Random r;
  /*
  std::vector<double> probs = sd.directMonteCarloDetection(
      params.graph(), params.realization(), params.simulations(), r);
      */
  std::vector<double> probs =
      sd.softMarginDetection(params.graph(), params.realization(),
                             params.simulations(), params.a(), r);

  for (double p : probs) {
    std::cout << p << " ";
  }
  std::cout << std::endl;

  return 0;
}
