#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "simul/simulator.h"
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
  SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  // DirectMCParams params = DirectMCParams::BenchmarkParams(1);
  //SourceDetectionParams params = SourceDetectionParams::LatticeCenter();

  SourceDetector sd;
  /*
  std::vector<double> probs = sd.directMonteCarloDetection(
      params.graph(), params.realization(), params.simulations());
      */
  std::vector<double> probs =
      sd.softMarginDetection(params.graph(), params.realization(),
                             params.simulations(), params.a());

  for (double p : probs) {
    std::cout << p << " ";
  }
  std::cout << std::endl;

  return 0;
}
