#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "simul/simulator.h"
#include "direct_mc_params.h"

#include <iostream>
#include <vector>

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::simul::Simulator;
using cmplx::DirectMCParams;

int main() {
  srand(time(NULL));
  DirectMCParams params = DirectMCParams::SupFig2Params();
  int simulations = 100000;

  SourceDetector sd;
  std::vector<double> probs = sd.directMonteCarloDetection(
      params.graph(), params.realization(), simulations, true);
  for (double p : probs) {
    std::cout << p << " ";
  }
  std::cout << std::endl;

  return 0;
}
