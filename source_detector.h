#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "simul/simulator.h"
#include "seq_mc_stat.h"

#include <vector>

namespace cmplx {
class SourceDetector {
 public:
  SourceDetector(const common::IGraph &g) : simulator_(g) {}

  // Return distribution of nodes being the source of SIR epidemic simulation
  // based on epidemic snapshot defined by sir_params.
  std::vector<double> directMonteCarloDetection(
      const common::Realization &realization, int no_simulations);

  int DMCSingleSourceSirSimulation(int source_id,
                                   const common::Realization &realization);

  // Return starting parameters for the epidemic that starts with a single
  // source defined by source_vertex and is capable of producing a snapshot
  // defined by
  // ending_params.
  common::SirParams paramsForSingleSource(
      int source_vertex, const common::Realization &realization);

  std::vector<double> softMarginDetection(
      const common::Realization &realization, int no_simulations, double a);

  double SMSingleSourceSirSimulation(int source_id,
                                     const common::Realization &realization);

  double likelihood(std::vector<double> fi, double a);

  std::vector<double> sequentialMCDetection(
      const common::Realization &realization);

 private:
  double w_(double x, double a) {
    return exp(-1.0 * (x - 1) * (x - 1) / (a * a));
  }

  // returns posterior probability P(R = r* | theta = v)
  double sequentialMCPosterior(int v, const common::Realization &realization);
  SeqSample forwardSeqSample(const SeqSample &seqSample);

  simul::Simulator simulator_;
};
}  // namespace cmplx
#endif  // CMPLX_SOURCE_DETECTOR_H
