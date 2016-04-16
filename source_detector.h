#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "simul/simulator.h"
#include "seq_sample.h"

#include <vector>
#include <set>

namespace cmplx {
enum ModelType {
  SIR,
  ISS
};

class SourceDetector {
 public:
  SourceDetector(const common::IGraph& g) : simulator_(g) {}

  // *******************DIRECT*************//
  // Return distribution of nodes being the source of SIR epidemic simulation
  // based on epidemic snapshot defined by sir_params.
  std::vector<double> directMonteCarloDetection(
      const common::Realization& realization, int no_simulations,
      ModelType model_type = ModelType::SIR);

  int DMCSingleSourceSimulation(int source_id,
                                const common::Realization& realization,
                                ModelType model_type = ModelType::SIR);

  // Return starting parameters for the epidemic that starts with a single
  // source defined by source_vertex and is capable of producing a snapshot
  // defined by
  // ending_params.
  common::SirParams paramsForSingleSource(
      int source_vertex, const common::Realization& realization);

  //*************SOFT ****************//

  std::vector<double> softMarginDetection(
      const common::Realization& realization, int no_simulations, double a,
      ModelType model_type = ModelType::SIR);

  double SMSingleSourceSimulation(int source_id,
                                  const common::Realization& realization,
                                  ModelType model_type = ModelType::SIR);

  double likelihood(std::vector<double> fi, double a);

  // ****************** SEQ ***************//
  std::vector<double> seqMonteCarloDetectionSIR(
      const common::Realization& realization);

  double seqPosterior(int v, const common::Realization& target_realization);

 private:
  double w_(double x, double a) {
    return exp(-1.0 * (x - 1) * (x - 1) / (a * a));
  }

  std::set<int> buildReachable(const common::BitArray& infected);
  std::pair<common::BitArray, common::BitArray> draw_g_1(
      double p, double q, const std::vector<int>& target_infected_idx,
      const common::BitArray& cur_inf, const common::BitArray& cur_rec);
  double g_1_cond(double p, double q,
                  const std::vector<int> target_infected_idx,
                  const common::BitArray& new_i,
                  const common::BitArray& new_rec,
                  const common::BitArray& old_i,
                  const common::BitArray& old_rec);
  double Pi_1_cond(double p, double q, const common::BitArray& new_i,
                   const common::BitArray& new_rec,
                   const common::BitArray& old_i,
                   const common::BitArray& old_rec);

  void printvc2(const std::vector<cmplx::SeqSample>& samples);

  simul::Simulator simulator_;
};
}  // namespace cmplx
#endif  // CMPLX_SOURCE_DETECTOR_H
