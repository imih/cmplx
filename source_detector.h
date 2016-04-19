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
  SourceDetector(const common::IGraph* g) : simulator_(g) {}

  // Return starting parameters for the epidemic that starts with a single
  // source defined by source_vertex and is capable of producing a snapshot
  // defined by
  // ending_params.
  common::SirParams paramsForSingleSource(
      int source_vertex, const common::Realization& realization);

 protected:
  simul::Simulator simulator_;
};

class DirectMonteCarloDetector : public SourceDetector {
 public:
  DirectMonteCarloDetector(const common::IGraph* g) : SourceDetector(g) {}

  // Return distribution of nodes being the source of SIR epidemic simulation
  // based on epidemic snapshot defined by sir_params.
  std::vector<double> directMonteCarloDetection(
      const common::Realization& realization, int no_simulations,
      ModelType model_type = ModelType::SIR);

  int DMCSingleSourceSimulation(int source_id,
                                const common::Realization& realization,
                                ModelType model_type = ModelType::SIR);
};

class SoftMarginDetector : public SourceDetector {
 public:
  SoftMarginDetector(const common::IGraph* g) : SourceDetector(g) {}

  std::vector<double> softMarginDetection(
      const common::Realization& realization, int no_simulations, double a,
      ModelType model_type = ModelType::SIR);

  double SMSingleSourceSimulation(int source_id,
                                  const common::Realization& realization,
                                  ModelType model_type = ModelType::SIR);

  double likelihood(std::vector<double> fi, double a);

 private:
  double w_(double x, double a) {
    return exp(-1.0 * (x - 1) * (x - 1) / (a * a));
  }
};

class SequentialMCDetector : public SourceDetector {
 public:
  SequentialMCDetector(const common::IGraph* g) : SourceDetector(g) {}

  std::vector<double> seqMonteCarloDetectionSIR(
      const common::Realization& realization, int sample_size);

  double seqPosterior(int v, int sample_size,
                      const common::Realization& target_realization);

 private:
  std::set<int> buildReachable(const common::BitArray& infected);

  struct NewSample {
    common::BitArray new_inf;
    common::BitArray new_rec;
    double new_g;
    double new_pi;
  };
  NewSample drawSample(int t, int tMax, double p, double q,
                       const std::vector<int>& target_infected_idx,
                       const common::BitArray& prev_inf,
                       const common::BitArray& prev_rec);

  double vc2(const std::vector<cmplx::SeqSample>& samples);
  double ESS(const std::vector<cmplx::SeqSample>& samples);
  void printvc2(const std::vector<cmplx::SeqSample>& samples);
};
}  // namespace cmplx
#endif  // CMPLX_SOURCE_DETECTOR_H
