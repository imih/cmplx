#ifndef CMPLX_SOURCE_DETECTOR_H
#define CMPLX_SOURCE_DETECTOR_H

#include "../common/igraph.h"
#include "../common/realization.h"
#include "../simul/simulator.h"
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

  /** Return starting parameters for the epidemic that starts with a single
   source defined by source_vertex and is capable of producing a snapshot
   defined by ending_params. It is assumed all nodes are susceptible.
   */
  common::Realization paramsForSingleSource(
      int source_vertex, const common::Realization& realization);

 protected:
  simul::Simulator simulator_;
};

class DirectMonteCarloDetector : public SourceDetector {
 public:
  DirectMonteCarloDetector(const common::IGraph* g) : SourceDetector(g) {}

  /** Return distribution of nodes being the source of SIR epidemic simulation
   based on epidemic snapshot defined by sir_params. */
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
    return exp(-1.0l * (x - 1) * (x - 1) / (a * a));
  }
};

enum ResamplingType {
  NONE,
  SIMPLE_RANDOM_SAMPLING,
  RESIDUAL_SAMPLING,
};

class SequentialMCDetector : public SourceDetector {
 public:
  SequentialMCDetector(const common::IGraph* g) : SourceDetector(g) {}

  virtual std::vector<double> seqMonteCarloDetectionSIR(
      const common::Realization& realization, int sample_size,
      ResamplingType resampling_type = ResamplingType::NONE);

  virtual double seqPosterior(int v, int sample_size,
                              const common::Realization& target_realization,
                              ResamplingType resampling_type,
                              bool maximize_hits);

  // Calculates unnormalized posterior probability based on the set of samples.
  virtual double posteriorFromSamples(
      const std::vector<SeqSample>& samples,
      const common::BitArray& target_realization);

 protected:
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
                       const common::BitArray& prev_rec,
                       bool maximize_hits = true);

  void printvc2(const std::vector<cmplx::SeqSample>& samples);

 private:
  std::vector<SeqSample> resampling(int t,
                                    const ResamplingType& resampling_type,
                                    const std::vector<SeqSample>& samplesOrg);

  double vc2(const std::vector<cmplx::SeqSample>& samples);
  double ESS(const std::vector<cmplx::SeqSample>& samples);
};

class SequentialSoftMCDetector : public SequentialMCDetector {
 public:
  SequentialSoftMCDetector(const common::IGraph* g) : SequentialMCDetector(g) {}

  std::vector<double> seqMonteCarloDetectionSIR(
      const common::Realization& realization, int sample_size,
      ResamplingType resampling_type);

  double posteriorFromSamples(const std::vector<SeqSample>& samples,
                              const common::BitArray& target_realization);

 private:
  double w_(double x, double a) {
    return exp(-1.0l * (x - 1) * (x - 1) / (a * a));
  }
};

class ConfigurationalBiasMCDetector : public SequentialMCDetector {
 public:
  ConfigurationalBiasMCDetector(const common::IGraph* g)
      : SequentialMCDetector(g) {}

  SeqSample drawFullSample(int v,
                           const common::Realization& target_realization);

  double seqPosterior(int v, int sample_size,
                      const common::Realization& target_realization,
                      ResamplingType resampling_type = ResamplingType::NONE,
                      bool maximize_hits = true);
};

}  // namespace cmplx
#endif  // CMPLX_SOURCE_DETECTOR_H
