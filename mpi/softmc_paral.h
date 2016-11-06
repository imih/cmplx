#ifndef SOFTMC_PARAL_H
#define SOFTMC_PARAL_H

#include "mpi_paral.h"
#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"

#include <vector>

namespace cmplx {

class SoftMCParal : public MpiParal {
 public:
  SoftMCParal();

  void benchmarkStepByStep(SourceDetectionParams* params, int benchmark_no,
                           ModelType model_type);

 private:
  std::vector<double> master(const SourceDetectionParams* params, bool end,
                             bool print);
  std::vector<double> softConvMaster(SourceDetectionParams* params, bool end);

  std::vector<double> convMaster(SourceDetectionParams* params);
  void worker(const SourceDetectionParams* params, ModelType model_type);
  void send_simul_end();
};

}  // namespace cmplx
#endif  // SOFTMC_PARAL_H
