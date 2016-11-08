#ifndef SOFTMC_PARAL_H
#define SOFTMC_PARAL_H

#include "omp_paral.h"
#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"

#include <vector>

namespace cmplx {

class OMPSoftMC : public OmpParal {
 public:
  OMPSoftMC() {}

  void benchmarkStepByStep(SourceDetectionParams* params, int benchmark_no);

 private:
  std::vector<double> master(const SourceDetectionParams* params);
  std::vector<double> softConvMaster(SourceDetectionParams* params);

  std::vector<double> convMaster(SourceDetectionParams* params);
  double work(const SourceDetectionParams* params, ModelType model_type,
              int node_id);
};

}  // namespace cmplx
#endif  // SOFTMC_PARAL_H
