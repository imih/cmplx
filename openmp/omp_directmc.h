#ifndef OPNMP_DIRECTMC_PARAL_H
#define OPNMP_DIRECTMC_PARAL_H

#include "omp_paral.h"
#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"

#include <vector>

namespace cmplx {

class OMPDirectMCParal : public OmpParal {
 public:
  OMPDirectMCParal() {}

  void benchmarkStepByStep(cmplx::SourceDetectionParams* params,
                           int benchmark_no);

 private:
  std::vector<double> master(const SourceDetectionParams* params);

  std::vector<double> convMaster(SourceDetectionParams* params);
  double work(const SourceDetectionParams* params, ModelType model_type,
              int source_id, int batch_size);
};

}  // namespace cmplx
#endif  // OPNMP_DIRECTMC_PARAL_H
