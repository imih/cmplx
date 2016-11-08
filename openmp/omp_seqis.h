#ifndef SEQIS_PARAL_H
#define SEQIS_PARAL_H

#include "omp_paral.h"
#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"

#include <vector>

namespace cmplx {

class OMPSeqIS : public OmpParal {
 public:
  OMPSeqIS() {}

  void benchmarkStepByStep(cmplx::SourceDetectionParams* params,
                           int benchmark_no);

 private:
  std::vector<double> master(const SourceDetectionParams* params);

  std::vector<double> convMaster(SourceDetectionParams* params);

  double work(const SourceDetectionParams* params, ModelType model_type,
              int source_id, int sample_size);
};

}  // namespace cmplx
#endif  // SEQIS_PARAL_H
