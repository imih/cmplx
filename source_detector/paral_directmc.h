#ifndef PARAL_DIRECTMC_H
#define PARAL_DIRECTMC_H

#include "source_detection_params.h"

#include <vector>
#include <string>

using std::vector;

namespace cmplx {
class ParalDirectMC {
 public:
  void benchmarkStepByStep(cmplx::SourceDetectionParams* params,
                           int benchmark_no);

 protected:
  std::vector<double> responseToProb(const std::vector<int>& events_resp,
                                     int vertices);

  std::vector<double> convMaster(SourceDetectionParams* params);

 private:
  virtual std::vector<double> master(
      const cmplx::SourceDetectionParams* params) = 0;
};
}

#endif  // PARAL_DIRECTMC_H
