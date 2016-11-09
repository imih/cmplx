#ifndef COMMON_MASTER_H
#define COMMON_MASTER_H

#include "../source_detector/source_detection_params.h"

#include <vector>
#include <string>

namespace cmplx {

class CommonMaster {
 public:
  CommonMaster() {}
  ~CommonMaster() = default;

  virtual std::vector<double> master(const SourceDetectionParams* params) = 0;
};

class CommonTraits {
 public:
  CommonTraits(std::vector<int> benchmarkStepByStepSims,
               std::string benchmarkStepByStepPrefix,
               std::vector<int> convMasterSims)
      : benchmarkStepByStepSims_(benchmarkStepByStepSims),
        benchmarkStepByStepPrefix_(benchmarkStepByStepPrefix),
        convMasterSims_(convMasterSims) {}

  const std::vector<int>& benchmarkStepByStepSims() {
    return benchmarkStepByStepSims_;
  }

  const std::string& benchmarkStepByStepPrefix() {
    return benchmarkStepByStepPrefix_;
  }

  const std::vector<int>& convMasterSims() { return convMasterSims_; }

 protected:
  const std::vector<int> benchmarkStepByStepSims_;

  const std::string benchmarkStepByStepPrefix_;

  const std::vector<int> convMasterSims_;
};

}  // namespace cmplx
#endif  // COMMON_MASTER_H
