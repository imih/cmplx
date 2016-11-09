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
  virtual const std::vector<int>& benchmarkStepByStepSims() = 0;

  virtual const std::string& benchmarkStepByStepPrefix() = 0;

  virtual const std::vector<int>& convMasterSims() = 0;
};

}  // namespace cmplx
#endif  // COMMON_MASTER_H
