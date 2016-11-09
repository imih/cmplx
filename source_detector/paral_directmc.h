#ifndef PARAL_DIRECTMC_H
#define PARAL_DIRECTMC_H

#include "source_detection_params.h"
#include "common_paral.h"

#include <vector>
#include <string>

using std::vector;

namespace cmplx {
class ParalDirectMC : public CommonTraits {
 public:
  ParalDirectMC() {}
  ~ParalDirectMC() = default;

  const vector<int>& benchmarkStepByStepSims() {
    return benchmarkStepByStepSims_;
  }

  const std::string& benchmarkStepByStepPrefix() {
    return benchmarkStepByStepPrefix_;
  }

  const vector<int>& convMasterSims() { return convMasterSims_; }

 private:
  const vector<int> benchmarkStepByStepSims_ = {(int)1e4, (int)1e5, (int)1e6,
                                                (int)1e7, (int)1e8, (int)1e9};
  const std::string benchmarkStepByStepPrefix_ = "DMCbench_SBS_";

  const vector<int> convMasterSims_ = {
      (int)1e6,     2 * (int)1e6, 4 * (int)1e6, (int)1e7,     2 * (int)1e7,
      4 * (int)1e7, (int)1e8,     2 * (int)1e8, 4 * (int)1e8, (int)1e9};
};
}

#endif  // PARAL_DIRECTMC_H
