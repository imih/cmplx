#ifndef OMP_SOFTMC_PARAL_H
#define OMP_SOFTMC_PARAL_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_softmc.h"
#include "../source_detector/common_master.h"

#include <vector>

namespace cmplx {

class OMPSoftMC : public CommonMaster {
 public:
  OMPSoftMC() {}
  ~OMPSoftMC() = default;

  std::vector<double> master(const SourceDetectionParams* params);

 private:
  double work(const SourceDetectionParams* params, ModelType model_type,
              int node_id);
};

}  // namespace cmplx
#endif  // OMP_SOFTMC_PARAL_H
