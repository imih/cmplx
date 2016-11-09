#ifndef OPNMP_DIRECTMC_PARAL_H
#define OPNMP_DIRECTMC_PARAL_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_directmc.h"
#include "../source_detector/common_master.h"

#include <vector>

namespace cmplx {

class OMPDirectMCParal : public CommonMaster {
 public:
  OMPDirectMCParal() {}
  ~OMPDirectMCParal() = default;

  std::vector<double> master(const SourceDetectionParams* params);

 private:
  double work(const SourceDetectionParams* params, ModelType model_type,
              int source_id, int batch_size);
};

}  // namespace cmplx
#endif  // OPNMP_DIRECTMC_PARAL_H
