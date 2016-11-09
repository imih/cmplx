#ifndef OMP_SEQISPARAL_H
#define OMP_SEQISPARAL_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_seqis.h"
#include "../source_detector/common_master.h"

#include <vector>

namespace cmplx {

class OMPSeqIS : public CommonMaster {
 public:
  OMPSeqIS() {}
  ~OMPSeqIS() = default;

  std::vector<double> master(const SourceDetectionParams* params);

 private:
  double work(const SourceDetectionParams* params, ModelType model_type,
              int source_id, int sample_size);
};

}  // namespace cmplx
#endif  // OMP_SEQISPARAL_H
