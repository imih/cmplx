#ifndef OMP_PARAL_H
#define OMP_PARAL_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"

namespace cmplx {

class OmpParal {
 public:
  OmpParal() {}

  void generateDistribution(SourceDetectionParams *params, ModelType model_type,
                            std::string &filename_prefix);

  void benchmark(SourceDetectionParams *params, int benchmark_no,
                 ModelType model_type, std::string filename_prefix);

 protected:
  double dabs(double x) {
    if (x < 0) return x * -1;
    return x;
  }

  virtual std::vector<double> convMaster(SourceDetectionParams *params) = 0;
};

}  // namespace cmplx
#endif  // OMP_PARAL_H
