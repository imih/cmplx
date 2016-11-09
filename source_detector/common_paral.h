#ifndef COMMON_PARAL_H
#define COMMON_PARAL_H

#include "source_detection_params.h"
#include "source_detector.h"
#include "common_master.h"

namespace cmplx {

class CommonParal {
 public:
  CommonParal(std::unique_ptr<CommonTraits> common_traits)
      : common_traits_(std::move(common_traits)) {}

  CommonParal(std::unique_ptr<CommonMaster> common_master,
              std::unique_ptr<CommonTraits> common_traits)
      : common_master_(std::move(common_master)),
        common_traits_(std::move(common_traits)) {}

  ~CommonParal() = default;

  void generateDistribution(SourceDetectionParams* params, ModelType model_type,
                            std::string& filename_prefix);

  void benchmark(SourceDetectionParams* params, int benchmark_no,
                 ModelType model_type, std::string filename_prefix);

  void benchmarkStepByStep(cmplx::SourceDetectionParams* params,
                           int benchmark_no, ModelType model_type);

  std::vector<double> convMaster(SourceDetectionParams* params,
                                 const std::vector<int>& sims);

 protected:
  void benchmarkStepByStep(std::string filename_prefix,
                           const std::vector<int>& sims,
                           cmplx::SourceDetectionParams* params,
                           int benchmark_no);

  double dabs(double x) {
    if (x < 0) return x * -1;
    return x;
  }

  std::vector<double> master(SourceDetectionParams* params) {
    return common_master_->master(params);
  }
  std::vector<double> convMaster(SourceDetectionParams* params) {
    return convMaster(params, common_traits_->convMasterSims());
  }

 private:
  std::unique_ptr<CommonMaster> common_master_;
  std::unique_ptr<CommonTraits> common_traits_;
};

}  // namespace cmplx
#endif  // COMMON_PARAL_H
