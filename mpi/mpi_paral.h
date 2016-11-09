#ifndef MPI_PARAL_H
#define MPI_PARAL_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/common_paral.h"

#include "mpi_master.h"

namespace cmplx {

class MpiParal : public CommonParal {
 public:
  MpiParal(std::unique_ptr<MpiMaster> mpi_master,
           std::unique_ptr<CommonTraits> common_traits);
  ~MpiParal() = default;

  void benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                           int benchmark_no, ModelType model_type);

  void generateDistribution(SourceDetectionParams *params, ModelType model_type,
                            std::string &filename_prefix);

  void benchmark(SourceDetectionParams *params, int benchmark_no,
                 ModelType model_type, std::string filename_prefix);

 protected:
  int rank_;
  int processes_;

  int nextV(int cur_v, const common::BitArray &realization);

  std::vector<double> master(SourceDetectionParams *params) {
    return mpi_master_->master(params);
  }

 private:
  std::unique_ptr<MpiMaster> mpi_master_;
};

}  // namespace cmplx
#endif  // MPI_PARAL_H
