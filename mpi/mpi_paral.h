#ifndef MPI_PARAL_H
#define MPI_PARAL_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"

namespace cmplx {

enum MessageType {
  SIMUL_PREREQUEST,
  SIMUL_REQUEST,
  SIMUL_RESPONSE,
  SIMUL_END,
  SIMUL_PARAMS
};

class MpiParal {
 public:
  MpiParal();

  void generateDistribution(SourceDetectionParams *params, ModelType model_type,
                            std::string &filename_prefix);

  void benchmark(SourceDetectionParams *params, int benchmark_no,
                 ModelType model_type, std::string filename_prefix);

 protected:
  int rank_;
  int processes_;

  double dabs(double x) {
    if (x < 0) return x * -1;
    return x;
  }

  int nextV(int cur_v, const common::BitArray &realization);

  virtual std::vector<double> convMaster(SourceDetectionParams *params) = 0;
  virtual void worker(const SourceDetectionParams *params,
                      ModelType model_type) = 0;
  virtual void send_simul_end() = 0;
};

}  // namespace cmplx
#endif  // MPI_PARAL_H
