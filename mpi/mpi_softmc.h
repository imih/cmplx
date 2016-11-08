#ifndef MPI_SOFTMC_H
#define MPI_SOFTMC_H

#include "mpi_paral.h"
#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_directmc.h"

#include <vector>

namespace cmplx {

class MPISoftMC : public MpiParal {
 public:
  MPISoftMC() {}

  void benchmarkStepByStep(SourceDetectionParams* params, int benchmark_no,
                           ModelType model_type);

 private:
  std::vector<double> master(const SourceDetectionParams* params, bool end,
                             bool print);
  std::vector<double> softConvMaster(SourceDetectionParams* params, bool end);

  std::vector<double> convMaster(SourceDetectionParams* params);
  void worker(const SourceDetectionParams* params, ModelType model_type);
  void send_simul_end();
};

}  // namespace cmplx
#endif  // MPI_SOFTMC_H
