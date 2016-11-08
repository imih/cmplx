#ifndef MPI_DIRECTMC_H
#define MPI_DIRECTMC_H

#include "mpi_paral.h"
#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_directmc.h"

#include <vector>

namespace cmplx {

class MPIDirectMC : public MpiParal, public ParalDirectMC {
 public:
  MPIDirectMC() {}

  void benchmarkStepByStep(cmplx::SourceDetectionParams* params,
                           int benchmark_no, ModelType model_type);

 private:
  std::vector<double> master(const SourceDetectionParams* params, bool end,
                             bool print);

  std::vector<double> master(const SourceDetectionParams* params) {
    return master(params, false, false);
  }

  std::vector<double> convMaster(SourceDetectionParams* params) {
    return ParalDirectMC::convMaster(params);
  }

  void worker(const SourceDetectionParams* params, ModelType model_type);
  void send_simul_end();
};

}  // namespace cmplx
#endif  // MPI_DIRECTMC_H
