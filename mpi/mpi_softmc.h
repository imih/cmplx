#ifndef MPI_SOFTMC_H
#define MPI_SOFTMC_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_softmc.h"

#include "mpi_master.h"
#include <vector>

namespace cmplx {

class MPISoftMC : public MpiMaster {
 public:
  MPISoftMC() {}
  ~MPISoftMC() = default;

  std::vector<double> master(const SourceDetectionParams* params, bool end,
                             bool print);

  std::vector<double> master(const SourceDetectionParams* params) {
    return master(params, false, false);
  }

  void worker(const SourceDetectionParams* params, ModelType model_type);
  void send_simul_end();
};

}  // namespace cmplx
#endif  // MPI_SOFTMC_H
