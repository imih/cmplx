#ifndef MPI_DIRECTMC_H
#define MPI_DIRECTMC_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/paral_directmc.h"

#include "mpi_master.h"

#include <vector>

namespace cmplx {

class MPIDirectMC : public ParalDirectMC, public MpiMaster {
 public:
  MPIDirectMC() {}
  ~MPIDirectMC() = default;

  std::vector<double> master(const SourceDetectionParams* params, bool end,
                             bool print);

  std::vector<double> master(const SourceDetectionParams* params) {
    return master(params, false, false);
  }

  void worker(const SourceDetectionParams* params, ModelType model_type);

  void send_simul_end();
};

}  // namespace cmplx
#endif  // MPI_DIRECTMC_H
