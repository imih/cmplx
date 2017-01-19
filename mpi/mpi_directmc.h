#ifndef MPI_DIRECTMC_H
#define MPI_DIRECTMC_H

#include "../source_detector/paral_directmc.h"
#include "../source_detector/source_detection_params.h"

#include "mpi_master.h"

#include <vector>

namespace cmplx {

class MPIDirectMC : public MpiMaster {
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

  void master_steve(const common::BitArray& realization, int simulations,
                    int simul_per_req);

};


  static void master_steve_wrapper(MPIDirectMC* direct_mc,
                                   const common::BitArray& realization,
                                   int simulation, int simul_per_req);

}  // namespace cmplx
#endif  // MPI_DIRECTMC_H
