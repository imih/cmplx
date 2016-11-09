#ifndef MPI_SEQIS_H
#define MPI_SEQIS_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_seqis.h"

#include "mpi_master.h"

#include <vector>

namespace cmplx {

class MPISeqIS : public MpiMaster {
 public:
  MPISeqIS() {}
  ~MPISeqIS() = default;

  std::vector<double> master(const SourceDetectionParams* params, bool end,
                             bool print);
  std::vector<double> master(const SourceDetectionParams* params) {
    return master(params, false, false);
  }

  void worker(const SourceDetectionParams* params, ModelType model_type);
  void send_simul_end();
};

}  // namespace cmplx
#endif  // MPI_SEQIS_H
