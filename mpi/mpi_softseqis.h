#ifndef MPI_SOFT_SEQIS_H
#define MPI_SOFT_SEQIS_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/paral_seqis.h"

#include "mpi_seqis.h"

#include <vector>

namespace cmplx {

class MPISoftSeqIS : public MPISeqIS {
 public:
  MPISoftSeqIS() {}
  ~MPISoftSeqIS() = default;

  void worker(const SourceDetectionParams* params, ModelType model_type);
};

}  // namespace cmplx
#endif  // MPI_SOFT_SEQIS_H
