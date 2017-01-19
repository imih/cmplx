#ifndef MPI_MASTER_H
#define MPI_MASTER_H

#include "../source_detector/source_detection_params.h"
#include "../source_detector/source_detector.h"
#include "../source_detector/common_master.h"

#include <vector>
#include <mpi.h>

namespace cmplx {

class MpiMaster : public CommonMaster {
 public:
  MpiMaster() {
    processes_ = MPI::COMM_WORLD.Get_size();
    rank_ = MPI::COMM_WORLD.Get_rank();
  }
  ~MpiMaster() = default;

  virtual std::vector<double> master(const SourceDetectionParams* params) = 0;
  virtual void worker(const SourceDetectionParams* params,
                      ModelType model_type) = 0;
  virtual void send_simul_end() = 0;

 protected:
  int processes() { return processes_; }
  int rank() { return rank_; }

  int nextV(int cur_v, const common::BitArray& realization) {
    int vertices = realization.bits_num();
    while ((cur_v < vertices) && (realization.bit(cur_v) == false)) cur_v++;
    if (cur_v >= vertices) cur_v = -1;
    return cur_v;
  }

 private:
  int processes_;
  int rank_;
};

}  // namespace cmplx
#endif  // MPI_MASTER_H
