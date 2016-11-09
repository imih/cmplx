#include "mpi_paral.h"

#include <mpi.h>
#include <cassert>

#include "mpi_common.h"

namespace cmplx {

namespace {

void share_params(SourceDetectionParams *params) {
  int rank = MPI::COMM_WORLD.Get_rank();
  int processes = MPI::COMM_WORLD.Get_size();
  if (rank == 0) {
    std::vector<int> r_pos = params->realization().positions();
    r_pos.push_back(-1);
    for (int v = 1; v < processes; ++v) {
      MPI::COMM_WORLD.Send(&r_pos[0], (int)r_pos.size(), MPI_INT, v,
                           MessageType::SIMUL_PARAMS);
    }
  } else {
    std::vector<int> r_pos;
    r_pos.resize(params->graph()->vertices() + 1);
    MPI::COMM_WORLD.Recv(&r_pos[0], (int)r_pos.size(), MPI_INT, 0,
                         MessageType::SIMUL_PARAMS);
    common::BitArray r_ba(params->graph()->vertices());
    for (int p : r_pos) {
      if (p == -1) break;
      r_ba.set(p, true);
    }
    params->setRealization(r_ba);
  }
}

}  // namespace

MpiParal::MpiParal(std::unique_ptr<MpiMaster> mpi_master,
                   std::unique_ptr<CommonTraits> common_traits)
    : CommonParal(std::move(common_traits)),
      mpi_master_(std::move(mpi_master)) {
  rank_ = MPI::COMM_WORLD.Get_rank();
  assert(rank_ > 0);
  processes_ = MPI::COMM_WORLD.Get_size();
}

void MpiParal::generateDistribution(SourceDetectionParams *params,
                                    ModelType model_type,
                                    std::string &filename_prefix) {
  MPI::COMM_WORLD.Barrier();
  share_params(params);
  MPI::COMM_WORLD.Barrier();
  // sleep(2);

  if (rank_ == 0) {
    generateDistribution(params, model_type, filename_prefix);
    mpi_master_->send_simul_end();
  } else {
    mpi_master_->worker(params, model_type);
  }
  exit(0);
}

void MpiParal::benchmark(SourceDetectionParams *params, int benchmark_no,
                         ModelType model_type, std::string filename_prefix) {
  if (rank_ == 0) {
    benchmark(params, benchmark_no, model_type, filename_prefix);
    mpi_master_->send_simul_end();
  } else {
    mpi_master_->worker(params, model_type);
  }
  exit(0);
}

void MpiParal::benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                                   int benchmark_no, ModelType model_type) {
  if (rank_ == 0) {
    benchmarkStepByStep(params, benchmark_no, model_type);
    mpi_master_->send_simul_end();
  } else {
    mpi_master_->worker(params, model_type);
  }
  exit(0);
}

}  // namespace cmplx
