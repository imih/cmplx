#include "mpi_paral.h"

#include <mpi.h>
#include <cassert>

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

void print_to_file(std::string filename, std::vector<double> P,
                   SourceDetectionParams *params) {
  FILE *f = fopen(filename.c_str(), "w+");
  fprintf(f, "%s\n", params->summary().c_str());
  fprintf(f, "s: %lld\n", params->simulations());
  for (int i = 0; i < (int)P.size(); ++i)
    fprintf(f, "%.10lf%c", P[i], i + 1 == (int)P.size() ? '\n' : ' ');
  fclose(f);
}

}  // namespace

MpiParal::MpiParal() {
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
    std::vector<double> P = convMaster(params);
    std::string filename = filename_prefix + params->summary();
    FILE *f = fopen(filename.c_str(), "a");
    fprintf(f, "s:%lld -%d ", params->simulations(), params->sourceID());
    for (int j = 0; j < (int)P.size(); ++j) {
      fprintf(f, "%.10lf%c", P[j], j == ((int)P.size() - 1) ? '\n' : ' ');
    }
    fclose(f);
    send_simul_end();
  } else {
    worker(params, model_type);
  }
  exit(0);
}

void MpiParal::benchmark(SourceDetectionParams *params, int benchmark_no,
                         ModelType model_type, std::string filename_prefix) {
  if (rank_ == 0) {
    std::vector<double> P = convMaster(params);
    std::string filename =
        filename_prefix + std::to_string(benchmark_no) + ".info";
    print_to_file(filename, P, params);
    send_simul_end();
  } else {
    worker(params, model_type);
  }
  exit(0);
}

int MpiParal::nextV(int cur_v, const common::BitArray &realization) {
  int vertices = realization.bits_num();
  while ((cur_v < vertices) && (realization.bit(cur_v) == false)) cur_v++;
  if (cur_v >= vertices) cur_v = -1;
  return cur_v;
}

}  // namespace cmplx
