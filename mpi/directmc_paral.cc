#include "directmc_paral.h"

#include "../common/bit_array.h"
#include "../common/igraph.h"
#include "../common/realization.h"

#include <mpi.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>

#include <cstdio>
#include <cstdlib>

#include <cassert>

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::RealizationRead;
using cmplx::SourceDetectionParams;

using std::string;
using std::vector;

namespace {
const int SIMUL_PER_REQ = 10000;

struct Message {
  int source_id;
  int event_outcome;
  int batch_size;
};

MPI::Datatype datatypeOfMessage() {
  int blockLen[3] = {1, 1, 1};
  MPI::Aint offsets[3] = {offsetof(Message, source_id),
                          offsetof(Message, event_outcome),
                          offsetof(Message, batch_size)};
  MPI::Datatype types[3] = {MPI::INT, MPI::INT, MPI::INT};
  return MPI::Datatype::Create_struct(3, blockLen, offsets, types);
}

typedef long long ll;
}  // anonymous

// TODO There's not parameter sharing in DMC since it was ran on grid only!

namespace cmplx {

DirectMCParal::DirectMCParal() : MpiParal() {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
}

void DirectMCParal::send_simul_end() {
  MPI::Datatype message_type = datatypeOfMessage();
  for (int v = 1; v < processes_; ++v) {
    Message end_message;
    MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                          MessageType::SIMUL_END);
  }
}

void DirectMCParal::benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                                        int benchmark_no,
                                        ModelType model_type) {
  if (rank_ == 0) {
    string filename = "DMCbench_SBS_" + std::to_string(benchmark_no) + ".info";
    FILE *f = fopen(filename.c_str(), "a+");
    vector<int> sims = {(int)1e4, (int)1e5, (int)1e6,
                        (int)1e7, (int)1e8, (int)1e9};
    for (int sim : sims) {
      fprintf(f, "s: %d\n", sim);
      params->setSimulations(sim);
      vector<double> p = master(params, false, false);
      for (int i = 0; i < (int)p.size(); ++i)
        fprintf(f, "%.10lf%c", p[i], i + 1 == (int)p.size() ? '\n' : ' ');
      fflush(f);
    }
    fclose(f);
    send_simul_end();
  } else {
    worker(params, model_type);
  }
  exit(0);
}

vector<double> DirectMCParal::master(const SourceDetectionParams *params,
                                     bool end, bool print) {
  MPI::Datatype message_type = datatypeOfMessage();
  assert(rank == 0);

  const long long simulations = params->simulations();

  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  const common::RealizationRead &snapshot = params->realization();

  long long cur_simul_count = 0;
  int cur_v = nextV(0, snapshot.realization());
  vector<int> events_resp(vertices, 0);
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();
  int SIMUL_PER_REQ = std::max(10000LL, simulations / 10000);
  assert(simulations % (int)SIMUL_PER_REQ == 0);
  while (jobs_remaining > 0) {
    for (int i = 0; i < processes_ - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_PREREQUEST)) {
        Message init_message;
        MPI::COMM_WORLD.Recv(&init_message, 1, message_type, i + 1,
                             MessageType::SIMUL_PREREQUEST);
        if (cur_simul_count < simulations) {
          cur_simul_count += (int)SIMUL_PER_REQ;
        } else {
          cur_simul_count = (int)SIMUL_PER_REQ;
          cur_v = nextV(cur_v + 1, snapshot.realization());
        }
        if (cur_v == -1) {
          // There's no more jobs
          if (end == true) {
            Message message = {-1, -1};
            MPI::COMM_WORLD.Send(&message, 1, message_type, i + 1,
                                 MessageType::SIMUL_END);
          } else {
            // do nothing
          }
        } else {
          Message m = {cur_v, 0, SIMUL_PER_REQ};
          MPI::COMM_WORLD.Send(&m, 1, message_type, i + 1,
                               MessageType::SIMUL_REQUEST);
        }
      }
    }

    for (int i = 0; i < processes_ - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_RESPONSE)) {
        Message received;
        MPI::COMM_WORLD.Recv(&received, 1, message_type, i + 1,
                             MessageType::SIMUL_RESPONSE);
        jobs_remaining -= (int)SIMUL_PER_REQ;
        events_resp[received.source_id] += received.event_outcome;
        if (!i) {
          printf("\r%.5f",
                 jobs_remaining * 100 /
                     ((double)simulations * snapshot.realization().bitCount()));
          fflush(stdout);
        }
      }
    }
  }

  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    sum += events_resp[v];
  }
  assert(sum > 0);

  printf("\r\r\n");
  vector<double> p;
  p.clear();
  for (int v = 0; v < vertices; ++v) {
    p.push_back(events_resp[v] / sum);
  }
  return p;
}

vector<double> DirectMCParal::convMaster(SourceDetectionParams *params) {
  MPI::Datatype message_type = datatypeOfMessage();

  vector<int> sims = {100 * (int)1e4,   200 * (int)1e4,   1000 * (int)1e4,
                      2000 * (int)1e4,  10000 * (int)1e4, 20000 * (int)1e4,
                      100000 * (int)1e4};

  /***************  */
  double c = 0.05;  //
  /**************   */
  int s0 = sims[0];
  params->setSimulations(s0);
  vector<double> p0 = master(params, false, false);

  int MAP0 = std::max_element(p0.begin(), p0.end()) - p0.begin();
  double pml0 = *std::max_element(p0.begin(), p0.end());
  int bits = params->realization().bitCount();
  for (int s_id = 1; s_id < (int)sims.size(); ++s_id) {
    int s1 = sims[s_id];
    printf("s: %d\n", s1);
    params->setSimulations(s1);
    vector<double> p1 = master(params, false, false);

    int MAP1 = std::max_element(p1.begin(), p1.end()) - p1.begin();
    double pml1 = *std::max_element(p1.begin(), p1.end());
    bool converge = true;
    if (MAP0 != MAP1) converge = false;
    if (isnan(pml1)) {
      converge = false;
      printf("NAN! %d \n", s1);
    }
    double delta = dabs(pml1 - pml0) / pml1;
    printf("c: %lf\n", delta);
    if (delta > c) converge = false;
    int pos = 0;
    for (int i = 0; i < (int)p1.size(); ++i) {
      if (dabs(p0[i] - p1[i]) > c) converge = false;
      if (p1[i] > 0) pos++;
    }
    if (pos == 0) converge = false;
    if (converge) {
      printf("Converged for %d\n", s0);
      send_simul_end();
      break;
    }
    pml0 = pml1;
    s0 = s1;
    p0 = p1;
    MAP0 = MAP1;
  }
  params->setSimulations(s0);
  return p0;
}

void DirectMCParal::worker(const SourceDetectionParams *params,
                           ModelType model_type) {
  MPI::Datatype message_type = datatypeOfMessage();

  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  const common::RealizationRead &snapshot = params->realization();

  // workers
  // Performs simulation on request.
  DirectMonteCarloDetector sd(graph);
  while (true) {
    Message message = {-1, 0, 0};
    MPI::COMM_WORLD.Send(&message, 1, message_type, 0 /* dest */,
                         MessageType::SIMUL_PREREQUEST);
    if (MPI::COMM_WORLD.Iprobe(0, MessageType::SIMUL_REQUEST)) {
      Message message_recv;
      MPI::COMM_WORLD.Recv(&message_recv, 1, message_type, 0 /* source */,
                           MessageType::SIMUL_REQUEST);
      int outcomes = 0;
      for (int t = 0; t < message_recv.batch_size; ++t) {
        common::RealizationRead sp0 = snapshot;
        outcomes += sd.DMCSingleSourceSimulation(message_recv.source_id, sp0,
                                                 model_type);
      }

      message_recv.event_outcome = outcomes;
      Message toSend = message_recv;
      MPI::COMM_WORLD.Isend(&toSend, 1, message_type, 0 /* dest */,
                            MessageType::SIMUL_RESPONSE);
    }

    if (MPI::COMM_WORLD.Iprobe(0, MessageType::SIMUL_END)) {
      Message m;
      MPI::COMM_WORLD.Recv(&m, 1, message_type, 0, MessageType::SIMUL_END);
      break;
    }
  }
}

}  // namespace cmplx
