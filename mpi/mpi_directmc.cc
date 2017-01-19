#include "mpi_directmc.h"

#include "../common/bit_array.h"
#include "../common/igraph.h"
#include "../common/realization.h"

#include "mpi_common.h"

#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>

#include <cassert>

#include <thread>

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

void MPIDirectMC::send_simul_end() {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int processes = MpiMaster::processes();
  for (int v = 1; v < processes; ++v) {
    Message end_message;
    MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                          MessageType::SIMUL_END);
  }
}

static void master_steve_wrapper(MPIDirectMC *direct_mc,
                                 const common::BitArray &realization,
                                 int simulation, int simul_per_req) {
  direct_mc->master_steve(realization, simulation, simul_per_req);
}

void MPIDirectMC::master_steve(const common::BitArray &realization,
                               int simulations, int simul_per_req) {
  long long cur_simul_count = 0;

  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int cur_v = nextV(0, realization);

  while (true) {
    Message init_message;
    MPI::COMM_WORLD.Recv(&init_message, 1, message_type, MPI_ANY_SOURCE,
                         MessageType::SIMUL_PREREQUEST);
    int worker_id = init_message.source_id;

    if (cur_simul_count < simulations) {
      cur_simul_count += simul_per_req;
    } else {
      cur_simul_count = simul_per_req;
      cur_v = nextV(cur_v + 1, realization);
    }
    if (cur_v == -1) {
      return;
    } else {
      Message m = {cur_v, 0, simul_per_req};
      MPI::COMM_WORLD.Send(&m, 1, message_type, worker_id,
                           MessageType::SIMUL_REQUEST);
    }
  }
}

vector<double> MPIDirectMC::master(const SourceDetectionParams *params,
                                   bool end, bool print) {

  const long long simulations = params->simulations();
  const common::RealizationRead &snapshot = params->realization();
  int simul_per_req = std::max(10000LL, simulations / 10000);
  assert(simulations % (int)simul_per_req == 0);

  MPIDirectMC *ja = this;
  std::thread steve(&master_steve_wrapper, ja, snapshot.realization(),
                    simulations, simul_per_req);

  const IGraph *graph = params->graph().get();
  int vertices = params->graph()->vertices();
  vector<int> events_resp(vertices, 0);

  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();

  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  while (jobs_remaining > 0) {
    Message received;
    MPI::COMM_WORLD.Recv(&received, 1, message_type, MPI_ANY_SOURCE,
                         MessageType::SIMUL_RESPONSE);
    jobs_remaining -= (int)SIMUL_PER_REQ;
    events_resp[received.source_id] += received.event_outcome;
    if (rand() % 100 == 0) {
      printf("\r%.5f",
             jobs_remaining * 100 /
                 ((double)simulations * snapshot.realization().bitCount()));
      fflush(stdout);
    }
  }

  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    sum += events_resp[v];
  }

  printf("\r\r\n");
  vector<double> p;
  p.clear();
  for (int v = 0; v < vertices; ++v) {
    p.push_back(events_resp[v] / sum);
  }
  steve.join();
  return p;
}

void MPIDirectMC::worker(const SourceDetectionParams *params,
                         ModelType model_type) {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  const IGraph *graph = params->graph().get();
  const common::RealizationRead &snapshot = params->realization();

  // workers
  // Performs simulation on request.
  auto sd = std::unique_ptr<DirectMonteCarloDetector>(
      new DirectMonteCarloDetector(graph));
  while (true) {
    Message message = {rank(), 0, 0};
    MPI::COMM_WORLD.Send(&message, 1, message_type, 0 /* dest */,
                         MessageType::SIMUL_PREREQUEST);
    if (MPI::COMM_WORLD.Iprobe(0, MessageType::SIMUL_REQUEST)) {
      Message message_recv;
      MPI::COMM_WORLD.Recv(&message_recv, 1, message_type, 0 /* source */,
                           MessageType::SIMUL_REQUEST);
      int outcomes = 0;
      for (int t = 0; t < message_recv.batch_size; ++t) {
        common::RealizationRead sp0 = snapshot;
        outcomes += sd->DMCSingleSourceSimulation(message_recv.source_id, sp0,
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
