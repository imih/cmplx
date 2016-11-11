#include "mpi_softmc.h"

#include "../common/bit_array.h"
#include "../common/igraph.h"
#include "../common/realization.h"

#include "mpi_common.h"

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
using std::vector;

namespace {
const int SIMUL_PER_REQ = 10000;

struct Message {
  int source_id;
  double event_outcome;
  double a;
  Message() : Message(-1, -1, -1) {}
  Message(int source_id, double event_outcome, double a)
      : source_id(source_id), event_outcome(event_outcome), a(a) {}
};

MPI::Datatype datatypeOfMessage() {
  int blockLen[3] = {1, 1, 1};
  MPI::Aint offsets[3] = {offsetof(Message, source_id),
                          offsetof(Message, event_outcome),
                          offsetof(Message, a)};
  MPI::Datatype types[3] = {MPI::INT, MPI::DOUBLE, MPI::DOUBLE};
  return MPI::Datatype::Create_struct(3, blockLen, offsets, types);
}

typedef long long ll;

std::vector<double> responseToProb(
    const std::vector<double> &events_resp_sum,
    const std::vector<long long> &events_resp_size, int vertices) {
  printf("\r\n");
  /*****/
  std::vector<double> P;
  P.clear();
  double sum = 0;

  int first_size = 0;
  bool normalize = false;
  for (int v = 0; v < vertices; ++v) {
    if (!first_size) first_size = events_resp_size[0];
    if (first_size && events_resp_size[v] &&
        events_resp_size[v] != first_size) {
      std::cout << std::endl << events_resp_size[v] << " " << first_size
                << std::endl;
      puts("PROBLEM!");
      normalize = true;
      break;
    }
  }

  for (int v = 0; v < vertices; ++v) {
    P.push_back((normalize && events_resp_size[v])
                    ? (events_resp_sum[v] / events_resp_size[v])
                    : events_resp_sum[v]);
    sum += P.back();
    if (P.back() > 0.01) printf("%.2lf ", P.back());
  }
  printf("\n");
  for (int v = 0; v < vertices; ++v) {
    if (sum > 0) P[v] /= sum;
  }
  return P;
}
}  // anonymous

namespace cmplx {

void MPISoftMC::send_simul_end() {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int processes = MpiMaster::processes();
  for (int v = 1; v < processes; ++v) {
    Message end_message;
    MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                          MessageType::SIMUL_END);
  }
}

vector<double> MPISoftMC::master(const SourceDetectionParams *params, bool end,
                                 bool print) {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  double p = params->realization().p();
  double q = params->realization().q();
  const int simulations = params->simulations();
  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  const common::RealizationRead &snapshot = params->realization();

  // master process
  int cur_simul_count = 0;
  int cur_v = 0;
  while ((cur_v < vertices) && (snapshot.realization().bit(cur_v) == false))
    cur_v++;
  vector<double> events_resp_sum(vertices, 0);
  vector<long long> events_resp_size(vertices, 0);
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();
  assert(jobs_remaining % (int)SIMUL_PER_REQ == 0);

  int processes = MpiMaster::processes();
  while (jobs_remaining > 0) {
    for (int i = 0; i < processes - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_PREREQUEST)) {
        Message init_message;
        MPI::COMM_WORLD.Recv(&init_message, 1, message_type, i + 1,
                             MessageType::SIMUL_PREREQUEST);
        if (cur_simul_count < simulations) {
          cur_simul_count += (int)SIMUL_PER_REQ;
        } else {
          cur_simul_count = (int)SIMUL_PER_REQ;
          cur_v++;
          while ((cur_v < vertices) &&
                 (snapshot.realization().bit(cur_v) == false)) {
            cur_v++;
          }
          if (cur_v >= vertices) {
            cur_v = -1;
          }
        }
        if (cur_v == -1) {
          if (end) {
            // There's no more jobs
            Message message;
            MPI::COMM_WORLD.Isend(&message, 1, message_type, i + 1,
                                  MessageType::SIMUL_END);
          }
        } else {
          Message m = Message(cur_v, 0, params->a());
          MPI::COMM_WORLD.Isend(&m, 1, message_type, i + 1,
                                MessageType::SIMUL_REQUEST);
        }
      }
    }

    for (int i = 0; i < processes - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_RESPONSE)) {
        Message received;
        MPI::COMM_WORLD.Recv(&received, 1, message_type, i + 1,
                             MessageType::SIMUL_RESPONSE);
        jobs_remaining -= (int)SIMUL_PER_REQ;
        events_resp_sum[received.source_id] += received.event_outcome;
        events_resp_size[received.source_id]++;
        if (!i) {
          printf("\r%.5f",
                 jobs_remaining * 100 /
                     ((double)simulations * snapshot.realization().bitCount()));
          fflush(stdout);
        }
      }
    }
  }

  return responseToProb(events_resp_sum, events_resp_size, vertices);
}

void MPISoftMC::worker(const SourceDetectionParams *params,
                       ModelType model_type) {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int rank = MPI::COMM_WORLD.Get_rank();

  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  common::RealizationRead snapshot = params->realization();

  // workers
  // Performs simulation on request.
  auto sd = std::unique_ptr<SoftMarginDetector>(new SoftMarginDetector(graph));

  while (true) {
    Message message;
    MPI::COMM_WORLD.Isend(&message, 1, message_type, 0 /* dest */,
                          MessageType::SIMUL_PREREQUEST);
    if (MPI::COMM_WORLD.Iprobe(0, MessageType::SIMUL_REQUEST)) {
      Message message_recv;
      MPI::COMM_WORLD.Recv(&message_recv, 1, message_type, 0 /* source */,
                           MessageType::SIMUL_REQUEST);

      /***/
      vector<double> fi;
      fi.clear();
      for (int t = 0; t < (int)SIMUL_PER_REQ; ++t) {
        RealizationRead sp0 = snapshot;
        fi.push_back(sd->SMSingleSourceSimulation(message_recv.source_id, sp0,
                                                  model_type));
      }
      message_recv.event_outcome = sd->likelihood(fi, message_recv.a);
      /****/

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
