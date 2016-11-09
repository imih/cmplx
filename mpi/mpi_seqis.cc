#include "mpi_seqis.h"

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
  int sample_size;
  Message() : Message(-1, -1, -1) {}
  Message(int source_id, double event_outcome, int sample_size)
      : source_id(source_id),
        event_outcome(event_outcome),
        sample_size(sample_size) {}
};

MPI::Datatype datatypeOfMessage() {
  int blockLen[3] = {1, 1, 1};
  MPI::Aint offsets[3] = {offsetof(Message, source_id),
                          offsetof(Message, event_outcome),
                          offsetof(Message, sample_size)};
  MPI::Datatype types[3] = {MPI::INT, MPI::DOUBLE, MPI::INT};
  return MPI::Datatype::Create_struct(3, blockLen, offsets, types);
}

std::vector<double> responseToProb(const std::vector<double> &events_resp,
                                   int vertices) {
  printf("\r\n");
  /*****/
  std::vector<double> P = events_resp;

  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    sum += P[v];
  }

  for (int v = 0; v < vertices; ++v) {
    if (sum > 0) P[v] /= sum;
  }
  return P;
}
}  // anonymous

namespace cmplx {

void MPISeqIS::send_simul_end() {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int processes = MpiMaster::processes();
  for (int v = 1; v < processes; ++v) {
    Message end_message;
    MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                          MessageType::SIMUL_END);
  }
}

vector<double> MPISeqIS::master(const SourceDetectionParams *params, bool end,
                                bool print) {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  double p = params->realization().p();
  double q = params->realization().q();
  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  const RealizationRead &snapshot = params->realization();

  // master process
  int cur_v = nextV(0, snapshot.realization());

  vector<double> events_resp(vertices, 0);
  long long jobs_remaining = 1LL * snapshot.realization().bitCount();
  int processes = MpiMaster::processes();
  while (jobs_remaining > 0) {
    for (int i = 0; i < processes - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_PREREQUEST)) {
        Message init_message;
        MPI::COMM_WORLD.Recv(&init_message, 1, message_type, i + 1,
                             MessageType::SIMUL_PREREQUEST);
        if (cur_v == -1) {
          if (end) {
            // There's no more jobs
            Message message;
            MPI::COMM_WORLD.Isend(&message, 1, message_type, i + 1,
                                  MessageType::SIMUL_END);
          }
        } else {
          Message m = Message(cur_v, 0, params->simulations());
          MPI::COMM_WORLD.Isend(&m, 1, message_type, i + 1,
                                MessageType::SIMUL_REQUEST);
          cur_v = nextV(cur_v + 1, snapshot.realization());
        }
      }
    }

    for (int i = 0; i < processes - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_RESPONSE)) {
        Message received;
        MPI::COMM_WORLD.Recv(&received, 1, message_type, i + 1,
                             MessageType::SIMUL_RESPONSE);
        jobs_remaining--;
        events_resp[received.source_id] = received.event_outcome;
        if (!i) {
          printf("\r%.5f", jobs_remaining * 100 /
                               ((double)snapshot.realization().bitCount()));
          fflush(stdout);
        }
      }
    }
  }

  return responseToProb(events_resp, vertices);
}

void MPISeqIS::worker(const SourceDetectionParams *params,
                      ModelType model_type) {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  RealizationRead snapshot = params->realization();

  // Performs simulation on request.
  SequentialMCDetector sd(graph);
  // SequentialSoftMCDetector sd(graph);

  while (true) {
    Message message;
    MPI::COMM_WORLD.Isend(&message, 1, message_type, 0 /* dest */,
                          MessageType::SIMUL_PREREQUEST);
    if (MPI::COMM_WORLD.Iprobe(0, MessageType::SIMUL_REQUEST)) {
      Message message_recv;
      MPI::COMM_WORLD.Recv(&message_recv, 1, message_type, 0 /* source */,
                           MessageType::SIMUL_REQUEST);

      /***/
      int sample_size = message_recv.sample_size;
      double Pos =
          sd.seqPosterior(message_recv.source_id, sample_size, snapshot,
                          cmplx::ResamplingType::SIMPLE_RANDOM_SAMPLING,
                          true /* p = 1 @ T = 5*/);
      message_recv.event_outcome = Pos;
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
    sleep(1);
  }
}

}  // namespace cmplx
