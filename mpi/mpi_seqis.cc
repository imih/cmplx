#include "mpi_seqis.h"

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
}  // anonymous

namespace cmplx {

MPISeqIS::MPISeqIS() : MpiParal() {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
}

void MPISeqIS::benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                                   int benchmark_no, ModelType model_type) {
  MPI::Datatype message_type = datatypeOfMessage();

  if (rank_ == 0) {
    std::string filename =
        "SEQSoftbench_SBS_" + std::to_string(benchmark_no) + ".info";
    FILE *f = fopen(filename.c_str(), "w+");
    std::vector<int> sims = {(int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6};
    for (int sim : sims) {
      fprintf(f, "s: %d\n", sim);
      params->setSimulations(sim);
      vector<double> p = master(params, false, false);
      for (int i = 0; i < (int)p.size(); ++i)
        fprintf(f, "%.10lf%c", p[i], i + 1 == (int)p.size() ? '\n' : ' ');
      fflush(f);
    }
    fclose(f);

    for (int v = 1; v < processes_; ++v) {
      Message end_message;
      MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                            MessageType::SIMUL_END);
    }
  } else {
    worker(params, model_type);
  }
  exit(0);
}

void MPISeqIS::send_simul_end() {
  MPI::Datatype message_type = datatypeOfMessage();
  for (int v = 1; v < processes_; ++v) {
    Message end_message;
    MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                          MessageType::SIMUL_END);
  }
}

vector<double> MPISeqIS::convMaster(cmplx::SourceDetectionParams *params) {
  assert(rank_ == 0);
  std::vector<int> sims = {(int)1e4,      2 * (int)1e4,  4 * (int)1e4,
                           8 * (int)1e4,  10 * (int)1e4, 20 * (int)1e4,
                           40 * (int)1e4, 80 * (int)1e4, 100 * (int)1e4};
  int s0 = sims[0];
  printf("s0: %d\n", s0);

  params->setSimulations(s0);
  vector<double> p0 = master(params, false, true);
  double pMAP0 = *std::max_element(p0.begin(), p0.end());
  double c = 0.05;

  vector<double> res;
  res.clear();
  int bits = params->realization().realization().bitCount();
  int convergeG = 0;
  int s_id = 1;
  for (int s = 1; s < (int)sims.size(); ++s) {
    int s1 = sims[s];
    printf("s2: %d\n", s1);
    params->setSimulations(s1);

    bool converge = true;
    vector<double> p1 = master(params, false, false);
    double pMAP1 = *std::max_element(p1.begin(), p1.end());
    double delta = dabs(pMAP1 - pMAP0) / pMAP1;
    printf("\rc: %lf ", delta);
    if (delta >= c) converge = false;
    int pos = 0;
    for (int j = 0; j < (int)p1.size(); ++j) {
      if (p1[j] > 0 && (dabs(p1[j] - p0[j]) > c)) converge = false;
      if (p1[j] > 0) {
        pos++;
      }
    }
    if (pos == 0) converge = false;
    printf("\n");
    if (converge)
      convergeG++;
    else
      convergeG = 0;
    if (convergeG > 0) {
      printf("Converged for n=%d\n", s0);
      res = p0;
      /*
      if (end) {
        send_simul_end();
      }
      */
      break;
    } else {
      printf("Not converged.\n");
    }

    s0 = s1;
    p0.clear();
    res = p1;
    p0.assign(p1.begin(), p1.end());
    pMAP0 = pMAP1;
  }
  params->setSimulations(s0);
  return res;
}

vector<double> MPISeqIS::master(const SourceDetectionParams *params, bool end,
                                bool print) {
  MPI::Datatype message_type = datatypeOfMessage();

  assert(rank_ == 0);

  double p = params->realization().p();
  double q = params->realization().q();
  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  const RealizationRead &snapshot = params->realization();

  // master process
  int cur_v = nextV(0, snapshot.realization());

  vector<double> events_resp(vertices, 0);
  long long jobs_remaining = 1LL * snapshot.realization().bitCount();
  while (jobs_remaining > 0) {
    for (int i = 0; i < processes_ - 1; ++i) {
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

    for (int i = 0; i < processes_ - 1; ++i) {
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

  printf("\r\n");
  /*****/
  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    sum += events_resp[v];
  }

  for (int v = 0; v < vertices; ++v) {
    if (sum > 0) events_resp[v] /= sum;
  }
  return events_resp;
}

void MPISeqIS::worker(const SourceDetectionParams *params,
                      ModelType model_type) {
  MPI::Datatype message_type = datatypeOfMessage();
  assert(rank_ > 0);

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
