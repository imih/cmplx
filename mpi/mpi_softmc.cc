#include "mpi_softmc.h"

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
}  // anonymous

namespace cmplx {

void MPISoftMC::benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                                    int benchmark_no, ModelType model_type) {

  if (rank_ == 0) {
    std::string filename =
        "SMbench_SBS_" + std::to_string(benchmark_no) + ".info";
    FILE *f = fopen(filename.c_str(), "a+");
    std::vector<int> sims = {(int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6};

    double a = pow(2, -5);
    params->setA(a);
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

void MPISoftMC::send_simul_end() {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  for (int v = 1; v < processes_; ++v) {
    Message end_message;
    MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                          MessageType::SIMUL_END);
  }
}

vector<double> MPISoftMC::master(const SourceDetectionParams *params, bool end,
                                 bool print) {
  assert(rank_ == 0);
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

    for (int i = 0; i < processes_ - 1; ++i) {
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

std::vector<double> MPISoftMC::convMaster(SourceDetectionParams *params) {
  std::vector<int> sims = {(int)1e4,       2 * (int)1e4,   4 * (int)1e4,
                           10 * (int)1e4,  20 * (int)1e4,  40 * (int)1e4,
                           100 * (int)1e4, 200 * (int)1e4, 400 * (int)1e4,
                           1000 * (int)1e4};

  int convergeG = 0;
  double c = 0.05;
  double a = pow(2, -5);
  int s0 = sims[0];
  printf("s0: %d\n", s0);
  params->setSimulations(s0);
  params->setA(a);
  vector<double> p0 = master(params, false, false);
  int MAP0 = std::max_element(p0.begin(), p0.end()) - p0.begin();
  double pMAP0 = *std::max_element(p0.begin(), p0.end());
  vector<double> P;
  P.clear();
  for (int s_id = 1; s_id < (int)sims.size(); ++s_id) {
    int s1 = sims[s_id];
    params->setSimulations(s1);
    printf("\ns: %d\n", s1);
    vector<double> p1 = master(params, false, false);
    int MAP1 = std::max_element(p1.begin(), p1.end()) - p1.begin();
    double pMAP1 = *std::max_element(p1.begin(), p1.end());
    double delta = dabs(pMAP1 - pMAP0) / pMAP1;
    printf("c: %lf\n", delta);
    double converge = true;
    if (MAP0 != MAP1) converge = false;
    if (delta > c) converge = false;
    int pos = 0;
    for (int j = 0; j < (int)p1.size(); ++j) {
      if (dabs(p1[j] - p0[j]) > c) converge = false;
      if (p1[j] > 0) pos++;
    }
    if (pos == 0) converge = false;
    if (converge) {
      convergeG++;
      printf("Converged for n=%d\n", s0);
      P = p0;
      params->setSimulations(s0);
      if (convergeG > 0) break;
    } else {
      convergeG = 0;
      printf("Not converged.\n");
    }
    s0 = s1;
    p0.assign(p1.begin(), p1.end());
    pMAP0 = pMAP1;
    MAP0 = MAP1;
    P = p1;
  }
  return P;
}

vector<double> MPISoftMC::softConvMaster(cmplx::SourceDetectionParams *params,
                                         bool end) {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int convergeG = 0;
  std::vector<int> sims = {(int)1e4,       2 * (int)1e4,   4 * (int)1e4,
                           10 * (int)1e4,  20 * (int)1e4,  40 * (int)1e4,
                           100 * (int)1e4, 200 * (int)1e4, 400 * (int)1e4,
                           800 * (int)1e4};

  double c = 0.05;
  const int MAXA = 9;
  int s0 = sims[0];
  printf("s0: %d\n", s0);
  vector<double> a(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    a[i] = 1.0 / (double)(1 << i);
  }
  vector<double> p0[MAXA + 1];
  vector<double> pMAP0(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    params->setSimulations(s0);
    params->setA(a[i]);
    printf("a[i]: %lf\n", a[i]);
    p0[i] = master(params, false, false);
    pMAP0[i] = *std::max_element(p0[i].begin(), p0[i].end());
  }

  vector<int> convergeGlobal(MAXA + 1, 0);
  vector<double> P;
  int bits = params->realization().realization().bitCount();
  for (int s = 1; s < (int)sims.size(); ++s) {
    int s1 = sims[s];
    printf("\ns: %d\n", s1);
    params->setSimulations(s1);
    vector<double> p1[MAXA + 1];
    vector<double> pMAP1(MAXA + 1, 0);

    for (int i = MAXA; i >= 3; --i) {
      printf("\ns: %d a: %.10lf\n", s1, a[i]);
      params->setA(a[i]);
      double converge = true;
      p1[i] = master(params, false, false);
      pMAP1[i] = *std::max_element(p1[i].begin(), p1[i].end());
      double delta = dabs(pMAP1[i] - pMAP0[i]) / pMAP1[i];
      printf("c: %lf\n", delta);
      if (delta > c) converge = false;
      int pos = 0;
      for (int j = 0; j < (int)p1[i].size(); ++j) {
        if (dabs(p1[i][j] - p0[i][j]) > c) converge = false;
        if (p1[i][j] > 0) pos++;
      }
      if (pos == 0) converge = false;
      if (converge) {
        convergeGlobal[i]++;
        printf("Converged for n=%d a=%lf\n", s0, a[i]);
        if (convergeGlobal[i] > 0) break;
      } else {
        convergeGlobal[i] = 0;
        printf("Not converged.\n");
      }
    }

    bool done = false;
    for (int i = MAXA; i >= 3; --i) {
      if (convergeGlobal[i] > 0) {
        P = p0[i];
        params->setSimulations(s0);
        params->setA(a[i]);
        if (end) {
          send_simul_end();
        }
        done = true;
        break;
      }
    }

    s0 = s1;
    for (int i = 3; i <= MAXA; ++i) p0[i] = p1[i];
    pMAP0 = pMAP1;
    if (done) break;
    P = p1[MAXA];
  }
  return P;
}

void MPISoftMC::worker(const SourceDetectionParams *params,
                       ModelType model_type) {
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int rank = MPI::COMM_WORLD.Get_rank();
  assert(rank > 0);

  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  common::RealizationRead snapshot = params->realization();

  // workers
  // Performs simulation on request.
  SoftMarginDetector sd(graph);

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
        fi.push_back(sd.SMSingleSourceSimulation(message_recv.source_id, sp0,
                                                 model_type));
      }
      message_recv.event_outcome = sd.likelihood(fi, message_recv.a);
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
