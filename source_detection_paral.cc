#include "source_detection_paral.h"

#include <mpi.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>

#include <cstdio>
#include <cstdlib>

#include <vector>

#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::SourceDetectionParams;
using std::vector;

namespace {
const int SIMUL_PER_REQ = 10000;

enum MessageType {
  SIMUL_PREREQUEST,
  SIMUL_REQUEST,
  SIMUL_RESPONSE,
  SIMUL_END,
  SIMUL_PARAMS
};

double dabs(double x) {
  if (x < 0) return x * -1;
  return x;
}
}  // namespace

namespace DMC {

struct Message {
  int source_id;
  int event_outcome;
};

MPI::Datatype datatypeOfMessage() {
  int blockLen[2] = {1, 1};
  MPI::Aint offsets[2] = {offsetof(Message, source_id),
                          offsetof(Message, event_outcome)};
  MPI::Datatype types[2] = {MPI::INT, MPI::INT};
  return MPI::Datatype::Create_struct(2, blockLen, offsets, types);
}
}  // namespace DMC

namespace cmplx {

void DirectMCSimulParalWorker(const SourceDetectionParams &, ModelType);
vector<double> DirectMCSimulParalMaster(const SourceDetectionParams &, bool,
                                        bool);

void DirectMCSimulParal(const SourceDetectionParams &params,
                        ModelType model_type) {
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank == 0) {
    DirectMCSimulParalMaster(params, true, true);
  } else {
    DirectMCSimulParalWorker(params, model_type);
  }
}

void DirectMCSimulParalConv(const SourceDetectionParams &params,
                            ModelType model_type) {
  int rank = MPI::COMM_WORLD.Get_rank();
  int processes = MPI::COMM_WORLD.Get_size();
  if (rank != 0) {
    DirectMCSimulParalWorker(params, model_type);
  } else {
    using namespace DMC;
    MPI::Datatype message_type = datatypeOfMessage();
    message_type.Commit();

    SourceDetectionParams params0(params);
    /***************  */
    double c = 0.05;  //
    /**************   */
    int s0 = SIMUL_PER_REQ;
    double pml0 = 0;
    params0.setSimulations(s0);
    vector<double> p0 = DirectMCSimulParalMaster(params0, false, false);
    pml0 = *std::max_element(p0.begin(), p0.end());
    while (true) {
      int s1 = 2 * s0;
      printf("s: %d\n", s1);
      params0.setSimulations(s1);
      vector<double> p1 = DirectMCSimulParalMaster(params0, false, false);

      double pml1 = *std::max_element(p0.begin(), p0.end());
      if (isnan(pml1)) {
        printf("NAN! %d \n", s1);
      }
      bool converge = true;
      double delta = dabs(pml1 - pml0) / pml1;
      printf("c: %lf\n", delta);
      if (delta >= c) converge = false;
      for (int i = 0; i < (int)p1.size(); ++i) {
        if (dabs(p0[i] - p1[i]) >= c) converge = false;
      }

      if (converge) {
        for (int v = 1; v < processes; ++v) {
          Message end_message;
          MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                                MessageType::SIMUL_END);
        }
        break;
      }
      pml0 = pml1;
      s0 = s1;
      p0 = p1;
    }
  }
}

void DirectMCSimulParalWorker(const SourceDetectionParams &params,
                              ModelType model_type) {
  using namespace DMC;
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  assert(rank > 0);

  int vertices = params.graph().vertices();
  const IGraph &graph = params.graph();
  const Realization &snapshot = params.realization();

  // workers
  // Performs simulation on request.
  SourceDetector sd(graph);
  while (true) {
    Message message = {-1, 0};
    MPI::COMM_WORLD.Send(&message, 1, message_type, 0 /* dest */,
                         MessageType::SIMUL_PREREQUEST);
    if (MPI::COMM_WORLD.Iprobe(0, MessageType::SIMUL_REQUEST)) {
      Message message_recv;
      MPI::COMM_WORLD.Recv(&message_recv, 1, message_type, 0 /* source */,
                           MessageType::SIMUL_REQUEST);
      int outcomes = 0;
      for (int t = 0; t < SIMUL_PER_REQ; ++t) {
        Realization sp0 = snapshot;
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

vector<double> DirectMCSimulParalMaster(const SourceDetectionParams &params,
                                        bool end = true, bool print = true) {
  using namespace DMC;
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  assert(rank == 0);

  const int simulations = params.simulations();

  int vertices = params.graph().vertices();
  const IGraph &graph = params.graph();
  const Realization &snapshot = params.realization();
  int cur_simul_count = 0;
  int cur_v = 0;
  while ((cur_v < vertices) && (snapshot.realization().bit(cur_v) == false))
    cur_v++;
  vector<int> events_resp(vertices, 0);
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();
  assert(simulations % SIMUL_PER_REQ == 0);
  while (jobs_remaining > 0) {
    for (int i = 0; i < processes - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_PREREQUEST)) {
        Message init_message;
        MPI::COMM_WORLD.Recv(&init_message, 1, message_type, i + 1,
                             MessageType::SIMUL_PREREQUEST);
        if (cur_simul_count < simulations) {
          cur_simul_count += SIMUL_PER_REQ;
        } else {
          cur_simul_count = SIMUL_PER_REQ;
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
          // There's no more jobs
          if (end == true) {
            Message message = {-1, -1};
            MPI::COMM_WORLD.Send(&message, 1, message_type, i + 1,
                                 MessageType::SIMUL_END);
          } else {
            // do nothing
            // MPI::COMM_WORLD.Send(&message, 1, message_type, i + 1,
            // MessageType::SIMUL_REQUEST);
          }
        } else {
          Message m = {cur_v, 0};
          MPI::COMM_WORLD.Send(&m, 1, message_type, i + 1,
                               MessageType::SIMUL_REQUEST);
        }
      }
    }

    for (int i = 0; i < processes - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_RESPONSE)) {
        Message received;
        MPI::COMM_WORLD.Recv(&received, 1, message_type, i + 1,
                             MessageType::SIMUL_RESPONSE);
        jobs_remaining -= SIMUL_PER_REQ;
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
  for (int v = 0; v < vertices; ++v) {
    if (print) printf("%.10f\n", events_resp[v] / sum);
    p.push_back(events_resp[v] / sum);
  }
  return p;
}

/******* SMP *******/
namespace SMP {
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
}  // namespace SMP

void SoftMarginSimulParalWorker(const SourceDetectionParams &, ModelType);
vector<double> SoftMarginSimulParalMaster(const SourceDetectionParams &, bool,
                                          bool);

void SoftMarginParal(const SourceDetectionParams &params,
                     ModelType model_type) {
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank == 0) {
    SoftMarginSimulParalMaster(params, true, true);
  } else {
    // workers
    SoftMarginSimulParalWorker(params, model_type);
  }
}

vector<double> SoftMarginParalConvMaster(
    const cmplx::SourceDetectionParams &params, bool end) {
  using namespace SMP;
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  SourceDetectionParams params0(params);
  double c = 0.05;
  // SourceDetectionParams& params;
  const int MAXA = 15;
  int s0 = SIMUL_PER_REQ;
  printf("s0: %d\n", s0);
  vector<double> a(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    a[i] = 1.0 / (double)(1 << i);
  }
  vector<double> p0[MAXA + 1];
  vector<double> pMAP0(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    params0.setSimulations(s0);
    params0.setA(a[i]);
    printf("a[i]: %lf\n", a[i]);
    p0[i] = SoftMarginSimulParalMaster(params0, false, false);
    pMAP0[i] = *std::max_element(p0[i].begin(), p0[i].end());
  }

  vector<int> convergeGlobal(MAXA + 1, 0);
  vector<double> res;
  int bits = params.realization().realization().bitCount();
  while (true) {
    int s1 = 2 * s0;
    printf("s: %d\n", s1);
    params0.setSimulations(s1);
    vector<double> p1[MAXA + 1];
    vector<double> pMAP1(MAXA + 1, 0);

    for (int i = 15; i >= 3; --i) {
      printf("s: %d a: %.10lf\n", s1, a[i]);
      params0.setA(a[i]);
      double converge = true;
      p1[i] = SoftMarginSimulParalMaster(params0, false, false);
      pMAP1[i] = *std::max_element(p1[i].begin(), p1[i].end());
      double delta = dabs(pMAP1[i] - pMAP0[i]) / pMAP1[i];
      printf("c: %lf\n", delta);
      if (delta >= c) converge = false;
      int pos = 0;
      for (int j = 0; j < (int)p1[i].size(); ++j) {
        if (dabs(p1[i][j] - p0[i][j]) >= c) converge = false;
        if (p1[i][j] > 0) pos++;
      }
      if (s1 > 1000000) converge = true;
      if (pos != bits) converge = false;
      if (converge) {
        convergeGlobal[i]++;
        printf("Converged for n=%d a=%lf\n", s1, a[i]);
        if (convergeGlobal[i]) break;
      } else {
        convergeGlobal[i] = 0;
        printf("Not converged.\n");
      }
    }

    bool done = false;
    for (int i = 15; i >= 3; --i) {
      if (convergeGlobal[i]) {
        res = p1[i];
        MPI::COMM_WORLD.Get_size();
        if (end) {
          int processes = MPI::COMM_WORLD.Get_size();
          for (int v = 1; v < processes; ++v) {
            Message end_message;
            MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                                  MessageType::SIMUL_END);
          }
        }
        done = true;
        break;
      }
    }

    s0 = s1;
    for (int i = 1; i <= 15; ++i) p0[i] = p1[i];
    pMAP0 = pMAP1;
    if (done) break;
  }
  return res;
}

void SoftMarginParalConv(const SourceDetectionParams &params,
                         ModelType model_type) {
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank != 0) {
    SoftMarginSimulParalWorker(params, model_type);
  } else {
    SoftMarginParalConvMaster(params);
  }
}

vector<double> SoftMarginSimulParalMaster(const SourceDetectionParams &params,
                                          bool end = true, bool print = true) {
  using namespace SMP;
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  assert(rank == 0);

  double p = params.realization().p();
  double q = params.realization().q();
  const int simulations = params.simulations();
  int vertices = params.graph().vertices();
  const IGraph &graph = params.graph();
  const Realization &snapshot = params.realization();

  // master process
  int cur_simul_count = 0;
  int cur_v = 0;
  while ((cur_v < vertices) && (snapshot.realization().bit(cur_v) == false))
    cur_v++;
  vector<vector<double> > events_resp(vertices, vector<double>());
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();
  while (jobs_remaining > 0) {
    for (int i = 0; i < processes - 1; ++i) {
      if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_PREREQUEST)) {
        Message init_message;
        MPI::COMM_WORLD.Recv(&init_message, 1, message_type, i + 1,
                             MessageType::SIMUL_PREREQUEST);
        if (cur_simul_count < simulations) {
          cur_simul_count += SIMUL_PER_REQ;
        } else {
          cur_simul_count = SIMUL_PER_REQ;
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
          Message m = Message(cur_v, 0, params.a());
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
        jobs_remaining -= SIMUL_PER_REQ;
        events_resp[received.source_id].push_back(received.event_outcome);
        if (!i) {
          printf("\r%.5f",
                 jobs_remaining * 100 /
                     ((double)simulations * snapshot.realization().bitCount()));
          fflush(stdout);
        }
      }
    }
  }

  // fprintf(stderr, "Simulations finished");
  printf("\r\n");
  /*****/
  std::vector<double> P;
  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    double P_v = 0;
    for (double d : events_resp[v]) {
      //        std::cout << d << std::endl;
      P_v += d;
    }
    if (events_resp[v].size()) P_v /= (int)events_resp[v].size();
    P.push_back(P_v);
    sum += P_v;
  }
  // assert(sum > 0);
  // fprintf(file, "\n\n%.10lf %.10lf\n\n", snapshot.p(), snapshot.q());
  if (print) params.realization().print();
  for (int v = 0; v < vertices; ++v) {
    // printf("%.10lf\n", P[v])
    if (sum > 0) P[v] /= sum;
    if (print) printf("%.10lf\n", P[v]);
    // fprintf(file, "%.10lf ", P[v]);
  }
  printf("\n");
  // fprintf(file, "\n");
  /******/
  //  fclose(file);
  return P;
}

void SoftMarginSimulParalWorker(const SourceDetectionParams &params,
                                ModelType model_type) {
  using namespace SMP;
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int rank = MPI::COMM_WORLD.Get_rank();
  assert(rank > 0);

  int vertices = params.graph().vertices();
  const IGraph &graph = params.graph();
  Realization snapshot = params.realization();

  // workers
  // Performs simulation on request.
  SourceDetector sd(graph);

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
      for (int t = 0; t < SIMUL_PER_REQ; ++t) {
        Realization sp0 = snapshot;
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

void GenerateSoftMarginDistributions(const SourceDetectionParams &params,
                                     int distributions, ModelType model_type) {
  int rank = MPI::COMM_WORLD.Get_rank();
  int processes = MPI::COMM_WORLD.Get_size();
  std::string filename = "distr_" + params.summary();
  FILE *f = fopen(filename.c_str(), "a");
  SourceDetectionParams params0 = params;

  for (int d = 0; d < distributions; ++d) {
    MPI::COMM_WORLD.Barrier();

    if (rank == 0) {
      MPI::COMM_WORLD.Barrier();
      SourceDetectionParams params_novi = SourceDetectionParams::ParamsFromGrid(
          params.realization().p(), params.realization().q(),
          (int)sqrt(params.realization().population_size()));
      vector<int> r_pos = params_novi.realization().realization().positions();
      for (int v = 1; v < processes; ++v) {
        MPI::COMM_WORLD.Send(&r_pos[0], (int)r_pos.size(), MPI_INT, v,
                             MessageType::SIMUL_PARAMS);
      }
      params0.setRealization(params_novi.realization().realization());
    } else {
      vector<int> r_pos;
      r_pos.resize(params.graph().vertices());
      MPI::COMM_WORLD.Recv(&r_pos[0], (int)r_pos.size(), MPI_INT, 0,
                           MessageType::SIMUL_PARAMS);
      BitArray r_ba(params.graph().vertices());
      for (int p : r_pos) r_ba.set(p, true);
      params0.setRealization(r_ba);
    }

    MPI::COMM_WORLD.Barrier();

    if (rank == 0) {
      std::vector<double> P = SoftMarginParalConvMaster(params0, true);

      for (int j = 0; j < (int)P.size(); ++j) {
        fprintf(f, "%.10lf%c", P[j], j == ((int)P.size() - 1) ? '\n' : ' ');
        printf("%.10lf\n", P[j]);
      }
      fflush(f);

      using namespace SMP;
      MPI::Datatype message_type = datatypeOfMessage();
      message_type.Commit();
      for (int v = 1; v < processes; ++v) {
        Message end_message;
        MPI::COMM_WORLD.Isend(&end_message, 1, message_type, v,
                              MessageType::SIMUL_END);
      }
    } else {
      SoftMarginSimulParalWorker(params0, model_type);
    }
  }
  fclose(f);
  exit(0);
}

// estimates the posterior probabilty of full match
void FEstimatorParal(const SourceDetectionParams &params,
                     ModelType model_type) {
  using namespace DMC;
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  double p = params.realization().p();
  double q = params.realization().q();
  const int simulations = params.simulations();
  int vertices = params.graph().vertices();
  const IGraph &graph = params.graph();
  const Realization &snapshot = params.realization();

  // TODO(iva) Make only one copy!
  // snapshot.print();

  if (rank == 0) {
    // master process
    int cur_simul_count = 0;
    int cur_v = 0;
    while ((cur_v < vertices) && (snapshot.realization().bit(cur_v) == false))
      cur_v++;
    vector<double> events_resp(vertices, 0);
    long long jobs_remaining =
        1LL * simulations * snapshot.realization().bitCount();
    while (jobs_remaining > 0) {
      for (int i = 0; i < processes - 1; ++i) {
        if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_PREREQUEST)) {
          Message init_message;
          MPI::COMM_WORLD.Recv(&init_message, 1, message_type, i + 1,
                               MessageType::SIMUL_PREREQUEST);
          if (cur_simul_count < simulations) {
            cur_simul_count += SIMUL_PER_REQ;
          } else {
            cur_simul_count = SIMUL_PER_REQ;
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
            puts("No more jobs");
            // There's no more jobs
            Message message = {-1, -1};
            MPI::COMM_WORLD.Isend(&message, 1, message_type, i + 1,
                                  MessageType::SIMUL_REQUEST);
          } else {
            Message m = {cur_v, 0};
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
          jobs_remaining -= SIMUL_PER_REQ;
          /***/
          double prev = events_resp[received.source_id];
          events_resp[received.source_id] +=
              received.event_outcome / (double)SIMUL_PER_REQ;
          events_resp[received.source_id] /= 2;
          // printf("%.10lf\n", (events_resp[received.source_id] - prev));
          /***/
          if (!i)
            printf("%.5f\n",
                   jobs_remaining * 100 / ((double)simulations *
                                           snapshot.realization().bitCount()));
        }
      }
    }
    for (int i = 0; i < processes - 1; ++i) {
      Message m = {-1, -1};
      MPI::COMM_WORLD.Isend(&m, 1, message_type, i + 1,
                            MessageType::SIMUL_REQUEST);
    }

    fprintf(stderr, "Simulations finished");
    for (int v = 0; v < graph.vertices(); ++v)
      printf("%.10lf\n", events_resp[v]);
    /******/
  } else {
    // workers
    // Performs simulation on request.
    SourceDetector sd(graph);
    while (true) {
      Message message = {-1, 0};
      MPI::COMM_WORLD.Send(&message, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_PREREQUEST);
      Message message_recv;
      MPI::COMM_WORLD.Recv(&message_recv, 1, message_type, 0 /* source */,
                           MessageType::SIMUL_REQUEST);
      if (message_recv.source_id == -1) {
        break;
      }

      int fi = 0;
      for (int t = 0; t < SIMUL_PER_REQ; ++t) {
        Realization sp0 = snapshot;
        double jac = sd.SMSingleSourceSimulation(message_recv.source_id, sp0,
                                                 model_type);
        if (jac == 1.0) fi++;
      }

      message_recv.event_outcome = fi;
      Message toSend = message_recv;
      MPI::COMM_WORLD.Send(&toSend, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_RESPONSE);
    }
  }
}

}  // namespace cmplx
