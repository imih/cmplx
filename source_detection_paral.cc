#include "source_detection_paral.h"

#include <mpi.h>
#include <unistd.h>

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
  SIMUL_RESPONSE
};
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

void DirectMCSimulParal(SourceDetectionParams params) {
  using namespace DMC;
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();
  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();

  const int simulations = params.simulations();

  int vertices = params.graph().vertices();
  const IGraph &graph = params.graph();
  const Realization &snapshot = params.realization();
  // TODO(iva) Make only one copy!

  if (rank == 0) {
    // master process
    int cur_simul_count = 0;
    int cur_v = 0;
    while ((cur_v < vertices) && (snapshot.realization().bit(cur_v) == false))
      cur_v++;
    vector<int> events_resp(vertices, 0);
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
          events_resp[received.source_id] += received.event_outcome;
          if (!i)
            printf("%.5f\n",
                   jobs_remaining * 100 / ((double)simulations *
                                           snapshot.realization().bitCount()));
        }
      }
    }

    fprintf(stderr, "Simulations finished");
    double sum = 0;
    for (int v = 0; v < vertices; ++v) {
      sum += events_resp[v];
    }

    for (int v = 0; v < vertices; ++v) {
      printf("%.10f\n", events_resp[v] / sum);
    }
    printf("\n");
    for (int i = 0; i < processes - 1; ++i) {
      Message m = {-1, -1};
      MPI::COMM_WORLD.Isend(&m, 1, message_type, i + 1,
                            MessageType::SIMUL_REQUEST);
    }
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
      int outcomes = 0;
      for (int t = 0; t < SIMUL_PER_REQ; ++t) {
        Realization sp0 = snapshot;
        outcomes +=
            sd.DMCSingleSourceSirSimulation(message_recv.source_id, sp0);
      }

      message_recv.event_outcome = outcomes;

      Message toSend = message_recv;
      MPI::COMM_WORLD.Send(&toSend, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_RESPONSE);
    }
  }
}

// estimates the posterior probabilty of full match
void FEstimatorParal(SourceDetectionParams params) {
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
          printf("%.10lf\n", (events_resp[received.source_id] - prev));
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
        double jac =
            sd.SMSingleSourceSirSimulation(message_recv.source_id, sp0);
        if (jac == 1.0) fi++;
      }

      message_recv.event_outcome = fi;
      Message toSend = message_recv;
      MPI::COMM_WORLD.Send(&toSend, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_RESPONSE);
    }
  }
}

namespace SMP {
struct Message {
  int source_id;
  double event_outcome;
};

MPI::Datatype datatypeOfMessage() {
  int blockLen[2] = {1, 1};
  MPI::Aint offsets[2] = {offsetof(Message, source_id),
                          offsetof(Message, event_outcome)};
  MPI::Datatype types[2] = {MPI::INT, MPI::DOUBLE};
  return MPI::Datatype::Create_struct(2, blockLen, offsets, types);
}
}  // namespace SMP

void SoftMarginParal(SourceDetectionParams params) {
  using namespace SMP;
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

  if (rank == 0) {
    std::string file_name = "distribution_sm-" + std::to_string((int)(10 * p)) +
                            "-" + std::to_string((int)(10 * q)) + "_grid" +
                            std::to_string((int)sqrt(vertices));
    FILE *file = fopen(file_name.c_str(), "a");
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
          events_resp[received.source_id].push_back(received.event_outcome);
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
    // fprintf(file, "\n\n%.10lf %.10lf\n\n", snapshot.p(), snapshot.q());
    for (int v = 0; v < vertices; ++v) {
      P[v] /= sum;
      printf("%.10lf\n", P[v]);
      fprintf(file, "%.10lf ", P[v]);
    }
    printf("\n");
    fprintf(file, "\n");
    /******/
    fclose(file);
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

      /***/
      vector<double> fi;
      for (int t = 0; t < SIMUL_PER_REQ; ++t) {
        Realization sp0 = snapshot;
        fi.push_back(
            sd.SMSingleSourceSirSimulation(message_recv.source_id, sp0));
      }

      message_recv.event_outcome = sd.likelihood(fi, params.a());
      /****/

      Message toSend = message_recv;
      MPI::COMM_WORLD.Send(&toSend, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_RESPONSE);
    }
  }
}

}  // namespace cmplx
