#include <mpi.h>
#include <unistd.h>

#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <numeric>

#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "common/random.h"
#include "source_detection_params.h"

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::SourceDetectionParams;
using std::vector;

const int SIMUL_PER_REQ = 10000;

struct Message {
  int source_id;
  double event_outcome;
};

enum MessageType { SIMUL_PREREQUEST, SIMUL_REQUEST, SIMUL_RESPONSE };

MPI::Datatype datatypeOfMessage() {
  int blockLen[2] = {1, 1};
  MPI::Aint offsets[2] = {offsetof(Message, source_id),
                          offsetof(Message, event_outcome)};
  MPI::Datatype types[2] = {MPI::INT, MPI::DOUBLE};
  return MPI::Datatype::Create_struct(2, blockLen, offsets, types);
}

// TODO unite paral for direct mc and soft

int main(int argc, char **argv) {
  // Paralelized
  MPI::Init(argc, argv);
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  MPI_Barrier(MPI::COMM_WORLD);

  //SourceDetectionParams params = SourceDetectionParams::BenchmarkParams(1);
  SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  double p = params.realization().p();
  double q = params.realization().q();
  const int simulations = params.simulations();
  int vertices = params.graph().vertices();
  const IGraph &graph = params.graph();
  const Realization &snapshot = params.realization();

  // TODO(iva) Make only one copy!
  // snapshot.print();

  if (rank == 0) {
    FILE *file = fopen("distributions-sm-106", "a");
    // master process
    int cur_simul_count = 0;
    int cur_v = 0;
    while ((cur_v < vertices) && (snapshot.realization().bit(cur_v) == false))
      cur_v++;
    vector<vector<double>> events_resp(vertices, vector<double>());
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
      if (events_resp[v].size())
        P_v /= (int)events_resp[v].size();
      P.push_back(P_v);
      sum += P_v;
    }
    fprintf(file, "\n\n%.10lf %.10lf\n\n", snapshot.p(), snapshot.q());
    for (int v = 0; v < vertices; ++v) {
      P[v] /= sum;
      printf("%.10lf\n", P[v]);
      fprintf(file, "%.10lf ", P[v]);
    }
    printf("\n");
    fprintf(file, "\n\n");
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
  MPI::Finalize();
  return 0;
}
