#include <mpi.h>
#include <unistd.h>

#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include <vector>

#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "simul/simulator.h"
#include "direct_mc_params.h"

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::simul::Simulator;
using cmplx::common::Realization;
using cmplx::DirectMCParams;
using std::vector;

const int SIMUL_PER_REQ = 10000;

struct Message {
  int source_id;
  int event_outcome;
};

enum MessageType { SIMUL_PREREQUEST, SIMUL_REQUEST, SIMUL_RESPONSE };

MPI::Datatype datatypeOfMessage() {
  int blockLen[2] = {1, 1};
  MPI::Aint offsets[2] = {offsetof(Message, source_id),
                          offsetof(Message, event_outcome)};
  MPI::Datatype types[2] = {MPI::INT, MPI::INT};
  return MPI::Datatype::Create_struct(2, blockLen, offsets, types);
}

// -s simulations_no -l lattice_size
// TODO(iva) paralelize for multiple realizations
int main(int argc, char **argv) {
  // Paralelized
  MPI::Init(argc, argv);
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();

  /*
  int command = 0;
  int simulations = 100000000;
  int lattice_size = 30;
  while ((command = getopt(argc, argv, "s:l:n:")) != -1) {
    switch (command) {
    case 's':
      simulations = atoi(optarg);
      break;
    case 'l':
      lattice_size = atoi(optarg);
      break;
    }
  }
  */

  srand(time(NULL));
  DirectMCParams params = DirectMCParams::SupFig2Params();
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
	    printf("%.10lf\n", 100 * cur_v / (double) vertices);
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
        }
      }
    }

    fprintf(stderr, "Simulations finished");
    for (int v = 0; v < vertices; ++v) {
      printf("%.8f ", events_resp[v] / (double)simulations);
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
    SourceDetector sd;
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
        outcomes += sd.SSSirSimulation(message_recv.source_id, graph, sp0);
      }

      message_recv.event_outcome = outcomes;

      Message toSend = message_recv;
      MPI::COMM_WORLD.Send(&toSend, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_RESPONSE);
    }
  }
  MPI::Finalize();
  return 0;
}
