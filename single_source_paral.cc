#include "source_detector.h"
#include "common/igraph.h"
#include "common/bit_array.h"
#include "common/sir_params.h"

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <unistd.h>
#include <vector>

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using std::vector;

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
int main(int argc, char **argv) {
  // Paralelized
  MPI::Init(argc, argv);
  MPI::Datatype message_type = datatypeOfMessage();
  message_type.Commit();

  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();

  int command = 0;
  int simulations = 10;
  int lattice_size = 3;
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

  // TODO Make only one copy!
  IGraph graph = IGraph::UndirectedLattice({lattice_size, lattice_size});
  int vertices = graph.vertices();
  SirParams sp(0.5 /* p*/, 0.5 /* q */, 100 /* T */,
               BitArray::ones(vertices) /* infected*/,
               BitArray::zeros(vertices) /* susceptible */);

  if (rank == 0) {
    // master process
    fprintf(stderr, "P0: Performing simulations with %d processes\n",
            processes);
    int cur_simul_count = 0;
    int cur_v = 0;
    vector<int> events_resp(vertices, 0);
    long long jobs_remaining = vertices * simulations;
    while (jobs_remaining > 0) {
      for (int i = 0; i < processes - 1; ++i) {
        if (MPI::COMM_WORLD.Iprobe(i + 1, MessageType::SIMUL_PREREQUEST)) {
          Message init_message;
          MPI::COMM_WORLD.Recv(&init_message, 1, message_type, i + 1,
                               MessageType::SIMUL_PREREQUEST);
          //   printf("P0 Received init request from %d\n", i + 1);
          if (cur_simul_count < simulations) {
            cur_simul_count++;
          } else {
            cur_simul_count = 1;
            cur_v++;
            if (cur_v == vertices) {
              cur_v = -1;
            }
          }
          if (cur_v == -1) {
            //    puts("No more jobs");
            // There's no more jobs
            Message message = {-1, -1};
            MPI::COMM_WORLD.Send(&message, 1, message_type, i + 1,
                                 MessageType::SIMUL_REQUEST);
          } else {
            //      printf("P0 Sending request for %d to %d\n", cur_v, i + 1);
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
          jobs_remaining--;
          //      printf("PO Received simulation of %d %d\n",
          //      received.source_id,
          //            received.event_outcome);
          events_resp[received.source_id] += received.event_outcome;
        }
      }
    }

    fprintf(stderr, "Simulations finished");
    for (int v = 0; v < vertices; ++v) {
      printf("\n%d %d\n", v, events_resp[v]);
    }
    for (int i = 0; i < processes - 1; ++i) {
      Message m = {-1, -1};
      MPI::COMM_WORLD.Send(&m, 1, message_type, i + 1,
                           MessageType::SIMUL_REQUEST);
    }
  } else {
    // workers
    // Performs simulation on request.
    SourceDetector sd;
    while (true) {
      Message message = {-1, 0};
      //   printf("P%d: Sending init request\n", rank);
      MPI::COMM_WORLD.Send(&message, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_PREREQUEST);
      Message message_recv;
      MPI::COMM_WORLD.Recv(&message_recv, 1, message_type, 0 /* source */,
                           MessageType::SIMUL_REQUEST);
      if (message_recv.source_id == -1) {
        break;
      }
      //    printf("P%d: Request for %d\n", rank, message_recv.source_id);
      SirParams sp0 = sp;
      int outcome = sd.SSSirSimulation(message_recv.source_id, graph, sp0);
      message_recv.event_outcome = outcome;

      Message toSend = message_recv;
      MPI::COMM_WORLD.Send(&toSend, 1, message_type, 0 /* dest */,
                           MessageType::SIMUL_RESPONSE);
      //     printf("P%d: Sending response %d for %d\n", rank, outcome,
      //            message_recv.source_id);
    }
    printf("P%d: Done.", rank);
  }
  MPI::Finalize();
  return 0;
}
