#include "source_detector.h"
#include "common/igraph.h"
#include "common/bit_array.h"
#include "common/sir_params.h"

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
MPI::Datatype message_type = MPI::INT.Create_contiguous(2);
enum MessageType { SIMUL_PREREQUEST, SIMUL_REQUEST, SIMUL_RESPONSE };

// -s simulations_no -l lattice_size
int main(int argc, char **argv) {
  // Paralelized
  MPI::Init(argc, argv);
  int processes = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();

  int command = 0;
  int simulations = 100;
  int lattice_size = 3;
  while ((command = getopt(argc, argv, "s:l:")) != -1) {
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
    fprintf(stderr, "Performing simulations with %d processes\n", processes);
    int cur_simul_count = 0;
    int cur_v = 0;
    vector<int> events_resp(0, vertices);
    long long jobs_remaining = vertices * simulations;
    Message init_message;
    MPI::Request init_request =
        MPI::COMM_WORLD.Irecv(&init_message, 1, message_type, MPI::ANY_SOURCE,
                              MessageType::SIMUL_PREREQUEST);
    Message done_message;
    MPI::Request done_request =
        MPI::COMM_WORLD.Irecv(&done_message, 1, message_type, MPI::ANY_SOURCE,
                              MessageType::SIMUL_RESPONSE);
    while (jobs_remaining > 0) {
      MPI::Status init_request_status;
      if (init_request.Test(init_request_status)) {
        int sender = init_request_status.Get_source();
        init_request = MPI::COMM_WORLD.Irecv(&init_message, 1, message_type,
                                             MPI::ANY_SOURCE,
                                             MessageType::SIMUL_PREREQUEST);
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
          // There's no more jobs
          Message message = {-1, -1};
          MPI::COMM_WORLD.Isend(&message, 1, message_type, sender,
                                MessageType::SIMUL_REQUEST);
        } else {
          Message m = {cur_v, 0};
          MPI::COMM_WORLD.Isend(&m, 1, message_type, sender,
                                MessageType::SIMUL_REQUEST);
        }
      }

      if (done_request.Test()) {
        cur_simul_count--;
        Message received = done_message;
        if (cur_simul_count) {
          done_request = MPI::COMM_WORLD.Irecv(&done_message, 1, message_type,
                                               MPI::ANY_SOURCE,
                                               MessageType::SIMUL_RESPONSE);
        }
        events_resp[received.source_id] += received.event_outcome;
      }
    }
    fprintf(stderr, "Simulations finished");
    for (int v = 0; v < vertices; ++v) {
      printf("%d %d\n", v, events_resp[v]);
    }
  } else {
    // workers
    // Performs simulation on request.
    SourceDetector sd;
    while (true) {
      Message message = {-1, 0};
      MPI::Request request =
          MPI::COMM_WORLD.Isend(&message, 1, message_type, 0 /* dest */,
                                MessageType::SIMUL_PREREQUEST);
      request.Wait();
      request = MPI::COMM_WORLD.Irecv(&message, 1, message_type, 0 /* source */,
                                      MessageType::SIMUL_REQUEST);
      request.Wait();
      if (message.source_id == -1) {
        break;
      }
      int outcome = sd.SSSirSimulation(message.source_id, graph, sp);
      message.event_outcome = outcome;
      request = MPI::COMM_WORLD.Isend(&message, 1, message_type, 0 /* dest */,
                                      MessageType::SIMUL_RESPONSE);
    }
  }
  MPI::Finalize();
  return 0;
}
