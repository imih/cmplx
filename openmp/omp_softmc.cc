#include "omp_softmc.h"

#include "../common/bit_array.h"
#include "../common/igraph.h"
#include "../common/realization.h"

#include <omp.h>

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

typedef long long ll;
}  // anonymous

namespace cmplx {

vector<double> OMPSoftMC::master(const SourceDetectionParams *params) {
  const int simulations = params->simulations();
  int vertices = params->graph()->vertices();
  const common::RealizationRead snapshot = params->realization();

  // master process
  vector<double> events_resp_sum(vertices, 0);
  vector<long long> events_resp_size(vertices, 0);
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();
  assert(jobs_remaining % (int)SIMUL_PER_REQ == 0);
  vector<int> activeNodes = snapshot.realization().positions();

  const IGraph graph = *params->graph().get();

  int j = 0;
  int node_id = 0;
  double event_outcome = 0;
  double a = params->a();
#pragma omp parallel for default(none)                                     \
    shared(events_resp_sum, events_resp_size, jobs_remaining, activeNodes, \
           graph, a) private(j, node_id, event_outcome)
  for (j = 0; j < jobs_remaining; j += SIMUL_PER_REQ) {
    node_id = activeNodes[jobs_remaining % simulations];
    event_outcome = work(graph, snapshot, a, ModelType::SIR, node_id);
    events_resp_sum[node_id] += event_outcome;
    events_resp_size[node_id]++;
  }

  return responseToProb(events_resp_sum, events_resp_size, vertices);
}

double OMPSoftMC::work(const IGraph &graph, const RealizationRead &snapshot,
                       double a, ModelType model_type, int node_id) {

  // workers
  // Performs simulation on request.
  auto sd = std::unique_ptr<SoftMarginDetector>(new SoftMarginDetector(&graph));

  /***/
  vector<double> fi;
  fi.clear();
  for (int t = 0; t < (int)SIMUL_PER_REQ; ++t) {
    RealizationRead sp0 = snapshot;
    fi.push_back(sd->SMSingleSourceSimulation(node_id, sp0, model_type));
  }
  return sd->likelihood(fi, a);
}

}  // namespace cmplx
