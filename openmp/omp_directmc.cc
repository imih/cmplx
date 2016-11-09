#include "omp_directmc.h"

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

using std::string;
using std::vector;

namespace {
const int SIMUL_PER_REQ = 10000;

typedef long long ll;
}  // anonymous

namespace cmplx {

vector<double> OMPDirectMCParal::master(const SourceDetectionParams *params) {
  int vertices = params->graph()->vertices();
  vector<int> events_resp(vertices, 0);

  const long long simulations = params->simulations();
  const common::RealizationRead &snapshot = params->realization();
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();

  int SIMUL_PER_REQ = std::max(10000LL, simulations / 10000);
  assert(simulations % (int)SIMUL_PER_REQ == 0);

  vector<int> activeNodes = snapshot.realization().positions();

  int j = 0;
  int node_id = 0;
  double event_outcome = 0;
#pragma omp parallel for default(none)                              \
    shared(events_resp, jobs_remaining, SIMUL_PER_REQ, activeNodes, \
           params) private(j, node_id, event_outcome)
  for (j = 0; j < jobs_remaining; j += SIMUL_PER_REQ) {
    node_id = activeNodes[jobs_remaining % simulations];
    event_outcome = work(params, ModelType::SIR, node_id, SIMUL_PER_REQ);
    events_resp[node_id] += event_outcome;
  }

  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    sum += events_resp[v];
  }

  printf("\r\r\n");
  vector<double> p;
  p.clear();
  for (int v = 0; v < vertices; ++v) {
    p.push_back(events_resp[v] / sum);
  }
  return p;
}

// TODO make thread safe!
double OMPDirectMCParal::work(const SourceDetectionParams *params,
                              ModelType model_type, int source_id,
                              int batch_size) {
  const IGraph *graph = params->graph().get();
  const common::RealizationRead &snapshot = params->realization();

  // workers
  // Performs simulation on request.
  DirectMonteCarloDetector sd(graph);
  int outcomes = 0;
  for (int t = 0; t < batch_size; ++t) {
    common::RealizationRead sp0 = snapshot;
    outcomes += sd.DMCSingleSourceSimulation(source_id, sp0, model_type);
  }
  return outcomes;
}

}  // namespace cmplx
