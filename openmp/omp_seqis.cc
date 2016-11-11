#include "omp_seqis.h"

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
typedef long long ll;

std::vector<double> responseToProb(const std::vector<double> &events_resp,
                                   int vertices) {
  printf("\r\n");
  /*****/
  std::vector<double> P = events_resp;

  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    sum += P[v];
  }

  for (int v = 0; v < vertices; ++v) {
    if (sum > 0) P[v] /= sum;
  }
  return P;
}

}  // anonymous

namespace cmplx {

vector<double> OMPSeqIS::master(const SourceDetectionParams *params) {
  int vertices = params->graph()->vertices();
  const RealizationRead snapshot = params->realization();

  vector<double> events_resp(vertices, 0);
  vector<int> activeNodes = snapshot.realization().positions();
  const int simulations = params->simulations();
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();

  int j = 0;
  double events_response = 0;

  const IGraph graph = *params->graph().get();

#pragma omp parallel for default(none)                       \
    shared(events_resp, jobs_remaining, activeNodes, graph), \
    private(j, params, events_response)
  for (j = 0; j < jobs_remaining; j += SIMUL_PER_REQ) {
    int node_id = activeNodes[jobs_remaining % simulations];
    events_response =
        work(graph, snapshot, ModelType::SIR, node_id, SIMUL_PER_REQ);
    events_resp[node_id] += events_response;
  }

  return responseToProb(events_resp, vertices);
}

double OMPSeqIS::work(const IGraph &graph, const RealizationRead &snapshot,
                      ModelType model_type, int source_id, int sample_size) {

  // Performs simulation on request.
  auto sd =
      std::unique_ptr<SequentialMCDetector>(new SequentialMCDetector(&graph));
  // SequentialSoftMCDetector sd(graph);

  return sd->seqPosterior(source_id, sample_size, snapshot,
                          cmplx::ResamplingType::NONE, true /* p = 1 @ T = 5*/);
}

}  // namespace cmplx
