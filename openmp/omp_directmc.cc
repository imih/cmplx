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

void OMPDirectMCParal::benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                                           int benchmark_no) {
  string filename = "DMCbench_SBS_" + std::to_string(benchmark_no) + ".info";
  FILE *f = fopen(filename.c_str(), "a+");
  vector<int> sims = {(int)1e4, (int)1e5, (int)1e6,
                      (int)1e7, (int)1e8, (int)1e9};
  for (int sim : sims) {
    fprintf(f, "s: %d\n", sim);
    params->setSimulations(sim);
    vector<double> p = master(params);
    for (int i = 0; i < (int)p.size(); ++i)
      fprintf(f, "%.10lf%c", p[i], i + 1 == (int)p.size() ? '\n' : ' ');
    fflush(f);
  }
  fclose(f);
}

vector<double> OMPDirectMCParal::master(const SourceDetectionParams *params) {
  const long long simulations = params->simulations();

  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  const common::RealizationRead &snapshot = params->realization();

  vector<int> events_resp(vertices, 0);
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
  assert(sum > 0);

  printf("\r\r\n");
  vector<double> p;
  p.clear();
  for (int v = 0; v < vertices; ++v) {
    p.push_back(events_resp[v] / sum);
  }
  return p;
}

vector<double> OMPDirectMCParal::convMaster(SourceDetectionParams *params) {
  vector<int> sims = {100 * (int)1e4,   200 * (int)1e4,   1000 * (int)1e4,
                      2000 * (int)1e4,  10000 * (int)1e4, 20000 * (int)1e4,
                      100000 * (int)1e4};

  /***************  */
  double c = 0.05;  //
  /**************   */
  int s0 = sims[0];
  params->setSimulations(s0);
  vector<double> p0 = master(params);

  int MAP0 = std::max_element(p0.begin(), p0.end()) - p0.begin();
  double pml0 = *std::max_element(p0.begin(), p0.end());
  for (int s_id = 1; s_id < (int)sims.size(); ++s_id) {
    int s1 = sims[s_id];
    printf("s: %d\n", s1);
    params->setSimulations(s1);
    vector<double> p1 = master(params);

    int MAP1 = std::max_element(p1.begin(), p1.end()) - p1.begin();
    double pml1 = *std::max_element(p1.begin(), p1.end());
    bool converge = true;
    if (MAP0 != MAP1) converge = false;
    if (std::isnan(pml1)) {
      converge = false;
      printf("NAN! %d \n", s1);
    }
    double delta = dabs(pml1 - pml0) / pml1;
    printf("c: %lf\n", delta);
    if (delta > c) converge = false;
    int pos = 0;
    for (int i = 0; i < (int)p1.size(); ++i) {
      if (dabs(p0[i] - p1[i]) > c) converge = false;
      if (p1[i] > 0) pos++;
    }
    if (pos == 0) converge = false;
    if (converge) {
      printf("Converged for %d\n", s0);
      break;
    }
    pml0 = pml1;
    s0 = s1;
    p0 = p1;
    MAP0 = MAP1;
  }
  params->setSimulations(s0);
  return p0;
}

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
