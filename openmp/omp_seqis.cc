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

}  // anonymous

namespace cmplx {

void OMPSeqIS::benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                                   int benchmark_no) {
  std::string filename =
      "SEQSoftbench_SBS_" + std::to_string(benchmark_no) + ".info";
  FILE *f = fopen(filename.c_str(), "w+");
  std::vector<int> sims = {(int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6};
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

vector<double> OMPSeqIS::convMaster(cmplx::SourceDetectionParams *params) {
  std::vector<int> sims = {(int)1e4,      2 * (int)1e4,  4 * (int)1e4,
                           8 * (int)1e4,  10 * (int)1e4, 20 * (int)1e4,
                           40 * (int)1e4, 80 * (int)1e4, 100 * (int)1e4};
  int s0 = sims[0];
  printf("s0: %d\n", s0);

  params->setSimulations(s0);
  vector<double> p0 = master(params);
  double pMAP0 = *std::max_element(p0.begin(), p0.end());
  double c = 0.05;

  vector<double> res;
  res.clear();
  int convergeG = 0;
  for (int s = 1; s < (int)sims.size(); ++s) {
    int s1 = sims[s];
    printf("s2: %d\n", s1);
    params->setSimulations(s1);

    bool converge = true;
    vector<double> p1 = master(params);
    double pMAP1 = *std::max_element(p1.begin(), p1.end());
    double delta = dabs(pMAP1 - pMAP0) / pMAP1;
    printf("\rc: %lf ", delta);
    if (delta >= c) converge = false;
    int pos = 0;
    for (int j = 0; j < (int)p1.size(); ++j) {
      if (p1[j] > 0 && (dabs(p1[j] - p0[j]) > c)) converge = false;
      if (p1[j] > 0) {
        pos++;
      }
    }
    if (pos == 0) converge = false;
    printf("\n");
    if (converge)
      convergeG++;
    else
      convergeG = 0;
    if (convergeG > 0) {
      printf("Converged for n=%d\n", s0);
      res = p0;
      /*
      if (end) {
        send_simul_end();
      }
      */
      break;
    } else {
      printf("Not converged.\n");
    }

    s0 = s1;
    p0.clear();
    res = p1;
    p0.assign(p1.begin(), p1.end());
    pMAP0 = pMAP1;
  }
  params->setSimulations(s0);
  return res;
}

vector<double> OMPSeqIS::master(const SourceDetectionParams *params) {
  double p = params->realization().p();
  double q = params->realization().q();
  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  const RealizationRead &snapshot = params->realization();

  vector<double> events_resp(vertices, 0);
  vector<int> activeNodes = snapshot.realization().positions();
  const int simulations = params->simulations();
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();

  int j = 0;
  int node_id = 0;
  double event_outcome = 0;
#pragma omp parallel for default(none)               \
    shared(events_resp, jobs_remaining, activeNodes, \
           params) private(j, node_id, event_outcome)
  for (j = 0; j < jobs_remaining; j += SIMUL_PER_REQ) {
    node_id = activeNodes[jobs_remaining % simulations];
    events_resp[node_id] +=
        work(params, ModelType::SIR, node_id, SIMUL_PER_REQ);
  }

  printf("\r\n");
  /*****/
  double sum = 0;
  for (int v = 0; v < vertices; ++v) {
    sum += events_resp[v];
  }

  for (int v = 0; v < vertices; ++v) {
    if (sum > 0) events_resp[v] /= sum;
  }
  return events_resp;
}

double OMPSeqIS::work(const SourceDetectionParams *params, ModelType model_type,
                      int source_id, int sample_size) {
  int vertices = params->graph()->vertices();
  const IGraph *graph = params->graph().get();
  RealizationRead snapshot = params->realization();

  // Performs simulation on request.
  SequentialMCDetector sd(graph);
  // SequentialSoftMCDetector sd(graph);

  /***/
  double Pos = sd.seqPosterior(source_id, sample_size, snapshot,
                               cmplx::ResamplingType::SIMPLE_RANDOM_SAMPLING,
                               true /* p = 1 @ T = 5*/);
  return Pos;
}

}  // namespace cmplx
