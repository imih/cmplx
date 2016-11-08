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

typedef long long ll;
}  // anonymous

namespace cmplx {

void OMPSoftMC::benchmarkStepByStep(cmplx::SourceDetectionParams *params,
                                    int benchmark_no) {
  std::string filename =
      "SMbench_SBS_" + std::to_string(benchmark_no) + ".info";
  FILE *f = fopen(filename.c_str(), "a+");
  std::vector<int> sims = {(int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6};

  double a = pow(2, -5);
  params->setA(a);
  for (int sim : sims) {
    fprintf(f, "s: %d\n", sim);
    params->setSimulations(sim);
    vector<double> p = master(params);
    for (int i = 0; i < (int)p.size(); ++i)
      fprintf(f, "%.10lf%c", p[i], i + 1 == (int)p.size() ? '\n' : ' ');
    fflush(f);
  }
  fclose(f);
  exit(0);
}

vector<double> OMPSoftMC::master(const SourceDetectionParams *params) {
  const int simulations = params->simulations();
  int vertices = params->graph()->vertices();
  const common::RealizationRead &snapshot = params->realization();

  // master process
  vector<double> events_resp_sum(vertices, 0);
  vector<long long> events_resp_size(vertices, 0);
  long long jobs_remaining =
      1LL * simulations * snapshot.realization().bitCount();
  assert(jobs_remaining % (int)SIMUL_PER_REQ == 0);
  vector<int> activeNodes = snapshot.realization().positions();

  int j = 0;
  int node_id = 0;
  double event_outcome = 0;
#pragma omp parallel for default(none)                                     \
    shared(events_resp_sum, events_resp_size, jobs_remaining, activeNodes, \
           params) private(j, node_id, event_outcome)
  for (j = 0; j < jobs_remaining; j += SIMUL_PER_REQ) {
    node_id = activeNodes[jobs_remaining % simulations];
    event_outcome = work(params, ModelType::SIR, node_id);
    events_resp_sum[node_id] += event_outcome;
    events_resp_size[node_id]++;
  }

  printf("\r\n");
  /*****/
  std::vector<double> P;
  P.clear();
  double sum = 0;

  int first_size = 0;
  bool normalize = false;
  // TODO possibly not needed
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

std::vector<double> OMPSoftMC::convMaster(SourceDetectionParams *params) {
  std::vector<int> sims = {(int)1e4,       2 * (int)1e4,   4 * (int)1e4,
                           10 * (int)1e4,  20 * (int)1e4,  40 * (int)1e4,
                           100 * (int)1e4, 200 * (int)1e4, 400 * (int)1e4,
                           1000 * (int)1e4};

  int convergeG = 0;
  double c = 0.05;
  double a = pow(2, -5);
  int s0 = sims[0];
  printf("s0: %d\n", s0);
  params->setSimulations(s0);
  params->setA(a);
  vector<double> p0 = master(params);
  int MAP0 = std::max_element(p0.begin(), p0.end()) - p0.begin();
  double pMAP0 = *std::max_element(p0.begin(), p0.end());
  vector<double> P;
  P.clear();
  for (int s_id = 1; s_id < (int)sims.size(); ++s_id) {
    int s1 = sims[s_id];
    params->setSimulations(s1);
    printf("\ns: %d\n", s1);
    vector<double> p1 = master(params);
    int MAP1 = std::max_element(p1.begin(), p1.end()) - p1.begin();
    double pMAP1 = *std::max_element(p1.begin(), p1.end());
    double delta = dabs(pMAP1 - pMAP0) / pMAP1;
    printf("c: %lf\n", delta);
    double converge = true;
    if (MAP0 != MAP1) converge = false;
    if (delta > c) converge = false;
    int pos = 0;
    for (int j = 0; j < (int)p1.size(); ++j) {
      if (dabs(p1[j] - p0[j]) > c) converge = false;
      if (p1[j] > 0) pos++;
    }
    if (pos == 0) converge = false;
    if (converge) {
      convergeG++;
      printf("Converged for n=%d\n", s0);
      P = p0;
      params->setSimulations(s0);
      if (convergeG > 0) break;
    } else {
      convergeG = 0;
      printf("Not converged.\n");
    }
    s0 = s1;
    p0.assign(p1.begin(), p1.end());
    pMAP0 = pMAP1;
    MAP0 = MAP1;
    P = p1;
  }
  return P;
}

vector<double> OMPSoftMC::softConvMaster(cmplx::SourceDetectionParams *params) {
  std::vector<int> sims = {(int)1e4,       2 * (int)1e4,   4 * (int)1e4,
                           10 * (int)1e4,  20 * (int)1e4,  40 * (int)1e4,
                           100 * (int)1e4, 200 * (int)1e4, 400 * (int)1e4,
                           800 * (int)1e4};

  double c = 0.05;
  const int MAXA = 9;
  int s0 = sims[0];
  printf("s0: %d\n", s0);
  vector<double> a(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    a[i] = 1.0 / (double)(1 << i);
  }
  vector<double> p0[MAXA + 1];
  vector<double> pMAP0(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    params->setSimulations(s0);
    params->setA(a[i]);
    printf("a[i]: %lf\n", a[i]);
    p0[i] = master(params);
    pMAP0[i] = *std::max_element(p0[i].begin(), p0[i].end());
  }

  vector<int> convergeGlobal(MAXA + 1, 0);
  vector<double> P;
  for (int s = 1; s < (int)sims.size(); ++s) {
    int s1 = sims[s];
    printf("\ns: %d\n", s1);
    params->setSimulations(s1);
    vector<double> p1[MAXA + 1];
    vector<double> pMAP1(MAXA + 1, 0);

    for (int i = MAXA; i >= 3; --i) {
      printf("\ns: %d a: %.10lf\n", s1, a[i]);
      params->setA(a[i]);
      double converge = true;
      p1[i] = master(params);
      pMAP1[i] = *std::max_element(p1[i].begin(), p1[i].end());
      double delta = dabs(pMAP1[i] - pMAP0[i]) / pMAP1[i];
      printf("c: %lf\n", delta);
      if (delta > c) converge = false;
      int pos = 0;
      for (int j = 0; j < (int)p1[i].size(); ++j) {
        if (dabs(p1[i][j] - p0[i][j]) > c) converge = false;
        if (p1[i][j] > 0) pos++;
      }
      if (pos == 0) converge = false;
      if (converge) {
        convergeGlobal[i]++;
        printf("Converged for n=%d a=%lf\n", s0, a[i]);
        if (convergeGlobal[i] > 0) break;
      } else {
        convergeGlobal[i] = 0;
        printf("Not converged.\n");
      }
    }

    bool done = false;
    for (int i = MAXA; i >= 3; --i) {
      if (convergeGlobal[i] > 0) {
        P = p0[i];
        params->setSimulations(s0);
        params->setA(a[i]);
        done = true;
        break;
      }
    }

    s0 = s1;
    for (int i = 3; i <= MAXA; ++i) p0[i] = p1[i];
    pMAP0 = pMAP1;
    if (done) break;
    P = p1[MAXA];
  }
  return P;
}

double OMPSoftMC::work(const SourceDetectionParams *params,
                       ModelType model_type, int node_id) {
  const IGraph *graph = params->graph().get();
  common::RealizationRead snapshot = params->realization();

  // workers
  // Performs simulation on request.
  SoftMarginDetector sd(graph);

  /***/
  vector<double> fi;
  fi.clear();
  for (int t = 0; t < (int)SIMUL_PER_REQ; ++t) {
    RealizationRead sp0 = snapshot;
    fi.push_back(sd.SMSingleSourceSimulation(node_id, sp0, model_type));
  }
  return sd.likelihood(fi, params->a());
}

}  // namespace cmplx
