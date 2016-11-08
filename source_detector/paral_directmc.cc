#include "paral_directmc.h"

#include <vector>
#include <string>
#include <algorithm>

namespace cmplx {

namespace {
double dabs(double x) {
  if (x < 0) return x * -1;
  return x;
}
}  // namespace

void ParalDirectMC::benchmarkStepByStep(cmplx::SourceDetectionParams* params,
                                        int benchmark_no) {
  std::string filename =
      "DMCbench_SBS_" + std::to_string(benchmark_no) + ".info";
  FILE* f = fopen(filename.c_str(), "a+");
  vector<int> sims = {(int)1e4, (int)1e5, (int)1e6,
                      (int)1e7, (int)1e8, (int)1e9};
  for (int sim : sims) {
    fprintf(f, "s: %d\n", sim);
    params->setSimulations(sim);
    std::vector<double> p = master(params);
    for (int i = 0; i < (int)p.size(); ++i)
      fprintf(f, "%.10lf%c", p[i], i + 1 == (int)p.size() ? '\n' : ' ');
    fflush(f);
  }
  fclose(f);
}

std::vector<double> ParalDirectMC::responseToProb(
    const std::vector<int>& events_resp, int vertices) {
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

vector<double> ParalDirectMC::convMaster(SourceDetectionParams* params) {
  vector<int> sims = {(int)1e6,     2 * (int)1e6, 4 * (int)1e6, (int)1e7,
                      2 * (int)1e7, 4 * (int)1e7, (int)1e8,     2 * (int)1e8,
                      4 * (int)1e8, (int)1e9};

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
    if (isnan(pml1)) {
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
}
