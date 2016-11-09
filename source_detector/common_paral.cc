#include "common_paral.h"

#include <cassert>
#include <algorithm>

namespace cmplx {

namespace {

void print_to_file(std::string filename, std::vector<double> P,
                   SourceDetectionParams *params) {
  FILE *f = fopen(filename.c_str(), "w+");
  fprintf(f, "%s\n", params->summary().c_str());
  fprintf(f, "s: %lld\n", params->simulations());
  for (int i = 0; i < (int)P.size(); ++i)
    fprintf(f, "%.10lf%c", P[i], i + 1 == (int)P.size() ? '\n' : ' ');
  fclose(f);
}

}  // namespace

void CommonParal::generateDistribution(SourceDetectionParams *params,
                                       ModelType model_type,
                                       std::string &filename_prefix) {
  std::vector<double> P = convMaster(params);
  std::string filename = filename_prefix + params->summary();
  FILE *f = fopen(filename.c_str(), "a");
  fprintf(f, "s:%lld -%d ", params->simulations(), params->sourceID());
  for (int j = 0; j < (int)P.size(); ++j) {
    fprintf(f, "%.10lf%c", P[j], j == ((int)P.size() - 1) ? '\n' : ' ');
  }
  fclose(f);
}

void CommonParal::benchmark(SourceDetectionParams *params, int benchmark_no,
                            ModelType model_type, std::string filename_prefix) {
  std::vector<double> P = convMaster(params);
  std::string filename =
      filename_prefix + std::to_string(benchmark_no) + ".info";
  print_to_file(filename, P, params);
}

void CommonParal::benchmarkStepByStep(std::string filename_prefix,
                                      const std::vector<int> &sims,
                                      cmplx::SourceDetectionParams *params,
                                      int benchmark_no) {
  std::string filename =
      filename_prefix + std::to_string(benchmark_no) + ".info";
  FILE *f = fopen(filename.c_str(), "a+");

  params->setA(pow(2, -5));
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

std::vector<double> CommonParal::convMaster(SourceDetectionParams *params,
                                            const std::vector<int> &sims) {
  /***************  */
  double c = 0.05;  //
  /**************   */
  int s0 = sims[0];
  params->setSimulations(s0);
  params->setA(pow(2, -5));
  std::vector<double> p0 = master(params);

  int MAP0 = std::max_element(p0.begin(), p0.end()) - p0.begin();
  double pml0 = *std::max_element(p0.begin(), p0.end());
  for (int s_id = 1; s_id < (int)sims.size(); ++s_id) {
    int s1 = sims[s_id];
    printf("s: %d\n", s1);
    params->setSimulations(s1);
    std::vector<double> p1 = master(params);

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

}  // namespace cmplx
