#include "omp_paral.h"

#include <cassert>

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

void OmpParal::generateDistribution(SourceDetectionParams *params,
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
  exit(0);
}

void OmpParal::benchmark(SourceDetectionParams *params, int benchmark_no,
                         ModelType model_type, std::string filename_prefix) {
  std::vector<double> P = convMaster(params);
  std::string filename =
      filename_prefix + std::to_string(benchmark_no) + ".info";
  print_to_file(filename, P, params);
  exit(0);
}

}  // namespace cmplx
