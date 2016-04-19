#include <mpi.h>
#include <ctime>

#include "source_detection_params.h"
#include "source_detection_paral.h"
#include "source_detector.h"

using cmplx::SourceDetectionParams;

// -type {-P 10p -Q 10q}
int main(int argc, char **argv) {
  // clock_t begin = std::clock();
  // Paralelized
  MPI::Init(argc, argv);
  bool seq = false;
  int P = 5, Q = 5;
  int n = 5;
  {
    int c;
    while ((c = getopt(argc, argv, "p:q:n:s")) != EOF) {
      switch (c) {
        case 'p':
          P = atoi(optarg);
          break;
        case 'q':
          Q = atoi(optarg);
          break;
        case 'n':
          n = atoi(optarg);
          break;
        case 's':
          seq = true;
          break;
      }
    }
  }

  // SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  // SourceDetectionParams params = SourceDetectionParams::BenchmarkParams(1);
  // cmplx::DirectMCSimulParalConv(params, cmplx::ModelType::SIR);
  // cmplx::DirectMCSimulParal(params);
  // cmplx::SoftMarginParal(params);
  // cmplx::SoftMarginParalConv(params);

  std::unique_ptr<SourceDetectionParams> params =
      SourceDetectionParams::ParamsFromGrid(P / 10.0, Q / 10.0, n);
  if (seq) {
    cmplx::GenerateSeqMonteCarloDistributions(params.get(), 1);
  } else {
    cmplx::GenerateSoftMarginDistributions(params.get(), 1);
  }

  MPI::Finalize();
  // clock_t end = clock();
  // printf("%lf sec\n", double(end - begin) / CLOCKS_PER_SEC);
  return 0;
}
