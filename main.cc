#include <mpi.h>
#include <ctime>
#include <string>
#include "source_detection_params.h"
#include "source_detection_paral.h"
#include "source_detector.h"

using cmplx::SourceDetectionParams;

// -n bench_no
// [no flag] - DirectMC
// -m - SoftMargin
// -s SoftMargin SIS
int main(int argc, char **argv) {
  // Paralelized
  MPI::Init(argc, argv);
  bool seq = false;
  bool sm = false;
  int n = 0;
  {
    int c;
    while ((c = getopt(argc, argv, "n:sm")) != EOF) {
      switch (c) {
        case 'n':
          n = atoi(optarg);
          break;
        case 's':
          seq = true;
          break;
        case 'm':
          sm = true;
          break;
      }
    }
  }

  auto params = SourceDetectionParams::BenchmarkParams(n);
  if (!seq && !sm) {
    cmplx::DirectMCBenchmark(params.get(), n);
  } else if (sm) {
    cmplx::SoftMarginBenchmarkConv(params.get(), n);
    // cmplx::SoftMarginBenchmarkStepByStep(params.get(), n);
  } else {
    // SeqMonteCarloBenchmarkStepByStep(params.get(), n);
    cmplx::SeqMonteCarloBenchmark(params.get(), n);
  }

  MPI::Finalize();
  return 0;
}
