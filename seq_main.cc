#include <mpi.h>
#include <ctime>
#include <string>
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
  } else if(sm) {
   cmplx::SoftMarginBenchmarkConv(params.get(), n);
  } else {
    cmplx::SeqMonteCarloBenchmark(params.get(), n);
    //SeqMonteCarloBenchmarkStepByStep(params.get(), n);
  }

  MPI::Finalize();
  // clock_t end = clock();
  // printf("%lf sec\n", double(end - begin) / CLOCKS_PER_SEC);
  return 0;
}
