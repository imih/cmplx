#include <mpi.h>
#include <ctime>

#include "source_detection_params.h"
#include "source_detection_paral.h"

using cmplx::SourceDetectionParams;

// -type {-P 10p -Q 10q}
int main(int argc, char **argv) {
 // clock_t begin = std::clock();
  // Paralelized
  MPI::Init(argc, argv);

  SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  // SourceDetectionParams params = SourceDetectionParams::BenchmarkParams(1);
  //cmplx::DirectMCSimulParalConv(params);
  //cmplx::DirectMCSimulParal(params);
  //cmplx::SoftMarginParal(params);
  cmplx::SoftMarginParalConv(params);

  MPI::Finalize();
 // clock_t end = clock();
 // printf("%lf sec\n", double(end - begin) / CLOCKS_PER_SEC);
  return 0;
}
