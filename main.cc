#include <mpi.h>
#include <ctime>

#include "source_detection_params.h"
#include "source_detection_paral.h"

using cmplx::SourceDetectionParams;

/*
 *

  int P = 0, Q = 0;
  {
    int c;
    while((c = getopt(argc, argv, "p:q:")) != EOF) {
      switch (c) {
        case 'p':
          P = atoi(optarg);
          break;
        case 'q':
          Q = atoi(optarg);
          break;
    }
  }
  }

  SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  //SourceDetectionParams::ParamsFromGrid(P / 10.0, Q / 10.0);
 *
 * */

// -type {-P 10p -Q 10q}
int main(int argc, char **argv) {
  clock_t begin = std::clock();
  // Paralelized
  MPI::Init(argc, argv);

  /******TODO ********** DirectMC Paral *********/
  // SourceDetectionParams params = SourceDetectionParams::SupFig2Params();
  SourceDetectionParams params = SourceDetectionParams::BenchmarkParams(1);
  cmplx::DirectMCSimulParal(params);

  MPI::Finalize();
  clock_t end = clock();
  printf("%lf sec\n", double(end - begin) / CLOCKS_PER_SEC);
  return 0;
}
