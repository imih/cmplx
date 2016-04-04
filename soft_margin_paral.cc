#include <mpi.h>
#include <unistd.h>

#include <ctime>

#include "source_detection_params.h"
#include "source_detection_paral.h"

using cmplx::SourceDetectionParams;

int main(int argc, char **argv) {
  clock_t begin = std::clock();

  // Paralelized
  MPI::Init(argc, argv);
  int P = 0, Q = 0;
  int c;
  while ((c = getopt(argc, argv, "p:q:")) != EOF) {
    switch (c) {
      case 'p':
        P = atoi(optarg);
        break;
      case 'q':
        Q = atoi(optarg);
        break;
    }
  }

  SourceDetectionParams params =
      SourceDetectionParams::ParamsFromGrid(P / 10.0, Q / 10.0);
  cmplx::SoftMarginParal(params);

  MPI::Finalize();
  clock_t end = clock();
  printf("%lf sec\n", double(end - begin) / CLOCKS_PER_SEC);
  return 0;
}
