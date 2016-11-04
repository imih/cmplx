#include "source_detector/source_detection_params.h"
#include "source_detector/source_detector.h"
#include "mpi/source_detection_paral.h"

#include <mpi.h>
#include <ctime>

using cmplx::SourceDetectionParams;

// -type {-P 10p -Q 10q}
int main(int argc, char **argv) {

  // Paralelized
  MPI::Init(argc, argv);
  bool seq = false;
  int P = 5, Q = 5;
  int n = 30;
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

  std::unique_ptr<SourceDetectionParams> params =
      SourceDetectionParams::ParamsFromGridISS(P / 10.0, Q / 10.0, n);
  cmplx::GenerateSoftMarginDistributions(params.get(), 1,
                                         cmplx::ModelType::ISS);

  MPI::Finalize();
  return 0;
}
