#include "source_detector/source_detection_params.h"
#include "source_detector/source_detector.h"
#include "mpi/mpi_source_detection.h"

#include <mpi.h>
#include <ctime>
#include <string>

using cmplx::SourceDetectionParams;

// -n bench_no
// [no flag] - DirectMC
// -m - SoftMargin
// -s SoftMargin SIS
int main(int argc, char** argv) {
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
  cmplx::MpiParal* mpi_paral;
  std::string filename_prefix = "";
  if (!seq && !sm) {
    mpi_paral = new cmplx::MPIDirectMC();
    filename_prefix += "DMC_";
  } else if (sm) {
    mpi_paral = new cmplx::MPISoftMC();
    filename_prefix += "SM_";
  } else {
    mpi_paral = new cmplx::MPISeqIS();
    filename_prefix += "Seq_";
  }

  mpi_paral->benchmark(params.get(), n, cmplx::ModelType::SIR, filename_prefix);
  MPI::Finalize();
  delete mpi_paral;
  return 0;
}