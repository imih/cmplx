#include "source_detector/source_detection_params.h"
#include "source_detector/source_detector.h"
#include "mpi/mpi_source_detection.h"
#include "mpi/mpi_paral.h"

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

  auto params = SourceDetectionParams::ParamsFromGML("WattsStrogatz30_0.5.gml",
                                                     30 * 30 / 2, 0.3, 0.3);
  std::string filename_prefix = "";
  std::unique_ptr<cmplx::MpiMaster> mpi_master;
  std::unique_ptr<cmplx::CommonTraits> common_traits;
  if (!seq && !sm) {
    mpi_master = std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPIDirectMC());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalDirectMC());
    filename_prefix += "DMCWS0.5_";
  } else if (sm) {
    mpi_master = std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPISoftMC());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalSoftMC());
    filename_prefix += "SMWS0.5_";
  } else {
    mpi_master = std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPISeqIS());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalSeqIS());
    filename_prefix += "SeqWS0.5_";
  }

  std::unique_ptr<cmplx::MpiParal> mpi_paral(
      new cmplx::MpiParal(std::move(mpi_master), std::move(common_traits)));
  mpi_paral->benchmark(params.get(), n, cmplx::ModelType::SIR, filename_prefix);
  MPI::Finalize();
  return 0;
}
