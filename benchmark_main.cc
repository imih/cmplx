#include "source_detector/source_detection_params.h"
#include "source_detector/source_detector.h"
#include "mpi/mpi_source_detection.h"
#include "mpi/mpi_paral.h"

#include <mpi.h>
#include <ctime>
#include <string>
#include <cassert>

using cmplx::SourceDetectionParams;

// -n bench_no
// -d -> DirectMC
// -s -> SequentialMC
// -m -> SoftMargin
// -f -> SoftSeqMC
// -b -> step_by_step mode
int main(int argc, char** argv) {
  // Paralelized
  int provided_level = MPI::Init_thread(argc, argv, MPI_THREAD_MULTIPLE);
  printf("Provided level of thread support: %d", provided_level);

  bool seq = false;
  bool sm = false;
  bool direct = false;
  bool softseq = false;
  bool sbs = false;
  int n = 0;
  {
    int c;
    while ((c = getopt(argc, argv, "n:dsmfb")) != EOF) {
      switch (c) {
        case 'n':
          n = atoi(optarg);
          break;
        case 'd':
          direct = true;
          break;
        case 's':
          seq = true;
          break;
        case 'm':
          sm = true;
          break;
        case 'f':
          softseq = true;
          break;
        case 'b':
          sbs = true;
          break;
      }
    }
  }

  auto params = SourceDetectionParams::BenchmarkParams(n);
  std::string filename_prefix = "";
  std::unique_ptr<cmplx::MpiMaster> mpi_master;
  std::unique_ptr<cmplx::CommonTraits> common_traits;
  if (direct) {
    mpi_master = std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPIDirectMC());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalDirectMC());
    filename_prefix += "DMC_";
  } else if (sm) {
    mpi_master = std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPISoftMC());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalSoftMC());
    filename_prefix += "SM_";
  } else if (seq) {
    mpi_master = std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPISeqIS());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalSeqIS());
    filename_prefix += "Seq_";
  } else if (softseq) {
    mpi_master = std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPISoftSeqIS());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalSoftSeqIS());
    filename_prefix += "SoftSeq_";
  }

  std::unique_ptr<cmplx::MpiParal> mpi_paral(
      new cmplx::MpiParal(std::move(mpi_master), std::move(common_traits)));
  assert(!(sbs && direct));
  if (sbs) {
    mpi_paral->benchmarkStepByStep(params.get(), n, cmplx::ModelType::SIR);
  } else {
    mpi_paral->benchmark(params.get(), n, cmplx::ModelType::SIR,
                         filename_prefix);
  }
  MPI::Finalize();
  return 0;
}
