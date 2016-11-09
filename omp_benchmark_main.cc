#include "source_detector/source_detection_params.h"
#include "source_detector/source_detector.h"
#include "openmp/omp_source_detection.h"

#include <ctime>
#include <string>

#include <omp.h>

using cmplx::SourceDetectionParams;

// -n bench_no
// [no flag] - DirectMC
// -m - SoftMargin
// -s SoftMargin SIS
int main(int argc, char** argv) {
  // Paralelized
  //
  //(void)omp_set_dynamic(1);

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
  std::string filename_prefix = "";

  std::unique_ptr<cmplx::CommonMaster> common_master;
  std::unique_ptr<cmplx::CommonTraits> common_traits;
  if (!seq && !sm) {
    common_master =
        std::unique_ptr<cmplx::CommonMaster>(new cmplx::OMPDirectMCParal());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalDirectMC());
    filename_prefix += "DMC_";
  } else if (sm) {
    common_master =
        std::unique_ptr<cmplx::CommonMaster>(new cmplx::OMPSoftMC());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalSoftMC());
    filename_prefix += "SM_";
  } else {
    common_master = std::unique_ptr<cmplx::CommonMaster>(new cmplx::OMPSeqIS());
    common_traits =
        std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalSeqIS());
    filename_prefix += "Seq_";
  }

  std::unique_ptr<cmplx::CommonParal> omp_paral(new cmplx::CommonParal(
      std::move(common_master), std::move(common_traits)));
  omp_paral->benchmark(params.get(), n, cmplx::ModelType::SIR, filename_prefix);
  return 0;
}
