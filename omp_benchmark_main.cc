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
  cmplx::OmpParal* omp_paral;
  std::string filename_prefix = "";
  if (!seq && !sm) {
    omp_paral = new cmplx::OMPDirectMCParal();
    filename_prefix += "DMC_";
  } else if (sm) {
    omp_paral = new cmplx::OMPSoftMC();
    filename_prefix += "SM_";
  } else {
    omp_paral = new cmplx::OMPSeqIS();
    filename_prefix += "Seq_";
  }

  omp_paral->benchmark(params.get(), n, cmplx::ModelType::SIR, filename_prefix);
  return 0;
}
