#include "source_detector/source_detection_params.h"
#include "source_detector/source_detector.h"
#include "mpi/mpi_source_detection.h"
#include "mpi/mpi_paral.h"

#include <mpi.h>
#include <ctime>
#include <string>

using cmplx::SourceDetectionParams;

void build_realizations() {
  std::string graph_path = "WattsStrogatz30_0.5.gml";
  int n = 30;
  int source_id = n * n / 2;
  for (int id = 1; id <= 160; ++id) {
    double p = 0.7;
    double q = 0.3;
    if (id <= 80) p = 0.3;
    if (id % 2 == 0) q = 0.7;
    auto params =
        SourceDetectionParams::ParamsFromGML(graph_path, source_id, p, q);
    std::string file_name =
        "realizations/realization_" + std::to_string(id) + ".txt";
    FILE* f = fopen(file_name.c_str(), "w");
    fprintf(f, "#source: %d\n", source_id);
    fprintf(f, "#p: %.lf\n#q: %.lf\n#T: 5\n#node_states:\n", p, q);
    for (int i = 0; i < 900; ++i)
      fprintf(f, "%d\n", params->realization().realization().bit(i) ? 1 : 0);
    fclose(f);
  }
}
// -n bench_no
// [no flag] - DirectMC
// -m - SoftMargin
// -s SoftMargin SIS
int main(int argc, char** argv) {
  // Paralelized
  std::string graph_path = "WattsStrogatz30_0.5.gml";
  // build_realizations();

  int id = 0;
  char c = 0;
  while ((c = getopt(argc, argv, "n:")) != EOF) {
    switch (c) {
      case 'n':
        id = atoi(optarg);
    }
  }

  MPI::Init(argc, argv);
  std::unique_ptr<cmplx::MpiMaster> mpi_master =
      std::unique_ptr<cmplx::MpiMaster>(new cmplx::MPIDirectMC());
  std::unique_ptr<cmplx::CommonTraits> common_traits =
      std::unique_ptr<cmplx::CommonTraits>(new cmplx::ParalDirectMC());
  std::unique_ptr<cmplx::MpiParal> mpi_paral(
      new cmplx::MpiParal(std::move(mpi_master), std::move(common_traits)));

  auto params =
      SourceDetectionParams::ParamsFromGraphRealization(graph_path, id);
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank == 0) {
    std::vector<double> P = mpi_paral->convMaster(params.get());
    std::string file_name =
        "solutions/inverse_solution_" + std::to_string(id) + ".txt";
    FILE* f = fopen(file_name.c_str(), "w");
    fprintf(f, "#Number of simulations: %lld\n", params->simulations());
    fprintf(f, "#Node source probabilities:\n");
    for (int i = 0; i < 900; ++i) {
      printf("%.10lf\n", P[i]);
    }
    fclose(f);
  } else {
    mpi_master->worker(params.get(), cmplx::ModelType::SIR);
  }

  MPI::Finalize();

  return 0;
}
