#include <mpi.h>
#include <ctime>
#include <string>
#include "source_detection_params.h"
#include "source_detection_paral.h"
#include "source_detector.h"
using cmplx::SourceDetectionParams;

// -type {-P 10p -Q 10q}
int main(int argc, char** argv) {
  // clock_t begin = std::clock();
  // Paralelized
  MPI::Init(argc, argv);
  bool seq = false;
  int P = 5, Q = 5;
  int n = 30;
  int g = 0;
  int v = 0;
  bool erdos = false;
  {
    int c;
    while ((c = getopt(argc, argv, "g:v:p:q:n:se")) != EOF) {
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
        case 'g':
          g = atoi(optarg);
          break;
        case 'v':
          v = atoi(optarg);
          break;
        case 'e':
          erdos = true;
          break;
      }
    }
  }

  std::string file_name =
      "/home/imiholic/graphs/barabasi2_100_" + std::to_string(g) + ".gml";
  if (erdos) {
    file_name = "/home/imiholic/graphs/erdos_renyi_100_0.01_" +
                std::to_string(g) + ".gml";
  }
  auto params = SourceDetectionParams::ParamsFromGML(
    file_name, v, P / 10.0, Q / 10.0);

/*
  auto params = SourceDetectionParams::ParamsFromGrid(
    P / 10.0, Q / 10.0, n);
*/
  int rank = MPI::COMM_WORLD.Get_rank();
  if (rank == 0) {
    std::string filename = "barabasi2_100_" + params->summary();
    FILE* f = fopen(filename.c_str(), "a");
    fprintf(f, "g: %d, bc: %d, ", g, params.get()->realization().realization().bitCount());
    fclose(f);
  }
  MPI::COMM_WORLD.Barrier();
  cmplx::GenerateSoftMarginDistributions(
    params.get(), 1);
  //cmplx::SeqMonteCarloBenchmark(params.get(), 1);
  MPI::Finalize();
  // clock_t end = clock();
  // printf("%lf sec\n", double(end - begin) / CLOCKS_PER_SEC);
  return 0;
}
