#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <unistd.h>
#include "source_detector.h"
#include "common/igraph.h"
#include "common/bit_array.h"
#include "common/sir_params.h"

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;

int main(int argc, char **argv) {
  IGraph g = IGraph::UndirectedLattice({3, 3});
  int n = g.vertices();
  SourceDetector sd;
  SirParams sp(0.5 /* p*/, 0.5 /* q */, 100 /* T */, BitArray::ones(n),
               BitArray::zeros(n));
  sd.directMonteCarloDetection(g, sp, 10);
  return 0;
}
