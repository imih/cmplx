#include "direct_mc_params.h"

#include <vector>

using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::Realization;

namespace cmplx {
DirectMCParams DirectMCParams::SupFig2Params() {
  int lattice_size1 = 5;
  int lattice_size2 = 4;
  IGraph graph = IGraph::UndirectedLattice({lattice_size1, lattice_size2});

  int vertices = graph.vertices();
  BitArray r = BitArray::zeros(vertices);
  std::vector<int> infected = {2, 6, 7, 8, 9, 12};
  for (int inf_id : infected) {
    r.set(inf_id, true);
  }
  double p = 0.2;
  double q = 0.3;
  int maxT = 5;
  Realization realization = Realization(p, q, maxT, r);

  int simulations = 1000000000;
  return DirectMCParams(graph, realization, simulations);
}

} // namespace cmplx

