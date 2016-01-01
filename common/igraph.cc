#include "igraph.h"

#include <cassert>

#include "ivector.h"

namespace cmplx {
namespace common {
IGraph IGraph::UndirectedLattice(const std::vector<int>& dimensions, int nei) {
  IVector<int> dim(dimensions);
  igraph_t graph;
  assert(!igraph_lattice(&graph, &dim.vector(), nei, IGRAPH_UNDIRECTED,
                         0 /* mutual */, 0 /* circular */));
  return IGraph(std::move(graph));
}

}  // namespace common
}  // namespace cmplx

